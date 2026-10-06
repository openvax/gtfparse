# GTF reader backend comparison

The 3.0.0 migration uses PyArrow for fixed-column parsing and pandas for attribute
processing and writing. Polars is an optional output adapter. This directory
contains the measurement harness and the evidence for that decision.

The follow-up [3.0.1 optimization report](optimization.md) compares the optimized
complete reader with both the 2.9.1 Polars and 3.0.0 Arrow/pandas releases.

## Reproduction

The recorded runs use Python 3.12.6, NumPy 2.5.3, PyArrow 25.0.1, Polars 1.44.2,
and pandas 2.3.3 or 3.0.6. Install the checkout with `pip install '.[dev,polars]'`
in isolated environments. Run one fresh process at a time with four threads:

```sh
POLARS_MAX_THREADS=4 python benchmarks/compare_readers.py sample.gtf --backend production
POLARS_MAX_THREADS=4 python benchmarks/compare_readers.py sample.gtf --backend pandas
POLARS_MAX_THREADS=4 python benchmarks/compare_readers.py sample.gtf --backend arrow
POLARS_MAX_THREADS=4 python benchmarks/compare_readers.py sample.gtf --backend arrow-extension
```

For the old production baseline, extract commit `3e0912c` (gtfparse 2.9.1) into
a separate directory and use the same harness with `--library-root /path/to/baseline`.
Repeat each combination three times, interleaving readers and reversing their
order for the second repetition. Imports and validation are outside wall/CPU
timing. Peak RSS is the process high-water mark, including imports, captured
before checksumming or constructing a reference frame. OS caches are not flushed.
The shared workstation and short runs limit precise throughput claims.

Each invocation emits JSON with wall/CPU time, stage timings, pre-validation
peak RSS, versions, dimensions, column order and a content checksum. Hashing
normalizes equivalent string/categorical representations and ordinary numeric
representations. The corpus coordinates are below float64's exact-integer limit;
separate regression tests verify exact parsing of larger int64 coordinates.
Use `--verify` for DataFrame equality against production on small inputs.
`--features exon`, `--usecols seqname start end gene_id transcript_id`, and
`--result-type polars` exercise selective reads and optional conversion.

## Corpus and equivalence

- [Ensembl release 114, human GRCh38](https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz)
- [GENCODE v48, human primary assembly](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.primary_assembly.annotation.gtf.gz)

The repeated sweep uses the first 250,000 annotation rows of each source with
its original header, in plain and gzip form. The full GENCODE gzip file provides
a separate scaling check on 4,119,244 rows. Original GTFs are not committed.
Input byte sizes, SHA-256 hashes, source URLs and selection rules are in
[results.json](results.json), alongside all 122 retained measurement records,
three-run medians/ranges and stage medians. Discarded production measurements
from before the final filtering changes are excluded.

All measured readers produce matching hashes, including the full-file pair.
The five real repository fixtures also pass 50 old/new API comparisons and
30 prototype comparisons on each pandas version (160 checks). The API checks
include raw attributes, filtering, selection, aliases, version casts, biotype
inference and converters; separate regressions cover streams, nullable fields,
comments, malformed inputs and progress/cancellation semantics.

## 3.0.0 migration results

The table shows three-run medians for equivalent pandas output. Complete load
time includes fixed-column reading, attribute expansion, assembly and transforms.

| pandas | Input (250k rows) | Old seconds | New seconds | Old peak MiB | New peak MiB |
| --- | --- | ---: | ---: | ---: | ---: |
| 2.3.3 | ensembl-250k.gtf | 2.51 | 2.75 | 560 | 614 |
| 2.3.3 | ensembl-250k.gtf.gz | 2.54 | 2.79 | 561 | 613 |
| 2.3.3 | gencode-250k.gtf | 2.08 | 2.42 | 567 | 581 |
| 2.3.3 | gencode-250k.gtf.gz | 2.16 | 2.42 | 568 | 581 |
| 3.0.6 | ensembl-250k.gtf | 2.28 | 2.69 | 552 | 601 |
| 3.0.6 | ensembl-250k.gtf.gz | 2.27 | 2.69 | 553 | 603 |
| 3.0.6 | gencode-250k.gtf | 2.01 | 2.41 | 560 | 597 |
| 3.0.6 | gencode-250k.gtf.gz | 2.03 | 2.44 | 560 | 592 |

Across these subsets, the new reader takes 1.10–1.20x the old total time and
1.03–1.10x the peak RSS. The single full GENCODE check with pandas 2.3.3 takes
38.88s / 6189 MiB for the baseline and 45.40s / 5348 MiB for version 3: matching
content, 1.17x time, and lower peak memory. This pair verifies scaling and
correctness; one observation per reader does not establish stable throughput.

The complete-load result meets #69's proposed 1.5x total-time threshold. The
callback `read` stages are different contracts for timing: version 3 includes
stream preprocessing and Arrow-to-pandas conversion; the old reader constructs
Polars frames and converts to pandas later. The reported read-stage times do
not establish the proposed 3x raw-CSV threshold. This migration is justified by
the measured complete-load cost and the user's preference to remove the required
Polars dependency, rather than a claim that Arrow is faster at raw CSV parsing.
Python attribute expansion remains the largest stage (roughly 1.7–1.8s on these
subsets). Eliminating it or changing its semantics is outside this change.

## Candidate evaluation and limits

The pandas C prototype is slower end to end on this corpus. The direct Arrow
and Arrow-extension prototypes have similar complete-load costs; extension
scalars introduce iteration/conversion overhead, and nullable dictionary casts
need an object bridge on pandas 2.3.3. These costs are included in the results.
The conventional Arrow-to-pandas path avoids that bridge and preserves the
established fixed-column and attribute output types.

The three prototypes are measurement tools. They support the paths, gzip,
filters, column selection and result types used here, but not the entire public
API. They omit buffers, custom converters, aliases, biotype inference, quote
cleanup, raw mode and complete progress/failure behavior. The pandas prototype
skips leading comments; the Arrow prototypes skip comment rows only on a field
count mismatch. The production reader implements whole-line comment filtering,
BOM and stream handling, exact field counts, declared numeric types and filtering
before pandas conversion so discarded null coordinates do not change dtypes.

Pandas reader behavior is documented in
[read_csv](https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html),
and Arrow's typed CSV options in its
[CSV guide](https://arrow.apache.org/docs/python/csv.html).
