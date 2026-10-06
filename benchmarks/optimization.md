# Complete-load optimization in 3.0.2

Version 3.0.2 retains the pandas default and optional Polars output introduced
in 3.0.0. It improves the shared Python attribute parser and the construction
of expanded pandas columns while keeping the supported reader behavior.

The benchmark candidate was labeled 3.0.1 when measured. The concurrent
documentation PR targets that version, so this independent optimization
ships as 3.0.2. Raw records retain their measured version labels;
result tables refer to the release containing that implementation.

## Changes

- Capture attribute names and double-quoted values directly in the compiled
  regex. Other quoted/unquoted values and pre-split sequences keep their
  permissive cleanup and custom quote handling. Avoid splitting and slicing
  every ordinary double-quoted pair and looking up each name twice.
- Convert raw attributes with missing entries mapped to None. For configured
  Arrow strings, expand and construct columns in batches of 250,000 rows, then
  retain compact Arrow chunks. Release each batch and temporary expanded list
  after conversion. Fill missing prefixes/gaps with empty strings and preserve
  first-seen column order, fixed-field collisions and global progress counts.
  This bounds temporary Python values for large files, addressing #91. Other
  string representations keep their compatible full-list construction.
- When pandas' configured inferred string dtype uses Arrow storage, construct
  the Arrow string array directly and wrap it in that same pandas dtype. This
  avoids an intermediate object-array inference pass. Other configured string
  representations retain their ordinary pandas construction behavior.

Column order, collisions with fixed fields, repeated and empty attributes,
custom quote/missing-value parameters, callbacks/cancellation, aliases, casts,
raw mode, filtering, gzip and caller-owned streams keep their previous behavior.
No process-wide pandas settings change. Numeric/categorical fixed-field reading
and the nine-field validation remain the same.

## Reproduction and method

Use the same benchmark interpreter/dependency versions for all release roots.
The recorded subsets use Python 3.12.6, pandas 2.3.3 or 3.0.6, NumPy 2.5.3,
Arrow 25.0.1 and optional Polars 1.44.2 with four threads. Extract the old package
sources into isolated directories, without changing their dependency versions:

```sh
mkdir -p /tmp/gtfparse-2.9.1 /tmp/gtfparse-3.0.0
git archive 3e0912c gtfparse | tar -x -C /tmp/gtfparse-2.9.1
git archive be449cb gtfparse | tar -x -C /tmp/gtfparse-3.0.0
python benchmarks/compare_releases.py \
  ensembl-250k.gtf ensembl-250k.gtf.gz gencode-250k.gtf gencode-250k.gtf.gz \
  --baseline 2.9.1-polars /tmp/gtfparse-2.9.1 \
  --baseline 3.0.0-arrow /tmp/gtfparse-3.0.0 \
  --python /path/to/pandas-2.3-env/bin/python \
  --python /path/to/pandas-3.0-env/bin/python \
  --repeats 3 --threads 4 --output comparisons.json
```

Install matching dependencies, including optional Polars for the historical
reader, in both environments. Each child starts a fresh process and runs the
complete production `read_gtf`; processes execute sequentially. Reader order
rotates and corpus order reverses between repetitions. Imports are outside the
reported load time. Peak RSS includes imported libraries and is captured before
content hashing. No OS-cache flushing or machine isolation is claimed.

The first 250,000 annotation rows and original header of Ensembl release 114
human GRCh38 and GENCODE v48 human primary assembly are the same inputs used in
[the migration benchmark](README.md#corpus-and-equivalence). Original URLs,
selection rules and hashes are retained there and in results.json. Complete
files provide separate scaling checks. Every run checks all column contents
and order against the other releases; different representations are normalized
only for hashing. Exact dtypes are checked separately against 3.0.0.

Use `--features exon --usecols seqname start end gene_id transcript_id` for
selective-load comparisons and `--result-type polars` for optional output.
The JSON records contain versions, content hashes, dimensions, stage times,
median/min/max and pre-validation peak RSS. A JSON-lines sidecar is flushed
on each successful read, so failed experiments retain their completed records.
Compare complete times; the old and current read-stage callbacks cover different
conversion work and do not establish a pure CSV-engine speed ratio.

## Results

The table contains three-run medians of complete loads with equivalent pandas
output. All 72 runs match content hashes and column order across releases and
pandas versions. Times and memory are fresh measurements from this sweep;
compare readers within this table rather than mixing earlier timings from a
shared workstation.

| pandas | Input (250k rows) | 2.9.1 Polars s | 3.0.0 Arrow s | 3.0.2 s | Polars MiB | 3.0.0 MiB | 3.0.2 MiB |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2.3.3 | ensembl-250k.gtf | 2.57 | 2.91 | 2.29 | 560 | 617 | 583 |
| 2.3.3 | ensembl-250k.gtf.gz | 2.86 | 4.58 | 2.36 | 561 | 620 | 584 |
| 2.3.3 | gencode-250k.gtf | 3.18 | 2.59 | 3.72 | 567 | 588 | 550 |
| 2.3.3 | gencode-250k.gtf.gz | 3.61 | 3.90 | 3.74 | 568 | 584 | 558 |
| 3.0.6 | ensembl-250k.gtf | 2.54 | 2.82 | 1.91 | 552 | 600 | 577 |
| 3.0.6 | ensembl-250k.gtf.gz | 2.35 | 2.80 | 1.88 | 553 | 600 | 574 |
| 3.0.6 | gencode-250k.gtf | 2.05 | 2.68 | 1.71 | 560 | 593 | 577 |
| 3.0.6 | gencode-250k.gtf.gz | 2.10 | 2.52 | 1.73 | 560 | 595 | 574 |

On pandas 3.0.6, the optimized reader takes 17–25% less time than the Polars
baseline and 31–36% less time than 3.0.0 across all four cases. On pandas 2.3.3,
the Ensembl medians improve, while GENCODE wall times are mixed and highly
variable (individual complete reads span about 2–5.6 seconds). All retained
CPU-time medians improve, but that does not establish a wall-time speedup in
those noisy cases. No runs are discarded. Subset peak RSS is 3–6% lower than
3.0.0 and ranges from about 3% lower to 4% higher than Polars.

The repeated corpus is a prefix of each annotation, so feature/attribute
distributions can differ from complete files. This demonstrates results for
gtfparse's complete pipeline, rather than a universal pandas-versus-Polars
engine comparison. Batched assembly bounds temporary Python attribute values
when pandas uses Arrow strings; it does not make the entire reader streaming.
Python/object string storage remains supported without this memory reduction.

The full-file checks use pandas 3.0.6 and the same versions/thread count, with
one fresh process per reader. All six complete loads match content and column
order. These single observations demonstrate scaling and correctness, and do
not establish repeated full-file throughput.

| Full gzip input | Rows | Polars s | 3.0.0 s | 3.0.2 s | Polars MiB | 3.0.0 MiB | 3.0.2 MiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Ensembl GRCh38.114 | 4,116,048 | 47.11 | 67.63 | 36.86 | 6504 | 7199 | 3639 |
| GENCODE v48 primary | 4,119,244 | 43.74 | 52.35 | 29.88 | 7091 | 6560 | 3698 |

Here the optimized reader takes 22–32% less time and 44–48% less peak RSS than
Polars. Compared with 3.0.0, it takes 43–45% less time and 44–49% less peak RSS.
The shared workstation was running other workloads; RSS and wall times vary
between sweeps. Compare the matched observations above, not runs from different
measurement windows. The benefit applies to configured Arrow string storage.

Selective reads and explicit Polars output match all contents in separate
one-pass comparisons; these checks do not establish repeated throughput for
those configurations. [optimization-results.json](optimization-results.json)
retains every final record, version, input hash, stage median and range.
The rejected unbatched full-file candidate is retained separately: its full
GENCODE peak RSS reached 5.2 GiB and 5.8 GiB on a repeat. That small-sample
optimization was replaced with bounded temporary attribute assembly before
this release. Its timings are excluded from the final result tables.

## Correctness and compatibility

The differential tests freeze the pre-optimization tokenization and compare
quoted/unquoted, custom quote delimiters, Unicode/whitespace, malformed fragments,
repeated/empty values, pre-split rows, missing-value sentinels and column
selection. Additional tests check configured Python/Arrow string storage, batch
boundaries, late-discovered columns, missing batches, fixed-field attribute
collisions, first-seen order and global progress counts. The batch regression
opts pandas 2 into inferred strings so both configured storage backends are
exercised there; default object-string behavior has separate coverage. Existing tests cover streams, null fields, exact
integer coordinates, comments/BOM, aliases, converters, progress/cancellation,
optional Polars, dictionary output and writer round trips.

The verification sweep compares 71 full API/dtype combinations on each pandas
version against 3.0.0, including all five real fixtures and configured string
storage. CI exercises Python 3.9+, minimum/current pandas and Arrow, the base
installation without Polars, optional Polars versions and packaged source tests.
