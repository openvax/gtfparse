[![Tests](https://github.com/openvax/gtfparse/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/gtfparse/actions/workflows/tests.yml)
[![Coverage Status](https://coveralls.io/repos/openvax/gtfparse/badge.svg?branch=master&service=github)](https://coveralls.io/github/openvax/gtfparse?branch=master)
<a href="https://pypi.python.org/pypi/gtfparse/">
    <img src="https://img.shields.io/pypi/v/gtfparse.svg?maxAge=1000" alt="PyPI" />
</a>

gtfparse
========
Parsing tools for GTF (gene transfer format) files.

## Version 3 migration

`read_gtf`, `parse_gtf`, and `parse_gtf_and_expand_attributes` now return pandas
DataFrames by default. GTF reading uses PyArrow; attribute processing and writing
use pandas. Polars is no longer a required dependency and is imported only for
explicit Polars conversion. Ordinary installations need only pandas and PyArrow.

```sh
pip install gtfparse
# Optional, for applications that still request Polars frames:
pip install 'gtfparse[polars]'
```

```python
from gtfparse import read_gtf, write_gtf

df = read_gtf("gene_annotations.gtf")  # pandas.DataFrame
write_gtf(df, "gene_annotations.gtf.gz")

# Requires the optional extra:
polars_df = read_gtf("gene_annotations.gtf", result_type="polars")
```

Existing pandas and dictionary callers retain their output shape, fixed-column
dtypes, attribute handling, filters, converters and progress callbacks. Applications
that used the previous default Polars output should add `result_type="polars"`
and install the extra, or use pandas methods. `write_gtf` accepts either frame
type. The legacy `parse_with_polars_lazy` helper remains an optional adapter:
it reads the file eagerly and returns a Polars LazyFrame, as before.

Empty or malformed GTF fields now raise `gtfparse.ParsingError`, including wrong
field counts and invalid numeric fields. Whole comment lines are ignored;
literal `#` inside attributes is preserved. Text and binary streams stay open
after parsing, and gzip is detected from its contents.

## Reporting parsing progress

Pass `progress_callback(stage, completed, total)` to `read_gtf` to connect
parsing to your own progress UI. The callback runs synchronously; exceptions
from it stop parsing and propagate to the caller. Without a callback, gtfparse
does not display a progress bar or add a progress-library dependency.

| Stage | Reports |
| --- | --- |
| `read` | `(0, None)` before loading/filtering the file; `(rows, rows)` after it succeeds |
| `attributes` | `(0, rows)`, then every 10,000 rows, then `(rows, rows)` |
| `convert` | `(0, None)` before output conversion/transforms; `(rows, rows)` after success |

Counts refer to rows retained after `features` filtering, not bytes or overall
percentages. `None` means the total is unknown: loading and conversion expose
only stage boundaries, so use an indeterminate indicator during those steps.
Attribute expansion is omitted when `expand_attribute_column=False`. A stage
that fails does not emit completion. An empty filtered result is supported;
its attributes stage emits a single `(0, 0)` event. Empty input files raise `gtfparse.ParsingError`.

For example, with the optional `tqdm` package installed:

```python
from gtfparse import read_gtf
from tqdm.auto import tqdm

with tqdm(unit="rows") as bar:
    def show_progress(stage, completed, total):
        if stage != bar.desc:
            bar.total = total
            bar.reset()
            bar.set_description_str(stage)
        bar.total = total
        bar.update(completed - bar.n)
        bar.refresh()

    df = read_gtf("gene_annotations.gtf", progress_callback=show_progress)
```

`parse_gtf` reports the `read` stage, `parse_gtf_and_expand_attributes` reports
`read` and `attributes`, `parse_gtf_pandas` reports `read` and `convert`, and
`expand_attribute_strings` reports only `attributes` using the same callback.

Quoted attribute values preserve spaces, semicolons, and apostrophes. Raw
attributes (`expand_attribute_column=False`) retain their original quote marks.

# Example usage

## Parsing all rows of a GTF file into a Pandas DataFrame

```python
from gtfparse import read_gtf

# returns GTF with essential columns such as "feature", "seqname", "start", "end"
# alongside the names of any optional keys which appeared in the attribute column
df = read_gtf("gene_annotations.gtf")

# filter DataFrame to gene entries on chrY
df_genes = df[df["feature"] == "gene"]
df_genes_chrY = df_genes[df_genes["seqname"] == "Y"]
```


## Getting gene FPKM values from a StringTie GTF file

```python
from gtfparse import read_gtf

df = read_gtf(
    "Transcripts.gtf",
    column_converters={"FPKM": float})

gene_fpkms = {
    gene_name: fpkm
    for (gene_name, fpkm, feature)
    in zip(df["seqname"], df["FPKM"], df["feature"])
    if feature == "gene"
}
```
