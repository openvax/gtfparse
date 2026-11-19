# Version 3 migration

Read this if your code was written for gtfparse 2.x, whose functions returned
Polars DataFrames by default.

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
type.

Two 2.x compatibility wrappers are deprecated since 3.2.0. They emit a
`DeprecationWarning` and will be removed in 4.0:

- `gtfparse.read_gtf.parse_with_polars_lazy` reads the whole file eagerly
  despite its name and then returns a Polars LazyFrame. Use
  `read_gtf(..., result_type="polars").lazy()`, or
  `polars.from_pandas(parse_gtf(...)).lazy()` to keep the unexpanded
  `attribute` and `attribute_split` columns.
- `parse_gtf_pandas` returns the same pandas DataFrame as `parse_gtf`. Use
  `parse_gtf`.

Empty or malformed GTF fields now raise `gtfparse.ParsingError`, including wrong
field counts and invalid numeric fields. Whole comment lines are ignored;
literal `#` inside attributes is preserved. Text and binary streams stay open
after parsing, and gzip is detected from its contents.
