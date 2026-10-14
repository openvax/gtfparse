# Version 3 migration

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
