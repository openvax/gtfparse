[![Tests](https://github.com/openvax/gtfparse/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/gtfparse/actions/workflows/tests.yml)
[![Coverage Status](https://coveralls.io/repos/openvax/gtfparse/badge.svg?branch=master&service=github)](https://coveralls.io/github/openvax/gtfparse?branch=master)
<a href="https://pypi.python.org/pypi/gtfparse/">
    <img src="https://img.shields.io/pypi/v/gtfparse.svg?maxAge=1000" alt="PyPI" />
</a>

gtfparse
========
Parsing tools for GTF (gene transfer format) files.

# Example usage

## Parsing all rows of a GTF file into a DataFrame

`read_gtf` returns a Polars DataFrame. Pass `result_type="pandas"` to get a
pandas DataFrame instead.

```python
import polars
from gtfparse import read_gtf

# returns GTF with essential columns such as "feature", "seqname", "start", "end"
# alongside the names of any optional keys which appeared in the attribute column
df = read_gtf("gene_annotations.gtf")

# filter DataFrame to gene entries on chrY
df_genes = df.filter(polars.col("feature") == "gene")
df_genes_chrY = df_genes.filter(polars.col("seqname") == "Y")

# all rows (gene, transcripts, exons, ...) for one gene
df_tp53 = df.filter(polars.col("gene_id") == "ENSG00000141510")
```


## Getting transcript FPKM values from a StringTie GTF file

```python
import polars
from gtfparse import read_gtf

df = read_gtf(
    "Transcripts.gtf",
    column_converters={"FPKM": float})

transcripts = df.filter(polars.col("feature") == "transcript")
transcript_fpkms = dict(zip(transcripts["transcript_id"], transcripts["FPKM"]))
```


