# gtfparse

Read GTF annotations into a pandas DataFrame, with one row per feature and
attribute values such as gene and transcript IDs expanded into columns.

## Install

```sh
python -m pip install gtfparse
```

Python 3.9 or later is required. Version 3 uses pandas and PyArrow; Polars
output is optional and requires `gtfparse[polars]`.

## Read an annotation and find genes

This small artificial GTF runs without a download. Replace the stream with
your annotation filename when working with a real dataset.

```python
from io import StringIO
from gtfparse import read_gtf

annotation = StringIO(
    '1\texample\tgene\t101\t200\t.\t+\t.\tgene_id "g1"; gene_name "EXAMPLE";\n'
    '1\texample\texon\t101\t150\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
)
features = read_gtf(annotation)
genes = features.loc[features["feature"] == "gene"]
print(genes["gene_name"].tolist())
for row in features.itertuples():
    print(row.feature, row.seqname, row.start, row.end, row.gene_id)
```

```text
['EXAMPLE']
gene 1 101 200 g1
exon 1 101 150 g1
```

`feature` identifies the row type. `seqname` is the chromosome or contig name;
it is not the gene ID. The ninth GTF field supplies attributes such as
`gene_id`, `gene_name` and `transcript_id`.

GTF coordinates are one-based and inclusive: 101–150 spans 50 bases. Parsing
preserves this convention rather than converting it to BED coordinates. See
the [GENCODE format description](https://www.gencodegenes.org/pages/data_format.html).

## Use your files

This template requires an existing annotation file:

```python
from gtfparse import read_gtf

features = read_gtf("annotations.gtf.gz", features={"gene", "transcript"})
```

The `features` option filters row types while reading. Compressed inputs are
detected from their contents. Attribute expansion does not turn every row into
a gene or automatically supply missing parent features.

Continue with [filter and write annotations](guides/annotations.md), then
[progress callbacks](guides/progress.md) or [version 3 migration](guides/migration.md).
The [API reference](reference.md) documents the public parsing functions.
