# gtfparse

[![Tests](https://github.com/openvax/gtfparse/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/gtfparse/actions/workflows/tests.yml)
[![Coverage Status](https://coveralls.io/repos/openvax/gtfparse/badge.svg?branch=master&service=github)](https://coveralls.io/github/openvax/gtfparse?branch=master)
[![PyPI](https://img.shields.io/pypi/v/gtfparse.svg)](https://pypi.org/project/gtfparse/)

gtfparse reads GTF (Gene Transfer Format) annotation files, such as those from
Ensembl and GENCODE, into pandas DataFrames: one row per feature, with
attributes such as gene and transcript IDs expanded into columns. It can also
write a table back to GTF. In the [OpenVax](https://github.com/openvax) tools,
[PyEnsembl](https://github.com/openvax/pyensembl) uses it to read Ensembl and
custom annotations before indexing them.

## Install

```sh
python -m pip install gtfparse
```

Python 3.11 or later is required. gtfparse depends on pandas and PyArrow; for
optional Polars output, install `gtfparse[polars]`.

## Quick start

<a id="read-an-annotation-and-find-genes"></a><a id="example-usage"></a><a id="parsing-all-rows-of-a-gtf-file-into-a-pandas-dataframe"></a>This
small artificial GTF runs without a download. Replace the stream with your
annotation filename when working with a real dataset.

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
`gene_id`, `gene_name` and `transcript_id`. Coordinates stay one-based and
inclusive, as in the file.

## Documentation

The guides are Markdown files in
[`docs/`](https://github.com/openvax/gtfparse/tree/master/docs). In reading order:

1. <a id="getting-gene-fpkm-values-from-a-stringtie-gtf-file"></a>[Filter and write annotations](https://github.com/openvax/gtfparse/blob/master/docs/guides/annotations.md):
   select chromosomes and features, write GTF and read numeric attributes.
2. <a id="reporting-parsing-progress"></a>[Report parsing progress](https://github.com/openvax/gtfparse/blob/master/docs/guides/progress.md)
   with your own progress bar.
3. <a id="version-3-migration"></a>[Version 3 migration](https://github.com/openvax/gtfparse/blob/master/docs/guides/migration.md):
   pandas by default, and optional Polars output.
4. [API reference](https://github.com/openvax/gtfparse/blob/master/docs/reference.md)
   for every public function.

To build the documentation locally, run
`python -m pip install -r requirements-docs.txt` and `./docs.sh`.
`python scripts/check_docs_examples.py` checks the printed example output.
