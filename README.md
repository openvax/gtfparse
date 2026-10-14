[![Tests](https://github.com/openvax/gtfparse/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/gtfparse/actions/workflows/tests.yml)
[![Coverage Status](https://coveralls.io/repos/openvax/gtfparse/badge.svg?branch=master&service=github)](https://coveralls.io/github/openvax/gtfparse?branch=master)
<a href="https://pypi.python.org/pypi/gtfparse/">
    <img src="https://img.shields.io/pypi/v/gtfparse.svg?maxAge=1000" alt="PyPI" />
</a>

gtfparse
========
Parsing tools for GTF (gene transfer format) files.


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


<a id="version-3-migration"></a>
<a id="reporting-parsing-progress"></a>
<a id="example-usage"></a>
<a id="parsing-all-rows-of-a-gtf-file-into-a-pandas-dataframe"></a>
<a id="getting-gene-fpkm-values-from-a-stringtie-gtf-file"></a>

## Documentation

- [Filter features, select chromosomes, write GTF and read numeric attributes](docs/guides/annotations.md)
- [Progress callbacks](docs/guides/progress.md)
- [Version 3 migration and optional Polars output](docs/guides/migration.md)
- [Complete API reference](docs/reference.md)

Build with `python -m pip install -r requirements-docs.txt` and `./docs.sh`.
Run `python scripts/check_docs_examples.py` to verify the printed examples.
