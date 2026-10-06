# Filter and write annotations

## Select a chromosome and feature

This standalone artificial input has two gene rows on different chromosomes:

```python
from io import StringIO
from gtfparse import read_gtf

annotation = StringIO(
    '1\texample\tgene\t101\t200\t.\t+\t.\tgene_id "g1"; gene_name "FIRST";\n'
    'Y\texample\tgene\t301\t400\t.\t+\t.\tgene_id "g2"; gene_name "SECOND";\n'
)
genes = read_gtf(annotation, features={"gene"})
selected = genes.loc[genes["seqname"] == "Y"]
print(selected["gene_id"].tolist())
```

```text
['g2']
```

Match the contig spelling in your file: `Y` and `chrY` are different strings.
The parser does not infer reference assembly compatibility from those names.

## Write a compressed GTF

Continuing with `genes` above, write and reread a temporary gzip file:

```python
from pathlib import Path
from tempfile import TemporaryDirectory
from gtfparse import write_gtf

with TemporaryDirectory() as directory:
    path = Path(directory) / "genes.gtf.gz"
    write_gtf(genes, path)
    restored = read_gtf(path)
    print(restored["gene_id"].tolist())
```

```text
['g1', 'g2']
```

The writer reconstructs the attribute column. It is intended to preserve
annotation values, rather than the original file's exact spacing and ordering.

## Read numeric attributes

GTF attributes initially represent text. Supply `column_converters` when an
attribute has a numeric meaning in your input file. This template requires a
file with gene rows carrying an FPKM attribute:

```python
from gtfparse import read_gtf

features = read_gtf("expression.gtf", column_converters={"FPKM": float})
genes = features.loc[features["feature"] == "gene"]
gene_fpkms = dict(zip(genes["gene_id"], genes["FPKM"]))
```

Use a gene identifier as the key, rather than `seqname`, which would collapse
all genes on a chromosome into one entry. For transcript-level expression,
select transcript rows and key by `transcript_id`. Parsing does not aggregate
transcript expression into gene expression, and a file without gene rows
will not produce a gene-level table from this filter.
