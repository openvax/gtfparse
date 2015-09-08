# gtftools
Parsing tools for GTF (gene transfer format) files.

# Example usage

## Parsing all rows of a GTF file into a Pandas DataFrame

```python
from gtf_tools import read_gtf_as_dataframe

# returns GTF with essential columns such as "feature", "seqname", "start", "end"
# alongside the names of any optional keys which appeared in the attribute column
df = read_gtf_as_dataframe("gene_annotations.gtf")

# filter DataFrame to gene entries on chrY
df_genes = df[df["feature"] == "gene"]
df_genes_chrY = df_genes[df_genes["seqname"] == "Y"]
```


