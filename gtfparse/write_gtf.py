import polars
from pathlib import Path
import typing as t


COMMONS_COL = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame']


def write_gtf(df: polars.DataFrame, export_path: str | Path, headers:  t.List[str] = None):
    headers = headers or []
    with open(export_path, 'w') as f:
        for header in headers:
            f.write(f"{header}\n")
        for row in df.iter_rows(named=True):
            f.write(f"{commons_cols(row)}\t{custom_fields(row)}\n")


def commons_cols(row) -> str :
    return "\t".join([str(row[field] or '.') for field in COMMONS_COL])


def custom_fields(row) -> str:
    return "; ".join([f'{field} "{row[field]}"' for field in row.keys() if (field not in COMMONS_COL) and (row[field])])
