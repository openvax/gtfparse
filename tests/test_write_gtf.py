from gtfparse import read_gtf, write_gtf
from .data import data_path
from polars import DataFrame

REFSEQ_GTF_PATH = data_path("refseq.ucsc.small.gtf")


def test_write_gtf(tmp_path):
    expected_gtf = read_gtf(REFSEQ_GTF_PATH)
    write_gtf(expected_gtf, tmp_path/"dummy_gtf.gtf")
    created_gtf = read_gtf(str(tmp_path/"dummy_gtf.gtf"))
    assert  isinstance(created_gtf, DataFrame)
    assert expected_gtf == created_gtf

