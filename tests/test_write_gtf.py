from gtfparse import read_gtf, write_gtf
from .data import data_path
from polars import DataFrame

REFSEQ_GTF_PATH = data_path("refseq.ucsc.small.gtf")


def test_write_gtf(tmp_path):
    gtf_dict = read_gtf(REFSEQ_GTF_PATH)
    write_gtf(gtf_dict, tmp_path/"dummy_gtf.gtf")
    assert  isinstance(read_gtf(str(tmp_path/"dummy_gtf.gtf")), DataFrame)
