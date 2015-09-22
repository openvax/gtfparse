from gtfparse import read_gtf_as_dict, read_gtf_as_dataframe
import numpy as np
from data import data_path

REFSEQ_GTF_PATH = data_path("refseq.ucsc.small.gtf")

def _check_required_columns(gtf_dict):
    assert "feature" in gtf_dict, "Expected column named 'feature' in RefSeq GTF"
    assert "gene_id" in gtf_dict, "Expected column named 'gene_id' in RefSeq GTF"
    assert "transcript_id" in gtf_dict, "Expected column named 'transcript_id' in RefSeq GTF"
    features = set(gtf_dict["feature"])
    assert "exon" in features, "No exon features in GTF (available: %s)" % features
    assert "CDS" in features, "No CDS features in GTF (available: %s)" % features

def test_read_refseq_gtf_as_dict():
    gtf_dict = read_gtf_as_dict(REFSEQ_GTF_PATH)
    _check_required_columns(gtf_dict)

def test_read_refseq_gtf_as_dataframe():
    gtf_df = read_gtf_as_dataframe(REFSEQ_GTF_PATH)
    _check_required_columns(gtf_df)
