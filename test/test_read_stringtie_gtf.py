from gtfparse import read_gtf_as_dict, read_gtf_as_dataframe
import numpy as np
from data import data_path

B16_GTF_PATH = data_path("B16.stringtie.head.gtf")

def _check_required_columns(gtf_dict):
    assert "feature" in gtf_dict, "Expected column named 'feature' in StringTie GTF"
    assert "cov" in gtf_dict, "Expected column named 'cov' in StringTie GTF"
    assert "FPKM" in gtf_dict, "Expected column named 'FPKM' in StringTie GTF"
    features = set(gtf_dict["feature"])
    assert "exon" in features, "No exons in GTF (available: %s)" % features
    assert "transcript" in features, "No transcripts in GTF (available: %s)" % features

def _check_string_cov_and_FPKM(gtf_dict):
    for i, feature_name in enumerate(gtf_dict["feature"]):
        cov = gtf_dict["cov"][i]
        fpkm = gtf_dict["FPKM"][i]
        if feature_name == "exon":
            assert len(fpkm) == 0, \
                "Expected missing FPKM for exon, got %s" % (fpkm,)
            assert  len(cov) > 0 and float(cov) >= 0, \
                "Expected non-negative cov for exon, got %s" % (cov,)
        elif feature_name == "transcript":
            assert len(cov) and float(cov) >= 0, \
                "Expected non-negative cov for transcript, got %s" % (cov,)
            assert len(fpkm) > 0 and float(fpkm) >= 0, \
                "Expected non-negative FPKM for transcript, got %s" % (fpkm,)

def _check_float_cov_and_FPKM(gtf_dict):
    for i, feature_name in enumerate(gtf_dict["feature"]):
        cov = gtf_dict["cov"][i]
        fpkm = gtf_dict["FPKM"][i]
        assert isinstance(cov, float), \
                "Expected cov to be float but got %s : %s" % (cov, type(cov))
        if feature_name == "exon":
            assert cov >= 0, "Expected non-negative cov for exon, got %s" % (cov,)
        elif feature_name == "transcript":
            assert isinstance(fpkm, float), \
                "Expected FPKM to be float but got %s : %s" % (fpkm, type(fpkm))
            assert cov >= 0, "Expected non-negative cov for transcript, got %s" % (cov,)
            assert fpkm >= 0, "Expected non-negative FPKM for transcript, got %s" % (fpkm,)

def test_read_string_gtf_as_dict():
    gtf_dict = read_gtf_as_dict(B16_GTF_PATH)
    _check_required_columns(gtf_dict)
    _check_string_cov_and_FPKM(gtf_dict)

def test_read_stringtie_gtf_as_dataframe():
    gtf_df = read_gtf_as_dataframe(B16_GTF_PATH)
    _check_required_columns(gtf_df)
    _check_string_cov_and_FPKM(gtf_df)

def test_read_string_gtf_as_dict_float_values():
    gtf_dict = read_gtf_as_dict(B16_GTF_PATH,
        column_converters={"cov": float, "FPKM": float})
    _check_required_columns(gtf_dict)
    _check_float_cov_and_FPKM(gtf_dict)

def test_read_stringtie_gtf_as_dataframe_float_values():
    gtf_df = read_gtf_as_dataframe(B16_GTF_PATH,
        column_converters={"cov": float, "FPKM": float})
    _check_required_columns(gtf_df)
    _check_float_cov_and_FPKM(gtf_df)

