from gtftools import read_gtf_as_dict, read_gtf_as_dataframe

B16_GTF_PATH = "B16.stringtie.head.gtf"

def _check_required_columns(gtf_dict):
    assert "feature" in gtf_dict, "Expected column named 'feature' in StringTie GTF"
    assert "cov" in gtf_dict, "Expected column named 'cov' in StringTie GTF"
    assert "FPKM" in gtf_dict, "Expected column named 'FPKM' in StringTie GTF"
    features = set(gtf_dict["feature"])
    assert "exon" in features, "No exons in GTF (available: %s)" % features
    assert "transcript" in features, "No transcripts in GTF (available: %s)" % features

def _check_cov_and_FPKM(gtf_dict):
    for i, feature_name in enumerate(gtf_dict["feature"]):
        cov = gtf_dict["cov"][i]
        fpkm = gtf_dict["FPKM"][i]

        if feature_name == "exon":
            # coverage is a string, expect it to be non-empty for exons
            assert len(cov) > 0
            assert float(cov) >= 0
        elif feature_name == "transcript":
            assert len(fpkm) > 0
            assert float(fpkm) >= 0

def test_read_string_gtf_as_dict():
    gtf_dict = read_gtf_as_dict(B16_GTF_PATH)
    _check_required_columns(gtf_dict)
    _check_cov_and_FPKM(gtf_dict)

def test_read_stringtie_gtf_as_dataframe():
    gtf_df = read_gtf_as_dataframe(B16_GTF_PATH)
    _check_required_columns(gtf_df)
    _check_cov_and_FPKM(gtf_df)
