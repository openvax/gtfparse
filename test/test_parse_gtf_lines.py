import numpy as np
from nose.tools import eq_, assert_raises
from gtfparse import parse_gtf_lines, REQUIRED_COLUMNS, ParsingError

gtf_lines = """
# sample GTF data copied from:
# http://useast.ensembl.org/info/website/upload/gff.html?redirect=no
1\ttranscribed_unprocessed_pseudogene\tgene\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene";
1\tprocessed_transcript\ttranscript\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-002"; transcript_source "havana";
""".split("\n")

def test_parse_gtf_lines_with_expand_attributes():
    parsed_dict = parse_gtf_lines(gtf_lines, expand_attribute_column=True)
    # excluding 'attribute' column from required names
    expected_columns = REQUIRED_COLUMNS[:8] + [
        "gene_id",
        "gene_name",
        "gene_source",
        "gene_biotype",
        "transcript_id",
        "transcript_name",
        "transcript_source",
    ]
    # convert to list since Py3's dictionary keys are a distinct collection type
    eq_(list(parsed_dict.keys()), expected_columns)
    eq_(parsed_dict["seqname"], ["1", "1"])
    # convert to list for comparison since numerical columns may be NumPy arrays
    eq_(list(parsed_dict["start"]), [11869, 11869])
    eq_(list(parsed_dict["end"] ), [14409, 14409])
    # can't compare NaN with equality
    scores = list(parsed_dict["score"])
    assert np.isnan(scores).all() , "Unexpected scores: %s" % scores
    eq_(parsed_dict["gene_id"], ["ENSG00000223972", "ENSG00000223972"])
    eq_(parsed_dict["transcript_id"], ["", "ENST00000456328"])


def test_parse_gtf_lines_without_expand_attributes():
    parsed_dict = parse_gtf_lines(gtf_lines, expand_attribute_column=False)

    # convert to list since Py3's dictionary keys are a distinct collection type
    eq_(list(parsed_dict.keys()), REQUIRED_COLUMNS)
    eq_(parsed_dict["seqname"], ["1", "1"])
    # convert to list for comparison since numerical columns may be NumPy arrays
    eq_(list(parsed_dict["start"]), [11869, 11869])
    eq_(list(parsed_dict["end"] ), [14409, 14409])
    # can't compare NaN with equality
    scores = list(parsed_dict["score"])
    assert np.isnan(scores).all() , "Unexpected scores: %s" % scores
    assert len(parsed_dict["attribute"]) == 2

def test_parse_gtf_lines_error_too_many_fields():
    bad_gtf_lines = [line.replace(" ", "\t") for line in gtf_lines]
    # pylint: disable=no-value-for-parameter
    with assert_raises(ParsingError):
        parse_gtf_lines(bad_gtf_lines)

def test_parse_gtf_lines_error_too_few_fields():
    bad_gtf_lines = [line.replace("\t", " ") for line in gtf_lines]
    # pylint: disable=no-value-for-parameter
    with assert_raises(ParsingError):
        parse_gtf_lines(bad_gtf_lines)
