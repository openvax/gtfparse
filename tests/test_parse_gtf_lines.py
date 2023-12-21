from pytest import raises
from gtfparse import (
    parse_gtf,
    parse_gtf_and_expand_attributes,
    REQUIRED_COLUMNS,
    ParsingError
)
from io import StringIO 

gtf_text = """
# sample GTF data copied from:
# http://useast.ensembl.org/info/website/upload/gff.html?redirect=no
1\ttranscribed_unprocessed_pseudogene\tgene\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene";
1\tprocessed_transcript\ttranscript\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-002"; transcript_source "havana";
"""

def test_parse_gtf_lines_with_expand_attributes():
    df = parse_gtf_and_expand_attributes(StringIO(gtf_text))


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
    assert list(df.columns) ==  expected_columns
    assert list(df["seqname"]) == ["1", "1"]
    # convert to list for comparison since numerical columns may be NumPy arrays
    assert list(df["start"]) == [11869, 11869]
    assert list(df["end"]) == [14409, 14409]

    assert df["score"].is_null().all(), "Unexpected scores: %s" % (df["score"],)
    assert list(df["gene_id"]) == ["ENSG00000223972", "ENSG00000223972"]
    assert list(df["transcript_id"]) == ["", "ENST00000456328"]


def test_parse_gtf_lines_without_expand_attributes():
    df = parse_gtf(StringIO(gtf_text), split_attributes=False)

    # convert to list since Py3's dictionary keys are a distinct collection type
    assert list(df.columns) == REQUIRED_COLUMNS
    assert list(df["seqname"]) == ["1", "1"]
    # convert to list for comparison since numerical columns may be NumPy arrays
    assert list(df["start"]) == [11869, 11869]
    assert list(df["end"]) == [14409, 14409]
    assert df["score"].is_null().all(), "Unexpected scores: %s" % (df["score"],)
    assert len(df["attribute"]) == 2

def test_parse_gtf_lines_error_too_few_fields():
    bad_gtf_text = gtf_text.replace("\t", " ")
    # pylint: disable=no-value-for-parameter
    with raises(ParsingError):
        parse_gtf(StringIO(bad_gtf_text))
