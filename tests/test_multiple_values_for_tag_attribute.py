from io import StringIO
from gtfparse import parse_gtf_and_expand_attributes

# failing example from https://github.com/openvax/gtfparse/issues/2
GTF_TEXT = (
    """1\tprotein_coding\texon\t860260\t860328\t.\t+\t.\t"""
    """gene_id "ENSG00000187634"; transcript_id "ENST00000420190"; """
    """exon_number "1"; gene_name "SAMD11"; gene_source "ensembl_havana"; """
    """gene_biotype "protein_coding"; transcript_name "SAMD11-011"; """
    """transcript_source "havana"; exon_id "ENSE00001637883"; """
    """tag "cds_end_NF"; tag "mRNA_end_NF"; """
)

def test_parse_tag_attributes():
    parsed = parse_gtf_and_expand_attributes(StringIO(GTF_TEXT))
    tag_column = parsed["tag"]
    assert len(tag_column) == 1
    tags = tag_column[0]
    assert tags == 'cds_end_NF,mRNA_end_NF'

def test_parse_tag_attributes_with_usecols():
    parsed = parse_gtf_and_expand_attributes(
        StringIO(GTF_TEXT),
        restrict_attribute_columns=["tag"])
    tag_column = parsed["tag"]
    assert len(tag_column) == 1
    tags = tag_column[0]
    assert tags == 'cds_end_NF,mRNA_end_NF'

def test_parse_tag_attributes_with_usecols_other_column():
    parsed = parse_gtf_and_expand_attributes(
        StringIO(GTF_TEXT),
        restrict_attribute_columns=["exon_id"])

    assert "tag" not in parsed, "Expected 'tag' to get dropped but got %s" % (parsed,)
