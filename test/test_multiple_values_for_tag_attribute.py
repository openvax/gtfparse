from gtfparse import parse_gtf_lines
from nose.tools import eq_

# failing example from https://github.com/openvax/gtfparse/issues/2
GTF_LINES = [
    """1\tprotein_coding\texon\t860260\t860328\t.\t+\t.\t"""
    """gene_id "ENSG00000187634"; transcript_id "ENST00000420190"; """
    """exon_number "1"; gene_name "SAMD11"; gene_source "ensembl_havana"; """
    """gene_biotype "protein_coding"; transcript_name "SAMD11-011"; """
    """transcript_source "havana"; exon_id "ENSE00001637883"; """
    """tag "cds_end_NF"; tag "mRNA_end_NF"; """
]

def test_parse_tag_attributes():
    parsed = parse_gtf_lines(GTF_LINES)
    tag_column = parsed["tag"]
    eq_(len(tag_column), 1)
    tags = tag_column[0]
    eq_(tags, 'cds_end_NF,mRNA_end_NF')
