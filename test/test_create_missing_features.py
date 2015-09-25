
import tempfile
from six.moves import StringIO
from gtfparse import (
    create_missing_features, parse_gtf_lines, read_gtf_as_dataframe
)
import pandas

# two lines from the Ensembl 54 human GTF containing only a stop_codon and
# exon features, but from which gene and transcript information could be
# inferred
GTF_DATA = """
# seqname biotype feature start end score strand frame attribute
18\tprotein_coding\tstop_codon\t32630766\t32630768\t.\t-\t0\tgene_id "ENSG00000134779"; transcript_id "ENST00000334295"; exon_number "7"; gene_name "C18orf10"; transcript_name "C18orf10-201";
18\tprotein_coding\texon\t32663078\t32663157\t.\t+\t.\tgene_id "ENSG00000150477"; transcript_id "ENST00000383055"; exon_number "1"; gene_name "KIAA1328"; transcript_name "KIAA1328-202";
"""

GTF_LINES = GTF_DATA.split("\n")

GTF_DICT = parse_gtf_lines(GTF_LINES)
GTF_DATAFRAME = pandas.DataFrame(GTF_DICT)

def test_create_missing_features_identity():
    df_should_be_same = create_missing_features(GTF_DATAFRAME, {})
    assert len(GTF_DATAFRAME) == len(df_should_be_same), \
        "GTF DataFrames should be same size"

def _check_expanded_dataframe(df):
    assert "gene" in set(df["feature"]), \
        "Extended GTF should contain gene feature"
    assert "transcript"  in set(df["feature"]), \
        "Extended GTF should contain transcript feature"

    C18orf10_201_transcript_mask = (
        (df["feature"] == "transcript") &
        (df["transcript_name"] == "C18orf10-201"))
    assert len(df[C18orf10_201_transcript_mask]) == 1, \
        "Expected only 1 gene entry for C18orf10-201, got %s" % (
            df[C18orf10_201_transcript_mask],)
    transcript_seqname = df[C18orf10_201_transcript_mask].seqname.irow(0)
    assert (transcript_seqname == "18"), \
        "Wrong seqname for C18orf10-201: %s" % transcript_seqname
    transcript_start = df[C18orf10_201_transcript_mask].start.irow(0)
    assert (transcript_start == 32630766), \
        "Wrong start for C18orf10-201: %s" % transcript_start
    transcript_end = df[C18orf10_201_transcript_mask].end.irow(0)
    assert (transcript_end == 32630768), \
        "Wrong end for C18orf10-201: %s" % transcript_end
    transcript_strand = df[C18orf10_201_transcript_mask].strand.irow(0)
    assert (transcript_strand == "-"), \
        "Wrong strand for C18orf10-201: %s" % transcript_strand

    KIAA1328_gene_mask = (
        (df["feature"] == "gene") &
        (df["gene_name"] == "KIAA1328"))
    assert len(df[KIAA1328_gene_mask]) == 1, "Expected only 1 gene entry for KIAA1328"
    gene_seqname = df[KIAA1328_gene_mask].seqname.irow(0)
    assert (gene_seqname == "18"), \
        "Wrong seqname for KIAA1328: %s" % gene_seqname
    gene_start = df[KIAA1328_gene_mask].start.irow(0)
    assert (gene_start == 32663078), \
        "Wrong start for KIAA1328: %s" % (gene_start,)
    gene_end = df[KIAA1328_gene_mask].end.irow(0)
    assert (gene_end == 32663157), \
        "Wrong end for KIAA1328: %s" % (gene_end,)
    gene_strand = df[KIAA1328_gene_mask].strand.irow(0)
    assert (gene_strand == "+"), \
        "Wrong strand for KIAA1328: %s" % gene_strand



def test_create_missing_features():
    assert "gene" not in set(GTF_DATAFRAME["feature"]), \
        "Original GTF should not contain gene feature"
    assert "transcript" not in set(GTF_DATAFRAME["feature"]), \
        "Original GTF should not contain transcript feature"
    df_extra_features = create_missing_features(
        GTF_DATAFRAME,
        unique_keys={
            "gene": "gene_id",
            "transcript": "transcript_id"
        },
        extra_columns={
            "gene": {"gene_name"},
            "transcript": {"gene_id", "gene_name", "transcript_name"},
        })
    _check_expanded_dataframe(df_extra_features)
