from io import StringIO

from gtfparse import create_missing_features, parse_gtf_and_expand_attributes

# two lines from the Ensembl 54 human GTF containing only a stop_codon and
# exon features, but from which gene and transcript information could be
# inferred
GTF_TEXT = "\n".join(
    [
        "# seqname biotype feature start end score strand frame attribute",
        "".join(
            [
                """18\tprotein_coding\tstop_codon\t32630766\t32630768\t.\t-\t0\t""",
                """gene_id "ENSG00000134779"; transcript_id "ENST00000334295"; exon_number "7";"""
                """ gene_name "C18orf10";""",
                """ transcript_name "C18orf10-201";""",
            ]
        ),
        "".join(
            [
                """18\tprotein_coding\texon\t32663078\t32663157\t.\t+\t.\tgene_id "ENSG00000150477"; """,
                """transcript_id "ENST00000383055"; exon_number "1"; gene_name "KIAA1328"; """,
                """transcript_name "KIAA1328-202";""",
            ]
        ),
    ]
)


GTF_DATAFRAME = parse_gtf_and_expand_attributes(StringIO(GTF_TEXT))


def test_create_missing_features_identity():
    df_should_be_same = create_missing_features(GTF_DATAFRAME, {})
    assert len(GTF_DATAFRAME) == len(df_should_be_same), "GTF DataFrames should be same size"


def _check_expanded_dataframe(df):
    assert "gene" in set(df["feature"]), "Extended GTF should contain gene feature"
    assert "transcript" in set(df["feature"]), "Extended GTF should contain transcript feature"

    C18orf10_201_transcript_mask = (df["feature"] == "transcript") & (
        df["transcript_name"] == "C18orf10-201"
    )
    assert len(df[C18orf10_201_transcript_mask]) == 1, (
        "Expected only 1 gene entry for C18orf10-201, got %s" % (df[C18orf10_201_transcript_mask],)
    )
    transcript_seqname = df[C18orf10_201_transcript_mask].seqname.iloc[0]
    assert transcript_seqname == "18", "Wrong seqname for C18orf10-201: %s" % transcript_seqname
    transcript_start = df[C18orf10_201_transcript_mask].start.iloc[0]
    assert transcript_start == 32630766, "Wrong start for C18orf10-201: %s" % transcript_start
    transcript_end = df[C18orf10_201_transcript_mask].end.iloc[0]
    assert transcript_end == 32630768, "Wrong end for C18orf10-201: %s" % transcript_end
    transcript_strand = df[C18orf10_201_transcript_mask].strand.iloc[0]
    assert transcript_strand == "-", "Wrong strand for C18orf10-201: %s" % transcript_strand

    KIAA1328_gene_mask = (df["feature"] == "gene") & (df["gene_name"] == "KIAA1328")
    assert len(df[KIAA1328_gene_mask]) == 1, "Expected only 1 gene entry for KIAA1328"
    gene_seqname = df[KIAA1328_gene_mask].seqname.iloc[0]
    assert gene_seqname == "18", "Wrong seqname for KIAA1328: %s" % gene_seqname
    gene_start = df[KIAA1328_gene_mask].start.iloc[0]
    assert gene_start == 32663078, "Wrong start for KIAA1328: %s" % (gene_start,)
    gene_end = df[KIAA1328_gene_mask].end.iloc[0]
    assert gene_end == 32663157, "Wrong end for KIAA1328: %s" % (gene_end,)
    gene_strand = df[KIAA1328_gene_mask].strand.iloc[0]
    assert gene_strand == "+", "Wrong strand for KIAA1328: %s" % gene_strand


def test_create_missing_features():
    assert "gene" not in set(GTF_DATAFRAME["feature"]), (
        "Original GTF should not contain gene feature"
    )
    assert "transcript" not in set(GTF_DATAFRAME["feature"]), (
        "Original GTF should not contain transcript feature"
    )
    df_extra_features = create_missing_features(
        GTF_DATAFRAME,
        unique_keys={"gene": "gene_id", "transcript": "transcript_id"},
        extra_columns={
            "gene": {"gene_name"},
            "transcript": {"gene_id", "gene_name", "transcript_name"},
        },
    )
    _check_expanded_dataframe(df_extra_features)


# A gene row without a transcript_id, then exons of two transcripts. The gene
# row's empty transcript_id must never become a synthesized transcript.
GENE_AND_EXONS_GTF = "".join(
    "1\tsrc\t%s\t%d\t%d\t.\t+\t.\t%s\n" % row
    for row in [
        ("gene", 100, 500, 'gene_id "G1";'),
        ("exon", 100, 200, 'gene_id "G1"; transcript_id "T1";'),
        ("exon", 300, 500, 'gene_id "G1"; transcript_id "T1";'),
        ("exon", 150, 250, 'gene_id "G1"; transcript_id "T2";'),
    ]
)


def synthesized_transcript_spans(df):
    """Return sorted (transcript_id, start, end) tuples of df's transcript rows.

    Parameters
    ----------
    df : pandas.DataFrame
        Output of create_missing_features.

    Returns
    -------
    list of tuple
        One (transcript_id, start, end) tuple per transcript row, sorted.
    """
    transcripts = df[df["feature"] == "transcript"]
    return sorted(zip(transcripts["transcript_id"], transcripts["start"], transcripts["end"]))


def test_create_missing_features_on_filtered_frame():
    # Filtering leaves gaps in the index, which a positional row mask can't align to.
    df = parse_gtf_and_expand_attributes(StringIO(GENE_AND_EXONS_GTF))
    exons = df[df["feature"] == "exon"]
    assert list(exons.index) == [1, 2, 3]
    result = create_missing_features(exons, unique_keys={"transcript": "transcript_id"})
    assert synthesized_transcript_spans(result) == [("T1", 100, 500), ("T2", 150, 250)]


def test_create_missing_features_on_reverse_sorted_frame():
    # Reversing keeps every index label, so a positional mask would silently
    # select the wrong rows: the gene row instead of T2's exon.
    df = parse_gtf_and_expand_attributes(StringIO(GENE_AND_EXONS_GTF))
    reversed_df = df.iloc[::-1]
    assert list(reversed_df.index) == [3, 2, 1, 0]
    result = create_missing_features(reversed_df, unique_keys={"transcript": "transcript_id"})
    assert synthesized_transcript_spans(result) == [("T1", 100, 500), ("T2", 150, 250)]


def test_create_missing_features_skips_null_keys():
    df = parse_gtf_and_expand_attributes(StringIO(GENE_AND_EXONS_GTF))
    df.loc[df["transcript_id"] == "", "transcript_id"] = None
    assert df["transcript_id"].isna().sum() == 1
    result = create_missing_features(df, unique_keys={"transcript": "transcript_id"})
    assert synthesized_transcript_spans(result) == [("T1", 100, 500), ("T2", 150, 250)]
