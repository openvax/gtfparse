from io import StringIO

import pandas as pd
import pytest

from gtfparse import attach_transcript_completeness, read_gtf, transcript_completeness

MIXED_COMPLETENESS_GTF = "".join(
    "1\tsrc\t%s\t%d\t%d\t.\t+\t.\t%s\n" % row
    for row in [
        # Gene row without transcript_id
        ("gene", 100, 900, 'gene_id "G1";'),
        # Complete coding transcript (has CDS + start_codon + stop_codon)
        ("transcript", 100, 300, 'gene_id "G1"; transcript_id "T_COMPLETE";'),
        ("exon", 100, 300, 'gene_id "G1"; transcript_id "T_COMPLETE";'),
        ("CDS", 120, 270, 'gene_id "G1"; transcript_id "T_COMPLETE";'),
        ("start_codon", 120, 122, 'gene_id "G1"; transcript_id "T_COMPLETE";'),
        ("stop_codon", 271, 273, 'gene_id "G1"; transcript_id "T_COMPLETE";'),
        # Coding fragment missing stop_codon (like TAIR chr_patch_hapl_scaff)
        ("transcript", 400, 550, 'gene_id "G1"; transcript_id "T_NO_STOP";'),
        ("CDS", 410, 540, 'gene_id "G1"; transcript_id "T_NO_STOP";'),
        ("start_codon", 410, 412, 'gene_id "G1"; transcript_id "T_NO_STOP";'),
        # Coding fragment missing both start and stop codons
        ("CDS", 600, 700, 'gene_id "G1"; transcript_id "T_NO_CODONS";'),
        # Non-coding transcript (no CDS -> complete)
        ("transcript", 750, 850, 'gene_id "G1"; transcript_id "T_NONCODING";'),
        ("exon", 750, 850, 'gene_id "G1"; transcript_id "T_NONCODING";'),
    ]
)


def test_transcript_completeness_companion_dataframe():
    df = read_gtf(StringIO(MIXED_COMPLETENESS_GTF))
    summary = transcript_completeness(df)

    assert summary.index.name == "transcript_id"
    assert list(summary.index) == ["T_COMPLETE", "T_NO_STOP", "T_NO_CODONS", "T_NONCODING"]
    assert list(summary.columns) == [
        "has_cds",
        "has_start_codon",
        "has_stop_codon",
        "transcript_complete",
    ]

    assert summary.loc["T_COMPLETE"].to_dict() == {
        "has_cds": True,
        "has_start_codon": True,
        "has_stop_codon": True,
        "transcript_complete": True,
    }
    assert summary.loc["T_NO_STOP"].to_dict() == {
        "has_cds": True,
        "has_start_codon": True,
        "has_stop_codon": False,
        "transcript_complete": False,
    }
    assert summary.loc["T_NO_CODONS"].to_dict() == {
        "has_cds": True,
        "has_start_codon": False,
        "has_stop_codon": False,
        "transcript_complete": False,
    }
    assert summary.loc["T_NONCODING"].to_dict() == {
        "has_cds": False,
        "has_start_codon": False,
        "has_stop_codon": False,
        "transcript_complete": True,
    }


def test_read_gtf_flag_incomplete_transcripts():
    default_df = read_gtf(StringIO(MIXED_COMPLETENESS_GTF))
    assert "transcript_complete" not in default_df.columns
    assert "has_start_codon" not in default_df.columns
    assert "has_stop_codon" not in default_df.columns

    flagged = read_gtf(StringIO(MIXED_COMPLETENESS_GTF), flag_incomplete_transcripts=True)
    assert {"has_start_codon", "has_stop_codon", "transcript_complete"}.issubset(flagged.columns)

    # Gene row without transcript_id has NA completeness flags
    gene_row = flagged.loc[flagged["feature"] == "gene"].iloc[0]
    assert pd.isna(gene_row["has_start_codon"])
    assert pd.isna(gene_row["has_stop_codon"])
    assert pd.isna(gene_row["transcript_complete"])

    # Every row for T_COMPLETE is complete
    complete_rows = flagged.loc[flagged["transcript_id"] == "T_COMPLETE"]
    assert complete_rows["has_start_codon"].all()
    assert complete_rows["has_stop_codon"].all()
    assert complete_rows["transcript_complete"].all()

    # Every row for T_NO_STOP has start_codon=True, stop_codon=False, complete=False
    no_stop_rows = flagged.loc[flagged["transcript_id"] == "T_NO_STOP"]
    assert no_stop_rows["has_start_codon"].all()
    assert not no_stop_rows["has_stop_codon"].any()
    assert not no_stop_rows["transcript_complete"].any()

    # Non-coding transcript is complete even without start/stop codons
    noncoding_rows = flagged.loc[flagged["transcript_id"] == "T_NONCODING"]
    assert not noncoding_rows["has_start_codon"].any()
    assert not noncoding_rows["has_stop_codon"].any()
    assert noncoding_rows["transcript_complete"].all()


def test_read_gtf_flag_incomplete_transcripts_with_usecols():
    flagged = read_gtf(
        StringIO(MIXED_COMPLETENESS_GTF),
        usecols=["seqname", "gene_id"],
        flag_incomplete_transcripts=True,
    )
    assert list(flagged.columns) == [
        "seqname",
        "gene_id",
        "has_start_codon",
        "has_stop_codon",
        "transcript_complete",
    ]
    assert int((flagged["transcript_complete"] == False).sum()) == 4  # noqa: E712


def test_read_gtf_flag_incomplete_transcripts_requires_expanded_attributes():
    with pytest.raises(ValueError, match="expand_attribute_column=True"):
        read_gtf(
            StringIO(MIXED_COMPLETENESS_GTF),
            expand_attribute_column=False,
            flag_incomplete_transcripts=True,
        )


def test_transcript_completeness_empty_and_no_transcript_id():
    gene_only_gtf = '1\tsrc\tgene\t100\t200\t.\t+\t.\tgene_id "G1";\n'
    df = read_gtf(StringIO(gene_only_gtf))
    summary = transcript_completeness(df)
    assert len(summary) == 0
    assert list(summary.columns) == [
        "has_cds",
        "has_start_codon",
        "has_stop_codon",
        "transcript_complete",
    ]
    attached = attach_transcript_completeness(df)
    assert pd.isna(attached["transcript_complete"].iloc[0])
