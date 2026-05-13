"""Tests for GENCODE attribute aliases (#63) and version-column casting (#64)."""

import logging
import os
import tempfile

from gtfparse import GENCODE_BIOTYPE_ALIASES, INTEGER_VERSION_COLUMNS, read_gtf

from .data import data_path

GENCODE_GTF_PATH = data_path("gencode.head.gtf")


# -----------------------------
# attribute_aliases (#63)
# -----------------------------


def test_no_aliases_by_default_leaves_columns_alone():
    df = read_gtf(GENCODE_GTF_PATH, result_type="pandas")
    # gene_type / transcript_type come through as-is, gene_biotype / transcript_biotype absent
    assert "gene_type" in df.columns
    assert "transcript_type" in df.columns
    assert "gene_biotype" not in df.columns
    assert "transcript_biotype" not in df.columns


def test_attribute_aliases_renames_gencode_to_ensembl():
    df = read_gtf(
        GENCODE_GTF_PATH,
        attribute_aliases=GENCODE_BIOTYPE_ALIASES,
        result_type="pandas",
    )
    assert "gene_biotype" in df.columns
    assert "transcript_biotype" in df.columns
    assert "gene_type" not in df.columns
    assert "transcript_type" not in df.columns
    # values are preserved
    gene_rows = df[df["feature"] == "gene"]
    assert set(gene_rows["gene_biotype"]) == {
        "transcribed_unprocessed_pseudogene",
        "protein_coding",
    }


def test_attribute_aliases_with_polars_result_type():
    df = read_gtf(
        GENCODE_GTF_PATH,
        attribute_aliases=GENCODE_BIOTYPE_ALIASES,
        result_type="polars",
    )
    assert "gene_biotype" in df.columns
    assert "gene_type" not in df.columns


def test_attribute_aliases_when_canonical_already_present_prefers_canonical(caplog):
    """If both the alias and canonical column exist, drop the alias and warn."""
    # Hand-build a GTF where gene_biotype is already populated AND gene_type is also present
    contents = (
        "chr1\tHAVANA\tgene\t1\t100\t.\t+\t.\t"
        'gene_id "ENSGX"; gene_type "alias_value"; gene_biotype "canonical_value";\n'
    )
    fd, path = tempfile.mkstemp(suffix=".gtf")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(contents)
        with caplog.at_level(logging.WARNING):
            df = read_gtf(
                path,
                attribute_aliases={"gene_type": "gene_biotype"},
                result_type="pandas",
            )
        assert "gene_biotype" in df.columns
        assert "gene_type" not in df.columns
        # canonical wins
        assert list(df["gene_biotype"]) == ["canonical_value"]
        # warning was logged
        assert any("Both alias column" in rec.message for rec in caplog.records)
    finally:
        os.unlink(path)


def test_attribute_aliases_with_missing_alias_is_noop():
    """Aliases that don't appear in the file are silently skipped."""
    df = read_gtf(
        GENCODE_GTF_PATH,
        attribute_aliases={"never_in_file": "gene_biotype"},
        result_type="pandas",
    )
    # gene_biotype did not exist and the alias did not exist -> still no gene_biotype
    assert "gene_biotype" not in df.columns


def test_attribute_aliases_with_empty_dict_is_noop():
    df = read_gtf(GENCODE_GTF_PATH, attribute_aliases={}, result_type="pandas")
    assert "gene_type" in df.columns


# -----------------------------
# cast_version_columns (#64)
# -----------------------------


def test_version_columns_cast_to_int_by_default():
    df = read_gtf(GENCODE_GTF_PATH, result_type="pandas")
    # all four version columns should be present and Int64-typed (nullable)
    for column_name in INTEGER_VERSION_COLUMNS:
        assert column_name in df.columns, column_name
        assert str(df[column_name].dtype) == "Int64", (
            f"{column_name} dtype is {df[column_name].dtype}, expected Int64"
        )
    # spot-check actual values are integers, not strings
    gene_rows = df[df["feature"] == "gene"].reset_index(drop=True)
    assert gene_rows.loc[0, "gene_version"] == 5
    assert gene_rows.loc[1, "gene_version"] == 6


def test_version_columns_handle_missing_values():
    """Rows where a version attribute is absent should produce pd.NA, not raise."""
    df = read_gtf(GENCODE_GTF_PATH, result_type="pandas")
    # gene rows don't have transcript_version / protein_version / exon_version
    gene_rows = df[df["feature"] == "gene"]
    assert gene_rows["transcript_version"].isna().all()
    assert gene_rows["protein_version"].isna().all()
    assert gene_rows["exon_version"].isna().all()
    # gene_version IS populated on gene rows
    assert gene_rows["gene_version"].notna().all()


def test_cast_version_columns_false_keeps_strings():
    df = read_gtf(GENCODE_GTF_PATH, cast_version_columns=False, result_type="pandas")
    # column is still object/string when opted out
    assert df["gene_version"].dtype == object
    gene_rows = df[df["feature"] == "gene"].reset_index(drop=True)
    assert gene_rows.loc[0, "gene_version"] == "5"


def test_version_columns_present_in_polars_result():
    df = read_gtf(GENCODE_GTF_PATH, result_type="polars")
    for column_name in INTEGER_VERSION_COLUMNS:
        assert column_name in df.columns


def test_version_columns_missing_when_attribute_absent():
    """When the underlying GTF doesn't carry version attributes, the columns
    simply don't exist — casting must not raise."""
    # The classic Ensembl release-75 fixture has no *_version attributes
    df = read_gtf(data_path("ensembl_grch37.head.gtf"), result_type="pandas")
    for column_name in INTEGER_VERSION_COLUMNS:
        assert column_name not in df.columns


# -----------------------------
# Interaction between #63 and #64
# -----------------------------


def test_aliases_and_version_casting_work_together():
    df = read_gtf(
        GENCODE_GTF_PATH,
        attribute_aliases=GENCODE_BIOTYPE_ALIASES,
        result_type="pandas",
    )
    assert "gene_biotype" in df.columns
    assert "transcript_biotype" in df.columns
    assert str(df["gene_version"].dtype) == "Int64"
    transcript_rows = df[df["feature"] == "transcript"].reset_index(drop=True)
    assert transcript_rows.loc[0, "transcript_biotype"] == "processed_transcript"
    assert transcript_rows.loc[0, "transcript_version"] == 2


def test_aliases_visible_to_infer_biotype_column():
    """attribute_aliases is applied before infer_biotype_column, so an aliased
    gene_biotype must prevent the inference path from overwriting it."""
    df = read_gtf(
        GENCODE_GTF_PATH,
        attribute_aliases=GENCODE_BIOTYPE_ALIASES,
        infer_biotype_column=True,
        result_type="pandas",
    )
    # Should keep the GENCODE-derived biotypes, not infer from "source"
    gene_rows = df[df["feature"] == "gene"]
    assert "protein_coding" in set(gene_rows["gene_biotype"])
    assert "transcribed_unprocessed_pseudogene" in set(gene_rows["gene_biotype"])
    # source column is HAVANA, definitely not a biotype value
    assert "HAVANA" not in set(gene_rows["gene_biotype"])


# -----------------------------
# Sanity: constants are stable
# -----------------------------


def test_gencode_biotype_aliases_constant():
    assert GENCODE_BIOTYPE_ALIASES == {
        "gene_type": "gene_biotype",
        "transcript_type": "transcript_biotype",
    }


def test_integer_version_columns_constant():
    assert set(INTEGER_VERSION_COLUMNS) == {
        "gene_version",
        "transcript_version",
        "protein_version",
        "exon_version",
    }
