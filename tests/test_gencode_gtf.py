"""Tests for GENCODE attribute aliases (#63) and version-column casting (#64)."""

import gzip
import logging
import os
import tempfile

from gtfparse import GENCODE_BIOTYPE_ALIASES, INTEGER_VERSION_COLUMNS, read_gtf

from .data import data_path

GENCODE_GTF_PATH = data_path("gencode.head.gtf")
GENCODE_REAL_GTF_PATH = data_path("gencode.real.head.gtf")


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
    import pandas as pd

    df = read_gtf(GENCODE_GTF_PATH, cast_version_columns=False, result_type="pandas")
    # column is still string-typed when opted out — accept either the legacy
    # object dtype or pandas' newer StringDtype, depending on the installed
    # pandas version.
    assert pd.api.types.is_string_dtype(df["gene_version"])
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


# -----------------------------
# Explicit-None defaults
# -----------------------------


def test_attribute_aliases_none_is_explicit_noop():
    """The kwarg default is None; passing it explicitly must match
    the implicit behavior (no rename, no error)."""
    df_default = read_gtf(GENCODE_GTF_PATH, result_type="pandas")
    df_explicit = read_gtf(GENCODE_GTF_PATH, attribute_aliases=None, result_type="pandas")
    assert list(df_default.columns) == list(df_explicit.columns)
    assert "gene_type" in df_explicit.columns
    assert "gene_biotype" not in df_explicit.columns


# -----------------------------
# result_type="dict" with new kwargs
# -----------------------------


def test_attribute_aliases_with_dict_result_type():
    result = read_gtf(
        GENCODE_GTF_PATH,
        attribute_aliases=GENCODE_BIOTYPE_ALIASES,
        result_type="dict",
    )
    assert isinstance(result, dict)
    assert "gene_biotype" in result
    assert "transcript_biotype" in result
    assert "gene_type" not in result
    assert "transcript_type" not in result


def test_version_columns_with_dict_result_type():
    result = read_gtf(GENCODE_GTF_PATH, result_type="dict")
    assert "gene_version" in result
    # values come back as ints when present, pd.NA when missing
    values = list(result["gene_version"].values())
    non_null = [v for v in values if v is not None and not _is_na(v)]
    assert non_null and all(isinstance(v, int) for v in non_null), non_null


def _is_na(v):
    # pd.NA isn't equal-comparable; use bool() guarded form
    try:
        return v != v  # NaN-style: pd.NA also returns NA from != itself
    except TypeError:
        return True


# -----------------------------
# Multi-alias and ordering edge cases
# -----------------------------


def test_multiple_aliases_mapping_to_same_canonical_first_wins(caplog):
    """Two aliases that both target the same canonical: the first in
    iteration order is renamed; the second is treated as a conflict
    and dropped. Prevents pandas producing two columns with the same
    name."""
    contents = (
        "chr1\tHAVANA\tgene\t1\t100\t.\t+\t.\t"
        'gene_id "ENSGX"; gene_type "alias_a"; gene_kind "alias_b";\n'
    )
    fd, path = tempfile.mkstemp(suffix=".gtf")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(contents)
        with caplog.at_level(logging.WARNING):
            df = read_gtf(
                path,
                # both aliases map to gene_biotype; gene_type comes first
                attribute_aliases={
                    "gene_type": "gene_biotype",
                    "gene_kind": "gene_biotype",
                },
                result_type="pandas",
            )
        # exactly one gene_biotype column, and the first alias's value wins
        assert list(df.columns).count("gene_biotype") == 1
        assert "gene_type" not in df.columns
        assert "gene_kind" not in df.columns
        assert list(df["gene_biotype"]) == ["alias_a"]
        # the collision was warned about
        assert any("Both alias column" in rec.message for rec in caplog.records)
    finally:
        os.unlink(path)


def test_aliases_iteration_order_matters_for_chains():
    """When canonical of one rule is the alias of another, behavior
    depends on iteration order against the original column set. This
    test pins the current behavior: each rule is decided against the
    pre-rename state, so chains aren't followed. Documents the contract."""
    contents = 'chr1\tHAVANA\tgene\t1\t100\t.\t+\t.\tgene_id "ENSGX"; a "1"; b "2";\n'
    fd, path = tempfile.mkstemp(suffix=".gtf")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(contents)
        # {"a": "b"} with b already in the file should drop a (both
        # present, canonical wins). NOT chain to {"a": "b", "b": "c"}.
        df = read_gtf(
            path,
            attribute_aliases={"a": "b"},
            result_type="pandas",
        )
        assert "a" not in df.columns
        assert "b" in df.columns
        assert list(df["b"]) == ["2"]
    finally:
        os.unlink(path)


def test_aliases_preserves_other_attribute_columns():
    """Renaming one attribute must not perturb the others."""
    df = read_gtf(
        GENCODE_GTF_PATH,
        attribute_aliases=GENCODE_BIOTYPE_ALIASES,
        result_type="pandas",
    )
    # spot-check several unrelated columns survive
    for col in ("gene_id", "transcript_id", "gene_name", "transcript_name"):
        assert col in df.columns, col


# -----------------------------
# usecols interaction
# -----------------------------


def test_aliases_with_usecols_filtering_canonical_column():
    """usecols is applied AFTER aliasing, so requesting the canonical
    column name returns the data even when the GTF only had the alias."""
    df = read_gtf(
        GENCODE_GTF_PATH,
        attribute_aliases=GENCODE_BIOTYPE_ALIASES,
        usecols=["gene_biotype"],
        result_type="pandas",
    )
    assert list(df.columns) == ["gene_biotype"]
    assert "protein_coding" in set(df["gene_biotype"])


# -----------------------------
# cast_version_columns edge cases
# -----------------------------


def test_version_zero_is_preserved_not_treated_as_missing():
    """0 must round-trip as Int64(0), not become pd.NA — the
    .replace("", None) call only targets empty strings."""
    contents = 'chr1\tHAVANA\tgene\t1\t100\t.\t+\t.\tgene_id "ENSGX"; gene_version "0";\n'
    fd, path = tempfile.mkstemp(suffix=".gtf")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(contents)
        df = read_gtf(path, result_type="pandas")
        assert df.loc[0, "gene_version"] == 0
        assert not df["gene_version"].isna().iloc[0]
    finally:
        os.unlink(path)


def test_version_with_corrupted_non_integer_value_coerces_to_na(caplog):
    """A malformed version string (e.g. 'v3' instead of '3') is
    coerced to pd.NA rather than raising. Documents the lenient
    behavior — corrupted-GTF data quality bugs become missing
    values, not exceptions."""
    contents = 'chr1\tHAVANA\tgene\t1\t100\t.\t+\t.\tgene_id "ENSGX"; gene_version "v3";\n'
    fd, path = tempfile.mkstemp(suffix=".gtf")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(contents)
        df = read_gtf(path, result_type="pandas")
        assert df["gene_version"].isna().iloc[0]
    finally:
        os.unlink(path)


def test_version_cast_is_idempotent_via_double_read():
    """Reading the same file twice with the cast on must produce the
    same dtypes — sanity check that the cast doesn't accumulate
    side effects via the global polars string cache or similar."""
    df1 = read_gtf(GENCODE_GTF_PATH, result_type="pandas")
    df2 = read_gtf(GENCODE_GTF_PATH, result_type="pandas")
    for col in INTEGER_VERSION_COLUMNS:
        assert str(df1[col].dtype) == str(df2[col].dtype) == "Int64"


# -----------------------------
# Real GENCODE-format fixture (versions in IDs, not in attributes)
# -----------------------------


def test_real_gencode_format_alias_rename():
    """A GENCODE GTF in its real on-disk shape — versioned IDs
    embedded in gene_id / transcript_id, GENCODE-only fields like
    level / hgnc_id / havana_gene / tag — should normalize to
    Ensembl-style biotype column names cleanly."""
    df = read_gtf(
        GENCODE_REAL_GTF_PATH,
        attribute_aliases=GENCODE_BIOTYPE_ALIASES,
        result_type="pandas",
    )
    assert "gene_biotype" in df.columns
    assert "transcript_biotype" in df.columns
    assert "gene_type" not in df.columns
    gene_rows = df[df["feature"] == "gene"]
    biotypes = set(gene_rows["gene_biotype"])
    assert biotypes == {"transcribed_unprocessed_pseudogene", "protein_coding"}


def test_real_gencode_format_has_no_separate_version_attrs():
    """Real GENCODE doesn't carry separate *_version attribute
    fields — versions are baked into gene_id/transcript_id strings.
    cast_version_columns must be a graceful no-op (not create empty
    columns, not raise)."""
    df = read_gtf(GENCODE_REAL_GTF_PATH, result_type="pandas")
    for col in INTEGER_VERSION_COLUMNS:
        assert col not in df.columns, (
            f"Real GENCODE has no '{col}' attribute; "
            f"version casting should not have synthesized one."
        )


def test_real_gencode_format_preserves_versioned_id_strings():
    """The version baked into gene_id / transcript_id / protein_id /
    exon_id must come through verbatim — those strings ARE the
    canonical identifier in GENCODE."""
    df = read_gtf(GENCODE_REAL_GTF_PATH, result_type="pandas")
    gene_ids = set(df["gene_id"])
    assert "ENSG00000223972.5" in gene_ids
    assert "ENSG00000186092.7" in gene_ids
    assert "ENSG00000198888.2" in gene_ids
    # protein_id and exon_id likewise carry .N
    protein_rows = df[df["feature"] == "CDS"]
    assert "ENSP00000493376.2" in set(protein_rows["protein_id"])


def test_real_gencode_format_preserves_gencode_only_fields():
    """GENCODE-specific attribute fields (level, hgnc_id,
    havana_gene, havana_transcript, tag, transcript_support_level)
    must survive parsing — they're useful downstream signal."""
    df = read_gtf(GENCODE_REAL_GTF_PATH, result_type="pandas")
    for col in (
        "level",
        "hgnc_id",
        "havana_gene",
        "tag",
        "transcript_support_level",
        "havana_transcript",
    ):
        assert col in df.columns, col
    # spot-check values
    transcript_rows = df[df["feature"] == "transcript"].reset_index(drop=True)
    assert "MANE_Select" in set(df["tag"])
    assert "HGNC:14825" in set(df["hgnc_id"])
    # transcript_support_level is "1" (string — we deliberately don't
    # auto-cast it since GENCODE allows non-numeric values like "NA")
    assert "1" in set(transcript_rows["transcript_support_level"])


def test_real_gencode_format_uses_chr_prefixed_seqnames():
    """GENCODE uses 'chr1' / 'chrM' while Ensembl uses bare '1' /
    'MT'. Confirming the seqname comes through unchanged so
    downstream tooling can normalize if needed."""
    df = read_gtf(GENCODE_REAL_GTF_PATH, result_type="pandas")
    seqnames = set(df["seqname"])
    assert seqnames <= {"chr1", "chrM"}
    assert "chr1" in seqnames
    assert "chrM" in seqnames


def test_real_gencode_format_via_gzip_compressed_input():
    """End-to-end: gzip a real-format GENCODE GTF, point read_gtf at
    the .gz, confirm aliasing + parsing both still work."""
    with open(GENCODE_REAL_GTF_PATH, "rb") as src:
        raw = src.read()
    fd, gz_path = tempfile.mkstemp(suffix=".gtf.gz")
    os.close(fd)
    try:
        with gzip.open(gz_path, "wb") as out:
            out.write(raw)
        df = read_gtf(
            gz_path,
            attribute_aliases=GENCODE_BIOTYPE_ALIASES,
            result_type="pandas",
        )
        assert "gene_biotype" in df.columns
        assert "ENSG00000186092.7" in set(df["gene_id"])
    finally:
        os.unlink(gz_path)


def test_real_gencode_format_with_infer_biotype_column_is_safe():
    """infer_biotype_column would normally fire if 'protein_coding'
    appears in the 'source' column. Real GENCODE's source is
    'HAVANA' / 'ENSEMBL', not 'protein_coding', so inference
    correctly stays silent. Regression guard against accidentally
    triggering it on GENCODE input."""
    df = read_gtf(
        GENCODE_REAL_GTF_PATH,
        attribute_aliases=GENCODE_BIOTYPE_ALIASES,
        infer_biotype_column=True,
        result_type="pandas",
    )
    # The aliased gene_biotype should be the authoritative one
    gene_rows = df[df["feature"] == "gene"]
    assert set(gene_rows["gene_biotype"]) == {
        "transcribed_unprocessed_pseudogene",
        "protein_coding",
    }


# -----------------------------
# Ensembl regression checks
# -----------------------------


def test_old_ensembl_fixture_still_uses_source_for_biotype_inference():
    """Ensembl release-75 puts biotype into the 'source' column. With
    infer_biotype_column=True the fix should still populate
    gene_biotype / transcript_biotype from source. Make sure the new
    kwargs don't break that legacy path."""
    df = read_gtf(
        data_path("ensembl_grch37.head.gtf"),
        infer_biotype_column=True,
        result_type="pandas",
    )
    assert "gene_biotype" in df.columns
    assert "transcript_biotype" in df.columns
    # The release-75 fixture has pseudogenes and protein_coding genes
    assert "protein_coding" in set(df["gene_biotype"])


def test_ensembl_attribute_aliases_no_op_when_aliases_absent():
    """Passing GENCODE_BIOTYPE_ALIASES against a pure Ensembl file
    that doesn't have gene_type / transcript_type must do nothing
    — no spurious columns, no exceptions."""
    df = read_gtf(
        data_path("ensembl_grch37.head.gtf"),
        attribute_aliases=GENCODE_BIOTYPE_ALIASES,
        result_type="pandas",
    )
    # gene_biotype exists in this fixture; aliasing didn't fight with it
    assert "gene_biotype" in df.columns
    assert "gene_type" not in df.columns
    assert "transcript_type" not in df.columns
