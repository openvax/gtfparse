"""Regression tests for #56: read_gtf(expand_attribute_column=False)
used to raise NameError because the else branch referenced `result_df`
before it had been assigned.
"""

import pandas as pd

from gtfparse import read_gtf

from .data import data_path

GTF_PATH = data_path("ensembl_grch37.head.gtf")


def test_expand_attribute_column_false_returns_raw_attribute_pandas():
    df = read_gtf(GTF_PATH, expand_attribute_column=False, result_type="pandas")
    assert isinstance(df, pd.DataFrame)
    # raw attribute column is preserved verbatim
    assert "attribute" in df.columns
    # none of the per-key attribute columns are produced
    assert "gene_name" not in df.columns
    assert "transcript_id" not in df.columns
    # the helper 'attribute_split' column from parse_gtf is also suppressed
    assert "attribute_split" not in df.columns
    # core GTF columns are present and populated
    for col in ("seqname", "source", "feature", "start", "end", "strand"):
        assert col in df.columns
    assert len(df) > 0
    # spot-check that the raw attribute string carries the original key/value form
    assert any("gene_id" in val for val in df["attribute"].astype(str))


def test_expand_attribute_column_false_returns_polars():
    df = read_gtf(GTF_PATH, expand_attribute_column=False, result_type="polars")
    # polars dataframe — has columns attribute but no per-key columns
    assert "attribute" in df.columns
    assert "gene_name" not in df.columns
    assert "attribute_split" not in df.columns


def test_expand_attribute_column_false_returns_dict():
    result = read_gtf(GTF_PATH, expand_attribute_column=False, result_type="dict")
    assert isinstance(result, dict)
    assert "attribute" in result
    assert "gene_name" not in result


def test_expand_attribute_column_false_with_features_filter():
    """The features filter must still apply when not expanding."""
    df = read_gtf(
        GTF_PATH,
        expand_attribute_column=False,
        features={"gene"},
        result_type="pandas",
    )
    assert set(df["feature"]) == {"gene"}


def test_expand_attribute_column_false_skips_alias_and_version_logic():
    """When attribute columns aren't expanded, attribute_aliases has
    nothing to rename and cast_version_columns has nothing to cast.
    Neither should raise — both must be graceful no-ops on the raw
    'attribute'-column-only frame."""
    df = read_gtf(
        GTF_PATH,
        expand_attribute_column=False,
        attribute_aliases={"gene_type": "gene_biotype"},
        cast_version_columns=True,
        result_type="pandas",
    )
    # alias source wasn't in columns → no rename happened → no canonical added
    assert "gene_biotype" not in df.columns
    # version columns weren't present → no cast → still nothing
    assert "gene_version" not in df.columns
