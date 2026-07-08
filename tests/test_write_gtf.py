import polars
import pytest
from polars.testing import assert_frame_equal

from gtfparse import read_gtf, write_gtf

from .data import data_path

REFSEQ_GTF_PATH = data_path("refseq.ucsc.small.gtf")
ENSEMBL_GTF_PATH = data_path("ensembl_grch37.head.gtf")


def _assert_round_trips(source_path, tmp_path, expand_attribute_column=True):
    """read -> write -> read should recover an equivalent DataFrame."""
    original = read_gtf(source_path, expand_attribute_column=expand_attribute_column)
    out_path = tmp_path / "round_trip.gtf"
    write_gtf(original, out_path)
    recovered = read_gtf(str(out_path), expand_attribute_column=expand_attribute_column)
    # categorical_as_str avoids spurious mismatches from category orderings
    assert_frame_equal(original, recovered, categorical_as_str=True)


def test_round_trip_expanded_refseq(tmp_path):
    _assert_round_trips(REFSEQ_GTF_PATH, tmp_path, expand_attribute_column=True)


def test_round_trip_expanded_ensembl(tmp_path):
    # Ensembl file exercises many attribute columns and null score values
    _assert_round_trips(ENSEMBL_GTF_PATH, tmp_path, expand_attribute_column=True)


def test_round_trip_unexpanded(tmp_path):
    # A raw, unexpanded 'attribute' column should be written verbatim
    _assert_round_trips(REFSEQ_GTF_PATH, tmp_path, expand_attribute_column=False)


def test_round_trip_from_pandas(tmp_path):
    # write_gtf should also accept a pandas DataFrame (result_type="pandas")
    original_polars = read_gtf(ENSEMBL_GTF_PATH)
    original_pandas = read_gtf(ENSEMBL_GTF_PATH, result_type="pandas")
    out_path = tmp_path / "from_pandas.gtf"
    write_gtf(original_pandas, out_path)
    recovered = read_gtf(str(out_path))
    assert_frame_equal(original_polars, recovered, categorical_as_str=True)


def test_falsy_attribute_values_are_preserved(tmp_path):
    """Attributes whose value is 0 or "" must be written; only None is omitted."""
    df = polars.DataFrame(
        {
            "seqname": ["chr1"],
            "source": ["test"],
            "feature": ["gene"],
            "start": [1],
            "end": [100],
            "score": [None],
            "strand": ["+"],
            "frame": [None],
            "gene_id": ["G1"],
            "zero_attr": ["0"],
            "empty_attr": [""],
            "missing_attr": [None],
        }
    )
    out_path = tmp_path / "falsy.gtf"
    write_gtf(df, out_path)
    line = out_path.read_text().strip()
    assert 'zero_attr "0"' in line
    assert 'empty_attr ""' in line
    # a None-valued attribute is omitted entirely
    assert "missing_attr" not in line


def test_missing_value_fixed_columns_use_dot(tmp_path):
    """None values in the fixed columns are serialized as '.'."""
    df = polars.DataFrame(
        {
            "seqname": ["chr1"],
            "source": ["test"],
            "feature": ["gene"],
            "start": [1],
            "end": [100],
            "score": [None],
            "strand": ["+"],
            "frame": [None],
            "gene_id": ["G1"],
        }
    )
    out_path = tmp_path / "dots.gtf"
    write_gtf(df, out_path)
    fields = out_path.read_text().strip().split("\t")
    assert fields[5] == "."  # score
    assert fields[7] == "."  # frame


def test_header_lines_are_written(tmp_path):
    df = read_gtf(REFSEQ_GTF_PATH)
    out_path = tmp_path / "with_header.gtf"
    header = ["##description: test", "##provider: gtfparse"]
    write_gtf(df, out_path, header_lines=header)
    lines = out_path.read_text().splitlines()
    assert lines[0] == "##description: test"
    assert lines[1] == "##provider: gtfparse"
    # header comment lines are ignored by read_gtf, so it still round-trips
    recovered = read_gtf(str(out_path))
    assert_frame_equal(df, recovered, categorical_as_str=True)


def test_missing_required_column_raises(tmp_path):
    df = polars.DataFrame({"seqname": ["chr1"], "gene_id": ["G1"]})
    with pytest.raises(ValueError, match="missing required GTF column"):
        write_gtf(df, tmp_path / "bad.gtf")
