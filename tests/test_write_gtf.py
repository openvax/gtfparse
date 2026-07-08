import gzip

import polars
import pytest
from polars.testing import assert_frame_equal

from gtfparse import read_gtf, write_gtf

from .data import data_path

# A spread of real GTF flavors: RefSeq (minimal attrs), Ensembl (many attrs +
# null scores + *_version columns), GENCODE (mixed feature types so different
# rows carry different attribute sets), and StringTie.
FIXTURES = [
    "refseq.ucsc.small.gtf",
    "ensembl_grch37.head.gtf",
    "gencode.head.gtf",
    "gencode.real.head.gtf",
    "B16.stringtie.head.gtf",
]


def _minimal_df(**extra_columns):
    """Build a one-row DataFrame with the fixed GTF columns plus extras.

    Each extra column is passed as a single-element list, e.g.
    ``_minimal_df(gene_id=["G1"])``.
    """
    row = {
        "seqname": ["chr1"],
        "source": ["test"],
        "feature": ["gene"],
        "start": [1],
        "end": [100],
        "score": [None],
        "strand": ["+"],
        "frame": [None],
    }
    row.update(extra_columns)
    return polars.DataFrame(row)


# ---------------------------------------------------------------------------
# write(read(...)): start from a GTF file, read it, write it back.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("fixture", FIXTURES)
@pytest.mark.parametrize("expand", [True, False])
def test_write_read_recovers_parsed_frame(fixture, expand, tmp_path):
    """read -> write -> read reproduces the parsed DataFrame exactly (column
    order included), for both expanded and unexpanded attribute modes."""
    original = read_gtf(data_path(fixture), expand_attribute_column=expand)
    out_path = tmp_path / "out.gtf"
    write_gtf(original, out_path)
    recovered = read_gtf(str(out_path), expand_attribute_column=expand)
    assert_frame_equal(original, recovered, categorical_as_str=True)


@pytest.mark.parametrize("fixture", FIXTURES)
def test_written_output_is_idempotent(fixture, tmp_path):
    """write(read(...)) is a fixed point: writing what we read, reading it, and
    writing it again yields byte-identical output."""
    df = read_gtf(data_path(fixture))
    first = tmp_path / "first.gtf"
    second = tmp_path / "second.gtf"
    write_gtf(df, first)
    write_gtf(read_gtf(str(first)), second)
    assert first.read_bytes() == second.read_bytes()


# ---------------------------------------------------------------------------
# read(write(...)): start from an in-memory DataFrame, write it, read it back.
# ---------------------------------------------------------------------------


def test_read_write_recovers_dataframe(tmp_path):
    """A DataFrame written out and read back preserves its attribute values."""
    df = _minimal_df(gene_id=["ENSG1"], gene_name=["DDX11L1"], gene_version=["5"])
    out_path = tmp_path / "df.gtf"
    write_gtf(df, out_path)
    recovered = read_gtf(str(out_path))
    assert recovered["gene_id"].to_list() == ["ENSG1"]
    assert recovered["gene_name"].to_list() == ["DDX11L1"]
    # writing the recovered frame is itself a fixed point
    again = tmp_path / "df2.gtf"
    write_gtf(recovered, again)
    assert_frame_equal(recovered, read_gtf(str(again)), categorical_as_str=True)


def test_round_trip_from_pandas(tmp_path):
    """write_gtf accepts a pandas DataFrame (read_gtf(result_type='pandas'))."""
    polars_df = read_gtf(data_path("ensembl_grch37.head.gtf"))
    pandas_df = read_gtf(data_path("ensembl_grch37.head.gtf"), result_type="pandas")
    out_path = tmp_path / "from_pandas.gtf"
    write_gtf(pandas_df, out_path)
    assert_frame_equal(polars_df, read_gtf(str(out_path)), categorical_as_str=True)


def test_gzip_output_round_trips(tmp_path):
    """A '.gz' path is gzip-compressed on disk and read back transparently."""
    df = read_gtf(data_path("ensembl_grch37.head.gtf"))
    out_path = tmp_path / "out.gtf.gz"
    write_gtf(df, out_path)
    # actually gzip-compressed on disk
    assert out_path.read_bytes()[:2] == b"\x1f\x8b"
    with gzip.open(out_path, "rt") as handle:
        assert "\t" in handle.readline()
    assert_frame_equal(df, read_gtf(str(out_path)), categorical_as_str=True)


def test_gzip_detection_is_case_insensitive(tmp_path):
    """An uppercase '.GZ' suffix is still gzip-compressed."""
    df = read_gtf(data_path("refseq.ucsc.small.gtf"))
    out_path = tmp_path / "out.GZ"
    write_gtf(df, out_path)
    assert out_path.read_bytes()[:2] == b"\x1f\x8b"


def test_empty_dataframe_writes_no_rows(tmp_path):
    """A zero-row DataFrame produces a file with only its header lines."""
    empty = read_gtf(data_path("refseq.ucsc.small.gtf")).clear()
    out_path = tmp_path / "empty.gtf"
    write_gtf(empty, out_path, header_lines=["##empty"])
    assert out_path.read_text() == "##empty\n"


def test_fixed_columns_only(tmp_path):
    """A DataFrame with only the fixed columns writes a valid 9-field line
    (empty attribute field) and reads back."""
    fixed_only = polars.DataFrame(
        {
            "seqname": ["chr1"],
            "source": ["test"],
            "feature": ["gene"],
            "start": [1],
            "end": [100],
            "score": [None],
            "strand": ["+"],
            "frame": [None],
        }
    )
    out_path = tmp_path / "fixed.gtf"
    write_gtf(fixed_only, out_path)
    # nine tab-separated fields, the last (attribute) empty
    assert out_path.read_text().strip("\n").split("\t") == [
        "chr1",
        "test",
        "gene",
        "1",
        "100",
        ".",
        "+",
        ".",
        "",
    ]
    # read_gtf accepts it (no attribute columns to expand)
    recovered = read_gtf(str(out_path), expand_attribute_column=False)
    assert recovered["seqname"].to_list() == ["chr1"]


# ---------------------------------------------------------------------------
# Attribute / missing-value semantics.
# ---------------------------------------------------------------------------


def test_nonempty_value_written_empty_and_null_omitted(tmp_path):
    """A non-empty value such as '0' is written; empty string and null are
    treated as absent and omitted (matching read_gtf's missing-value model)."""
    df = _minimal_df(
        gene_id=["G1"],
        zero_attr=["0"],
        empty_attr=[""],
        missing_attr=[None],
    )
    out_path = tmp_path / "attrs.gtf"
    write_gtf(df, out_path)
    line = out_path.read_text().strip()
    assert 'gene_id "G1";' in line
    assert 'zero_attr "0";' in line  # non-empty falsy value survives
    assert "empty_attr" not in line  # empty string omitted as absent
    assert "missing_attr" not in line  # null omitted as absent


def test_missing_value_fixed_columns_use_dot(tmp_path):
    """None in the fixed columns is serialized as '.'."""
    df = _minimal_df(gene_id=["G1"])  # score and frame are None
    out_path = tmp_path / "dots.gtf"
    write_gtf(df, out_path)
    fields = out_path.read_text().strip().split("\t")
    assert fields[5] == "."  # score
    assert fields[7] == "."  # frame


def test_structural_characters_are_not_round_trippable(tmp_path):
    """GTF has no escaping for '"' or ';'; read_gtf strips quotes and splits on
    ';'. Pin that such values cannot round-trip so a future change is noticed."""
    df = _minimal_df(gene_id=["A;B"], note=['say "hi"'])
    out_path = tmp_path / "special.gtf"
    write_gtf(df, out_path)
    recovered = read_gtf(str(out_path))
    # the semicolon split the value apart
    assert recovered["gene_id"].to_list() != ["A;B"]
    # the double quotes were stripped from the value
    assert '"' not in recovered["note"][0]


def test_header_lines_are_written(tmp_path):
    df = read_gtf(data_path("refseq.ucsc.small.gtf"))
    out_path = tmp_path / "with_header.gtf"
    write_gtf(df, out_path, header_lines=["##description: test", "##provider: gtfparse"])
    lines = out_path.read_text().splitlines()
    assert lines[0] == "##description: test"
    assert lines[1] == "##provider: gtfparse"
    # comment lines are ignored by read_gtf, so the data still round-trips
    assert_frame_equal(df, read_gtf(str(out_path)), categorical_as_str=True)


def test_missing_required_column_raises(tmp_path):
    df = polars.DataFrame({"seqname": ["chr1"], "gene_id": ["G1"]})
    with pytest.raises(ValueError, match="missing required GTF column"):
        write_gtf(df, tmp_path / "bad.gtf")
