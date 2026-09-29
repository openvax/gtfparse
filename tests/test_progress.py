import gzip
from io import BytesIO, StringIO

import pandas as pd
import polars
import pytest
from polars.testing import assert_frame_equal

from gtfparse import (
    expand_attribute_strings,
    parse_gtf,
    parse_gtf_and_expand_attributes,
    parse_gtf_pandas,
    read_gtf,
)

GTF = (
    '1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "G1"; gene_version "2";\n'
    '1\ttest\texon\t1\t50\t.\t+\t0\tgene_id "G1"; tag "basic";\n'
)


@pytest.mark.parametrize("result_type", ["polars", "pandas", "dict"])
@pytest.mark.parametrize("expand", [True, False])
@pytest.mark.parametrize("features,rows", [(None, 2), ({"gene"}, 1), ({"CDS"}, 0)])
def test_read_progress_preserves_results(result_type, expand, features, rows, capsys):
    events = []
    kwargs = {
        "result_type": result_type,
        "expand_attribute_column": expand,
        "features": features,
        "attribute_aliases": {"gene_id": "gene"},
        "usecols": ["feature", "start", "gene", "gene_version", "attribute"],
    }
    expected = read_gtf(StringIO(GTF), **kwargs)
    actual = read_gtf(
        StringIO(GTF), progress_callback=lambda *event: events.append(event), **kwargs
    )
    if result_type == "polars":
        assert_frame_equal(actual, expected)
    elif result_type == "pandas":
        pd.testing.assert_frame_equal(actual, expected)
    else:
        assert actual == expected
    expected_events = [("read", 0, None), ("read", rows, rows)]
    if expand:
        expected_events += [("attributes", 0, rows)]
        if rows:
            expected_events.append(("attributes", rows, rows))
    expected_events += [("convert", 0, None), ("convert", rows, rows)]
    assert events == expected_events
    assert capsys.readouterr() == ("", "")


@pytest.mark.parametrize("kind", ["path", "gzip", "bytes", "text"])
def test_progress_with_supported_inputs(kind, tmp_path):
    if kind == "path":
        source = tmp_path / "test.gtf"
        source.write_text(GTF)
    elif kind == "gzip":
        source = tmp_path / "test.gtf.gz"
        source.write_bytes(gzip.compress(GTF.encode()))
    elif kind == "bytes":
        source = BytesIO(GTF.encode())
    else:
        source = StringIO(GTF)
    events = []
    result = read_gtf(source, progress_callback=lambda *event: events.append(event))
    assert result.height == 2
    assert events[-1] == ("convert", 2, 2)


@pytest.mark.parametrize(
    "parser,stages",
    [
        (parse_gtf, ["read"]),
        (parse_gtf_and_expand_attributes, ["read", "attributes"]),
        (parse_gtf_pandas, ["read", "convert"]),
    ],
)
def test_lower_level_progress(parser, stages):
    events = []
    parser(StringIO(GTF), progress_callback=lambda *event: events.append(event))
    assert [event[0] for event in events[::2]] == stages
    assert all(event[1:] == (2, 2) for event in events[1::2])


@pytest.mark.parametrize("rows", [0, 1, 10_000, 20_003])
def test_attribute_progress_counts_rows_and_is_bounded(rows):
    events = []
    attributes = ['gene_id "G1";'] * rows
    result = expand_attribute_strings(
        attributes, progress_callback=lambda *event: events.append(event)
    )
    assert events[0] == ("attributes", 0, rows)
    assert events[-1] == ("attributes", rows, rows)
    counts = [event[1] for event in events]
    assert counts == sorted(set(counts))
    assert len(events) <= rows // 10_000 + 2
    assert result.get("gene_id", []) == ["G1"] * rows


@pytest.mark.parametrize("stage", ["read", "attributes", "convert"])
def test_callback_errors_propagate_without_completion(stage):
    events = []
    error = RuntimeError("cancel parsing")

    def callback(*event):
        events.append(event)
        if event[0] == stage:
            raise error

    source = StringIO(GTF)
    with pytest.raises(RuntimeError) as caught:
        read_gtf(source, progress_callback=callback)
    assert caught.value is error
    assert events[-1][0:2] == (stage, 0)
    assert not any(event[0] == stage and event[1] == 2 for event in events)
    if stage == "read":
        assert source.tell() == 0


def test_conversion_errors_do_not_report_completion():
    events = []
    with pytest.raises(ValueError):
        read_gtf(
            StringIO(GTF),
            column_converters={"gene_id": int},
            progress_callback=lambda *event: events.append(event),
        )
    assert events[-1] == ("convert", 0, None)


def test_callback_can_cancel_during_attribute_expansion():
    events = []

    def callback(*event):
        events.append(event)
        if event[:2] == ("attributes", 10_000):
            raise RuntimeError("cancel during expansion")

    with pytest.raises(RuntimeError, match="cancel during expansion"):
        read_gtf(StringIO(GTF * 10_000), progress_callback=callback)
    assert events[-1] == ("attributes", 10_000, 20_000)
    assert not any(event[0] == "convert" for event in events)


def test_empty_file_keeps_existing_error_without_completion():
    events = []
    with pytest.raises(polars.exceptions.NoDataError):
        read_gtf(StringIO(""), progress_callback=lambda *event: events.append(event))
    assert events == [("read", 0, None)]


@pytest.mark.parametrize("invalid", ["numpy", "Pandas", "", None, []])
def test_invalid_result_type_rejected_before_reading(invalid):
    source = StringIO(GTF)
    events = []
    with pytest.raises(ValueError, match="result_type"):
        read_gtf(source, result_type=invalid, progress_callback=lambda *event: events.append(event))
    assert source.tell() == 0
    assert events == []
