import gzip
import subprocess
import sys
from io import BytesIO, StringIO

import pandas as pd
import pytest

from gtfparse import ParsingError, parse_gtf, read_gtf, write_gtf
from gtfparse._read_csv import _GtfReader, _PrefixedReader

GTF = 'NA\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "G1"; note "café #1; ok";\n'


@pytest.mark.parametrize("size", [1, 7, 128, -1])
def test_arrow_input_adapter_respects_byte_read_sizes(size):
    payload = ("\ufeff#header\n" + GTF + "#comment\n" + GTF.rstrip("\n")).encode()
    with _GtfReader(BytesIO(payload)) as reader:
        assert reader.readable()
        assert reader.read(0) == b""
        chunks = []
        while chunk := reader.read(size):
            if size >= 0:
                assert len(chunk) <= size
            chunks.append(chunk)
        assert b"".join(chunks) == (GTF * 2).encode()


def test_gzip_prefix_adapter_handles_partial_and_complete_reads():
    reader = _PrefixedReader(BytesIO(b"rest"), b"prefix")
    assert reader.read(2) == b"pr"
    assert reader.read(5) == b"efixr"
    assert reader.read() == b"est"


def test_arrow_input_adapter_rejects_malformed_lines():
    with _GtfReader(StringIO(GTF + "bad\n")) as reader, pytest.raises(ParsingError):
        reader.read()


@pytest.mark.parametrize("binary", [False, True])
def test_comments_bom_crlf_and_literal_hash(binary):
    text = "\ufeff##header\n" + GTF + "#\tcomment\twith\teight\ttabs\t.\t.\t.\t.\n" + GTF
    text = text.replace("\n", "\r\n")
    source = BytesIO(text.encode()) if binary else StringIO(text)
    frame = read_gtf(source)
    assert isinstance(frame, pd.DataFrame)
    assert frame["seqname"].tolist() == ["NA", "NA"]
    assert frame["note"].tolist() == ["café #1; ok", "café #1; ok"]
    assert not source.closed


@pytest.mark.parametrize("compressed", [False, True])
def test_nonseekable_short_reads_across_utf8_and_line_boundaries(compressed):
    payload = ("\ufeff#header\n" + GTF.rstrip("\n")).encode()
    if compressed:
        payload = gzip.compress(payload)

    class Stream:
        def __init__(self):
            self.buffer = BytesIO(payload)

        def read(self, size=-1):
            return self.buffer.read(min(size, 1) if size >= 0 else size)

    frame = read_gtf(Stream())
    assert frame["note"].tolist() == ["café #1; ok"]


def test_gzip_detected_from_contents_and_caller_stream_stays_open(tmp_path):
    payload = gzip.compress(GTF.encode())
    path = tmp_path / "annotation.GZ"
    path.write_bytes(payload)
    stream = BytesIO(payload)
    pd.testing.assert_frame_equal(read_gtf(path), read_gtf(stream))
    assert not stream.closed


@pytest.mark.parametrize("text", ["", "#header\n", " \n\n"])
def test_empty_input_is_package_error(text):
    events = []
    with pytest.raises(ParsingError):
        read_gtf(StringIO(text), progress_callback=lambda *event: events.append(event))
    assert events == [("read", 0, None)]


@pytest.mark.parametrize("change", ["short", "extra", "bad_start", "negative_frame", "overflow"])
def test_invalid_rows_raise_package_error_without_completion(change):
    fields = GTF.rstrip("\n").split("\t")
    if change == "short":
        fields.pop()
    elif change == "extra":
        fields.append("extra")
    elif change == "bad_start":
        fields[3] = "abc"
    elif change == "negative_frame":
        fields[7] = "-1"
    else:
        fields[3] = "9223372036854775808"
    events = []
    with pytest.raises(ParsingError):
        read_gtf(StringIO("\t".join(fields)), progress_callback=lambda *e: events.append(e))
    assert events == [("read", 0, None)]


def test_fixed_column_types_missing_values_and_large_integers():
    text = (
        "1\ttest\tgene\t9007199254740993\t9007199254740995\t0.25\t+\t2\t.\n"
        '1\ttest\tgene\t1\t10\t.\t.\t.\tgene_id "G1";\n'
    )
    frame = read_gtf(StringIO(text))
    assert frame["start"].tolist() == [9007199254740993, 1]
    assert frame["start"].dtype == "int64"
    assert frame["score"].dtype == "float32"
    assert pd.isna(frame["score"][1])
    assert frame["frame"].dtype == "uint32"
    assert frame["frame"].tolist() == [2, 0]
    assert frame["strand"].dtype == "category"
    assert pd.isna(frame["strand"][1])
    assert frame["gene_id"].tolist() == ["", "G1"]


def test_split_attributes_keep_quoted_semicolons():
    frame = parse_gtf(StringIO(GTF))
    assert frame["attribute_split"][0] == ['gene_id "G1"', 'note "café #1; ok"']


def test_filtered_out_null_coordinates_do_not_change_integer_types():
    text = '1\ttest\texon\t.\t20\t.\t+\t.\tgene_id "G2";\n' + GTF
    frame = read_gtf(StringIO(text), features={"gene"})
    assert frame["start"].dtype == "int64"
    assert frame["start"].tolist() == [1]
    assert frame["end"].dtype == "int64"


def test_null_seqnames_are_dropped_before_attribute_expansion():
    frame = read_gtf(StringIO(GTF.replace("NA\t", ".\t") + GTF))
    assert frame["gene_id"].tolist() == ["G1"]
    assert frame.index.tolist() == [0]


def test_opt_in_quote_cleanup_precedes_feature_filter():
    frame = parse_gtf(
        StringIO(GTF.replace("\tgene\t", "\tgene;-x\t")),
        features={"gene-x"},
        fix_quotes_columns=["feature"],
    )
    assert frame["feature"].tolist() == ["gene-x"]


def test_invalid_writer_input_does_not_replace_destination(tmp_path):
    path = tmp_path / "existing.gtf"
    path.write_text("keep")
    with pytest.raises(TypeError, match="DataFrame"):
        write_gtf({}, path)
    assert path.read_text() == "keep"


def test_base_api_does_not_import_polars(tmp_path):
    code = """
import importlib.abc
import sys
from io import StringIO

class NoPolars(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "polars" or fullname.startswith("polars."):
            raise AssertionError("Base API attempted to import Polars")

sys.meta_path.insert(0, NoPolars())
from gtfparse import read_gtf, write_gtf, parse_gtf, parse_gtf_and_expand_attributes
text = '1\\ttest\\tgene\\t1\\t100\\t.\\t+\\t.\\tgene_id "G1";\\n'
frame = read_gtf(StringIO(text))
assert "polars" not in sys.modules
assert frame["gene_id"].tolist() == ["G1"]
parse_gtf(StringIO(text))
parse_gtf_and_expand_attributes(StringIO(text))
write_gtf(frame, sys.argv[1])
assert read_gtf(sys.argv[1])["gene_id"].tolist() == ["G1"]
assert "polars" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code, str(tmp_path / "out.gtf")], check=True)


def test_missing_polars_fails_before_reading(monkeypatch):
    monkeypatch.setitem(sys.modules, "polars", None)
    source = StringIO(GTF)
    events = []
    with pytest.raises(ImportError, match=r"gtfparse\[polars\]"):
        read_gtf(source, result_type="polars", progress_callback=lambda *e: events.append(e))
    assert source.tell() == 0
    assert events == []


def test_writer_handles_multiple_batches_and_nullable_attributes(tmp_path):
    frame = read_gtf(StringIO(GTF))
    frame = pd.concat([frame] * 10_003, ignore_index=True)
    frame["gene_version"] = pd.Series([1] + [pd.NA] * (len(frame) - 1), dtype="Int64")
    out = tmp_path / "out.gtf.gz"
    write_gtf(frame, out)
    pd.testing.assert_frame_equal(frame, read_gtf(out))
