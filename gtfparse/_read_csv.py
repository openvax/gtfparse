"""Read GTF fixed columns without importing an optional dataframe library."""

import gzip
import io
from contextlib import contextmanager
from os import PathLike

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pacsv

from .parsing_error import ParsingError

REQUIRED_COLUMNS = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute",
]
CATEGORICAL_COLUMNS = ("seqname", "source", "feature", "strand")


class _PrefixedReader:
    def __init__(self, stream, prefix):
        self.stream = stream
        self.prefix = prefix

    def read(self, size=-1):
        if size < 0:
            prefix, self.prefix = self.prefix, self.prefix[:0]
            return prefix + self.stream.read()
        prefix, self.prefix = self.prefix[:size], self.prefix[size:]
        return prefix + self.stream.read(size - len(prefix))


@contextmanager
def _uncompressed(stream):
    prefix = stream.read(2)
    while prefix and len(prefix) < 2:
        extra = stream.read(2 - len(prefix))
        if not extra:
            break
        prefix += extra
    prefixed = _PrefixedReader(stream, prefix)
    if prefix == b"\x1f\x8b":
        with gzip.GzipFile(fileobj=prefixed) as uncompressed:
            yield uncompressed
    else:
        yield prefixed


@contextmanager
def _open_input(source):
    if isinstance(source, (str, PathLike)):
        with open(source, "rb") as stream, _uncompressed(stream) as uncompressed:
            yield uncompressed
    else:
        with _uncompressed(source) as uncompressed:
            yield uncompressed


class _GtfReader(io.RawIOBase):
    """Filter whole comments and validate fields in bounded input chunks.

    Filtering whole lines preserves literal # inside attributes and skips
    comments even when they contain tabs. Validate fields before CSV parsing.
    """

    def __init__(self, stream):
        super().__init__()
        self.stream = stream
        self.pending = b""
        self.buffer = b""
        self.first = True
        self.finished = False

    def readable(self):
        return True

    def read(self, size=-1):
        while not self.finished and (size < 0 or len(self.buffer) < size):
            chunk = self.stream.read(max(size, 64 * 1024))
            self.finished = not chunk
            if isinstance(chunk, str):
                chunk = chunk.encode("utf-8")
            lines = (self.pending + chunk).split(b"\n")
            self.pending = b"" if self.finished else lines.pop()
            if self.first and lines:
                lines[0] = lines[0].removeprefix(b"\xef\xbb\xbf")
                self.first = False
            valid = []
            for line in lines:
                if not line.strip() or line.startswith(b"#"):
                    continue
                if line.count(b"\t") != 8:
                    raise ParsingError("Wrong number of columns: expected nine GTF fields")
                valid.append(line + b"\n")
            self.buffer += b"".join(valid)
        if size < 0:
            result, self.buffer = self.buffer, b""
        else:
            result, self.buffer = self.buffer[:size], self.buffer[size:]
        return result


def read_fixed_columns(source, features=None):
    types = dict.fromkeys(CATEGORICAL_COLUMNS, pa.dictionary(pa.int32(), pa.string()))
    types.update(
        start=pa.int64(),
        end=pa.int64(),
        score=pa.float32(),
        frame=pa.uint32(),
        attribute=pa.string(),
    )
    try:
        with _open_input(source) as stream, _GtfReader(stream) as cleaned:
            table = pacsv.read_csv(
                cleaned,
                read_options=pacsv.ReadOptions(column_names=REQUIRED_COLUMNS),
                parse_options=pacsv.ParseOptions(delimiter="\t", quote_char=False),
                convert_options=pacsv.ConvertOptions(
                    column_types=types, null_values=[".", ""], strings_can_be_null=True
                ),
            )
    except pa.ArrowInvalid as error:
        raise ParsingError("Invalid or empty GTF input: %s" % error) from error
    if table["seqname"].null_count:
        table = table.filter(pc.is_valid(table["seqname"]))
    if features is not None:
        values = pa.array(sorted(set(features)), type=pa.string())
        table = table.filter(pc.is_in(table["feature"], value_set=values))
    table = table.set_column(7, "frame", pc.fill_null(table["frame"], 0))
    return table.to_pandas()
