# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging
from collections.abc import Iterable
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union

import polars

if TYPE_CHECKING:
    import pandas

logger = logging.getLogger(__name__)

# The eight tab-separated columns that precede the attribute field of a GTF
# line, in the order they must be written. Any other column in a DataFrame is
# treated as an expanded attribute (see read_gtf's expand_attribute_column).
GTF_FIXED_COLUMNS = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
]

# GTF uses a single dot to denote a missing value in the fixed columns.
MISSING_VALUE = "."


def _format_fixed_value(value) -> str:
    """Render a fixed-column value, using '.' for missing values."""
    if value is None:
        return MISSING_VALUE
    return str(value)


def _format_attributes(row: dict, attribute_columns: list[str], raw_column: Optional[str]) -> str:
    """
    Build the GTF attribute field for a single row.

    If ``raw_column`` is given (an unexpanded 'attribute' column) its value is
    emitted verbatim. Otherwise each expanded attribute column is serialized as
    ``key "value";`` and the pairs are joined with a single space, matching the
    format produced by common GTF writers (and parsed by read_gtf).

    Attributes whose value is None are omitted, since that is how read_gtf
    represents a key that was absent on a given row. Note that a value must be
    *None* to be skipped -- falsy-but-present values such as 0 or the empty
    string are written out, so they survive a read/write round trip.
    """
    if raw_column is not None:
        raw = row[raw_column]
        return "" if raw is None else str(raw)

    parts = []
    for column_name in attribute_columns:
        value = row[column_name]
        if value is None:
            continue
        parts.append('%s "%s";' % (column_name, value))
    return " ".join(parts)


def write_gtf(
    df: Union[polars.DataFrame, "pandas.DataFrame"],
    path: Union[str, Path],
    header_lines: Optional[Iterable[str]] = None,
) -> None:
    """
    Write a DataFrame of genomic features back out to a GTF file.

    This is the inverse of :func:`read_gtf`. A DataFrame produced by
    ``read_gtf`` (whether the attribute column was expanded or not) can be
    written back out and re-read to recover an equivalent DataFrame.

    Parameters
    ----------
    df : polars.DataFrame or pandas.DataFrame
        Feature rows to write. Must contain the fixed GTF columns
        (seqname, source, feature, start, end, score, strand, frame).
        Any additional columns are written as attributes, except a column
        literally named 'attribute', which is treated as a pre-formatted
        attribute string and emitted verbatim.

    path : str or pathlib.Path
        Destination file path. Any existing file is overwritten.

    header_lines : iterable of str, optional
        Lines to write at the top of the file before any feature rows, e.g.
        ``["##description: example", "##provider: GENCODE"]``. Each is written
        verbatim on its own line, so include a leading '#' if you want it to be
        parsed back as a comment.
    """
    # Accept a pandas DataFrame too, since read_gtf(result_type="pandas")
    # returns one; convert to polars so the row iteration below is uniform.
    if not isinstance(df, polars.DataFrame):
        df = polars.from_pandas(df)

    columns = df.columns
    fixed_columns = [name for name in GTF_FIXED_COLUMNS if name in columns]
    missing_fixed = [name for name in GTF_FIXED_COLUMNS if name not in columns]
    if missing_fixed:
        raise ValueError(
            "DataFrame is missing required GTF column(s): %s" % ", ".join(missing_fixed)
        )

    # A column named 'attribute' is the raw, unexpanded attribute string;
    # everything else that isn't a fixed column is an expanded attribute.
    raw_column = "attribute" if "attribute" in columns else None
    attribute_columns = [
        name for name in columns if name not in GTF_FIXED_COLUMNS and name != "attribute"
    ]

    n_rows = 0
    with open(path, "w") as output_file:
        if header_lines is not None:
            for line in header_lines:
                output_file.write("%s\n" % line)
        for row in df.iter_rows(named=True):
            fixed_fields = [_format_fixed_value(row[name]) for name in fixed_columns]
            attribute_field = _format_attributes(row, attribute_columns, raw_column)
            output_file.write("%s\t%s\n" % ("\t".join(fixed_fields), attribute_field))
            n_rows += 1

    logger.info("Wrote %d GTF rows to %s", n_rows, path)
