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

import gzip
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

# Name of the raw, unexpanded attribute column produced by
# read_gtf(expand_attribute_column=False).
RAW_ATTRIBUTE_COLUMN = "attribute"


def _attribute_expr(attribute_columns: list[str]) -> polars.Expr:
    """
    Build a polars expression that renders the GTF attribute field for each row
    from a set of expanded attribute columns.

    Each column ``key`` becomes ``key "value";`` and the per-row pairs are
    joined with a single space, matching the format emitted by Ensembl/GENCODE
    and parsed by read_gtf.

    A pair is omitted for any row where the value is "absent". read_gtf marks an
    absent key with null (numeric columns, e.g. the ``*_version`` fields) or the
    empty string (string columns) and cannot distinguish absent from
    present-but-empty, so we treat both null and "" as absent. Values that are
    merely falsy but non-empty -- notably the string "0" -- are written and
    survive a round trip (this is the regression that the naive ``if value:``
    check in earlier drafts got wrong).
    """
    if not attribute_columns:
        return polars.lit("")
    pairs = []
    for name in attribute_columns:
        value = polars.col(name).cast(polars.String)
        pairs.append(
            polars.when(value.is_null() | (value == ""))
            .then(None)
            .otherwise(polars.format('{} "{}";', polars.lit(name), value))
        )
    # ignore_nulls drops absent keys; fill_null covers the all-absent row so the
    # surrounding line expression never collapses to null.
    return polars.concat_str(pairs, separator=" ", ignore_nulls=True).fill_null("")


def _line_series(df: polars.DataFrame) -> polars.Series:
    """
    Turn a DataFrame into a Series of fully-formatted GTF lines (without
    trailing newlines), building the whole thing with vectorized polars
    expressions rather than per-row Python.
    """
    columns = df.columns
    missing = [name for name in GTF_FIXED_COLUMNS if name not in columns]
    if missing:
        raise ValueError("DataFrame is missing required GTF column(s): %s" % ", ".join(missing))

    # Fixed columns are always written in canonical order; nulls become '.'.
    # Cast before fill_null so numeric columns accept the string sentinel.
    fixed = [
        polars.col(name).cast(polars.String).fill_null(MISSING_VALUE) for name in GTF_FIXED_COLUMNS
    ]

    if RAW_ATTRIBUTE_COLUMN in columns:
        # Unexpanded read: the attribute column is already a formatted string.
        attribute = polars.col(RAW_ATTRIBUTE_COLUMN).cast(polars.String).fill_null("")
    else:
        attribute_columns = [name for name in columns if name not in GTF_FIXED_COLUMNS]
        attribute = _attribute_expr(attribute_columns)

    # Always keep the 9th (attribute) field, even when empty, so every line has
    # the nine tab-separated columns read_gtf expects.
    line = polars.concat_str([*fixed, attribute], separator="\t")
    return df.select(line.alias("_gtf_line")).get_column("_gtf_line")


def write_gtf(
    df: Union[polars.DataFrame, "pandas.DataFrame"],
    path: Union[str, Path],
    header_lines: Optional[Iterable[str]] = None,
) -> None:
    """
    Write a DataFrame of genomic features back out to a GTF file.

    This is the inverse of :func:`read_gtf`. A DataFrame produced by
    ``read_gtf`` (whether the attribute column was expanded into one column per
    key or left as a raw ``attribute`` string) can be written back out and
    re-read to recover an equivalent DataFrame.

    Parameters
    ----------
    df : polars.DataFrame or pandas.DataFrame
        Feature rows to write. Must contain the fixed GTF columns
        (seqname, source, feature, start, end, score, strand, frame).
        Any additional column is written as an attribute, except a column
        literally named ``attribute``, which is treated as a pre-formatted
        attribute string and emitted verbatim.

    path : str or pathlib.Path
        Destination file path. Any existing file is overwritten. If the path
        ends in ``.gz`` the output is gzip-compressed (mirroring read_gtf,
        which transparently reads gzip-compressed GTFs).

    header_lines : iterable of str, optional
        Lines to write at the top of the file before any feature rows, e.g.
        ``["##description: example", "##provider: GENCODE"]``. Each is written
        verbatim on its own line, so include a leading ``#`` if you want it to
        be parsed back as a comment.

    Notes
    -----
    GTF has no escaping mechanism for the structural characters ``"`` and
    ``;`` inside attribute values (nor for the tab and newline that delimit
    columns and rows), and read_gtf strips quotes and splits on ``;`` when
    parsing. Values returned by read_gtf therefore never contain those
    characters, so any DataFrame obtained from read_gtf round-trips exactly. A
    DataFrame built by hand whose attribute values contain ``"``, ``;``, a tab,
    or a newline cannot be represented losslessly and will not round-trip.
    """
    # Accept a pandas DataFrame too, since read_gtf(result_type="pandas")
    # returns one; convert to polars so the formatting below is uniform.
    if not isinstance(df, polars.DataFrame):
        df = polars.from_pandas(df)

    lines = _line_series(df)

    open_file = gzip.open if str(path).lower().endswith(".gz") else open
    with open_file(path, "wt", encoding="utf-8", newline="\n") as output_file:
        if header_lines is not None:
            for header_line in header_lines:
                output_file.write("%s\n" % header_line)
        for line in lines:
            output_file.write("%s\n" % line)

    logger.info("Wrote %d GTF rows to %s", lines.len(), path)
