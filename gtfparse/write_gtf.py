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
import sys
from collections.abc import Iterable
from itertools import repeat
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union

import pandas as pd

if TYPE_CHECKING:
    import polars

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


def _lines(df):
    """Format bounded batches, preserving each column's string representation."""
    attribute_columns = [name for name in df.columns if name not in GTF_FIXED_COLUMNS]
    for start in range(0, len(df), 10_000):
        batch = df.iloc[start : start + 10_000]
        fixed = batch[GTF_FIXED_COLUMNS].astype("string").fillna(MISSING_VALUE)
        if RAW_ATTRIBUTE_COLUMN in df.columns:
            attributes = batch[RAW_ATTRIBUTE_COLUMN].astype("string").fillna("")
        elif not attribute_columns:
            attributes = repeat("", len(batch))
        else:
            expanded = batch[attribute_columns].astype("string").fillna("")
            attributes = (
                " ".join(
                    '%s "%s";' % (name, value)
                    for name, value in zip(attribute_columns, row)
                    if value != ""
                )
                for row in expanded.itertuples(index=False, name=None)
            )
        for row, attribute in zip(fixed.itertuples(index=False, name=None), attributes):
            yield "\t".join([*row, attribute]) + "\n"


def write_gtf(
    df: Union[pd.DataFrame, "polars.DataFrame"],
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
    df : pandas.DataFrame or polars.DataFrame
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
    Attribute values are enclosed in double quotes, so spaces, apostrophes,
    and semicolons round-trip. This writer does not escape embedded double
    quotes, tabs, or newlines; values containing those characters are outside
    its round-trip guarantees. Empty and missing attributes are both omitted.
    """
    if not isinstance(df, pd.DataFrame):
        polars = sys.modules.get("polars")
        if polars is None or not isinstance(df, polars.DataFrame):
            raise TypeError("df must be a pandas or Polars DataFrame")
        df = df.to_pandas()
    missing = [name for name in GTF_FIXED_COLUMNS if name not in df.columns]
    if missing:
        raise ValueError("DataFrame is missing required GTF column(s): %s" % ", ".join(missing))

    open_file = gzip.open if str(path).lower().endswith(".gz") else open
    with open_file(path, "wt", encoding="utf-8", newline="\n") as output_file:
        if header_lines is not None:
            for header_line in header_lines:
                output_file.write("%s\n" % header_line)
        output_file.writelines(_lines(df))

    logger.info("Wrote %d GTF rows to %s", len(df), path)
