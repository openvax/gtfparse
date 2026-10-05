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
import re
from collections import OrderedDict
from sys import intern

logger = logging.getLogger(__name__)

# Quoted values may contain semicolons. Keep this compatible with both Python's
# re module so the low-level split column agrees with
# attribute expansion. Unquoted values and single quotes are accepted as well.
ATTRIBUTE_PATTERN = r"""[^\s;]+[ \t]+(?:"[^"]*"|'[^']*'|[^;]+)"""
_ATTRIBUTE_PAIRS = re.compile(ATTRIBUTE_PATTERN)


def expand_attribute_strings(
    attribute_strings, quote_char="'", missing_value="", usecols=None, *, progress_callback=None
):
    """
    The last column of GTF has a variable number of key value pairs
    of the format: "key1 value1; key2 value2;"
    Parse these into a dictionary mapping each key onto a list of values,
    using missing_value for any row where the key was missing.

    Parameters
    ----------
    attribute_strings : sequence of str or sequences of str
        Raw attribute fields or already separated key/value pairs per row.
        None represents a row with no attributes.

    quote_char : str
        Additional surrounding quote character to remove from values.
        Standard double quotes are always recognized.

    missing_value : any
        If an attribute is missing from a row, give it this value.

    usecols : list of str or None
        If not None, then only expand columns included in this set,
        otherwise use all columns.

    progress_callback : callable, optional
        Called as ``callback("attributes", completed_rows, total_rows)`` at
        the start, every 10,000 rows, and at completion. An empty input emits
        one ``("attributes", 0, 0)`` event. Callback exceptions propagate.

    Returns OrderedDict of column->value list mappings, in the order they
    appeared in the attribute strings.
    """
    n = len(attribute_strings)

    extra_columns = {}
    column_order = []

    # While parsing millions of repeated strings (e.g. "gene_id" and "TP53"),
    # we can save a lot of memory by making sure there's only one string
    # object per unique column name. Cache interned names locally as well.
    column_interned_strings = {}

    if progress_callback is not None:
        progress_callback("attributes", 0, n)

    for i, kv_strings in enumerate(attribute_strings):
        if isinstance(kv_strings, str):
            kv_strings = _ATTRIBUTE_PAIRS.findall(kv_strings)
        elif kv_strings is None:
            kv_strings = ()
        for kv in kv_strings:
            # Split once: values such as transcript_support_level may include
            # spaces, and whitespace between a key and its value may vary.
            parts = kv.strip().split(None, 1)

            if len(parts) != 2:
                continue

            column_name, value = parts

            try:
                column_name = column_interned_strings[column_name]
            except KeyError:
                column_name = intern(column_name)
                column_interned_strings[column_name] = column_name

            if usecols is not None and column_name not in usecols:
                continue

            if len(value) >= 2 and value[0] in ('"', quote_char) and value[-1] == value[0]:
                value = value[1:-1]

            try:
                column = extra_columns[column_name]
                # if an attribute is used repeatedly then
                # keep track of all its values in a list
                old_value = column[i]
                if old_value is missing_value:
                    column[i] = value
                else:
                    column[i] = "%s,%s" % (old_value, value)
            except KeyError:
                column = [missing_value] * n
                column[i] = value
                extra_columns[column_name] = column
                column_order.append(column_name)

        completed = i + 1
        if progress_callback is not None and (completed % 10_000 == 0 or completed == n):
            progress_callback("attributes", completed, n)

    logger.info("Extracted GTF attributes: %s", column_order)
    return OrderedDict((column_name, extra_columns[column_name]) for column_name in column_order)
