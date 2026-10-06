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

logger = logging.getLogger(__name__)

# Quoted values may contain semicolons. Keep this compatible with both Python's
# re module so the low-level split column agrees with
# attribute expansion. Unquoted values and single quotes are accepted as well.
ATTRIBUTE_PATTERN = r"""[^\s;]+[ \t]+(?:"[^"]*"|'[^']*'|[^;]+)"""
# Capture double-quoted values without their delimiters. Other values retain
# their original quoting for custom quote_char and pre-split input support.
_ATTRIBUTE_PAIRS = re.compile(r"""([^\s;]+)[ \t]+(?:"([^"]*)"|('[^']*'|[^;]+))""")


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

    if progress_callback is not None:
        progress_callback("attributes", 0, n)

    for i, kv_strings in enumerate(attribute_strings):
        if isinstance(kv_strings, str):
            pairs = _ATTRIBUTE_PAIRS.findall(kv_strings)
        elif kv_strings is None:
            pairs = ()
        else:
            # Already separated pairs retain their permissive split behavior.
            pairs = (
                (parts[0], "", parts[1])
                for kv in kv_strings
                if len(parts := kv.strip().split(None, 1)) == 2
            )
        for column_name, double_quoted, other in pairs:
            if usecols is not None and column_name not in usecols:
                continue

            if other:
                value = other.strip()
                if not value:
                    continue
                if len(value) >= 2 and value[0] in ('"', quote_char) and value[-1] == value[0]:
                    value = value[1:-1]
            else:
                value = double_quoted

            column = extra_columns.get(column_name)
            if column is None:
                column = [missing_value] * n
                extra_columns[column_name] = column
                column[i] = value
            elif column[i] is missing_value:
                column[i] = value
            else:
                # Preserve repeated attributes and the missing-value identity
                # convention, including repeated empty quoted values.
                column[i] = "%s,%s" % (column[i], value)

        if progress_callback is not None:
            completed = i + 1
            if completed % 10_000 == 0 or completed == n:
                progress_callback("attributes", completed, n)

    logger.info("Extracted GTF attributes: %s", list(extra_columns))
    return OrderedDict(extra_columns)
