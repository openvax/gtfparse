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

from .attribute_parsing import expand_attribute_strings
from .create_missing_features import create_missing_features
from .parsing_error import ParsingError
from .read_gtf import (
    GENCODE_BIOTYPE_ALIASES,
    INTEGER_VERSION_COLUMNS,
    REQUIRED_COLUMNS,
    parse_gtf,
    parse_gtf_and_expand_attributes,
    parse_gtf_pandas,
    read_gtf,
)

# Per the Python logging HOWTO ("Configuring Logging for a Library"), attach a
# no-op handler to the package logger so that importing gtfparse neither emits
# log output nor reconfigures the root logger; applications opt in to gtfparse
# logs through their own logging configuration.
logging.getLogger(__name__).addHandler(logging.NullHandler())

__version__ = "2.7.1"

__all__ = [
    "GENCODE_BIOTYPE_ALIASES",
    "INTEGER_VERSION_COLUMNS",
    "REQUIRED_COLUMNS",
    "ParsingError",
    "__version__",
    "create_missing_features",
    "expand_attribute_strings",
    "parse_gtf",
    "parse_gtf_and_expand_attributes",
    "parse_gtf_pandas",
    "read_gtf",
]
