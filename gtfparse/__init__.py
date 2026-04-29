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

from .attribute_parsing import expand_attribute_strings
from .create_missing_features import create_missing_features
from .parsing_error import ParsingError
from .write_gtf import write_gtf
from .read_gtf import (
    REQUIRED_COLUMNS,
    parse_gtf,
    parse_gtf_and_expand_attributes,
    parse_gtf_pandas,
    read_gtf,
)

__version__ = "2.6.3"

__all__ = [
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
