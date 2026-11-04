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

"""
The nine tab-separated columns of a GTF file, shared by the reader and writer.

    seqname   - name of the chromosome or scaffold; chromosome names
                without a 'chr' in Ensembl (but sometimes with a 'chr'
                elsewhere)
    source    - name of the program that generated this feature, or
                the data source (database or project name)
    feature   - feature type name.
                Features currently in Ensembl GTFs:
                    gene
                    transcript
                    exon
                    CDS
                    Selenocysteine
                    start_codon
                    stop_codon
                    UTR
                Older Ensembl releases may be missing some of these features.
    start     - start position of the feature, with sequence numbering
                starting at 1.
    end       - end position of the feature, with sequence numbering
                starting at 1.
    score     - a floating point value indicating the score of a feature
    strand    - defined as + (forward) or - (reverse).
    frame     - one of '0', '1' or '2'. Frame indicates the number of base pairs
                before you encounter a full codon. '0' indicates the feature
                begins with a whole codon. '1' indicates there is an extra
                base (the 3rd base of the prior codon) at the start of this feature.
                '2' indicates there are two extra bases (2nd and 3rd base of the
                prior exon) before the first codon. All values are given with
                relation to the 5' end.
    attribute - a semicolon-separated list of tag-value pairs (separated by a space),
                providing additional information about each feature. A key can be
                repeated multiple times.

(from ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/README)
"""

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

# Name of the ninth, raw attribute column, which read_gtf keeps unexpanded
# when called with expand_attribute_column=False.
RAW_ATTRIBUTE_COLUMN = "attribute"

# All nine GTF columns in file order, as parse_gtf returns them.
REQUIRED_COLUMNS = [*GTF_FIXED_COLUMNS, RAW_ATTRIBUTE_COLUMN]

# Fixed columns with few distinct values, which the reader stores as categoricals.
CATEGORICAL_COLUMNS = ("seqname", "source", "feature", "strand")

# GTF uses a single dot to denote a missing value in the fixed columns.
MISSING_VALUE = "."
