# Copyright (c) 2015. Mount Sinai School of Medicine
#
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

from __future__ import print_function, division, absolute_import

"""
Columns of a GTF file:

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
    score     - a floating point value indiciating the score of a feature
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