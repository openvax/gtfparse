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
import logging
from collections import OrderedDict

from six.moves import intern
import numpy as np

from .util import memory_usage
from .attribute_parsing import expand_attribute_strings
from .parsing_error import ParsingError

def parse_gtf_lines(lines, expand_attribute_column=True):
    """
    Parse the lines of GTF file into a dictionary mapping
    column names to sequences of values.

    Parameters
    ----------
    lines : sequence or generator
        Any iterable object which contains strings representing the
        lines of a GTF file.

    expand_attribute_column : bool
        Replace strings of semi-colon separated key-value values in the
        'attribute' column with one column per distinct key, with a list of
        values for each row (using None for rows where key didn't occur).
    """
    logging.debug("Memory usage before GTF parsing: %0.4f MB" % memory_usage())
    seqname_values = []
    source_values = []
    feature_values = []
    start_values = []
    end_values = []
    score_values = []
    strand_values = []
    frame_values = []
    attribute_values = []

    for i, line in enumerate(lines):
        if not line or line[0] == "#":
            continue
        fields = line.split("\t", 9)
        if len(fields) != 9:
            raise ParsingError("Wrong number of fields %d (expected 9)" % len(fields))
        # GTF columns:
        # 1) seqname: str ("1", "X", "chrX", etc...)
        # 2) source : str
        #      Different versions of GTF use second column as of:
        #      (a) gene biotype
        #      (b) transcript biotype
        #      (c) the annotation source
        #      See: https://www.biostars.org/p/120306/#120321
        # 3) feature : str ("gene", "transcript", &c)
        # 4) start : int
        # 5) end : int
        # 6) score : float or "."
        # 7) strand : "+", "-", or "."
        # 8) frame : 0, 1, 2 or "."
        # 9) attribute : key-value pairs separated by semicolons
        # (see more complete description in docstring at top of file)
        seq, source, feature, start, end, score, strand, frame, attr = fields
        seqname_values.append(intern(str(seq)))
        source_values.append(intern(str(source)))
        feature_values.append(intern(str(feature)))

        # start/end values will be converted to integers collectively
        start_values.append(start)
        end_values.append(end)

        # scores will be converted to floats collectively
        score_values.append("nan" if score == "." else score)

        strand_values.append(intern(str(strand)))
        frame_values.append(intern(str(frame)))

        # Catch mistaken semicolons by replacing "xyz;" with "xyz"
        # Required to do this since the Ensembl GTF for Ensembl release 78 has
        # gene_name = "PRAMEF6;"
        # transcript_name = "PRAMEF6;-201"
        attribute_values.append(attr.replace(';\"', '\"').replace(";-", "-"))
    logging.debug("Memory usage after GTF parsing: %0.4f MB" % memory_usage())

    result_dict = OrderedDict([
        ("seqname", seqname_values),
        ("source", source_values),
        ("feature", feature_values),
        ("start", np.array(start_values, dtype=np.int64)),
        ("end", np.array(end_values, dtype=np.int64)),
        ("score", np.array(score_values, dtype=np.float32)),
        ("strand", strand_values),
        ("frame", frame_values)
    ])
    if expand_attribute_column:
        # remove references to the columns we converted since their data is
        # being duplicated and the increased memory usage of attribute expansion
        # might push us over e.g. the 3gb Travis limit
        del start_values
        del end_values
        del score_values
        result_dict.update(expand_attribute_strings(attribute_values))
    else:
        result_dict["attribute"] = attribute_values
    return result_dict
