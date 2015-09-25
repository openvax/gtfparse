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
from os.path import exists
import gzip


import numpy as np
import pandas as pd
from six.moves import intern

from .create_missing_features import create_missing_features
from .util import memory_usage
from .line_parsing import parse_gtf_lines


def read_gtf_as_dict(
        filename,
        expand_attribute_column=True,
        infer_biotype_column=False,
        column_converters={}):
    """
    Parse a GTF into a dictionary mapping column names to sequences of values.

    Parameters
    ----------
    filename : str
        Name of GTF file (may be gzip compressed)

    expand_attribute_column : bool
        Replace strings of semi-colon separated key-value values in the
        'attribute' column with one column per distinct key, with a list of
        values for each row (using None for rows where key didn't occur).

    infer_biotype_column : bool
        Due to the annoying ambiguity of the second GTF column across multiple
        Ensembl releases, figure out if an older GTF's source column is actually
        the gene_biotype or transcript_biotype.

    column_converters : dict, optional
        Dictionary mapping column names to conversion functions. Will replace
        empty strings with None and otherwise passes them to given conversion
        function.
    """
    if not exists(filename):
        raise ValueError("GTF file does not exist: %s" % filename)

    if filename.endswith("gz") or filename.endswith("gzip"):
        with gzip.open(filename, mode="rt") as f:
            result_dict = parse_gtf_lines(
                lines=f,
                expand_attribute_column=expand_attribute_column)
    else:
        with open(filename) as f:
            result_dict = parse_gtf_lines(
                lines=f,
                expand_attribute_column=expand_attribute_column)

    for column_name, column_type in list(column_converters.items()):
        result_dict[column_name] = [
            column_type(string_value) if len(string_value) else None
            for string_value
            in result_dict[column_name]
        ]
    # Hackishly infer whether the values in the 'source' column of this GTF
    # are actually representing a biotype by checking for the most common
    # gene_biotype and transcript_biotype value 'protein_coding'
    if infer_biotype_column and "protein_coding" in result_dict["source"]:
        # Disambiguate between the two biotypes by checking if
        # gene_biotype is already present in another column. If it is,
        # the 2nd column is the transcript_biotype (otherwise, it's the
        # gene_biotype)
        column_names = set(result_dict.keys())
        if "gene_biotype" not in column_names:
            result_dict["gene_biotype"] = result_dict["source"]
        if "transcript_biotype" not in column_names:
            result_dict["transcript_biotype"] = result_dict["source"]

    return result_dict

def read_gtf_as_dataframe(
        filename,
        expand_attribute_column=True,
        infer_biotype_column=False,
        column_converters={}):
    """
    Parse GTF and convert it to a DataFrame.

    Parameters
    ----------
    filename : str
        Name of GTF file (may be gzip compressed)

    expand_attribute_column : bool
        Replace strings of semi-colon separated key-value values in the
        'attribute' column with one column per distinct key, with a list of
        values for each row (using None for rows where key didn't occur).

    infer_biotype_column : bool
        Due to the annoying ambiguity of the second GTF column across multiple
        Ensembl releases, figure out if an older GTF's source column is actually
        the gene_biotype or transcript_biotype.

    column_converters : dict, optional
        Dictionary mapping column names to conversion functions. Will replace
        empty strings with None and otherwise passes them to given conversion
        function.
    """
    gtf_dict = read_gtf_as_dict(
        filename=filename,
        expand_attribute_column=expand_attribute_column,
        infer_biotype_column=infer_biotype_column,
        column_converters=column_converters)

    # add columns one at a time so we can remove potentially duplicated data
    # from the dictionary, saving on memory usage
    df = pd.DataFrame({})
    for column_name, column_values in list(gtf_dict.items()):
        df[column_name] = column_values
        del gtf_dict[column_name]

    logging.debug("Memory usage after DataFrame construction: %0.4f MB" % (
        memory_usage(),))

    return df
