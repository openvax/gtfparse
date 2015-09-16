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
from collections import OrderedDict

import numpy as np
import pandas as pd
from six.moves import intern

from .util import memory_usage
from .line_parsing import parse_gtf_lines


def pandas_series_is_biotype(series):
    """
    Hackishly infer whether a Pandas series is either from
    gene_biotype or transcript_biotype annotations by checking
    whether the 'protein_coding' biotype annotation is among its values.
    """
    return 'protein_coding' in series.values


def read_gtf_as_dict(filename, expand_attribute_column=True):
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
    """
    if not exists(filename):
        raise ValueError("GTF file does not exist: %s" % filename)

    if filename.endswith("gz") or filename.endswith("gzip"):
        with gzip.open(filename, mode="rt") as f:
            return parse_gtf_lines(
                lines=f,
                expand_attribute_column=expand_attribute_column)
    else:
        with open(filename) as f:
            return parse_gtf_lines(
                lines=f,
                expand_attribute_column=expand_attribute_column)

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

def read_gtf_as_dataframe(
        filename,
        expand_attribute_column=True,
        infer_second_column=True):
    """
    Parse GTF and convert it to a DataFrame.

    Parameters
    ----------
    filename : str

    expand_attribute_column : bool

    infer_second_column : bool
    """
    column_dict = read_gtf_as_dict(
        filename=filename,
        expand_attribute_column=expand_attribute_column)

    # add columns one at a time so we can remove potentially duplicated data
    # from the dictionary, saving on memory usage
    df = pd.DataFrame({})
    for column_name, column_values in column_dict.items():
        df[column_name] = column_values
        del column_dict[column_name]

    logging.debug("Memory usage after DataFrame construction: %0.4f MB" % (
        memory_usage(),))

    # very old GTF files use the second column to store the gene biotype
    # others use it to store the transcript biotype and
    # anything beyond release 77+ will use it to store
    # the source of the annotation (e.g. "havana")
    if infer_second_column and pandas_series_is_biotype(df['second_column']):
        # patch this later to either 'transcript_biotype' or 'gene_biotype'
        # depending on what other annotations are present
        column_name = 'biotype'
    else:
        column_name = 'source'
    df[column_name] = df["second_column"]
    del df["second_column"]
    return df


# In addition to the required 8 columns and IDs of genes & transcripts,
# there might also be annotations like 'transcript_biotype' but these aren't
# available for all species/releases.

# TODO: gene and transcript names, as well as biotypes, should have
# the option to be required.
REQUIRED_ATTRIBUTE_COLUMNS = [
    'gene_id',
    'transcript_id',
]

def _dataframe_from_groups(groups, feature):
    """
    Helper function used to construct a missing feature such as 'transcript'
    or 'gene'. For example, a sufficiently old release might only have
    'exon' entries, but these were tagged with which transcript_id and gene_id
    they're associated with. Grouping by transcript_id lets you reconstruct
    the feature='transcript' entries which are normally present in later
    releases.
    """
    start = groups.start.min()
    end = groups.end.max()
    strand = groups.strand.first()
    seqname = groups.seqname.first()

    # Include these columns when they're available. We also want to
    # include gene_id if we're grouping by transcript_id and vice versa.
    conditional_column_names = [
        "gene_name", "gene_biotype", "transcript_name",
        "transcript_biotype", "gene_id", "transcript_id"
    ]

    def pick_protein_id(candidates):
        for c in candidates:
            if c is not None and len(c) > 0:
                return c
        return None

    columns = [seqname, start, end, strand]

    if "protein_id" in groups.first().columns:
        protein_id = groups.protein_id.apply(pick_protein_id)
        columns.append(protein_id)

    for conditional_column_name in conditional_column_names:
        if conditional_column_name in groups.first().columns:
            column = groups[conditional_column_name].first()
            columns.append(column)

    df = pd.concat(columns, axis=1).reset_index()

    # score seems to be always '.' in Ensembl GTF files value,
    # not sure why it's even given
    df['score'] = '.'

    # frame values only make sense for CDS entries, but need this column
    # so we concatenate these rows with the rest of the Ensembl entries
    df['frame'] = '.'

    df['feature'] = feature
    return df

def reconstruct_gene_rows(df):
    gene_id_groups = df.groupby(['gene_id'])
    genes_df = _dataframe_from_groups(gene_id_groups, feature='gene')
    return pd.concat([df, genes_df], ignore_index=True)

def reconstruct_transcript_rows(df):
    transcript_id_groups = df.groupby(['transcript_id'])
    transcripts_df = _dataframe_from_groups(
        transcript_id_groups,
        feature='transcript'
    )
    return pd.concat([df, transcripts_df], ignore_index=True)

def load_gtf_as_dataframe(filename):
    # due to the annoying ambiguity of the second GTF column,
    # figure out if an older GTF's biotype is actually the gene_biotype
    # or transcript_biotype
    if 'biotype' in df.columns:
        assert 'transcript_biotype' not in df.columns, \
            "Inferred 2nd column as biotype but also found transcript_biotype"

        # Initially we could only figure out if the 2nd column was either
        # the source of the annotation or some kind of biotype (either
        # a gene_biotype or transcript_biotype).
        # Now we disambiguate between the two biotypes by checking if
        # gene_biotype is already present in another column. If it is,
        # the 2nd column is the transcript_biotype (otherwise, it's the
        # gene_biotype)
        if 'gene_biotype' in df.columns:
            rename_to = 'transcript_biotype'
        else:
            rename_to = 'gene_biotype'
        df.rename(columns={'biotype': rename_to}, inplace=True)

    for column_name in REQUIRED_ATTRIBUTE_COLUMNS:
        assert column_name in df.columns, \
            "Missing required column '%s', available: %s" % (
                column_name,
                list(sorted(df.columns)))

    # older Ensembl releases only had features:
    #   - exon
    #   - CDS
    #   - start_codon
    #   - stop_codon
    # (And this also applies to other non-Ensembl GTF files.)
    #
    # Might have to manually reconstruct gene & transcript entries
    # by grouping the gene_id and transcript_id columns of existing features

    distinct_features = df.feature.unique()

    if 'gene' not in distinct_features:
        logging.info("Creating entries for feature='gene'")
        df = reconstruct_gene_rows(df)

    if 'transcript' not in distinct_features:
        logging.info("Creating entries for feature='transcript'")
        df = reconstruct_transcript_rows(df)
    return df
