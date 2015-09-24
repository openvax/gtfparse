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

import pandas as pd

def create_missing_features(dataframe, feature_name_to_unique_key_dict):
    """
    Helper function used to construct a missing feature such as 'transcript'
    or 'gene'. Some GTF files only have 'exon' and 'CDS' entries, but have
    transcript_id and gene_id annotations which allow us to construct those
    missing features.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        Should contain at least the core GTF columns, such as "seqname",
        "start", and "end"

    feature_name_to_unique_key_dict : dict
        Mapping from feature names to the name of the column which should
        act as a unique key for that feature. Example: {"gene": "gene_id"}

    Returns original dataframe along with all extra rows created for missing
    features.
    """
    extra_dataframes = []
    all_column_names = list(dataframe.keys())
    # column names for which we have special logic for inferring
    # their values
    special_column_names = {"start", "end", "seqname", "featue"}
    missing_column_names = {
        column_name
        for column_name in all_column_names
        if column_name not in special_column_names
    }
    unique_features = set(dataframe["feature"])
    for (feature_name, groupby_key) in feature_name_to_unique_key_dict.items():
        if feature_name in unique_features:
            logging.warn(
                "Feature '%s' already exists in GTF data" % feature_name)
        logging.info("Creating rows for missing feature '%s'" % feature_name)
        groups = dataframe.groupby(groupby_key)
        # Each group corresponds to a unique feature entry for which the
        # other columns may or may not be uniquely defined. Start off by
        # assuming the values for every column are missing and fill them in
        # where possible.
        feature_values = OrderedDict([
            (column_name, [None] * len(groups))
            for column_name in all_column_names
        ])
        for i, (feature_id, group) in enumerate(groups):
            feature_values["start"][i] = group["start"].min()
            feature_values["end"][i] = group["end"].max()
            feature_values["seqname"][i] = group["seqname"].irow(0)
            feature_values["feature"][i] = feature_name
            for column_name in missing_column_names:
                # expect that all entries related to a reconstructed featue
                # are related and are thus within the same interval of
                # positions on the same chromosome
                unique_values = set(group[column_name].dropna())
                if len(unique_values) == 1:
                    feature_values[column_name][i] = unique_values.pop()

        extra_dataframes.append(pd.DataFrame(feature_values))
    return pd.concat([dataframe] + extra_dataframes, ignore_index=True)
