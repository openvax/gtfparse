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
from collections import OrderedDict

import pandas as pd

logger = logging.getLogger(__name__)


VALID_STRICT_MODES = ("best_effort", "warn", "raise")


def _non_empty_unique(series):
    """Return unique non-null, non-empty values from a pandas Series."""
    mask = series.notna() & (series != "")
    return series[mask].unique()


def create_missing_features(
    dataframe,
    unique_keys={},
    extra_columns={},
    missing_value=None,
    *,
    strict="best_effort",
    required_attributes_per_feature=None,
):
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

    unique_keys : dict
        Mapping from feature names to the name of the column which should
        act as a unique key for that feature. Example: {"gene": "gene_id"}

    extra_columns : dict
        By default the constructed feature row will include only the 8
        core columns and its unique key. Any other columns that should
        be included should be associated with the feature name in this
        dict.

    missing_value : any
        Which value to fill in for columns that we don't infer values for.

    strict : {"best_effort", "warn", "raise"}, optional
        Policy when a group lacks required attributes or has inconsistent
        parent metadata (conflicting ``seqname``/``strand`` or missing/conflicting
        ``extra_columns`` values):

        - ``"best_effort"`` (default): synthesize rows with available metadata,
          filling ambiguous or missing extra columns with ``missing_value``.
          Groups that fail ``required_attributes_per_feature`` are skipped.
        - ``"warn"``: log a warning and skip the synthesized row.
        - ``"raise"``: raise ``ValueError``.

    required_attributes_per_feature : dict of str -> iterable of str, optional
        Mapping from feature name to attribute columns that must have at least
        one non-null, non-empty value among a group's child rows in order to
        synthesize that parent row. Example: ``{"transcript": {"exon_id"}}``.

    Returns
    -------
    pandas.DataFrame
        Original rows (converted to pandas if necessary) together with the
        extra rows created for missing features.
    """
    if strict not in VALID_STRICT_MODES:
        raise ValueError("strict must be one of %r, got %r" % (VALID_STRICT_MODES, strict))

    if hasattr(dataframe, "to_pandas"):
        dataframe = dataframe.to_pandas()

    extra_dataframes = []

    existing_features = set(dataframe["feature"])
    existing_columns = set(dataframe.columns)
    required_map = required_attributes_per_feature or {}

    for feature_name, groupby_key in unique_keys.items():
        if feature_name in existing_features:
            logger.info("Feature '%s' already exists in GTF data", feature_name)
            continue
        if groupby_key not in existing_columns:
            msg = "Cannot synthesize feature '%s': key column '%s' does not exist in GTF" % (
                feature_name,
                groupby_key,
            )
            if strict == "warn":
                logger.warning(msg)
                continue
            raise ValueError(msg)

        logger.info("Creating rows for missing feature '%s'", feature_name)

        # User specifies which non-required columns we should try to infer
        # values for
        feature_columns = list(extra_columns.get(feature_name, []))
        missing_extra = [c for c in feature_columns if c not in existing_columns]
        if missing_extra:
            msg = "Column '%s' does not exist in GTF, columns = %s" % (
                missing_extra[0],
                existing_columns,
            )
            if strict == "warn":
                logger.warning("Skipping synthesized '%s' rows: %s", feature_name, msg)
                continue
            raise ValueError(msg)

        required_attrs = sorted(required_map.get(feature_name, ()))
        missing_required_cols = [c for c in required_attrs if c not in existing_columns]
        if missing_required_cols:
            msg = (
                "Cannot synthesize feature '%s': required attribute column(s) %s "
                "do not exist in GTF" % (feature_name, missing_required_cols)
            )
            if strict == "raise":
                raise ValueError(msg)
            if strict == "warn":
                logger.warning(msg)
            continue

        # don't include rows where the groupby key was missing; the mask keeps
        # the frame's own index, so filtered or re-sorted frames work (#104)
        keys = dataframe[groupby_key]
        has_key = keys.notna() & (keys != "")
        row_groups = dataframe[has_key].groupby(groupby_key)

        feature_values = OrderedDict([(column_name, []) for column_name in dataframe])

        for feature_id, group in row_groups:
            problems = []

            if required_attrs:
                missing_in_group = [
                    attr for attr in required_attrs if len(_non_empty_unique(group[attr])) == 0
                ]
                if missing_in_group:
                    problems.append("missing required attribute(s) %s" % (missing_in_group,))

            if strict != "best_effort":
                unique_seqnames = _non_empty_unique(group["seqname"])
                if len(unique_seqnames) != 1:
                    problems.append("inconsistent seqname values %r" % (list(unique_seqnames),))
                unique_strands = _non_empty_unique(group["strand"])
                if len(unique_strands) != 1:
                    problems.append("inconsistent strand values %r" % (list(unique_strands),))

            inferred_extra = {}
            for column_name in feature_columns:
                if strict == "best_effort":
                    unique_values = group[column_name].dropna().unique()
                else:
                    unique_values = _non_empty_unique(group[column_name])
                if len(unique_values) == 1:
                    inferred_extra[column_name] = unique_values[0]
                elif strict != "best_effort":
                    if len(unique_values) == 0:
                        problems.append("missing extra column '%s'" % (column_name,))
                    else:
                        problems.append(
                            "ambiguous extra column '%s' values %r"
                            % (column_name, list(unique_values))
                        )

            if problems:
                msg = "Cannot faithfully synthesize %s '%s': %s" % (
                    feature_name,
                    feature_id,
                    "; ".join(problems),
                )
                if strict == "raise":
                    raise ValueError(msg)
                if strict == "warn":
                    logger.warning(msg)
                continue

            # fill in the required columns by assuming that this feature
            # is the union of all intervals of other features that were
            # tagged with its unique ID (e.g. union of exons which had a
            # particular gene_id).
            row = dict.fromkeys(dataframe, missing_value)
            row["feature"] = feature_name
            row[groupby_key] = feature_id
            # set the source to 'gtfparse' to indicate that we made this
            # entry up from other data
            row["source"] = "gtfparse"
            row["start"] = group["start"].min()
            row["end"] = group["end"].max()
            # assume that seqname and strand are the same for all other
            # entries in the GTF which shared this unique ID
            row["seqname"] = group["seqname"].iat[0]
            row["strand"] = group["strand"].iat[0]
            row.update(inferred_extra)

            for column_name in dataframe:
                feature_values[column_name].append(row[column_name])

        if feature_values[groupby_key]:
            extra_dataframes.append(pd.DataFrame(feature_values))
    if not extra_dataframes:
        return dataframe.copy()
    return pd.concat([dataframe, *extra_dataframes], ignore_index=True)
