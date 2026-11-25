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

import pandas as pd

DEFAULT_REQUIRED_CODING_FEATURES = ("start_codon", "stop_codon")


def transcript_completeness(
    dataframe,
    *,
    required_features=DEFAULT_REQUIRED_CODING_FEATURES,
    coding_feature="CDS",
    transcript_id_column="transcript_id",
):
    """Compute per-transcript completeness of required coding features.

    A transcript with one or more ``coding_feature`` rows (default ``"CDS"``)
    is considered complete when every feature named in ``required_features``
    (default ``("start_codon", "stop_codon")``) also appears for that
    ``transcript_id``. Non-coding transcripts (those without any
    ``coding_feature`` row) are marked ``transcript_complete = True``.

    Parameters
    ----------
    dataframe : pandas.DataFrame or polars.DataFrame
        Parsed GTF table containing ``feature`` and ``transcript_id_column``.

    required_features : sequence of str, optional
        Feature names required whenever a transcript has ``coding_feature``
        rows. Defaults to ``("start_codon", "stop_codon")``.

    coding_feature : str, optional
        Feature name that identifies coding transcripts. Defaults to ``"CDS"``.

    transcript_id_column : str, optional
        Column containing transcript identifiers. Defaults to
        ``"transcript_id"``.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by non-empty ``transcript_id`` values with boolean
        columns ``has_cds`` (or ``has_<coding_feature.lower()>``),
        ``has_<feature>`` for each feature in ``required_features``, and
        ``transcript_complete``.
    """
    if hasattr(dataframe, "to_pandas"):
        dataframe = dataframe.to_pandas()

    if "feature" not in dataframe.columns:
        raise ValueError("DataFrame must contain a 'feature' column")

    required_features = tuple(required_features)
    coding_col = f"has_{coding_feature.lower()}"
    feature_cols = [(feat, f"has_{feat.lower()}") for feat in required_features]
    output_columns = [coding_col, *(col for _, col in feature_cols), "transcript_complete"]

    if transcript_id_column not in dataframe.columns or len(dataframe) == 0:
        empty = pd.DataFrame(
            {col: pd.Series(dtype=bool) for col in output_columns},
            index=pd.Index([], name=transcript_id_column),
        )
        return empty

    tids = dataframe[transcript_id_column]
    valid_mask = tids.notna() & (tids != "")
    if not valid_mask.any():
        empty = pd.DataFrame(
            {col: pd.Series(dtype=bool) for col in output_columns},
            index=pd.Index([], name=transcript_id_column),
        )
        return empty

    valid_tids = tids[valid_mask]
    valid_features = dataframe.loc[valid_mask, "feature"]
    unique_tids = pd.Index(valid_tids.unique(), name=transcript_id_column)

    def _has_feature(feature_name):
        matching = set(valid_tids[valid_features == feature_name].unique())
        return unique_tids.isin(matching)

    data = {coding_col: _has_feature(coding_feature)}
    all_required_present = pd.Series(True, index=unique_tids, dtype=bool)
    for feat, col_name in feature_cols:
        present = _has_feature(feat)
        data[col_name] = present
        all_required_present &= present

    data["transcript_complete"] = (~data[coding_col]) | all_required_present.to_numpy()
    return pd.DataFrame(data, index=unique_tids)


def attach_transcript_completeness(
    dataframe,
    *,
    required_features=DEFAULT_REQUIRED_CODING_FEATURES,
    coding_feature="CDS",
    transcript_id_column="transcript_id",
):
    """Attach per-transcript feature completeness columns to a GTF DataFrame.

    Adds nullable boolean columns ``has_<feature>`` for each feature in
    ``required_features`` (default ``has_start_codon`` and ``has_stop_codon``)
    plus ``transcript_complete``. Rows without a ``transcript_id`` (such as
    ``gene`` rows) receive ``pd.NA``.

    Parameters
    ----------
    dataframe : pandas.DataFrame or polars.DataFrame
        Parsed GTF table containing ``feature`` and ``transcript_id_column``.

    required_features : sequence of str, optional
        Feature names required whenever a transcript has ``coding_feature``
        rows. Defaults to ``("start_codon", "stop_codon")``.

    coding_feature : str, optional
        Feature name that identifies coding transcripts. Defaults to ``"CDS"``.

    transcript_id_column : str, optional
        Column containing transcript identifiers. Defaults to
        ``"transcript_id"``.

    Returns
    -------
    pandas.DataFrame
        Copy of ``dataframe`` with the completeness columns attached.
    """
    if hasattr(dataframe, "to_pandas"):
        dataframe = dataframe.to_pandas()

    required_features = tuple(required_features)
    summary = transcript_completeness(
        dataframe,
        required_features=required_features,
        coding_feature=coding_feature,
        transcript_id_column=transcript_id_column,
    )
    attach_columns = [*(f"has_{feat.lower()}" for feat in required_features), "transcript_complete"]

    out = dataframe.copy()
    if transcript_id_column not in out.columns:
        for col_name in attach_columns:
            out[col_name] = pd.Series(pd.NA, index=out.index, dtype="boolean")
        return out

    tids = out[transcript_id_column]
    valid_mask = tids.notna() & (tids != "")
    for col_name in attach_columns:
        series = pd.Series(pd.NA, index=out.index, dtype="boolean")
        if valid_mask.any():
            series.loc[valid_mask] = tids.loc[valid_mask].map(summary[col_name]).astype("boolean")
        out[col_name] = series
    return out
