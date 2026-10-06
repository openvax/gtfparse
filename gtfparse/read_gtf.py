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
from os.path import exists

import pandas as pd
import pyarrow as pa

from ._polars import require_polars
from ._read_csv import REQUIRED_COLUMNS as REQUIRED_COLUMNS
from ._read_csv import read_fixed_columns
from .attribute_parsing import ATTRIBUTE_PATTERN, expand_attribute_strings

logger = logging.getLogger(__name__)


# GENCODE GTFs use *_type where Ensembl GTFs use *_biotype. Pass this
# (or a superset) as `attribute_aliases` to read_gtf to normalize a
# GENCODE-format GTF onto the Ensembl column names that downstream
# tools like pyensembl expect.
GENCODE_BIOTYPE_ALIASES = {
    "gene_type": "gene_biotype",
    "transcript_type": "transcript_biotype",
}


# Ensembl-style attribute columns that are always integer-valued when
# present. read_gtf casts these from string to pandas nullable Int64 by
# default; pass cast_version_columns=False to keep them as strings.
INTEGER_VERSION_COLUMNS = (
    "gene_version",
    "transcript_version",
    "protein_version",
    "exon_version",
)


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


def parse_with_polars_lazy(
    filepath_or_buffer, split_attributes=True, features=None, fix_quotes_columns=None
):
    """Read eagerly and return an optional Polars LazyFrame.

    As in earlier releases, file I/O happens before this function returns.
    Install ``gtfparse[polars]`` to use this compatibility adapter.
    """
    polars = require_polars()
    return polars.from_pandas(
        parse_gtf(filepath_or_buffer, split_attributes, features, fix_quotes_columns)
    ).lazy()


def parse_gtf(
    filepath_or_buffer,
    split_attributes=True,
    features=None,
    fix_quotes_columns=None,
    *,
    progress_callback=None,
):
    """Read the fixed GTF columns, preserving raw attribute text.

    If split_attributes is True, add an attribute_split list column with
    quoted key/value pairs kept intact. Legacy semicolon cleanup is available
    by explicitly setting fix_quotes_columns; it is disabled by default.
    progress_callback follows read_gtf's contract and emits only the read stage.
    """
    if progress_callback is not None:
        progress_callback("read", 0, None)
    # Apply the normal filter before conversion so discarded null coordinates
    # do not change the retained columns' integer dtypes. Opt-in quote cleanup
    # can change feature names, so preserve its order ahead of filtering.
    read_features = None if fix_quotes_columns else features
    df = read_fixed_columns(filepath_or_buffer, features=read_features)
    for name in fix_quotes_columns or ():
        df[name] = (
            df[name]
            .str.replace(';"', '"', n=1, regex=False)
            .str.replace(";-", "-", n=1, regex=False)
        )
    if features is not None and fix_quotes_columns:
        df = df.loc[df["feature"].isin(set(features))].reset_index(drop=True)
    if split_attributes:
        df["attribute_split"] = df["attribute"].str.findall(ATTRIBUTE_PATTERN)
    if progress_callback is not None:
        progress_callback("read", len(df), len(df))
    return df


def parse_gtf_pandas(*args, **kwargs):
    """Like parse_gtf, returning pandas and reporting read/convert progress."""
    df = parse_gtf(*args, **kwargs)
    progress_callback = kwargs.get("progress_callback")
    if progress_callback is not None:
        progress_callback("convert", 0, None)
    result = df
    if progress_callback is not None:
        progress_callback("convert", len(result), len(result))
    return result


_ATTRIBUTE_BATCH_SIZE = 250_000


def _expand_attributes_as_arrow(attributes, usecols, progress_callback):
    """Keep temporary Python values bounded while retaining Arrow column chunks."""
    total = len(attributes)
    columns = {}
    string_type = pa.large_string()
    empty = pa.scalar("", type=string_type)
    if progress_callback is not None:
        progress_callback("attributes", 0, total)

    def report(stage, completed, _total):
        absolute = start + completed
        if completed and (absolute % 10_000 == 0 or absolute == total):
            progress_callback(stage, absolute, total)

    for start in range(0, total, _ATTRIBUTE_BATCH_SIZE):
        stop = min(start + _ATTRIBUTE_BATCH_SIZE, total)
        values = attributes.iloc[start:stop].to_numpy(dtype=object, na_value=None)
        expanded = expand_attribute_strings(
            values,
            usecols=usecols,
            progress_callback=report if progress_callback is not None else None,
        )
        del values
        missing = columns.keys() - expanded.keys()
        if missing:
            blank = pa.repeat(empty, stop - start)
            for name in missing:
                columns[name].append(blank)
        prefix = None
        for name in list(expanded):
            if name not in columns:
                if start and prefix is None:
                    prefix = pa.repeat(empty, start)
                columns[name] = [prefix] if start else []
            columns[name].append(pa.array(expanded.pop(name), type=string_type))
    return {name: pa.chunked_array(chunks, type=string_type) for name, chunks in columns.items()}


def parse_gtf_and_expand_attributes(
    filepath_or_buffer, restrict_attribute_columns=None, features=None, *, progress_callback=None
):
    """
    Parse a GTF into a pandas DataFrame and then expand
    the 'attribute' column into multiple columns. This expansion happens
    by replacing strings of semi-colon separated key-value values in the
    'attribute' column with one column per distinct key, with a list of
    values for each row (using empty strings for rows where the key didn't occur).

    Parameters
    ----------
    filepath_or_buffer : str or buffer object

    restrict_attribute_columns : list/set of str or None
        If given, then only use these attribute columns.

    features : set or None
        Ignore entries which don't correspond to one of the supplied features

    progress_callback : callable, optional
        Follows read_gtf's callback contract, emitting read and attributes stages.
    """
    df = parse_gtf(
        filepath_or_buffer=filepath_or_buffer,
        features=features,
        split_attributes=False,
        progress_callback=progress_callback,
    )
    if type(restrict_attribute_columns) is str:
        restrict_attribute_columns = {restrict_attribute_columns}
    elif restrict_attribute_columns:
        restrict_attribute_columns = set(restrict_attribute_columns)
    attributes = df.pop("attribute")
    string_dtype = pd.Series([""]).dtype
    arrow_strings = isinstance(string_dtype, pd.StringDtype) and string_dtype.storage == "pyarrow"
    if arrow_strings:
        expanded = _expand_attributes_as_arrow(
            attributes, restrict_attribute_columns, progress_callback
        )
    else:
        expanded = expand_attribute_strings(
            attributes.to_numpy(dtype=object, na_value=None),
            usecols=restrict_attribute_columns,
            progress_callback=progress_callback,
        )
    # Release raw strings before assembling the final pandas columns.
    del attributes
    # Assigning a colliding attribute replaces its fixed column, matching the
    # previous reader, while preserving the fixed columns' original positions.
    for name in list(expanded):
        values = expanded.pop(name)
        if arrow_strings:
            values = pd.array(values, dtype=string_dtype)
        df[name] = values
    return df


def _apply_attribute_aliases(result_df, attribute_aliases):
    """
    Rename alias attribute columns onto canonical names in-place.

    For each (alias -> canonical) pair, in iteration order:
      * if only the alias is present, rename it to the canonical name.
      * if both are present, drop the alias and warn (canonical wins).
      * if neither is present, do nothing.

    When two aliases target the same canonical (e.g. both ``gene_type``
    and a hypothetical ``gene_kind`` map to ``gene_biotype``), the first
    rename in iteration order wins; subsequent aliases targeting an
    already-renamed canonical are treated as collisions, dropped, and
    warned about.
    """
    if not attribute_aliases:
        return result_df
    columns_present = set(result_df.columns)
    rename_map = {}
    drop_aliases = []
    for alias, canonical in attribute_aliases.items():
        if alias not in columns_present:
            continue
        if canonical in columns_present:
            logger.warning(
                "Both alias column '%s' and canonical column '%s' are present; "
                "dropping alias and keeping canonical values.",
                alias,
                canonical,
            )
            drop_aliases.append(alias)
        else:
            rename_map[alias] = canonical
            # Reflect the rename in the running column set so a later
            # alias mapping to the same canonical sees the collision
            # instead of silently producing a duplicate-named column.
            columns_present.discard(alias)
            columns_present.add(canonical)
    if drop_aliases:
        result_df = result_df.drop(columns=drop_aliases)
    if rename_map:
        result_df = result_df.rename(columns=rename_map)
    return result_df


def _cast_version_columns(result_df, version_columns=INTEGER_VERSION_COLUMNS):
    """
    Cast known Ensembl *_version attribute columns from strings to
    pandas nullable Int64 in-place. Missing/empty values become pd.NA.
    """
    for column_name in version_columns:
        if column_name not in result_df.columns:
            continue
        result_df[column_name] = pd.to_numeric(
            result_df[column_name].replace("", None), errors="coerce"
        ).astype("Int64")
    return result_df


def read_gtf(
    filepath_or_buffer,
    expand_attribute_column=True,
    infer_biotype_column=False,
    column_converters={},
    column_cast_types={},
    usecols=None,
    features=None,
    result_type="pandas",
    attribute_aliases=None,
    cast_version_columns=True,
    *,
    progress_callback=None,
):
    """
    Parse a GTF into a Polars DataFrame, pandas DataFrame, or dictionary.

    Parameters
    ----------
    filepath_or_buffer : str or buffer object
        Path to GTF file (may be gzip compressed) or buffer object
        such as StringIO

    expand_attribute_column : bool
        Replace strings of semi-colon separated key-value values in the
        'attribute' column with one column per distinct key, with a list of
        values for each row (using empty strings for rows where the key didn't occur).

    infer_biotype_column : bool
        Due to the annoying ambiguity of the second GTF column across multiple
        Ensembl releases, figure out if an older GTF's source column is actually
        the gene_biotype or transcript_biotype.

    column_converters : dict, optional
        Dictionary mapping column names to conversion functions. Will replace
        empty strings with None and otherwise passes them to given conversion
        function.

    column_cast_types : dict, optional
        Dictionary mapping column names to dtypes. Will cast columns to given
        pandas-compatible types.

    usecols : list of str or None
        Restrict which columns are loaded to the give set. If None, then
        load all columns.

    features : set of str or None
        Drop rows which aren't one of the features in the supplied set

    result_type : One of 'pandas', 'dict', or 'polars' (requires the optional extra)
        Return a pandas DataFrame by default. Polars output requires
        ``pip install 'gtfparse[polars]'``; dictionary output is also supported.

    attribute_aliases : dict of str -> str, optional
        Maps alias attribute names onto canonical ones. After attributes
        are expanded into columns, each alias column is renamed to its
        canonical name when the canonical column is absent. If both are
        present the alias is dropped and a warning is logged. Pass
        `GENCODE_BIOTYPE_ALIASES` to normalize a GENCODE GTF's
        `gene_type`/`transcript_type` onto Ensembl's
        `gene_biotype`/`transcript_biotype`.

    cast_version_columns : bool
        When True (default), cast the well-known integer version
        attribute columns (`gene_version`, `transcript_version`,
        `protein_version`, `exon_version`) from strings to pandas
        nullable Int64 when present. Set to False to keep them as
        strings.

    progress_callback : callable, optional
        Called synchronously as ``callback(stage, completed, total)``. Counts
        are rows retained after feature filtering. Stages occur in order:
        ``read``, ``attributes`` (if expanded), and ``convert``. Read and convert
        emit ``(stage, 0, None)`` before work and ``(stage, rows, rows)`` after
        success; None means the total is unknown, so show an indeterminate
        indicator. Attribute expansion reports a known total at the start,
        every 10,000 rows, and at completion (one event for zero rows).
        Exceptions, including those raised by the callback, propagate without
        reporting completion of the failed stage. Default: no callback or UI.
    """
    if result_type not in ("polars", "pandas", "dict"):
        raise ValueError("result_type must be one of 'polars', 'pandas', or 'dict'")
    if result_type == "polars":
        require_polars()
    if type(filepath_or_buffer) is str and not exists(filepath_or_buffer):
        raise ValueError("GTF file does not exist: %s" % filepath_or_buffer)

    # If usecols asks for a canonical column that's only present in the
    # GTF under an alias name, expand the parse-time column filter to
    # also pull the alias through — otherwise it gets dropped at parse
    # time before _apply_attribute_aliases can see it. The end-of-function
    # usecols filter still narrows the result down to the canonical name.
    parse_usecols = usecols
    if usecols is not None and attribute_aliases:
        usecols_set = set(usecols)
        parse_usecols = set(usecols_set)
        for alias, canonical in attribute_aliases.items():
            if canonical in usecols_set:
                parse_usecols.add(alias)

    if expand_attribute_column:
        result_df = parse_gtf_and_expand_attributes(
            filepath_or_buffer,
            restrict_attribute_columns=parse_usecols,
            features=features,
            progress_callback=progress_callback,
        )
    else:
        # When the caller opts out of attribute expansion they want the raw
        # 'attribute' column verbatim — no need to also produce the
        # 'attribute_split' helper that parse_gtf adds by default.
        result_df = parse_gtf(
            filepath_or_buffer,
            features=features,
            split_attributes=False,
            progress_callback=progress_callback,
        )

    if progress_callback is not None:
        progress_callback("convert", 0, None)
    if column_converters or column_cast_types:

        def wrap_to_always_accept_none(f):
            def wrapped_fn(x):
                if x is None or x == "":
                    return None
                else:
                    return f(x)

            return wrapped_fn

        column_names = set(column_converters.keys()).union(column_cast_types.keys())
        for column_name in column_names:
            if column_name in column_converters:
                column_fn = wrap_to_always_accept_none(column_converters[column_name])
                result_df[column_name] = result_df[column_name].apply(column_fn)

            if column_name in column_cast_types:
                column_type = column_cast_types[column_name]
                result_df[column_name] = result_df[column_name].astype(column_type)

    # Rename alias attribute columns onto their canonical names. Done before
    # infer_biotype_column so an aliased gene_biotype/transcript_biotype is
    # visible to the inference logic.
    result_df = _apply_attribute_aliases(result_df, attribute_aliases)

    # Cast Ensembl *_version columns from strings to nullable integers so
    # downstream consumers (e.g. pyensembl) don't have to int(...) themselves.
    if cast_version_columns:
        result_df = _cast_version_columns(result_df)

    # Hackishly infer whether the values in the 'source' column of this GTF
    # are actually representing a biotype by checking for the most common
    # gene_biotype and transcript_biotype value 'protein_coding'
    if infer_biotype_column:
        unique_source_values = set(result_df["source"])
        if "protein_coding" in unique_source_values:
            column_names = set(result_df.columns)
            # Disambiguate between the two biotypes by checking if
            # gene_biotype is already present in another column. If it is,
            # the 2nd column is the transcript_biotype (otherwise, it's the
            # gene_biotype)
            if "gene_biotype" not in column_names:
                logger.info("Using column 'source' to replace missing 'gene_biotype'")
                result_df["gene_biotype"] = result_df["source"]
            if "transcript_biotype" not in column_names:
                logger.info("Using column 'source' to replace missing 'transcript_biotype'")
                result_df["transcript_biotype"] = result_df["source"]

    if usecols is not None:
        column_names = set(result_df.columns)
        valid_columns = [c for c in usecols if c in column_names]
        result_df = result_df[valid_columns]

    if result_type == "pandas":
        result = result_df
    elif result_type == "polars":
        result = require_polars().from_pandas(result_df)
    else:
        result = result_df.to_dict()
    if progress_callback is not None:
        progress_callback("convert", len(result_df), len(result_df))
    return result
