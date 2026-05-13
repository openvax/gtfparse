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
import polars

from .attribute_parsing import expand_attribute_strings
from .parsing_error import ParsingError

logging.basicConfig(level=logging.INFO)
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


DEFAULT_COLUMN_DTYPES = {
    "seqname": polars.Categorical,
    "source": polars.Categorical,
    "start": polars.Int64,
    "end": polars.Int64,
    "score": polars.Float32,
    "feature": polars.Categorical,
    "strand": polars.Categorical,
    "frame": polars.UInt32,
}


def parse_with_polars_lazy(
    filepath_or_buffer, split_attributes=True, features=None, fix_quotes_columns=["attribute"]
):
    # use a global string cache so that all strings get intern'd into
    # a single numbering system
    polars.enable_string_cache()
    kwargs = {
        "has_header": False,
        "separator": "\t",
        "comment_prefix": "#",
        "null_values": ".",
        "schema_overrides": DEFAULT_COLUMN_DTYPES,
    }
    try:
        df = polars.read_csv(filepath_or_buffer, new_columns=REQUIRED_COLUMNS, **kwargs).lazy()
    except polars.exceptions.ShapeError as err:
        raise ParsingError("Wrong number of columns") from err

    # Drop empty lines that may appear as all-null rows
    df = df.filter(polars.col("seqname").is_not_null())

    df = df.with_columns(
        [polars.col("frame").fill_null(0), polars.col("attribute").str.replace_all('"', "'")]
    )

    for fix_quotes_column in fix_quotes_columns:
        # Catch mistaken semicolons by replacing "xyz;" with "xyz"
        # Required to do this since the Ensembl GTF for Ensembl
        # release 78 has mistakes such as:
        #   gene_name = "PRAMEF6;" transcript_name = "PRAMEF6;-201"
        df = df.with_columns(
            [polars.col(fix_quotes_column).str.replace(';"', '"').str.replace(";-", "-")]
        )

    if features is not None:
        features = sorted(set(features))
        df = df.filter(polars.col("feature").is_in(features))

    if split_attributes:
        df = df.with_columns([polars.col("attribute").str.split(";").alias("attribute_split")])
    return df


def parse_gtf(
    filepath_or_buffer, split_attributes=True, features=None, fix_quotes_columns=["attribute"]
):
    df_lazy = parse_with_polars_lazy(
        filepath_or_buffer=filepath_or_buffer,
        split_attributes=split_attributes,
        features=features,
        fix_quotes_columns=fix_quotes_columns,
    )
    return df_lazy.collect()


def parse_gtf_pandas(*args, **kwargs):
    return parse_gtf(*args, **kwargs).to_pandas()


def parse_gtf_and_expand_attributes(
    filepath_or_buffer, restrict_attribute_columns=None, features=None
):
    """
    Parse lines into column->values dictionary and then expand
    the 'attribute' column into multiple columns. This expansion happens
    by replacing strings of semi-colon separated key-value values in the
    'attribute' column with one column per distinct key, with a list of
    values for each row (using None for rows where key didn't occur).

    Parameters
    ----------
    filepath_or_buffer : str or buffer object

    chunksize : int

    restrict_attribute_columns : list/set of str or None
        If given, then only use these attribute columns.

    features : set or None
        Ignore entries which don't correspond to one of the supplied features
    """
    df = parse_gtf(filepath_or_buffer=filepath_or_buffer, features=features, split_attributes=True)
    if type(restrict_attribute_columns) is str:
        restrict_attribute_columns = {restrict_attribute_columns}
    elif restrict_attribute_columns:
        restrict_attribute_columns = set(restrict_attribute_columns)
    df.drop_in_place("attribute")
    attribute_pairs = df.drop_in_place("attribute_split")
    return df.with_columns(
        [
            polars.Series(k, vs)
            for (k, vs) in expand_attribute_strings(attribute_pairs).items()
            if restrict_attribute_columns is None or k in restrict_attribute_columns
        ]
    )


def _apply_attribute_aliases(result_df, attribute_aliases):
    """
    Rename alias attribute columns onto canonical names in-place.

    For each (alias -> canonical) pair:
      * if only the alias is present, rename it to the canonical name.
      * if both are present, drop the alias and warn (canonical wins).
      * if neither is present, do nothing.
    """
    if not attribute_aliases:
        return result_df
    existing = set(result_df.columns)
    rename_map = {}
    drop_aliases = []
    for alias, canonical in attribute_aliases.items():
        if alias not in existing:
            continue
        if canonical in existing:
            logger.warning(
                "Both alias column '%s' and canonical column '%s' are present; "
                "dropping alias and keeping canonical values.",
                alias,
                canonical,
            )
            drop_aliases.append(alias)
        else:
            rename_map[alias] = canonical
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
    result_type="polars",
    attribute_aliases=None,
    cast_version_columns=True,
):
    """
    Parse a GTF into a dictionary mapping column names to sequences of values.

    Parameters
    ----------
    filepath_or_buffer : str or buffer object
        Path to GTF file (may be gzip compressed) or buffer object
        such as StringIO

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

    column_cast_types : dict, optional
        Dictionary mapping column names to dtypes. Will cast columns to given
        Polars types.

    usecols : list of str or None
        Restrict which columns are loaded to the give set. If None, then
        load all columns.

    features : set of str or None
        Drop rows which aren't one of the features in the supplied set

    result_type : One of 'polars', 'pandas', or 'dict'
        Default behavior is to return a Polars DataFrame, but will convert to
        Pandas DataFrame or dictionary if specified.

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
    """
    if type(filepath_or_buffer) is str and not exists(filepath_or_buffer):
        raise ValueError("GTF file does not exist: %s" % filepath_or_buffer)

    if expand_attribute_column:
        result_df = parse_gtf_and_expand_attributes(
            filepath_or_buffer, restrict_attribute_columns=usecols, features=features
        )
    else:
        result_df = parse_gtf(result_df, features=features)

    # converting back to pandas here because Polars bugs manifest
    # as `pyo3_runtime.PanicException: assertion `left == right` failed: impl error`
    # and are generally insane to chase down
    result_df = result_df.to_pandas()
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
                logging.info("Using column 'source' to replace missing 'gene_biotype'")
                result_df["gene_biotype"] = result_df["source"]
            if "transcript_biotype" not in column_names:
                logging.info("Using column 'source' to replace missing 'transcript_biotype'")
                result_df["transcript_biotype"] = result_df["source"]

    if usecols is not None:
        column_names = set(result_df.columns)
        valid_columns = [c for c in usecols if c in column_names]
        result_df = result_df[valid_columns]

    result = None
    if result_type == "pandas":
        result = result_df
    elif result_type == "polars":
        result = polars.from_pandas(result_df)
    elif result_type == "dict":
        result = result_df.to_dict()
    return result
