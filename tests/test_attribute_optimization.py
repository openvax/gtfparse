"""Differential coverage for the 3.0.0 attribute parser's permissive grammar."""

import random
import re
from collections import OrderedDict
from io import StringIO

import pytest

from gtfparse import expand_attribute_strings, read_gtf

# Freeze the pre-optimization tokenization independently of the captured regex.
LEGACY_PATTERN = re.compile(r"""[^\s;]+[ \t]+(?:"[^"]*"|'[^']*'|[^;]+)""")


def legacy_expand(rows, quote_char, missing_value, usecols):
    columns = OrderedDict()
    for i, row in enumerate(rows):
        pairs = LEGACY_PATTERN.findall(row) if isinstance(row, str) else row or ()
        for pair in pairs:
            parts = pair.strip().split(None, 1)
            if len(parts) != 2:
                continue
            name, value = parts
            if usecols is not None and name not in usecols:
                continue
            if len(value) >= 2 and value[0] in ('"', quote_char) and value[-1] == value[0]:
                value = value[1:-1]
            if name not in columns:
                columns[name] = [missing_value] * len(rows)
            previous = columns[name][i]
            columns[name][i] = value if previous is missing_value else "%s,%s" % (previous, value)
    return columns


@pytest.mark.parametrize("quote_char", ["'", '"', "~", ""])
@pytest.mark.parametrize("missing_value", ["", None, "missing"])
@pytest.mark.parametrize("usecols", [None, set(), {"a", "tag", "gene_id"}])
def test_optimized_expansion_matches_legacy_grammar(quote_char, missing_value, usecols):
    rows = [
        None,
        "",
        'gene_id ""; gene_id "G1"; tag "a; b"; tag "";',
        "a '  padded  '; a '~custom~'; a ~custom~;",
        'a "unterminated; b value; a "close"trailing;',
        "a    ; b\t\t; a \u00a0 trimmed \r\n; malformed;",
        "a \"é\n\t; quote'\"; a 'single \" quote';",
        ["broken", 'a ""', 'a "next"', "gene_id\t'G1'", " a    padded  "],
        [],
    ]
    rng = random.Random(85)
    alphabet = "ab_ .;\t\n\r\"'~é\u00a0"
    rows += ["".join(rng.choices(alphabet, k=rng.randrange(100))) for _ in range(500)]
    expected = legacy_expand(rows, quote_char, missing_value, usecols)
    actual = expand_attribute_strings(
        rows, quote_char=quote_char, missing_value=missing_value, usecols=usecols
    )
    assert isinstance(actual, OrderedDict)
    assert list(actual) == list(expected)
    assert actual == expected


def test_expanded_attribute_collisions_preserve_position_and_missing_values():
    text = (
        '1\ttest\tgene\t1\t10\t.\t+\t.\tgene_id "G1"; start "7"; note "a";\n'
        '1\ttest\tgene\t11\t20\t.\t+\t.\tgene_id "G2"; note "b";\n'
    )
    frame = read_gtf(StringIO(text))
    assert list(frame) == [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "gene_id",
        "note",
    ]
    assert frame["start"].tolist() == ["7", ""]
    assert frame["note"].tolist() == ["a", "b"]


@pytest.mark.parametrize("storage", ["python", "pyarrow"])
def test_expanded_columns_follow_pandas_string_configuration(storage):
    import pandas as pd

    text = (
        '1\ttest\tgene\t1\t10\t.\t+\t.\tgene_id "G1"; note "café; value";\n'
        '1\ttest\tgene\t11\t20\t.\t+\t.\tgene_id "G2";\n'
    )
    with pd.option_context("mode.string_storage", storage):
        frame = read_gtf(StringIO(text))
        pd.testing.assert_series_equal(frame["gene_id"], pd.Series(["G1", "G2"], name="gene_id"))
        pd.testing.assert_series_equal(frame["note"], pd.Series(["café; value", ""], name="note"))
        assert pd.get_option("mode.string_storage") == storage


@pytest.mark.parametrize("storage", ["python", "pyarrow"])
def test_batched_expansion_preserves_late_columns_gaps_collisions_and_progress(
    storage, monkeypatch
):
    import importlib

    import pandas as pd

    reader = importlib.import_module("gtfparse.read_gtf")
    monkeypatch.setattr(reader, "_ATTRIBUTE_BATCH_SIZE", 10_000)
    attrs = ['gene_id "G1"; tag "first"; tag "second";'] * 10_000
    attrs += ['gene_id "G2"; note "middle; value";'] * 10_000
    attrs += ['gene_id "G3"; start "42"; late "last";'] * 3
    text = "".join("1\ttest\tgene\t1\t100\t.\t+\t.\t%s\n" % a for a in attrs)
    events = []
    with pd.option_context("mode.string_storage", storage):
        frame = read_gtf(StringIO(text), progress_callback=lambda *event: events.append(event))
    assert list(frame)[-4:] == ["gene_id", "tag", "note", "late"]
    assert list(frame).index("start") == 3
    assert frame["tag"].tolist() == ["first,second"] * 10_000 + [""] * 10_003
    assert frame["note"].tolist() == [""] * 10_000 + ["middle; value"] * 10_000 + [""] * 3
    assert frame["late"].tolist() == [""] * 20_000 + ["last"] * 3
    assert frame["start"].tolist() == [""] * 20_000 + ["42"] * 3
    assert [e for e in events if e[0] == "attributes"] == [
        ("attributes", 0, 20_003),
        ("attributes", 10_000, 20_003),
        ("attributes", 20_000, 20_003),
        ("attributes", 20_003, 20_003),
    ]
