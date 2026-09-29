from io import StringIO

import pytest

from gtfparse import expand_attribute_strings, parse_gtf, read_gtf


def gtf_line(attributes):
    return "1\ttest\tgene\t1\t100\t.\t+\t.\t%s\n" % attributes


@pytest.mark.parametrize(
    "value",
    [
        "two words",
        "1 (assigned to previous version 5)",
        "part1; part2",
        "A;B",
        "PRAMEF6;-201",
        "O'Brien",
        "  padded  ",
        "",
    ],
)
def test_quoted_values_survive_all_parsing_paths(value):
    attributes = 'description "%s"; gene_id "G1";' % value
    expected = {"description": [value], "gene_id": ["G1"]}
    assert expand_attribute_strings([attributes]) == expected
    parsed = parse_gtf(StringIO(gtf_line(attributes)))
    assert parsed["attribute"][0] == attributes
    assert expand_attribute_strings(parsed["attribute_split"]) == expected
    result = read_gtf(StringIO(gtf_line(attributes)))
    assert result["description"].to_list() == [value]
    assert result["gene_id"].to_list() == ["G1"]


def test_single_quotes_unquoted_values_and_repeated_tags():
    attributes = "gene_id 'G1'; note 'one; two'; count 0; tag basic; tag more words;"
    expected = {
        "gene_id": ["G1"],
        "note": ["one; two"],
        "count": ["0"],
        "tag": ["basic,more words"],
    }
    assert expand_attribute_strings([attributes]) == expected
    assert read_gtf(StringIO(gtf_line(attributes)))["tag"].to_list() == ["basic,more words"]


def test_whitespace_missing_attributes_and_restricted_columns():
    attributes = [None, "", '; malformed; gene_id    "G1"; note "x; y";', 'gene_id "G2";']
    assert expand_attribute_strings(attributes, usecols={"gene_id"}) == {
        "gene_id": ["", "", "G1", "G2"]
    }
    assert expand_attribute_strings(['gene_id\t"G1";']) == {"gene_id": ["G1"]}


def test_raw_attributes_are_preserved():
    attributes = 'gene_id "G1"; note "O\'Brien; part2;-suffix";'
    result = read_gtf(StringIO(gtf_line(attributes)), expand_attribute_column=False)
    assert result["attribute"].to_list() == [attributes]


def test_missing_attribute_column_values():
    text = gtf_line(".") + gtf_line('gene_id "G1";')
    result = read_gtf(StringIO(text))
    assert result["gene_id"].to_list() == ["", "G1"]


def test_legacy_semicolon_cleanup_is_explicit():
    attributes = 'gene_name "PRAMEF6;"; transcript_name "PRAMEF6;-201";'
    result = parse_gtf(StringIO(gtf_line(attributes)), fix_quotes_columns=["attribute"])
    assert result["attribute"][0] == 'gene_name "PRAMEF6"; transcript_name "PRAMEF6-201";'


def test_pre_split_attributes_keep_missing_and_repeated_values():
    rows = [["broken", 'gene_id "G1"', 'tag "first; value"', 'tag "second"'], []]
    assert expand_attribute_strings(rows) == {
        "gene_id": ["G1", ""],
        "tag": ["first; value,second", ""],
    }
