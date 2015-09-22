from gtfparse import expand_attribute_strings
from nose.tools import eq_

def test_attributes_in_quotes():
    attributes = [
        "gene_id \"ENSG001\"; tag \"bogotron\"; version \"1\";",
        "gene_id \"ENSG002\"; tag \"wolfpuppy\"; version \"2\";"
    ]
    parsed_dict = expand_attribute_strings(attributes)
    eq_(list(sorted(parsed_dict.keys())), ["gene_id", "tag", "version"])
    eq_(parsed_dict["gene_id"], ["ENSG001", "ENSG002"])
    eq_(parsed_dict["tag"], ["bogotron", "wolfpuppy"])
    eq_(parsed_dict["version"], ["1", "2"])


def test_attributes_without_quotes():
    attributes = [
        "gene_id ENSG001; tag bogotron; version 1;",
        "gene_id ENSG002; tag wolfpuppy; version 2"
    ]
    parsed_dict = expand_attribute_strings(attributes)
    eq_(list(sorted(parsed_dict.keys())), ["gene_id", "tag", "version"])
    eq_(parsed_dict["gene_id"], ["ENSG001", "ENSG002"])
    eq_(parsed_dict["tag"], ["bogotron", "wolfpuppy"])
    eq_(parsed_dict["version"], ["1", "2"])


def test_optional_attributes():
    attributes = [
        "gene_id ENSG001; sometimes-present bogotron;",
        "gene_id ENSG002;",
        "gene_id ENSG003; sometimes-present wolfpuppy;",
    ]
    parsed_dict = expand_attribute_strings(attributes)
    eq_(list(sorted(parsed_dict.keys())), ["gene_id", "sometimes-present"])
    eq_(parsed_dict["gene_id"], ["ENSG001", "ENSG002", "ENSG003"])
    eq_(parsed_dict["sometimes-present"], ["bogotron", "", "wolfpuppy"])
