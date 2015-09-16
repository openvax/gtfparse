from gtftools import expand_attribute_strings
from nose.tools import eq_

def test_attributes_in_quotes():
    attributes = [
        "gene_id \"ENSG001\"; tag \"bogotron\";",
        "gene_id \"ENSG002\"; tag \"wolfpuppy\";"
    ]
    parsed_dict = expand_attribute_strings(attributes)
    eq_(list(sorted(parsed_dict.keys())), ["gene_id", "tag"])
    eq_(parsed_dict["gene_id"], ["ENSG001", "ENSG002"])
    eq_(parsed_dict["tag"], ["bogotron", "wolfpuppy"])