from gtfparse import expand_attribute_strings

def test_attributes_in_quotes():
    attributes = [
        "gene_id \"ENSG001\"; tag \"bogotron\"; version \"1\";",
        "gene_id \"ENSG002\"; tag \"wolfpuppy\"; version \"2\";"
    ]
    parsed_dict = expand_attribute_strings(attributes, quote_char='"')
    assert list(sorted(parsed_dict.keys())), ["gene_id", "tag", "version"]
    assert parsed_dict["gene_id"] == ["ENSG001", "ENSG002"]
    assert parsed_dict["tag"] == ["bogotron", "wolfpuppy"]
    assert parsed_dict["version"] == ["1", "2"]


def test_attributes_without_quotes():
    attributes = [
        "gene_id ENSG001; tag bogotron; version 1;",
        "gene_id ENSG002; tag wolfpuppy; version 2"
    ]
    parsed_dict = expand_attribute_strings(attributes)
    assert list(sorted(parsed_dict.keys())) == ["gene_id", "tag", "version"]
    assert parsed_dict["gene_id"] == ["ENSG001", "ENSG002"]
    assert parsed_dict["tag"] == ["bogotron", "wolfpuppy"]
    assert parsed_dict["version"] == ["1", "2"]


def test_optional_attributes():
    attributes = [
        "gene_id ENSG001; sometimes-present bogotron;",
        "gene_id ENSG002;",
        "gene_id ENSG003; sometimes-present wolfpuppy;",
    ]
    parsed_dict = expand_attribute_strings(attributes)
    assert list(sorted(parsed_dict.keys())) ==  ["gene_id", "sometimes-present"]
    assert parsed_dict["gene_id"] ==  ["ENSG001", "ENSG002", "ENSG003"]
    assert parsed_dict["sometimes-present"] ==  ["bogotron", "", "wolfpuppy"]
