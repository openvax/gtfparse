import importlib

import gtfparse
from gtfparse import columns


def test_required_columns_are_fixed_columns_then_attribute():
    assert [*columns.GTF_FIXED_COLUMNS, "attribute"] == columns.REQUIRED_COLUMNS
    assert columns.RAW_ATTRIBUTE_COLUMN == "attribute"
    assert set(columns.CATEGORICAL_COLUMNS) <= set(columns.GTF_FIXED_COLUMNS)


def test_column_constants_remain_importable_from_their_old_modules():
    # gtfparse.read_gtf and gtfparse.write_gtf are functions on the package,
    # so look the modules up by name.
    read_module = importlib.import_module("gtfparse.read_gtf")
    write_module = importlib.import_module("gtfparse.write_gtf")
    reader = importlib.import_module("gtfparse._read_csv")
    assert gtfparse.REQUIRED_COLUMNS is columns.REQUIRED_COLUMNS
    assert read_module.REQUIRED_COLUMNS is columns.REQUIRED_COLUMNS
    assert reader.REQUIRED_COLUMNS is columns.REQUIRED_COLUMNS
    assert reader.CATEGORICAL_COLUMNS is columns.CATEGORICAL_COLUMNS
    assert write_module.GTF_FIXED_COLUMNS is columns.GTF_FIXED_COLUMNS
    assert write_module.MISSING_VALUE == columns.MISSING_VALUE == "."
    assert write_module.RAW_ATTRIBUTE_COLUMN == columns.RAW_ATTRIBUTE_COLUMN
