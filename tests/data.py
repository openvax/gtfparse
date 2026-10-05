import os
from importlib.util import find_spec

import pytest

optional_polars = pytest.mark.skipif(find_spec("polars") is None, reason="Requires polars extra")


def data_path(name):
    """
    Return the absolute path to a fixture relative to tests/data.
    """
    return os.path.join(os.path.dirname(__file__), "data", name)
