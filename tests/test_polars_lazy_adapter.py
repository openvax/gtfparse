"""parse_with_polars_lazy behaviour that does not need Polars installed."""

import sys
import warnings

import pandas as pd
import pytest

from gtfparse import parse_gtf
from gtfparse.read_gtf import parse_with_polars_lazy

from .data import data_path

GTF_PATH = data_path("ensembl_grch37.head.gtf")


class RecordingPolars:
    """Stand-in for the polars module that records the frame it is handed."""

    class Frame:
        def __init__(self, pandas_frame):
            self.pandas_frame = pandas_frame

        def lazy(self):
            return ("lazy", self.pandas_frame)

    @classmethod
    def from_pandas(cls, pandas_frame):
        return cls.Frame(pandas_frame)


def test_missing_polars_raises_install_hint_without_deprecation_warning(monkeypatch):
    monkeypatch.setitem(sys.modules, "polars", None)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with pytest.raises(ImportError, match=r"gtfparse\[polars\]"):
            parse_with_polars_lazy(GTF_PATH)
    assert [w for w in caught if issubclass(w.category, DeprecationWarning)] == []


def test_adapter_warns_and_wraps_the_parse_gtf_frame(monkeypatch):
    monkeypatch.setattr(sys.modules["gtfparse.read_gtf"], "require_polars", lambda: RecordingPolars)
    with pytest.warns(DeprecationWarning, match="parse_with_polars_lazy is deprecated"):
        kind, frame = parse_with_polars_lazy(GTF_PATH)
    assert kind == "lazy"
    pd.testing.assert_frame_equal(frame, parse_gtf(GTF_PATH))
