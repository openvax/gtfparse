from io import StringIO

import pytest

from gtfparse import read_gtf, write_gtf
from gtfparse.read_gtf import parse_with_polars_lazy

polars = pytest.importorskip("polars")

GTF_TEXT = (
    '1\tensembl\tgene\t10\t20\t.\t+\t.\tgene_id "g1";\n'
    '2\thavana\texon\t30\t40\t.\t-\t.\tgene_id "g2";\n'
)

pytestmark = pytest.mark.filterwarnings("error::DeprecationWarning")


@pytest.mark.parametrize("result_type", ["pandas", "polars"])
def test_repeated_reads_preserve_categoricals(result_type):
    for _ in range(2):
        df = read_gtf(StringIO(GTF_TEXT), result_type=result_type, features={"exon"})
        assert list(df["seqname"]) == ["2"]
        assert list(df["feature"]) == ["exon"]
        assert list(df["gene_id"]) == ["g2"]
        if result_type == "polars":
            assert df["seqname"].dtype == polars.Categorical
        else:
            assert df["seqname"].dtype.name == "category"


def test_categoricals_from_separate_lazy_reads_can_join():
    left = parse_with_polars_lazy(StringIO(GTF_TEXT))
    reversed_text = "\n".join(reversed(GTF_TEXT.splitlines())) + "\n"
    right = parse_with_polars_lazy(StringIO(reversed_text))
    keys = ["seqname", "source", "feature", "strand"]
    joined = left.join(right, on=keys).sort("start").collect()

    assert joined.height == 2
    assert joined["start"].to_list() == [10, 30]
    assert joined["start_right"].to_list() == [10, 30]
    for key in keys:
        assert joined[key].dtype == polars.Categorical


@pytest.mark.skipif(not hasattr(polars, "Categories"), reason="Requires modern Polars")
def test_modern_reads_do_not_toggle_global_string_cache(monkeypatch):
    def fail_on_cache_toggle():
        pytest.fail("Reading a GTF must not toggle the global string cache on modern Polars")

    monkeypatch.setattr(polars, "enable_string_cache", fail_on_cache_toggle)
    monkeypatch.setattr(polars, "disable_string_cache", fail_on_cache_toggle)
    assert read_gtf(StringIO(GTF_TEXT), result_type="polars").height == 2


def test_optional_polars_writer_round_trip(tmp_path):
    original = read_gtf(StringIO(GTF_TEXT), result_type="polars")
    path = tmp_path / "out.gtf.gz"
    write_gtf(original, path)
    recovered = read_gtf(path, result_type="polars")
    from polars.testing import assert_frame_equal

    assert_frame_equal(original, recovered, categorical_as_str=True)
