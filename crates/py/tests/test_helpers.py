import json

import pytest

import classdiam


@pytest.fixture(scope="module")
def n6() -> classdiam.RunResult:
    return classdiam.run(n=6, unions=["2"], engine="exact")[0]


def test_support_matches_layers(n6, golden):
    for layer in golden["results"]["layers"]:
        expected = [golden["partition_order"]["partitions_reduced"][i] for i in layer["support"]]
        assert n6.support(layer["r"]) == expected
    with pytest.raises(ValueError):
        n6.support(n6.stop_radius + 1)
    with pytest.raises(ValueError):
        n6.support(-1)


def test_distance_of(n6):
    assert n6.distance_of([2]) == 1
    assert n6.distance_of([6]) == 5
    assert n6.distance_of([]) == 0  # identity
    assert n6.distance_of([2, 1, 1, 1, 1]) == 1  # fixed points tolerated
    with pytest.raises(ValueError):
        n6.distance_of([7])


def test_accessor_properties(n6, golden):
    assert n6.n == 6
    assert n6.label == "g2"
    assert n6.stop_radius == 5
    assert n6.partitions == golden["partition_order"]["partitions_reduced"]
    assert n6.first_hit_even == golden["results"]["first_hit_even"]
    assert n6.first_hit_odd == golden["results"]["first_hit_odd"]
    assert len(n6.layers) == n6.stop_radius + 1
    assert "RunResult" in repr(n6)


def test_to_dataframe(n6):
    pd = pytest.importorskip("pandas")
    df = n6.to_dataframe()
    assert list(df.columns) == [
        "partition",
        "distance",
        "first_hit_even",
        "first_hit_odd",
        "sign",
        "class_size",
    ]
    assert len(df) == 11
    assert df["distance"].max() == n6.diameter
    assert df["class_size"].sum() == 720  # 6!
    assert df.loc[df["partition"] == (2,), "distance"].item() == 1


def test_runset_metadata(tmp_path):
    out = tmp_path / "results"
    res = classdiam.run(n=6, unions=["2"], out_dir=out, engine="exact")
    assert res.run_dir == str(out)
    assert res.run_id
    assert isinstance(res, list) and len(res) == 1


def test_load_results_rejects_future_format(tmp_path):
    doc = {"format": "classdiam/result", "format_version": 2}
    (tmp_path / "n99_gX.json").write_text(json.dumps(doc), encoding="utf-8")
    with pytest.raises(ValueError, match="format_version"):
        classdiam.load_results(tmp_path)


def test_load_results_ignores_foreign_json(tmp_path):
    (tmp_path / "notes.json").write_text(json.dumps({"hello": 1}), encoding="utf-8")
    assert list(classdiam.load_results(tmp_path)) == []


def test_version_attribute():
    assert isinstance(classdiam.__version__, str) and classdiam.__version__
