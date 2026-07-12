import json

import pytest

import classdiam
from conftest import strip_volatile


def test_exact_run_matches_golden(golden):
    res = classdiam.run(n=6, unions=["2"], engine="exact")
    assert len(res) == 1
    assert strip_volatile(res[0].raw) == strip_volatile(golden)


def test_modular_results_match_golden(golden):
    res = classdiam.run(n=6, unions=["2"])  # default engine = modular
    assert res[0].raw["results"] == golden["results"]
    assert res[0].raw["arithmetic"]["mode"] == "modular+certified"
    assert res[0].raw["certification"] is not None


def test_out_dir_files_and_load_results(tmp_path):
    out = tmp_path / "results"
    res = classdiam.run(n=6, unions=["2"], out_dir=out, engine="exact")
    assert (out / "n06_g2.json").exists()
    assert (out / "manifest.json").exists()
    on_disk = json.loads((out / "n06_g2.json").read_text(encoding="utf-8"))
    assert on_disk == res[0].raw

    loaded = classdiam.load_results(out)
    assert len(loaded) == 1
    assert loaded[0].raw == res[0].raw
    assert loaded.run_id == res.run_id
    assert loaded.manifest["status"] == "completed"


def test_in_memory_run_writes_nothing(tmp_path):
    before = sorted(tmp_path.iterdir())
    classdiam.run(n=6, unions=["2"])
    assert sorted(tmp_path.iterdir()) == before


def test_union_spellings_agree():
    a = classdiam.run(n=7, unions=["3+2,2"], engine="exact")
    b = classdiam.run(n=7, unions=[[[3], [2, 2]]], engine="exact")
    assert a[0].distances == b[0].distances
    # a flat int list is a single class
    c = classdiam.run(n=7, unions=[[3, 2]], engine="exact")
    d = classdiam.run(n=7, unions=["3,2"], engine="exact")
    assert c[0].distances == d[0].distances


def test_n_spellings_agree():
    a = classdiam.run(n="6,7", unions=["2"], engine="exact")
    b = classdiam.run(n=[6, 7], unions=["2"], engine="exact")
    assert [r.n for r in a] == [r.n for r in b] == [6, 7]
    assert [r.diameter for r in a] == [r.diameter for r in b] == [5, 6]


def test_oversized_union_is_skipped():
    res = classdiam.run(n=4, unions=["2", "7"], engine="exact")
    assert len(res) == 1  # the 7-cycle job is skipped, not fatal
    assert res[0].label == "g2"


def test_deadline_without_out_dir_is_rejected():
    with pytest.raises(ValueError, match="out_dir"):
        classdiam.run(n=6, unions=["2"], deadline_s=1)


def test_bad_specs_raise_value_error():
    with pytest.raises(ValueError):
        classdiam.run(n="six", unions=["2"])
    with pytest.raises(ValueError):
        classdiam.run(n=6, unions=["0"])
    with pytest.raises(ValueError):
        classdiam.run(n=6, unions=["2"], engine="bogus")
    with pytest.raises(ValueError):
        classdiam.run(n=6, unions=[])


def test_progress_callback_events():
    events = []
    classdiam.run(n=9, unions=["2"], progress=events.append, progress_every_ms=0)
    assert events, "modular runs must emit progress"
    assert [e["radius"] for e in events] == list(range(1, len(events) + 1))
    keys = {"n", "job_name", "radius", "new", "support", "reachable", "job_index", "job_total"}
    assert keys <= set(events[0])
    assert events[0]["job_name"] == "n09_g2"


def test_progress_callback_error_propagates():
    def boom(event):
        raise RuntimeError("stop right there")

    with pytest.raises(RuntimeError, match="stop right there"):
        classdiam.run(n=9, unions=["2"], progress=boom, progress_every_ms=0)


def test_threads_kwarg():
    res = classdiam.run(n=8, unions=["2"], threads=1)
    assert res[0].diameter == 7
