import subprocess

import pytest

import classdiam


def _suspend(tmp_path):
    """Run n=9 transpositions with an already-expired deadline."""
    out = tmp_path / "results"
    with pytest.raises(classdiam.Suspended) as excinfo:
        classdiam.run(n=9, unions=["2"], deadline_s=0, out_dir=out)
    return out, excinfo.value


def test_deadline_suspends_then_resume_completes(tmp_path):
    out, suspended = _suspend(tmp_path)
    assert suspended.reason == "deadline"
    assert suspended.run_dir == str(out)
    assert suspended.checkpoints, "suspension must write checkpoints"
    assert (out / "checkpoints" / "n09_g2.ckpt").exists()
    assert list(suspended.completed) == []

    resumed = classdiam.resume(out)
    reference = classdiam.run(n=9, unions=["2"])
    assert resumed[0].distances == reference[0].distances
    assert resumed[0].diameter == reference[0].diameter
    assert resumed[0].raw["run"]["resumed_from_checkpoint"] is True
    assert not (out / "checkpoints" / "n09_g2.ckpt").exists()

    loaded = classdiam.load_results(out)
    assert loaded.manifest["status"] == "completed"
    assert len(loaded) == 1


def test_resume_deadline_can_expire_again(tmp_path):
    out, _ = _suspend(tmp_path)
    with pytest.raises(classdiam.Suspended) as excinfo:
        classdiam.resume(out, deadline_s=0)
    assert excinfo.value.reason == "deadline"
    # still resumable after that
    assert classdiam.resume(out)[0].diameter == 8


def test_resume_rejects_non_run_dir(tmp_path):
    with pytest.raises(ValueError):
        classdiam.resume(tmp_path / "nowhere")


def test_cli_resumes_python_suspended_dir(tmp_path, cli_bin):
    out, _ = _suspend(tmp_path)
    proc = subprocess.run(
        [cli_bin, "resume", str(out)], capture_output=True, text=True, timeout=600
    )
    assert proc.returncode == 0, proc.stderr
    loaded = classdiam.load_results(out)
    assert loaded[0].diameter == 8
    assert loaded.manifest["status"] == "completed"


def test_python_resumes_cli_suspended_dir(tmp_path, cli_bin):
    out = tmp_path / "cli_run"
    proc = subprocess.run(
        [cli_bin, "run", "-n", "9", "-u", "2", "--deadline", "0", "-o", str(out)],
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert proc.returncode == 75, proc.stderr  # suspended
    resumed = classdiam.resume(out)
    assert resumed[0].diameter == 8
    assert resumed[0].raw["run"]["resumed_from_checkpoint"] is True
