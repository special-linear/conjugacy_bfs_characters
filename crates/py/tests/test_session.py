import pytest

import classdiam


def test_session_matches_run_and_reuses_tables():
    session = classdiam.Session(6)
    a = session.run_union("2")
    reference = classdiam.run(n=6, unions=["2"])[0]
    assert a.distances == reference.distances
    assert a.diameter == reference.diameter

    b = session.run_union([[3], [2, 2]])  # second union on the same tables
    assert b.label == "g3+2.2"
    assert b.reachable_count == 6

    assert session.n == 6
    assert session.class_count == 11
    assert repr(session)


def test_session_partition_order_matches_documents():
    session = classdiam.Session(6)
    result = session.run_union("2")
    assert session.partition_order == result.partitions
    assert session.partition_order[0] == [6]
    assert session.partition_order[-1] == []


def test_session_exact_engine_and_threads():
    session = classdiam.Session(7, engine="exact", threads=1)
    assert session.run_union("2").diameter == 6


def test_session_label_override():
    session = classdiam.Session(6)
    assert session.run_union("2", label="transpositions").label == "transpositions"


def test_session_deadline_requires_checkpoint_dir():
    session = classdiam.Session(9)
    with pytest.raises(ValueError, match="checkpoint_dir"):
        session.run_union("2", deadline_s=1)


def test_session_suspension_is_resumable(tmp_path):
    session = classdiam.Session(9)
    with pytest.raises(classdiam.Suspended) as excinfo:
        session.run_union(
            "2",
            deadline_s=0,
            checkpoint_dir=tmp_path / "checkpoints",
            result_dir=tmp_path,
        )
    assert excinfo.value.reason == "deadline"
    assert (tmp_path / "checkpoints" / "n09_g2.ckpt").exists()

    # the parent dir is a valid run dir for resume (manifest is optional)
    resumed = classdiam.resume(tmp_path)
    assert resumed[0].diameter == 8
    assert (tmp_path / "n09_g2.json").exists()


def test_session_result_dir_writes_document(tmp_path):
    session = classdiam.Session(6)
    result = session.run_union("2", result_dir=tmp_path)
    assert (tmp_path / "n06_g2.json").exists()
    on_disk = classdiam.load_results(tmp_path)
    assert on_disk[0].raw == result.raw
