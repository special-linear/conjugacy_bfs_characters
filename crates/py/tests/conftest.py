import json
import os
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]
GOLDEN_PATH = REPO_ROOT / "fixtures" / "golden" / "n06_g2.json"

# Fields the golden-file convention treats as volatile, plus the
# run-configuration hash (the golden file carries a placeholder there).
_VOLATILE = ("run", "timings_s", "tool", "config_hash_blake3")


def strip_volatile(doc: dict) -> dict:
    doc = {k: v for k, v in doc.items() if k not in _VOLATILE}
    if "engine" in doc:
        doc["engine"] = {k: v for k, v in doc["engine"].items() if k != "threads"}
    return doc


@pytest.fixture(scope="session")
def golden() -> dict:
    return json.loads(GOLDEN_PATH.read_text(encoding="utf-8"))


@pytest.fixture(scope="session")
def cli_bin() -> str:
    """Path to the classdiam CLI for interop tests (set CLASSDIAM_CLI_BIN)."""
    path = os.environ.get("CLASSDIAM_CLI_BIN")
    if not path or not Path(path).exists():
        pytest.skip("CLASSDIAM_CLI_BIN not set; skipping CLI interop test")
    return path
