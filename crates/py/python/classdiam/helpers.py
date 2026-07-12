"""Pure-Python layer over ``classdiam._core``.

The native module returns schema-v1 result documents (the exact JSON the
CLI writes) as strings; this layer parses them into dicts and wraps them
in small conveniences. ``RunResult.raw`` is always the canonical document.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Callable, Iterable, List, Optional, Sequence, Union

from . import _core

ProgressCallback = Callable[[dict], None]
UnionSpec = Union[str, Sequence[Sequence[int]], Sequence[int]]

_RESULT_FORMAT = "classdiam/result"
_SUPPORTED_FORMAT_VERSION = 1


class Suspended(Exception):
    """A run stopped early with resumable state.

    Raised when the wall-clock deadline expired or the run was interrupted
    (Ctrl-C / notebook interrupt). Completed jobs are in :attr:`completed`;
    continue later with ``classdiam.resume(run_dir)``.

    Attributes:
        run_dir: directory holding checkpoints/manifest (``None`` for
            in-memory runs — their suspension state is dropped).
        checkpoints: written checkpoint file paths.
        completed: :class:`RunSet` of jobs that did finish.
        reason: ``"deadline"`` or ``"interrupt"``.
    """

    def __init__(
        self,
        message: str,
        *,
        run_dir: Optional[str],
        checkpoints: List[str],
        completed: "RunSet",
        reason: str,
    ) -> None:
        super().__init__(message)
        self.run_dir = run_dir
        self.checkpoints = checkpoints
        self.completed = completed
        self.reason = reason


class RunResult:
    """One ``(n, union)`` result. ``raw`` is the full schema-v1 document."""

    def __init__(self, raw: dict) -> None:
        self.raw = raw

    # -- identity ---------------------------------------------------------
    @property
    def n(self) -> int:
        return self.raw["n"]

    @property
    def label(self) -> str:
        return self.raw["generators"]["label"]

    # -- headline numbers -------------------------------------------------
    @property
    def diameter(self) -> int:
        """Diameter of the identity component (max finite distance)."""
        return self.raw["results"]["diameter_identity_component"]

    @property
    def stop_radius(self) -> int:
        return self.raw["results"]["stopping"]["stop_radius"]

    @property
    def reachable_count(self) -> int:
        return self.raw["results"]["reachable_count"]

    # -- per-class arrays (canonical partition order) ----------------------
    @property
    def partitions(self) -> List[List[int]]:
        """Canonical partition order, reduced form (parts >= 2 only;
        identity class is ``[]``). All indexed arrays follow this order."""
        return self.raw["partition_order"]["partitions_reduced"]

    @property
    def distances(self) -> List[int]:
        """Exact distance per canonical index; ``-1`` = unreachable."""
        return self.raw["results"]["distance"]

    @property
    def first_hit_even(self) -> List[int]:
        return self.raw["results"]["first_hit_even"]

    @property
    def first_hit_odd(self) -> List[int]:
        return self.raw["results"]["first_hit_odd"]

    @property
    def layers(self) -> List[dict]:
        return self.raw["results"]["layers"]

    def support(self, r: int) -> List[List[int]]:
        """Cycle types with a length-``r`` factorization (exact-length
        support of radius ``r``; NOT the same as distance == r)."""
        layers = self.layers
        if not 0 <= r < len(layers):
            raise ValueError(f"r={r} outside recorded layers 0..{len(layers) - 1}")
        partitions = self.partitions
        return [partitions[i] for i in layers[r]["support"]]

    def distance_of(self, cycle_type: Iterable[int]) -> int:
        """Distance of a conjugacy class given by its cycle type (fixed
        points optional): ``distance_of([2])``, ``distance_of([3, 2])``."""
        key = sorted((int(p) for p in cycle_type if int(p) >= 2), reverse=True)
        for i, parts in enumerate(self.partitions):
            if parts == key:
                return self.distances[i]
        raise ValueError(f"cycle type {list(cycle_type)!r} is not a class of S_{self.n}")

    def to_dataframe(self) -> Any:
        """One row per conjugacy class, in canonical order (needs pandas)."""
        try:
            import pandas as pd
        except ImportError as e:  # pragma: no cover
            raise ImportError(
                "RunResult.to_dataframe() needs pandas — pip install pandas"
            ) from e
        results = self.raw["results"]
        class_data = self.raw["class_data"]
        return pd.DataFrame(
            {
                "partition": [tuple(p) for p in self.partitions],
                "distance": results["distance"],
                "first_hit_even": results["first_hit_even"],
                "first_hit_odd": results["first_hit_odd"],
                "sign": class_data["sign"],
                "class_size": [int(s) for s in class_data["class_size"]],
            }
        )

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"<RunResult n={self.n} {self.label} diameter={self.diameter} "
            f"reachable={self.reachable_count}/{len(self.partitions)}>"
        )


class RunSet(list):
    """A list of :class:`RunResult` plus run-level metadata."""

    def __init__(
        self,
        results: Iterable[RunResult] = (),
        *,
        run_dir: Optional[str] = None,
        run_id: Optional[str] = None,
        manifest: Optional[dict] = None,
    ) -> None:
        super().__init__(results)
        self.run_dir = run_dir
        self.run_id = run_id
        self.manifest = manifest

    def __repr__(self) -> str:  # pragma: no cover
        inner = ", ".join(repr(r) for r in self)
        return f"RunSet([{inner}])"


def _normalize_ns(n: Union[int, str, Sequence[int]]) -> List[int]:
    if isinstance(n, str):
        ns = _core.parse_n_spec(n)
    elif isinstance(n, int):
        ns = [n]
    else:
        ns = [int(x) for x in n]
    if not ns:
        raise ValueError("no n values given")
    return ns


def _normalize_union(union: UnionSpec) -> List[List[int]]:
    """One union -> list of classes (each a list of parts >= 2)."""
    if isinstance(union, str):
        return _core.parse_union_spec(union)
    union = list(union)
    if union and all(isinstance(p, int) for p in union):
        # a bare cycle type like [3, 2] means the single class (3,2)
        return [[int(p) for p in union]]
    return [[int(p) for p in cls] for cls in union]


def _normalize_unions(unions: Union[UnionSpec, Sequence[UnionSpec]]) -> List[List[List[int]]]:
    if isinstance(unions, str):
        return [_normalize_union(unions)]
    normalized = [_normalize_union(u) for u in unions]
    if not normalized:
        raise ValueError("no generating unions given")
    return normalized


def _finish(outcome: dict) -> RunSet:
    results = RunSet(
        (RunResult(job["document"]) for job in outcome["jobs"] if job["status"] == "done"),
        run_dir=outcome.get("out_dir"),
        run_id=outcome.get("run_id"),
    )
    if outcome["any_suspended"] or outcome["interrupted"]:
        reason = "interrupt" if outcome["interrupted"] else "deadline"
        suspended = [j for j in outcome["jobs"] if j["status"] == "suspended"]
        checkpoints = [j["checkpoint"] for j in suspended if j.get("checkpoint")]
        run_dir = outcome.get("out_dir")
        hint = f"; resume with classdiam.resume({run_dir!r})" if run_dir else ""
        raise Suspended(
            f"{len(suspended)} job(s) suspended ({reason}){hint}",
            run_dir=run_dir,
            checkpoints=checkpoints,
            completed=results,
            reason=reason,
        )
    return results


def run(
    n: Union[int, str, Sequence[int]],
    unions: Union[UnionSpec, Sequence[UnionSpec]],
    *,
    engine: str = "modular",
    primes: int = 3,
    deadline_s: Optional[float] = None,
    allow_identity: bool = False,
    out_dir: Optional[Union[str, os.PathLike]] = None,
    progress: Optional[ProgressCallback] = None,
    progress_every_ms: int = 250,
    threads: Optional[int] = None,
) -> RunSet:
    """Compute distances/diameters for each ``n`` and each generating union.

    Args:
        n: ``12``, ``"6..=12"``, ``"6,8,10"``, or a sequence of ints.
        unions: one union or a list of unions. Each union is a CLI-style
            string (classes joined ``+``, parts joined ``,``: ``"2"``,
            ``"3+2,2"``) or a list of classes (``[[3], [2, 2]]``); a flat
            list of ints is a single class (``[3, 2]``).
        engine: ``"modular"`` (production, certified) or ``"exact"``
            (big-integer oracle; small n).
        primes: resident screening primes for the modular engine.
        deadline_s: wall-clock budget in seconds; requires ``out_dir``.
            On expiry, unfinished jobs suspend into checkpoints and a
            :class:`Suspended` is raised.
        allow_identity: permit the identity class as a generator.
        out_dir: run directory (results + manifest + checkpoints). With
            ``None`` everything stays in memory and nothing is written.
        progress: callable receiving one dict per committed radius
            (keys: n, job_name, radius, new, support, reachable,
            job_index, job_total).
        progress_every_ms: progress/interrupt polling throttle.
        threads: worker threads (default: all cores).

    Returns:
        A :class:`RunSet` of :class:`RunResult`, one per completed job,
        in batch order.

    Raises:
        Suspended: deadline expired or run interrupted (Ctrl-C).
        ValueError: malformed specs or options.
        ClassdiamError: internal engine failure.
    """
    if deadline_s is not None and out_dir is None:
        raise ValueError(
            "deadline_s requires out_dir (suspension writes checkpoints there)"
        )
    cfg = {
        "ns": _normalize_ns(n),
        "unions": _normalize_unions(unions),
        "engine": engine,
        "primes": primes,
        "deadline_s": deadline_s,
        "allow_identity": allow_identity,
        "out_dir": os.fspath(out_dir) if out_dir is not None else None,
    }
    outcome = json.loads(
        _core.run_batch(json.dumps(cfg), progress, threads, progress_every_ms)
    )
    return _finish(outcome)


def resume(
    run_dir: Union[str, os.PathLike],
    *,
    deadline_s: Optional[float] = None,
    progress: Optional[ProgressCallback] = None,
    progress_every_ms: int = 250,
    threads: Optional[int] = None,
) -> RunSet:
    """Resume every suspended job in a run directory (CLI- or Python-made).

    Returns only the jobs resumed in this call; use :func:`load_results`
    for the full directory. Raises :class:`Suspended` if the new deadline
    expires again.
    """
    outcome = json.loads(
        _core.resume_batch(
            os.fspath(run_dir), deadline_s, progress, threads, progress_every_ms
        )
    )
    return _finish(outcome)


def load_results(path: Union[str, os.PathLike]) -> RunSet:
    """Load all result documents from a run directory."""
    p = Path(path)
    manifest = None
    manifest_path = p / "manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    results = []
    for file in sorted(p.glob("*.json")):
        if file.name == "manifest.json":
            continue
        doc = json.loads(file.read_text(encoding="utf-8"))
        if doc.get("format") != _RESULT_FORMAT:
            continue
        if doc.get("format_version") != _SUPPORTED_FORMAT_VERSION:
            raise ValueError(
                f"{file}: unsupported format_version {doc.get('format_version')!r}"
            )
        results.append(RunResult(doc))
    return RunSet(
        results,
        run_dir=str(p),
        run_id=(manifest or {}).get("run_id"),
        manifest=manifest,
    )


class Session:
    """Per-``n`` reusable state: build the character tables once, run many
    unions. Cheaper than repeated :func:`run` calls at the same ``n``."""

    def __init__(
        self,
        n: int,
        *,
        engine: str = "modular",
        primes: int = 3,
        threads: Optional[int] = None,
    ) -> None:
        self._native = _core.Session(n, engine, primes, threads)

    @property
    def n(self) -> int:
        return self._native.n

    @property
    def class_count(self) -> int:
        """p(n) — number of conjugacy classes."""
        return self._native.class_count

    @property
    def partition_order(self) -> List[List[int]]:
        """Canonical partition order (reduced form), shared by all results."""
        return json.loads(self._native.partition_order_json())

    def run_union(
        self,
        classes: UnionSpec,
        *,
        label: Optional[str] = None,
        allow_identity: bool = False,
        deadline_s: Optional[float] = None,
        checkpoint_dir: Optional[Union[str, os.PathLike]] = None,
        result_dir: Optional[Union[str, os.PathLike]] = None,
        progress: Optional[ProgressCallback] = None,
        progress_every_ms: int = 250,
    ) -> RunResult:
        """Run one generating union; returns its :class:`RunResult`.

        ``deadline_s`` requires ``checkpoint_dir`` (that is where the
        resumable checkpoint goes; resume it with :func:`resume` on the
        parent directory, or the CLI).
        """
        if deadline_s is not None and checkpoint_dir is None:
            raise ValueError(
                "deadline_s requires checkpoint_dir for resumable suspension"
            )
        opts = {
            "label": label,
            "allow_identity": allow_identity,
            "deadline_s": deadline_s,
            "checkpoint_dir": os.fspath(checkpoint_dir) if checkpoint_dir else None,
            "result_dir": os.fspath(result_dir) if result_dir else None,
        }
        outcome = json.loads(
            self._native.run_union(
                _normalize_union(classes), json.dumps(opts), progress, progress_every_ms
            )
        )
        return _finish(outcome)[0]

    def __repr__(self) -> str:  # pragma: no cover
        return f"<classdiam.Session n={self.n} classes={self.class_count}>"
