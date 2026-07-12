"""classdiam — conjugacy-class BFS distances and diameters of symmetric-group
Cayley graphs, computed exactly via characters.

Quick start (Kaggle-style)::

    import classdiam
    res = classdiam.run(n="6..=20", unions=["2"])   # transpositions
    res[0].diameter
    res[0].to_dataframe()                            # needs pandas

Long runs: pass ``deadline_s=`` and ``out_dir=``; on expiry (or Ctrl-C) a
:class:`classdiam.Suspended` exception carries the resumable state, and
``classdiam.resume(out_dir)`` continues in a later session. Run directories
are fully interchangeable with the ``classdiam`` CLI.
"""

from ._core import ClassdiamError, __version__
from .helpers import (
    RunResult,
    RunSet,
    Session,
    Suspended,
    load_results,
    resume,
    run,
)

__all__ = [
    "ClassdiamError",
    "RunResult",
    "RunSet",
    "Session",
    "Suspended",
    "__version__",
    "load_results",
    "resume",
    "run",
]
