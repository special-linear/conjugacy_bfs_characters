#!/usr/bin/env python
"""Generate committed character-table fixtures for the Rust test suite.

SymPy (1.14) has no ready-made S_n character table, so this script implements
an INDEPENDENT Murnaghan-Nakayama evaluator (hook-cell/rim-path enumeration on
the Young diagram -- deliberately a different algorithm than the Rust abacus /
beta-set method) and refuses to write any fixture that does not pass the
layered certificate (design doc 03 section 4.1):

  1. identity column equals hook-length degrees (independent closed form);
  2. transposition column: chi * C(n,2) == f * (sum of contents);
  3. full exact row orthogonality X diag(|C|) X^T = n! I;
  4. transpose-sign relation chi^{rho'} = sgn(nu) chi^rho, and
     self-transpose rows vanish on odd classes;
  5. trivial row == 1, sign row == sgn.

SymPy is used for an independent partition-enumeration cross-check (count and
content of the canonical order).

Outputs (all JSON, committed to the repo):
  fixtures/characters/char_n{NN}.json       full tables, n <= --full-max
  fixtures/characters/spot_n{NN}.json       spot values, --spot-min..--spot-max
  fixtures/degrees/deg_n{NN}.json           degrees/class data, n <= --deg-max

Usage:  python tools/gen_fixtures.py [--full-max 14] [--spot-min 17]
        [--spot-max 22] [--spot-count 200] [--deg-max 30] [--out fixtures]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import random
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

SCRIPT_VERSION = "1"
ORDER_CONVENTION = "lex_desc_full_parts_v1"


# ----------------------------------------------------------------------------
# partitions, canonical order, class data
# ----------------------------------------------------------------------------

def partitions_desc(n: int) -> list[tuple[int, ...]]:
    """All partitions of n as descending part tuples, lex-descending order
    ([n] first, [1,...,1] last) -- the project's canonical order."""
    out: list[tuple[int, ...]] = []

    def rec(remaining: int, max_part: int, prefix: list[int]) -> None:
        if remaining == 0:
            out.append(tuple(prefix))
            return
        for k in range(min(max_part, remaining), 0, -1):
            prefix.append(k)
            rec(remaining - k, k, prefix)
            prefix.pop()

    rec(n, n if n else 1, [])
    return out


def transpose(p: tuple[int, ...]) -> tuple[int, ...]:
    if not p:
        return ()
    return tuple(sum(1 for x in p if x >= j) for j in range(1, p[0] + 1))


def z_value(p: tuple[int, ...]) -> int:
    z = 1
    for part in set(p):
        m = p.count(part)
        z *= part**m * math.factorial(m)
    return z


def sign(p: tuple[int, ...]) -> int:
    return -1 if (sum(p) - len(p)) % 2 else 1


def hook_degree(p: tuple[int, ...]) -> int:
    n = sum(p)
    t = transpose(p)
    hooks = 1
    for i, row in enumerate(p):
        for j in range(row):
            hooks *= row - j + t[j] - i - 1
    d = Fraction(math.factorial(n), hooks)
    assert d.denominator == 1
    return int(d)


def content_sum(p: tuple[int, ...]) -> int:
    return sum(j - i for i, row in enumerate(p) for j in range(row))


# ----------------------------------------------------------------------------
# independent Murnaghan-Nakayama via hook cells (rim-path removal)
# ----------------------------------------------------------------------------

def strip_removals(p: tuple[int, ...], length: int):
    """Yield (smaller_partition, height_sign) for every removable border
    strip of the given length: strips correspond to cells (i, j) with hook
    length == `length`; removal follows the rim path, and the sign is
    (-1)^(rows spanned - 1)."""
    t = transpose(p)
    for i, row in enumerate(p):  # 0-based
        for j in range(row):
            hook = row - j + t[j] - i - 1
            if hook != length:
                continue
            r = t[j]  # 1-based index of the last row of column j+1
            new_parts = list(p[: i])
            for k in range(i, r - 1):
                new_parts.append(p[k + 1] - 1)
            new_parts.append(j)  # row r keeps its first j cells
            new_parts.extend(p[r:])  # rows below the strip are untouched
            new_parts = [x for x in new_parts if x > 0]
            assert sum(new_parts) == sum(p) - length, \
                f"strip removal size bug at {p} ({i},{j})"
            assert all(
                new_parts[k] >= new_parts[k + 1] for k in range(len(new_parts) - 1)
            ), f"bad strip removal from {p} at ({i},{j})"
            yield tuple(new_parts), (-1) ** (r - 1 - i)


@lru_cache(maxsize=None)
def chi(rho: tuple[int, ...], nu: tuple[int, ...]) -> int:
    """chi^rho(nu) by Murnaghan-Nakayama, removing the largest cycle first."""
    assert sum(rho) == sum(nu)
    if not nu:
        return 1
    length, rest = nu[0], nu[1:]
    return sum(s * chi(smaller, rest) for smaller, s in strip_removals(rho, length))


# ----------------------------------------------------------------------------
# certificates
# ----------------------------------------------------------------------------

def certify_full_table(n: int, parts: list[tuple[int, ...]],
                       table: list[list[int]], degrees: list[int],
                       full_orthogonality: bool) -> list[str]:
    """table[nu_idx][rho_idx]. Raises on any failure; returns passed list."""
    q = len(parts)
    passed = []

    identity_idx = q - 1
    assert parts[identity_idx] == tuple([1] * n)
    assert table[identity_idx] == degrees, "identity column != hook degrees"
    passed.append("hook_degrees")

    if n >= 2:
        transposition = tuple([2] + [1] * (n - 2))
        ti = parts.index(transposition)
        choose2 = n * (n - 1) // 2
        for r, rho in enumerate(parts):
            assert table[ti][r] * choose2 == degrees[r] * content_sum(rho), \
                f"content formula fails at rho={rho}"
        passed.append("content_column")

    class_sizes = [math.factorial(n) // z_value(p) for p in parts]
    pairs = ((a, b) for a in range(q) for b in range(a, q))
    if not full_orthogonality:
        rng = random.Random(20260712)
        pairs = {(a, a) for a in range(q)}
        while len(pairs) < 3 * q:
            pairs.add(tuple(sorted((rng.randrange(q), rng.randrange(q)))))
        pairs = sorted(pairs)
    for a, b in pairs:
        s = sum(cs * table[v][a] * table[v][b] for v, cs in enumerate(class_sizes))
        expected = math.factorial(n) if a == b else 0
        assert s == expected, f"orthogonality fails at rows {a},{b}"
    passed.append("orthogonality" if full_orthogonality else "orthogonality_sampled")

    tmap = [parts.index(transpose(p)) for p in parts]
    for v, nu in enumerate(parts):
        sg = sign(nu)
        for r in range(q):
            assert table[v][tmap[r]] == sg * table[v][r], "transpose-sign fails"
    for r in range(q):
        if tmap[r] == r:
            for v, nu in enumerate(parts):
                if sign(nu) == -1:
                    assert table[v][r] == 0, "self-transpose row not vanishing"
    passed.append("transpose_sign")

    for v, nu in enumerate(parts):
        assert table[v][0] == 1, "trivial row"
        assert table[v][q - 1] == sign(nu), "sign row"
    passed.append("trivial_and_sign_rows")

    return passed


def sympy_cross_check(n: int, parts: list[tuple[int, ...]]) -> None:
    """Independent enumeration cross-check via SymPy's partition iterator."""
    from sympy.utilities.iterables import partitions as sympy_partitions

    sympy_set = set()
    for d in sympy_partitions(n) if n else [{}]:
        flat = []
        for part, mult in d.items():
            flat.extend([part] * mult)
        sympy_set.add(tuple(sorted(flat, reverse=True)))
    ours = set(parts)
    assert ours == sympy_set, f"partition enumeration disagrees with SymPy at n={n}"


# ----------------------------------------------------------------------------
# payload hashing (documented; reimplemented by the Rust consumer)
# ----------------------------------------------------------------------------

def compact(obj) -> str:
    return json.dumps(obj, separators=(",", ":"))


def payload_sha256(*segments: str) -> str:
    h = hashlib.sha256()
    h.update("|".join(segments).encode())
    return h.hexdigest()


def meta(seed: int | None = None) -> dict:
    import sympy

    m = {
        "script": "gen_fixtures.py",
        "script_version": SCRIPT_VERSION,
        "sympy_version": sympy.__version__,
        "python_version": platform.python_version(),
    }
    if seed is not None:
        m["seed"] = seed
    return m


# ----------------------------------------------------------------------------
# emitters
# ----------------------------------------------------------------------------

def write_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {path}")


def emit_full_table(n: int, out: Path) -> None:
    parts = partitions_desc(n)
    sympy_cross_check(n, parts)
    degrees = [hook_degree(p) for p in parts]
    table = [[chi(rho, nu) for rho in parts] for nu in parts]  # [nu][rho]
    certificates = certify_full_table(n, parts, table, degrees,
                                      full_orthogonality=(n <= 14))
    parts_json = [list(p) for p in parts]
    write_json(out / "characters" / f"char_n{n:02}.json", {
        "schema_version": "fixture.v1",
        "kind": "full_character_table",
        "generator": meta(),
        "n": n,
        "partition_order_convention": ORDER_CONVENTION,
        "row_col_convention": "table[nu_index][rho_index]",
        "partitions": parts_json,
        "degrees": [str(d) for d in degrees],
        "table": table,
        "certificates_passed": certificates,
        "payload_sha256": payload_sha256(compact(parts_json), compact(table)),
    })


def emit_spot_values(n: int, count: int, out: Path) -> None:
    parts = partitions_desc(n)
    sympy_cross_check(n, parts)
    q = len(parts)
    degrees = {p: hook_degree(p) for p in parts}
    rng = random.Random(1_000_000 + n)

    columns: set[tuple[int, ...]] = {
        tuple([1] * n),                      # identity: values = degrees
        tuple([2] + [1] * (n - 2)),          # transpositions: content formula
        tuple([2, 2] + [1] * (n - 4)),
        tuple([n]),                          # n-cycle: hooks only
    }
    structured_pairs = [(rho, nu) for nu in sorted(columns, reverse=True) for rho in parts]

    hard_rhos = [p for p in parts if p == transpose(p)]          # self-transpose
    hard_rhos += [p for p in parts if len(p) >= 3 and len(set(p)) == 1]  # rectangles
    random_pairs = set()
    while len(random_pairs) < count:
        random_pairs.add((parts[rng.randrange(q)], parts[rng.randrange(q)]))
    for rho in hard_rhos:
        for _ in range(3):
            random_pairs.add((rho, parts[rng.randrange(q)]))

    spots = []
    for rho, nu in structured_pairs + sorted(random_pairs):
        value = chi(rho, nu)
        # certify structured columns by their closed forms
        if nu == tuple([1] * n):
            assert value == degrees[rho]
        elif nu == tuple([2] + [1] * (n - 2)):
            assert value * (n * (n - 1) // 2) == degrees[rho] * content_sum(rho)
        elif nu == tuple([n]):
            is_hook = len(rho) == 1 or all(x == 1 for x in rho[1:])
            assert (value != 0) == is_hook
            if is_hook:
                assert value == (-1) ** (len(rho) - 1)
        spots.append([parts.index(rho), parts.index(nu), str(value)])

    parts_json = [list(p) for p in parts]
    write_json(out / "characters" / f"spot_n{n:02}.json", {
        "schema_version": "fixture.v1",
        "kind": "spot_character_values",
        "generator": meta(seed=1_000_000 + n),
        "n": n,
        "partition_order_convention": ORDER_CONVENTION,
        "row_col_convention": "spots are [rho_index, nu_index, value]",
        "partitions": parts_json,
        "spots": spots,
        "certificates_passed": ["hook_degrees", "content_column", "n_cycle_hooks"],
        "payload_sha256": payload_sha256(compact(parts_json), compact(spots)),
    })


def emit_degrees(n: int, out: Path) -> None:
    parts = partitions_desc(n)
    sympy_cross_check(n, parts)
    degrees = [hook_degree(p) for p in parts]
    class_sizes = [math.factorial(n) // z_value(p) for p in parts]
    assert sum(d * d for d in degrees) == math.factorial(n)
    assert sum(class_sizes) == math.factorial(n)
    tmap = [parts.index(transpose(p)) for p in parts]
    for i, t in enumerate(tmap):
        assert degrees[i] == degrees[t]
    parts_json = [list(p) for p in parts]
    degrees_s = [str(d) for d in degrees]
    sizes_s = [str(c) for c in class_sizes]
    z_s = [str(z_value(p)) for p in parts]
    signs = [sign(p) for p in parts]
    write_json(out / "degrees" / f"deg_n{n:02}.json", {
        "schema_version": "fixture.v1",
        "kind": "degrees_and_class_data",
        "generator": meta(),
        "n": n,
        "partition_order_convention": ORDER_CONVENTION,
        "partitions": parts_json,
        "degrees": degrees_s,
        "class_sizes": sizes_s,
        "z_values": z_s,
        "signs": signs,
        "transpose_map": tmap,
        "certificates_passed": ["degree_squares_sum", "class_size_sum",
                                "transpose_degree_symmetry"],
        "payload_sha256": payload_sha256(
            compact(parts_json), compact(degrees_s), compact(sizes_s),
            compact(z_s), compact(signs), compact(tmap)),
    })


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--full-max", type=int, default=14)
    ap.add_argument("--spot-min", type=int, default=17)
    ap.add_argument("--spot-max", type=int, default=22)
    ap.add_argument("--spot-count", type=int, default=200)
    ap.add_argument("--deg-max", type=int, default=30)
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent.parent / "fixtures")
    args = ap.parse_args()

    for n in range(1, args.full_max + 1):
        emit_full_table(n, args.out)
    for n in range(args.spot_min, args.spot_max + 1):
        emit_spot_values(n, args.spot_count, args.out)
        chi.cache_clear()
    for n in range(1, args.deg_max + 1):
        emit_degrees(n, args.out)
    print("all fixtures written and certified")
    return 0


if __name__ == "__main__":
    sys.exit(main())
