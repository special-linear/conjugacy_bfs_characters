#!/usr/bin/env python
"""One-time search for adversarial modular tuples (design doc 03 section 5.1).

Finds (n, union, r, nu) where the exact coefficient a_r(nu) is POSITIVE but
divisible by every prime in a small injected prime set (e.g. {11, 13}) — so a
modular engine running on exactly those primes sees all-zero residues and MUST
route the entry through exact certification instead of declaring it zero.
Committed as fixtures/adversarial_v1.json; consumed by the P2 modular-engine
tests (a false-zero layer must not stop the engine, spec section 23 F4/F9).

Tuples where r equals the type's first positive radius in its parity chain
are flagged `masks_first_hit` — those would silently corrupt distances if
screening were trusted.

Usage: python tools/find_adversarial.py [--max-n 9] [--max-r 12] [--out fixtures]
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_fixtures as gf  # noqa: E402  (independent MN + partition machinery)

PRIME_SETS = [[11], [13], [11, 13]]


def coefficients(parts, table, degrees, theta_power, n):
    """a(nu) = (1/n!) sum_rho f_rho chi^rho(nu) theta_power[rho], exact."""
    fact = math.factorial(n)
    out = []
    for v in range(len(parts)):
        numerator = sum(
            degrees[r] * table[v][r] * theta_power[r] for r in range(len(parts))
        )
        a = Fraction(numerator, fact)
        assert a.denominator == 1, "coefficient not integral"
        a = int(a)
        assert a >= 0
        out.append(a)
    return out


def search(max_n: int, max_r: int):
    found = []
    for n in range(6, max_n + 1):
        parts = gf.partitions_desc(n)
        q = len(parts)
        degrees = [gf.hook_degree(p) for p in parts]
        class_sizes = [math.factorial(n) // gf.z_value(p) for p in parts]
        table = [[gf.chi(rho, nu) for rho in parts] for nu in parts]  # [nu][rho]

        # unions: all single classes + all pairs (small n keeps this cheap)
        identity = q - 1
        singles = [[c] for c in range(q) if c != identity]
        pairs = [
            [a, b]
            for i, a in enumerate(range(q))
            if a != identity
            for b in list(range(q))[i + 1:]
            if b != identity
        ]
        for classes in singles + pairs:
            templates = [[x for x in parts[c] if x >= 2] for c in classes]
            theta = []
            for r in range(q):
                w = sum(
                    Fraction(class_sizes[c] * table[c][r], degrees[r]) for c in classes
                )
                assert w.denominator == 1
                theta.append(int(w))

            power = [1] * q
            first_hit = {}  # (nu, parity) -> first positive radius
            for radius in range(1, max_r + 1):
                power = [p * t for p, t in zip(power, theta)]
                coeffs = coefficients(parts, table, degrees, power, n)
                if all(c == 0 for c in coeffs):
                    break  # nilpotent impossible; safety
                for nu, a in enumerate(coeffs):
                    if a <= 0:
                        continue
                    key = (nu, radius % 2)
                    masks = key not in first_hit
                    if masks:
                        first_hit[key] = radius
                    for prime_set in PRIME_SETS:
                        if all(a % p == 0 for p in prime_set):
                            found.append({
                                "n": n,
                                "union_templates": templates,
                                "union_class_indices": classes,
                                "r": radius,
                                "nu_index": nu,
                                "nu": list(parts[nu]),
                                "primes": prime_set,
                                "a_r": str(a),
                                "masks_first_hit": masks,
                            })
    return found


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-n", type=int, default=9)
    ap.add_argument("--max-r", type=int, default=12)
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).resolve().parent.parent / "fixtures")
    args = ap.parse_args()

    found = search(args.max_n, args.max_r)
    # keep the catalog small but pointed: prefer first-hit-masking tuples and
    # both-prime tuples, then trim per (n, prime-set) group
    found.sort(key=lambda t: (not t["masks_first_hit"], -len(t["primes"]),
                              t["n"], t["r"], len(t["a_r"])))
    kept, seen = [], set()
    for t in found:
        key = (t["n"], tuple(t["primes"]), t["masks_first_hit"])
        if sum(1 for k in seen if k == key) >= 6:
            continue
        seen.add(key) if key not in seen else None
        count = sum(
            1 for u in kept
            if (u["n"], tuple(u["primes"]), u["masks_first_hit"]) == key
        )
        if count < 6:
            kept.append(t)

    masking = sum(1 for t in kept if t["masks_first_hit"])
    both = sum(1 for t in kept if len(t["primes"]) == 2)
    print(f"found {len(found)} tuples, keeping {len(kept)} "
          f"({masking} first-hit-masking, {both} with both primes)")
    assert kept, "no adversarial tuples found — widen the search"

    payload = {
        "schema_version": "fixture.v1",
        "kind": "adversarial_modular_tuples",
        "generator": gf.meta(),
        "search": {"max_n": args.max_n, "max_r": args.max_r,
                   "prime_sets": PRIME_SETS},
        "note": "a_r(nu) > 0 but divisible by every listed prime: a modular "
                "engine on exactly these primes sees all-zero residues and "
                "must certify, not stop (spec 23 F4/F9)",
        "tuples": kept,
    }
    out_path = args.out / "adversarial_v1.json"
    out_path.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
