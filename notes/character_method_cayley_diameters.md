# Computing Diameters of Conjugacy-Invariant Cayley Graphs of Symmetric Groups via Characters

## 1. Purpose and scope

This document specifies an exact computational approach for determining distances and diameters in Cayley graphs of symmetric groups when the generating set is a union of conjugacy classes.

The principal use cases are:

1. A fixed conjugacy class `C_lambda` is used as the generating set, and one wants the supports of
   \[
   C_\lambda,\ C_\lambda^2,\ C_\lambda^3,\ldots
   \]
   and the distance of every permutation from the identity.
2. Several conjugacy classes are investigated independently for the same `S_n`.
3. Many unions of a selected family of conjugacy classes are investigated.
4. Weighted unions or mixed products are studied after the same preprocessing.

The central idea is that multiplication by a conjugacy-class sum is diagonal in the irreducible-character basis. Therefore repeated products can be computed without constructing the full Cayley graph and, in the one-source setting, without constructing the full class-multiplication matrix.

The output of the main algorithm is:

- the exact distance from the identity to every reachable cycle type;
- the diameter of the identity component;
- optionally, the exact number of words of each length giving each element of a given cycle type;
- optionally, the same data for many generating classes or unions in one batch.

This document is written as an implementation specification. It emphasizes normalization, exactness, stopping criteria, memory layout, validation, and failure modes.

---

## 2. Mathematical setting

Let `S_n` be the symmetric group. Conjugacy classes are indexed by partitions of `n`.

Write a partition as

\[
\lambda=(1^{m_1}2^{m_2}\cdots),
\]

where `m_i` is the number of cycles of length `i`.

Let `C_lambda` be the conjugacy class of cycle type `lambda`, and let

\[
K_\lambda=\sum_{g\in C_\lambda}g
\]

be its class sum in the integral group algebra `Z[S_n]`.

For a union of pairwise distinct conjugacy classes

\[
U=C_{\lambda_1}\cup\cdots\cup C_{\lambda_t},
\]

define

\[
K_U=K_{\lambda_1}+\cdots+K_{\lambda_t}.
\]

More generally, for nonnegative integer weights `w_j`, define

\[
K_w=\sum_j w_jK_{\lambda_j}.
\]

The unweighted union corresponds to `w_j in {0,1}`.

Every conjugacy class of `S_n` is closed under inversion, since a permutation and its
inverse have the same cycle type. Hence the Cayley graphs considered here are undirected
unless the framework is deliberately extended to different groups or non-class generators.

### 2.1 Why cycle types are sufficient

If the generating set `U` is invariant under conjugation by `S_n`, then conjugation is an automorphism of the Cayley graph. Therefore the distance from the identity is constant on every conjugacy class.

Consequently, to compute the diameter it is enough to determine the distance of each cycle type. The quotient by cycle type is exact for one-source distances; it is not a heuristic compression.

### 2.2 The graph may be disconnected

A union of conjugacy classes need not generate all of `S_n`.

- If every generator is even, the identity component is contained in `A_n`.
- If at least one generator is odd, the generated subgroup may be `S_n`, but small-degree exceptions and degenerate generating sets must be handled.
- If the generating set is empty or consists only of the identity, the identity component is trivial.

The algorithm always computes the identity component. The software must distinguish:

- `diameter_of_identity_component`, which is always meaningful for a finite generating set;
- `diameter_of_cayley_graph_on_Sn`, which should be reported as disconnected or infinite unless the generators generate `S_n`.

For `n >= 5`, normal-subgroup structure gives useful shortcuts, but the implementation should not rely on these shortcuts for correctness. Reachability computed by exact powers is authoritative.

---

## 3. Class sizes, signs, and indexing data

For

\[
\lambda=(1^{m_1}2^{m_2}\cdots),
\]

define

\[
z_\lambda=\prod_{i\ge1} i^{m_i}m_i!.
\]

Then

\[
|C_\lambda|=\frac{n!}{z_\lambda}.
\]

The sign of the class is

\[
\operatorname{sgn}(\lambda)=(-1)^{n-\ell(\lambda)},
\]

where `ell(lambda)` is the number of parts.

The implementation must choose one canonical order of the partitions of `n` and use it consistently for:

- conjugacy classes;
- irreducible characters;
- rows and columns of the character table;
- distance arrays;
- union incidence matrices;
- serialized cache files.

A recommended canonical order is reverse lexicographic order or the native order supplied by the chosen character-table backend. Never assume that two libraries use the same order.

Every exported cache should include the explicit partition list and a versioned metadata record.

---

## 4. Character-theoretic diagonalization

Irreducible characters of `S_n` are also indexed by partitions `rho` of `n`.

Let

\[
f_\rho=\chi^\rho(1)
\]

be the degree of the irreducible character `chi^rho`.

The class sum `K_lambda` acts on the irreducible representation indexed by `rho` as the scalar

\[
\omega_\rho(\lambda)
=
\frac{|C_\lambda|\chi^\rho(\lambda)}{f_\rho}.
\]

For symmetric groups, these central character values are integers. The implementation should nevertheless compute them by exact division and assert divisibility:

```text
numerator = class_size[lambda] * character[rho, lambda]
assert numerator % degree[rho] == 0
omega[rho, lambda] = numerator // degree[rho]
```

A failed divisibility assertion almost always indicates one of the following:

- row/column indexing mismatch;
- wrong class-size normalization;
- a character-table orientation error;
- corrupted or approximate character values.

### 4.1 Expansion of a power

The primitive central idempotent corresponding to `rho` is

\[
e_\rho
=
\frac{f_\rho}{n!}
\sum_{\nu\vdash n}\chi^\rho(\nu)K_\nu,
\]

because symmetric-group characters are real-valued.

Hence

\[
K_\lambda^r
=
\sum_\rho \omega_\rho(\lambda)^r e_\rho.
\]

Therefore the coefficient of `K_nu` is

\[
a_r^{(\lambda)}(\nu)
=
\frac{1}{n!}
\sum_{\rho\vdash n}
 f_\rho\chi^\rho(\nu)
 \omega_\rho(\lambda)^r.
\]

This coefficient has a direct combinatorial meaning:

> `a_r^(lambda)(nu)` is the number of length-`r` words in elements of `C_lambda` whose product equals any fixed permutation of type `nu`.

In particular,

\[
\nu\in C_\lambda^r
\quad\Longleftrightarrow\quad
 a_r^{(\lambda)}(\nu)>0.
\]

**Do not Booleanize before taking powers.** The character basis diagonalizes the
integer class-sum multiplication operator, whose entries count factorizations. It does
not diagonalize the Boolean support matrix under ordinary arithmetic. The correct order
is:

1. take spectral powers in the integer class algebra;
2. apply the exact inverse character transform;
3. test the resulting nonnegative integer coefficients for positivity.

Replacing eigenvalues or intermediate coefficients by `0/1` values before the inverse
transform gives the wrong reachability relation.

### 4.2 Expansion for a union

For

\[
U=C_{\lambda_1}\cup\cdots\cup C_{\lambda_t},
\]

the eigenvalue on `rho` is

\[
\theta_\rho(U)
=
\sum_{j=1}^t\omega_\rho(\lambda_j).
\]

Then

\[
a_r^{(U)}(\nu)
=
\frac{1}{n!}
\sum_{\rho\vdash n}
 f_\rho\chi^\rho(\nu)
 \theta_\rho(U)^r.
\]

This automatically includes every mixed word in the chosen classes.

For a weighted union `w`, replace the eigenvalue by

\[
\theta_\rho(w)=\sum_j w_j\omega_\rho(\lambda_j).
\]

If all weights are positive, the support is the same as for the corresponding unweighted union. Negative weights must not be used for reachability because cancellation destroys the counting interpretation.

---

## 5. Distance and diameter from exact powers

For a generating set `U`, define

\[
d_U(\nu)=\min\{r\ge0:a_r^{(U)}(\nu)>0\}.
\]

Then `d_U(nu)` is exactly the Cayley-graph distance from the identity to every permutation of cycle type `nu`.

The diameter of the identity component is

\[
\operatorname{diam}(U)=\max_{\nu:\ d_U(\nu)<\infty}d_U(\nu).
\]

This maximum of distances from the identity is the graph diameter because the identity
component is a Cayley graph of the generated subgroup and is vertex-transitive. Moreover,
the subgroup generated by an `S_n`-conjugacy-invariant set is normal in `S_n`, so every
`S_n` conjugacy class is either wholly contained in the identity component or disjoint
from it.

### 5.1 Exact-length supports are not BFS frontiers

The support of `K_U^r` is the set of vertices reachable in exactly `r` steps. It generally contains vertices seen at smaller radii as well.

Therefore maintain:

- `support_r`: classes with positive coefficient at exact length `r`;
- `visited`: classes seen at any length at most `r`;
- `new_r = support_r & ~visited_before_r`.

Assign distance `r` only to `new_r`.

### 5.2 Correct stopping criterion

For a fixed union, once `new_r` is empty, no new class can appear at a later radius.

Reason: every class in `support_r` was already reached earlier; all of its generator-neighbors were therefore reachable no later than radius `r`. Thus radius `r+1` cannot introduce a new class.

This stopping rule is valid only when `support_r` has been determined exactly. It is not valid if some zero/nonzero decisions remain unresolved after modular screening.

### 5.3 Identity in the generating set

The identity class should normally be excluded from a Cayley generating set.

If it is included:

- graph distances are unchanged;
- word counts change;
- exact-length supports become padded by shorter words;
- the spectral formulas remain correct.

The API should either reject the identity class by default or require an explicit option such as `allow_identity_generator=True`.

---

## 6. Shared preprocessing for several classes

Fix `n` and let `q=p(n)` be the number of partitions of `n`.

Suppose the selected base classes are

\[
\lambda_1,\ldots,\lambda_t.
\]

Construct the central-eigenvalue matrix

\[
\Omega=(\omega_{\rho j})_{\rho,j},
\qquad
\omega_{\rho j}=\omega_\rho(\lambda_j).
\]

The expensive data shared by all base classes and all of their unions is:

1. the partition index;
2. the degrees `f_rho`;
3. the class sizes;
4. the relevant character values `chi^rho(nu)`;
5. the matrix `Omega`;
6. modular reductions and block layouts;
7. parity and transpose-partition metadata.

A union is represented by a column vector `b in {0,1}^t`, and its spectrum is

\[
\theta(U)=\Omega b.
\]

For `m` unions, collect their incidence vectors as columns of a `t x m` matrix `B`. Then

\[
\Theta=\Omega B
\]

contains all union spectra.

At radius `r`, define

\[
W_r(\rho,j)=f_\rho\Theta(\rho,j)^r.
\]

This reuse is specific to a fixed value of `n`. When `n` changes, the partition index,
character table, degrees, class sizes, and transform all change. Only the software,
recursion caches of a custom evaluator, and more advanced stable-center constructions
can be reused across different degrees.

If `X` is the character table with rows `rho` and columns `nu`, the numerator matrix is

\[
N_r=X^{\mathsf T}W_r,
\]

and the desired coefficients are

\[
A_r=N_r/n!.
\]

The division is entrywise and exact.

### Important implementation recommendation

Do not materialize the rational transform

\[
F_{\nu\rho}=f_\rho\chi^\rho(\nu)/n!.
\]

Instead compute the integer numerator `N_r`, then divide by `n!`. This avoids rational matrices and makes exactness checks straightforward.

---

## 7. Do we need the full character table?

Not always.

For one class `lambda`, only irreducibles with

\[
\chi^\rho(\lambda)\ne0
\]

contribute for positive radii.

For several base classes, define the safe active set

\[
R=\bigcup_j\{\rho:\chi^\rho(\lambda_j)\ne0\}.
\]

It is enough to retain the restricted rectangle

\[
X_R=(\chi^\rho(\nu))_{\rho\in R,\ \nu\vdash n}.
\]

For a particular union, additional cancellation may make

\[
\theta_\rho(U)=0
\]

for some rows in `R`; those rows may be dropped for that union **for every positive power**.

**Radius-zero exception:** zero-eigenvalue rows still contribute to the decomposition of
`K_U^0=1`. Therefore an implementation using restricted rows should initialize the
identity distance directly and begin spectral transforms at `r=1`. It must not expect a
restricted-row transform at `r=0` to reconstruct the identity class sum.

### 7.1 When restriction helps

Restriction can be dramatic for classes whose character columns are sparse. For example, an `n`-cycle has nonzero values only on hook representations.

For short cycles with many fixed points, the character column may be much denser. Do not assume that small support implies a sparse character column.

### 7.2 Streaming instead of storing

If the restricted table is still too large, process target cycle types in blocks:

1. generate or load a block of character columns;
2. multiply the block by all active union power vectors;
3. classify coefficients;
4. discard the block.

For many unions, each streamed block should be reused across the entire union batch before it is discarded.

---

## 8. Core exact algorithm

### 8.1 Reference algorithm for one union

Inputs:

- `n`;
- ordered partition list `P`;
- exact degrees `degree[rho]`;
- exact character data `X[rho, nu]` for active rows;
- exact central eigenvalue vector `theta[rho]`;
- `factorial_n`.

State:

- `distance[nu] = -1` initially;
- `distance[id_type] = 0`;
- `visited` bitset containing only the identity type;
- `power[rho] = 1`, representing `theta[rho]^0`.

Loop for `r = 1,2,...`:

```text
for rho in active_rows:
    power[rho] *= theta[rho]
    weighted[rho] = degree[rho] * power[rho]

numerator[nu] = sum_rho X[rho, nu] * weighted[rho]

for nu:
    assert numerator[nu] % factorial_n == 0
    coefficient[nu] = numerator[nu] // factorial_n
    assert coefficient[nu] >= 0

support = bitset(coefficient[nu] > 0)
new = support & ~visited

for nu in new:
    distance[nu] = r

if new is empty:
    stop

visited |= new
```

Return:

- `distance`;
- `diameter = max(distance[nu] for distance[nu] >= 0)`;
- `reachable_types = visited`.

### 8.2 Batched algorithm for many unions

For `m` unions, use matrices with union columns.

Maintain:

- `Theta[rho, j]`;
- `Power[rho, j]`, initialized to `1`;
- independent `visited[j]`, `distance[:,j]`, and `active_union[j]` flags.

At radius `r`:

```text
Power[:, active_unions] *= Theta[:, active_unions]
W[:, active_unions] = degree[:, None] * Power[:, active_unions]
N[:, active_unions] = X.T @ W[:, active_unions]
A[:, active_unions] = exact_divide(N, factorial_n)
```

Then update every union independently. A union can leave the batch once its `new` set is empty.

To avoid wasting work on finished unions, compact active union columns periodically.

---

## 9. Strong validation invariants

A production implementation should run these checks in debug mode and on every small test case.

### 9.1 Radius zero

With the **full** irreducible row set, `r=0` must give

\[
K^0=1=K_{(1^n)}.
\]

Thus the coefficient vector must be `1` at the identity type and `0` elsewhere.

If zero-eigenvalue rows have been removed, this test is no longer valid on the restricted
transform: the reduction is valid only for `r >= 1`. In that implementation, initialize
the identity type manually at distance zero and use the radius-zero identity test only
against the full character table in backend-validation tests.

### 9.2 Radius one

For a single class `lambda`, the coefficient vector at `r=1` must be the indicator of `lambda`.

For an unweighted union, it must be the indicator of the included classes.

For a weighted union, it must equal the input weights.

This is one of the best tests for row/column orientation and normalization.

### 9.3 Total-word count

For every radius,

\[
\sum_{\nu\vdash n}|C_\nu|a_r^{(U)}(\nu)=|U|^r.
\]

For a weighted union, replace `|U|` by

\[
\sum_j w_j|C_{\lambda_j}|.
\]

This detects many silent arithmetic errors.

### 9.4 Exact divisibility

Every numerator must be divisible by `n!`.

Never silently round or truncate.

### 9.5 Nonnegativity

Every coefficient for a nonnegative generating combination must be a nonnegative integer.

A negative result indicates a bug or arithmetic overflow.

### 9.6 Parity

If all generators in a union have the same sign `epsilon`, then at radius `r` only classes of sign `epsilon^r` can occur.

If the union contains both even and odd generators, this simple filter is invalid.

### 9.7 Small-degree direct verification

For small `n`, enumerate the full group and compare:

- exact supports of powers;
- distances by ordinary BFS;
- word counts for several radii;
- diameters.

This should be part of the automated test suite.

---

## 10. Character-table backends

There are three realistic implementation routes.

### 10.1 GAP or Sage as the authoritative backend

Use GAP's generic symmetric-group character tables or Sage's interfaces.

Advantages:

- reliable exact values;
- low implementation risk;
- useful for reference and validation;
- class and character metadata are already available.

Caveats:

- object overhead can be substantial;
- exporting very large tables may dominate memory;
- native class order must be recorded explicitly;
- high-level routines may materialize more data than needed.

A practical design is to use GAP to generate exact rows or blocks and serialize them into a compact binary format consumed by a faster implementation in Rust, C++, or Python/NumPy with modular backends.

### 10.2 Custom Murnaghan-Nakayama evaluator

A custom evaluator can compute `chi^rho(nu)` without constructing the entire table.

This route is useful for streaming and restricted-row computation, but it is easy to implement incorrectly.

Warnings:

- rim-hook enumeration is subtle;
- signs are `(-1)^(height-1)`;
- repeated cycle lengths require correct recursive multiplicity handling;
- partition normalization and removal order must be consistent;
- naive recursion recomputes enormous subproblems;
- memoization keys must include both the current partition and the remaining cycle multiset;
- integer overflow is possible even when final values are moderate.

A custom evaluator must be cross-checked against GAP for many random pairs `(rho, nu)` before being trusted.

### 10.3 Hybrid backend

Recommended for larger experiments:

1. GAP supplies partition order, degrees, class sizes, and spot checks.
2. A custom or exported character backend provides blocks.
3. The power engine performs batched exact or modular transforms.

---

## 11. Symmetries and reductions

### 11.1 Transposed irreducible partitions

Let `rho'` be the transpose partition. Then

\[
\chi^{\rho'}(\nu)=\operatorname{sgn}(\nu)\chi^\rho(\nu).
\]

Also

\[
f_{\rho'}=f_\rho.
\]

For a single class `lambda`,

\[
\omega_{\rho'}(\lambda)
=
\operatorname{sgn}(\lambda)\omega_\rho(\lambda).
\]

This can nearly halve stored character rows, apart from self-transpose partitions.

However, the paired-row reconstruction must be implemented carefully, especially for unions containing both even and odd classes. For a mixed-parity union,

\[
\theta_{\rho'}(U)
\ne \pm\theta_\rho(U)
\]

in general, because the even and odd contributions transform differently.

A safe decomposition is

\[
\theta_\rho(U)=\theta_\rho^{\mathrm{even}}+\theta_\rho^{\mathrm{odd}},
\]

\[
\theta_{\rho'}(U)=\theta_\rho^{\mathrm{even}}-\theta_\rho^{\mathrm{odd}}.
\]

Do not apply the simple sign-pair formula to mixed-parity unions.

### 11.2 Self-transpose rows

For self-transpose `rho`, the character vanishes on odd conjugacy classes. This is a useful check and may remove rows for odd-only generators.

### 11.3 Parity filtering of target columns

If all generators have one sign, skip all target classes of impossible parity at a given radius. This can almost halve the target work.

For mixed parity, no such single-parity filter is available.

---

## 12. Exact arithmetic strategies

The defining sums contain large cancellations. Floating-point arithmetic must not be used to decide whether a coefficient is zero.

### 12.1 Pure arbitrary-precision integers

This is the simplest correct implementation.

Advantages:

- easy correctness story;
- exact divisibility and nonnegativity checks;
- straightforward debugging.

Disadvantages:

- big-integer matrix multiplication is expensive;
- coefficients grow approximately exponentially with the radius;
- generic dense linear algebra libraries usually do not optimize big integers.

Use this as the reference implementation.

### 12.2 Modular screening with exact fallback

Choose primes `p` such that `p` does not divide `n!`; taking `p > n` is sufficient.

For each prime, compute

\[
a_r(\nu)\bmod p.
\]

A nonzero residue proves the exact nonnegative integer is positive.

A zero residue does **not** prove the coefficient is zero. A positive coefficient may be divisible by the chosen prime.

Recommended workflow:

1. evaluate all candidates modulo several primes;
2. mark any nonzero residue as definitely positive;
3. collect entries that are zero modulo every prime;
4. compute those entries exactly;
5. only then update `new` and test the stopping criterion.

Never label a coefficient zero merely because it vanished modulo one or several primes.

### 12.3 CRT certification

If the product `M` of used primes exceeds a rigorous upper bound for a coefficient, then residues determine the coefficient uniquely.

The elementary bound

\[
0\le a_r(\nu)\le |U|^r
\]

is always valid but quickly becomes enormous. Full CRT certification may therefore require many primes.

CRT is most useful when:

- radii are small;
- generator sets are small;
- sharper coefficient bounds are available;
- exact reconstruction of all counts is desired.

For diameter-only computation, modular screening plus selective exact fallback is usually better.

### 12.4 Modular division by `n!`

If working modulo `p`, compute

\[
a_r(\nu)\equiv N_r(\nu)(n!)^{-1}\pmod p.
\]

This requires `p` not to divide `n!`.

Alternatively, if central eigenvalues and all arithmetic are constructed over the integers first, one may reduce the exact quotient, but that defeats the purpose of modular acceleration.

---

## 13. Machine arithmetic and GPU caveats

Dense modular transforms are attractive for SIMD, multicore CPUs, and GPUs, but naive implementations overflow.

### 13.1 Integer overflow

Suppose residues are stored in signed 64-bit integers. A dot product of length `R` involving values near `p` can overflow long before the final modular reduction.

Do not write

```text
sum += a * b
```

for an entire long row unless a proved bound guarantees safety.

Use one of:

- small primes and blockwise reduction;
- 128-bit accumulators on CPUs;
- Montgomery or Barrett reduction;
- residue-number batching with bounded block lengths;
- limb decomposition.

### 13.2 Floating-point GEMM is not automatically exact

Using `float64` matrix multiplication for modular arithmetic is only safe if every intermediate integer is proved to stay below `2^53` before reduction. This bound must include the accumulation length.

Do not assume that integral inputs imply exact integral outputs.

### 13.3 GPU recommendation

Implement and validate the CPU exact version first. Then add a modular batched backend.

A good GPU design uses:

- several small primes;
- tiled matrix multiplication;
- frequent modular reduction;
- target-column blocks;
- many unions in the batch, so the character block is reused.

GPU acceleration is least effective when only one union and very few target types remain unresolved.

---

## 14. Memory layout and blocking

Let `R` be the number of retained character rows and `q=p(n)`.

The restricted character rectangle has `R*q` exact integers. It may dominate memory.

Recommended layouts:

- store character data row-major if generating by irreducible row;
- store column blocks if the primary operation is `X_block.T @ W`;
- for many unions, choose blocks large enough to amortize loading but small enough for cache or GPU memory;
- store modular copies separately by prime;
- avoid Python objects per entry.

Suggested block API:

```text
for target_block in character_source.iter_target_blocks():
    # X_block shape: (R, block_size)
    numerator_block = X_block.T @ W
    classify(numerator_block)
```

The character source may be:

- an in-memory array;
- a memory-mapped file;
- a compressed row store;
- an on-demand Murnaghan-Nakayama generator;
- a GAP subprocess/export.

Do not store coefficient matrices for every radius. For diameter computation, keep only:

- current spectral powers;
- current support;
- visited bitsets;
- distance arrays.

---

## 15. Repeated eigenvalues and recurrence optimization

For a fixed union `U`, many rows may share the same eigenvalue `theta_rho(U)`.

Let `E_U` be the set of distinct eigenvalues. Group rows with equal eigenvalue and define
the **integer numerator blocks**

\[
H_{\alpha}(\nu)
=
\sum_{\rho:\theta_\rho(U)=\alpha}
 f_\rho\chi^\rho(\nu).
\]

Then

\[
a_r^{(U)}(\nu)
=
\frac{1}{n!}
\sum_{\alpha\in E_U}H_\alpha(\nu)\alpha^r.
\]

Keep `H_alpha` integral and perform one exact division by `n!` at the end. The
individually normalized quantities `H_alpha/n!` need not be integers.

This can reduce per-radius cost when the number of distinct eigenvalues is much smaller than the number of active rows.

Tradeoff:

- grouping is union-specific;
- it is less reusable when many unions are batched;
- constructing all grouped transforms may cost more than it saves for short diameters.

Use it when one union requires many radii or when profiling shows a large eigenvalue multiplicity.

Each sequence `a_r^(U)(nu)` satisfies a linear recurrence whose characteristic polynomial divides

\[
\prod_{\alpha\in E_U}(x-\alpha).
\]

This is useful for isolated large powers, but ordinary iterative powers are usually simpler for finding first occurrence and diameter.

---

## 16. Many unions of a base family

Suppose central eigenvalue columns are known for base classes

\[
\lambda_1,\ldots,\lambda_t.
\]

Then any union is specified by a bit vector `b` and requires only

\[
\theta=\Omega b.
\]

No new character values are needed.

For many unions, batch them as columns of `B`:

\[
\Theta=\Omega B.
\]

### 16.1 Exponential output caveat

There are `2^t` possible unions. Character theory removes repeated preprocessing but does not remove the exponential number of requested outputs.

If all unions are requested, use:

- chunked union batches;
- symmetry or monotonicity filters where mathematically justified;
- early stopping per union;
- deduplication of identical input bit vectors;
- compressed distance output.

### 16.2 Monotonicity of distances under inclusion

If `U subseteq V`, then

\[
d_V(g)\le d_U(g)
\]

for every `g` reachable under `U`.

Therefore

\[
\operatorname{diam}(V)
\]

is not necessarily monotone in a naive way when the generated subgroup changes, but on a fixed common vertex set additional generators cannot increase distances.

This can provide bounds and pruning, but it does not by itself determine diameters of all unions.

### 16.3 Mixed products

After `Omega` and the shared transform are known, arbitrary mixed products

\[
K_{\lambda_1}^{r_1}\cdots K_{\lambda_t}^{r_t}
\]

use the spectral vector

\[
\prod_j\omega_\rho(\lambda_j)^{r_j}.
\]

This machinery is optional for diameter computation but useful for verifying class-product identities.

---

## 17. Reachability, subgroup interpretation, and special cases

### 17.1 Odd and even generators

- Odd-only generating sets produce a bipartite Cayley graph on the generated subgroup.
- Even-only generating sets remain inside `A_n`.
- Mixed-parity generating sets have no simple fixed parity at radius `r`.

### 17.2 Split classes in `A_n`

A crucial warning:

The method described here indexes states by `S_n` cycle types and assumes every generator is a union of full `S_n` conjugacy classes.

Some `S_n` classes split into two conjugacy classes in `A_n`. If the intended generating set is only one half of a split `A_n` class, then:

- cycle type is no longer a sufficient state label;
- the `S_n` character table is insufficient;
- one must use the `A_n` character table and distinguish split classes.

If the generating set contains the entire `S_n` class, even when it lies in `A_n`, the present method remains valid.

### 17.3 Non-conjugacy-invariant generators

If the generator set is not a union of full conjugacy classes, distance need not be constant on cycle types. The entire reduction fails.

Do not apply this machinery to arbitrary selected permutations.

### 17.4 Identity-only or empty generators

Handle explicitly. Spectral formulas work, but generic loop logic may otherwise produce confusing stopping behavior.

---

## 18. What the method does not provide

The character method gives exact reachability, distances, diameters, and word counts. It does not directly give:

- a concrete shortest word for a target permutation;
- predecessor pointers in the full Cayley graph;
- a canonical factorization witnessing membership in `U^r`.

To obtain witnesses, one needs an additional constructive layer, for example:

1. determine the target distance by characters;
2. work backward through class products to identify possible predecessor cycle types;
3. solve a relative-position or double-coset problem to construct actual permutations;
4. continue recursively.

Do not promise explicit geodesics from the spectral computation alone.

---

## 19. Recommended software architecture

### 19.1 Data types

```text
PartitionIndex
    n
    partitions: list[Partition]
    partition_to_index: map
    identity_index
    sign[index]
    z_value[index]
    class_size[index]

CharacterData
    partition_index
    degrees[rho]
    active_rows
    character_source
    transpose_row[rho]

BaseClassSpectra
    base_class_indices
    omega[rho, base_class]

UnionBatch
    incidence[base_class, union]
    theta[rho, union]
    active_union_mask

DiameterState
    distance[type, union]
    visited_bitsets[union]
    current_radius
    spectral_power[rho, union]
```

### 19.2 Main modules

1. `partitions`
   - partition generation;
   - signs;
   - `z_lambda`;
   - class sizes;
   - transpose partitions.
2. `characters`
   - backend abstraction;
   - degrees;
   - exact value blocks;
   - GAP import or Murnaghan-Nakayama.
3. `spectra`
   - central eigenvalues;
   - union spectra;
   - active-row reduction.
4. `transform`
   - exact integer backend;
   - modular backend;
   - block multiplication;
   - exact fallback.
5. `diameter`
   - radius iteration;
   - support classification;
   - distance updates;
   - stopping logic.
6. `validation`
   - radius-zero and radius-one checks;
   - word-count sum;
   - parity;
   - direct small-`n` comparison.
7. `cache`
   - partition-order metadata;
   - character blocks;
   - spectra;
   - checksums and versioning.

### 19.3 Cache warnings

A cached character array without its explicit partition order is unsafe.

Cache metadata should include:

- `n`;
- partition list hash;
- character backend and version;
- row and column conventions;
- integer encoding;
- active-row list;
- base-class list;
- checksum.

---

## 20. Implementation phases

### Phase 1: exact reference implementation

Target small and moderate `n`.

- use GAP/Sage character data;
- use arbitrary-precision integers;
- support one class and one union;
- implement all invariants;
- compare against direct BFS for small `n`.

### Phase 2: selected classes and union batching

- construct `Omega` for many base classes;
- support incidence matrices;
- batch unions;
- compact finished union columns;
- add parity filtering.

### Phase 3: restricted rows and streaming

- determine active-row union;
- stream target blocks;
- use memory-mapped character data;
- avoid materializing the full table when possible.

### Phase 4: modular acceleration

- implement several safe primes;
- use modular screening;
- exact fallback for unresolved zeros;
- validate every result against the reference backend.

### Phase 5: multicore/GPU optimization

- tiled modular matrix multiplication;
- blockwise reduction with proved overflow bounds;
- many-union batching;
- performance profiling.

### Phase 6: optional advanced features

- repeated-eigenvalue grouping;
- linear recurrences;
- mixed-product queries;
- witness reconstruction support;
- stable-center methods across varying `n`.

---

## 21. Pseudocode for a robust batched implementation

```text
function prepare(n, base_cycle_types, character_backend):
    P = build_partition_index(n)
    q = len(P.partitions)
    fact = factorial(n)

    degrees = character_backend.degrees(P)

    base_idx = [P.index(lambda) for lambda in base_cycle_types]
    class_sizes = P.class_size

    # First obtain only the character columns of the base classes.
    chi_base = character_backend.values(
        rows=P.partitions,
        columns=base_cycle_types,
    )

    omega = zeros_exact(q, len(base_idx))
    active_rows = empty_set()

    for rho in range(q):
        for j, lambda_idx in enumerate(base_idx):
            num = class_sizes[lambda_idx] * chi_base[rho, j]
            assert num % degrees[rho] == 0
            omega[rho, j] = num // degrees[rho]
            if omega[rho, j] != 0:
                active_rows.add(rho)

    character_source = character_backend.restricted_source(
        rows=sorted(active_rows),
        target_columns=P.partitions,
    )

    return PreparedData(P, fact, degrees, omega, active_rows, character_source)
```

```text
function compute_union_diameters(prepared, incidence_matrix):
    P = prepared.partition_index
    rows = prepared.active_rows
    fact = prepared.factorial_n

    Theta = prepared.omega[rows, :] @ incidence_matrix

    # Optional per-union row cancellation.
    # Keep a common row set for simple batching; specialize only if worthwhile.

    m = number_of_columns(incidence_matrix)
    distance = fill(-1, shape=(len(P), m))
    visited = [empty_bitset(len(P)) for j in range(m)]

    for j in range(m):
        distance[P.identity_index, j] = 0
        visited[j].add(P.identity_index)

    Power = ones_exact(shape=(len(rows), m))
    active_unions = bitset(all union columns)
    r = 0

    while active_unions not empty:
        r += 1
        Power[:, active_unions] *= Theta[:, active_unions]
        W = prepared.degrees[rows, None] * Power[:, active_unions]

        # Results assembled target block by target block.
        support = [empty_bitset(len(P)) for j in active_unions]

        for block in prepared.character_source.iter_target_blocks():
            # block.values shape: (len(rows), block_size)
            N = block.values.T @ W

            for local_union, j in enumerate(active_unions):
                for local_target, nu in enumerate(block.target_indices):
                    value = N[local_target, local_union]
                    assert value % fact == 0
                    coeff = value // fact
                    assert coeff >= 0
                    if coeff > 0:
                        support[j].add(nu)

        finished = []
        for j in active_unions:
            new = support[j] - visited[j]
            if new.empty():
                finished.append(j)
            else:
                for nu in new:
                    distance[nu, j] = r
                visited[j] |= new

        active_unions.remove_all(finished)

    diameters = []
    for j in range(m):
        diameters.append(max(distance[:, j]))

    return distance, diameters, visited
```

For a modular backend, the block classification stage must not finalize a zero until exact fallback or rigorous CRT certification has been completed.

---

## 22. Test plan

### 22.1 Unit tests

- partition generation and sums;
- `z_lambda` and class sizes;
- hook-length degrees;
- signs;
- transpose partitions;
- central-eigenvalue divisibility;
- radius-zero transform;
- radius-one transform;
- exact word-count sum.

### 22.2 Character tests

Compare random values against GAP for a grid of `n`, `rho`, and `nu`.

Include difficult cases:

- repeated cycle lengths;
- long rim hooks;
- self-transpose partitions;
- odd classes where self-transpose rows should vanish;
- identity class, where values equal degrees.

### 22.3 End-to-end tests

For small `n`:

1. enumerate `S_n`;
2. build each requested generating union explicitly;
3. run ordinary BFS;
4. aggregate distances by cycle type;
5. compare with the character method.

Also compare exact word counts at radii `0,1,2,3`.

### 22.4 Union tests

For two classes `A` and `B`, verify explicitly that

\[
(K_A+K_B)^2=K_A^2+2K_AK_B+K_B^2.
\]

At the support level, compare against direct mixed-word enumeration.

### 22.5 Modular tests

- compare every modular residue against exact computation;
- include coefficients divisible by one or more chosen primes;
- verify that zero residues are never treated as certified zero without fallback;
- test overflow bounds for every kernel configuration.

---

## 23. Common implementation failures

### Failure 1: wrong normalization

Symptom: radius one is not the input class indicator.

Likely cause: missing factor `f_rho`, incorrect `n!`, or transposed character table.

### Failure 2: confusing class coefficients with total words into a class

The coefficient of `K_nu` counts words reaching one fixed element of class `nu`.

The total number of words whose product lies somewhere in `C_nu` is

\[
|C_\nu|a_r(\nu).
\]

Do not multiply by the class size twice.

### Failure 3: floating-point zero tests

Small numerical values after cancellation are not evidence of exact zero.

Use exact or certified modular arithmetic.

### Failure 4: stopping after an unresolved modular-zero layer

A layer is complete only after every candidate zero has been certified.

### Failure 5: assuming parity filtering for mixed unions

If a union contains both even and odd classes, both target parities may occur at the same radius.

### Failure 6: applying the `S_n` quotient to one split `A_n` class

Cycle type no longer distinguishes the generator orbit. Use `A_n` data.

### Failure 7: losing partition order in serialization

Character values attached to the wrong partitions can still satisfy superficial size checks and produce meaningless results.

### Failure 8: integer overflow in a modular dot product

The final residue may look plausible. Overflow can remain silent. Prove accumulator bounds or reduce in blocks.

### Failure 9: assuming no new types after a modularly apparent empty layer

Only exact emptiness justifies termination.

### Failure 10: Booleanizing the spectral calculation

The character transform applies to integer class-sum multiplication. Converting
intermediate spectral data or class coefficients to `0/1` before completing the exact
inverse transform is incorrect.

### Failure 11: expecting geodesic witnesses

The output is distance by type, not an explicit shortest factorization.

---

## 24. Practical recommendations

For an initial coding agent implementation:

1. Use GAP or Sage for exact character data.
2. Implement the integer numerator formula, followed by exact division by `n!`.
3. Validate radius zero, radius one, total word counts, divisibility, and nonnegativity.
4. Compute distances from exact supports without building the class-multiplication matrix.
5. Add several base classes and represent unions by an incidence matrix.
6. Batch unions at each radius.
7. Store distances and supports, not all coefficients at all radii.
8. Add restricted active rows and target blocking only after the reference version passes exhaustive small-`n` tests.
9. Add modular acceleration only with exact fallback.
10. Treat split `A_n` classes, non-conjugacy-invariant generators, and explicit witnesses as separate extensions.

The most reusable computational object for fixed `n` is the character transform together with the central-eigenvalue columns of the selected base classes. Once these are available, a new union requires only a column sum in spectral coordinates, followed by the same repeated inverse transforms.

---

## 25. Summary of the central formulas

For a class `C_lambda`:

\[
\omega_\rho(\lambda)
=
\frac{|C_\lambda|\chi^\rho(\lambda)}{f_\rho}.
\]

For a union `U`:

\[
\theta_\rho(U)
=
\sum_{C_\lambda\subseteq U}\omega_\rho(\lambda).
\]

The number of length-`r` words giving any fixed element of type `nu` is

\[
a_r^{(U)}(\nu)
=
\frac{1}{n!}
\sum_{\rho\vdash n}
 f_\rho\chi^\rho(\nu)\theta_\rho(U)^r.
\]

Reachability and distance are

\[
\nu\in U^r
\iff
 a_r^{(U)}(\nu)>0,
\]

\[
d_U(\nu)=\min\{r:a_r^{(U)}(\nu)>0\}.
\]

The identity-component diameter is

\[
\operatorname{diam}(U)
=
\max_{\nu:d_U(\nu)<\infty}d_U(\nu).
\]

For many unions with incidence matrix `B`:

\[
\Theta=\Omega B,
\]

and at each radius

\[
N_r=X^{\mathsf T}
\left(\operatorname{diag}(f_\rho)\Theta^{\circ r}\right),
\qquad
A_r=N_r/n!.
\]

All zero tests used for reachability must be exact or rigorously certified.
