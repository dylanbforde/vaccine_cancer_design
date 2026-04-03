## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2026-02-18 - Pandas df.apply(axis=1) vs itertuples
**Learning:** `df.apply(axis=1)` introduces severe performance overhead for row-wise operations. Benchmarking showed a 2.3x speedup when converting rows to dictionaries via `dict(zip(df.columns, row))` paired with `df.itertuples(index=False, name=None)`.
**Action:** Always avoid `apply(axis=1)`. Use the `itertuples` + `dict(zip)` pattern for complex row-wise processing to maintain performance and prevent namedtuple key-mangling bugs.
