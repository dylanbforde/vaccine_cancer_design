## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2026-10-24 - Pandas `apply(axis=1)` Overhead
**Learning:** `DataFrame.apply(axis=1)` creates a pandas Series for every row, which is extremely slow (4.5s for 100k rows). Using list comprehension with `itertuples(index=False)` avoids this overhead and is significantly faster (0.7s, ~6.4x speedup).
**Action:** For row-wise operations that cannot be fully vectorized, prefer `[func(row) for row in df.itertuples(index=False)]` over `df.apply(func, axis=1)`. Ensure helper functions accept attribute access (named tuples) instead of dictionary access.
