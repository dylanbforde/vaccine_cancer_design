## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2026-03-05 - Pandas DataFrame Row Iteration Optimization
**Learning:** Using `apply(axis=1)` creates a new `pd.Series` object for every single row, causing tremendous overhead, especially on large datasets. `itertuples(index=False)` avoids this by yielding lightweight `NamedTuple` objects, proving significantly faster (e.g., ~40% faster on large DataFrame operations).
**Action:** When row-wise iteration is strictly necessary and vectorization is impossible, use a list comprehension with `itertuples(index=False)` and use `getattr()` in the applied function to handle the `NamedTuple` instead of dictionary-style access.
