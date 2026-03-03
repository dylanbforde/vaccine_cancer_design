## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2026-02-17 - Pandas Apply vs Itertuples
**Learning:** `df.apply(axis=1)` creates a massive overhead by instantiating a new `pd.Series` object for every single row. Benchmarking showed that using `[func(row) for row in df.itertuples(index=False)]` combined with a refactored target function (using `getattr` instead of dictionary indexing to handle NamedTuples) was up to 6x faster.
**Action:** Replace `df.apply(axis=1)` with a list comprehension over `df.itertuples(index=False)` and use `getattr` inside the row-processing function for dramatic speedups in large datasets.
