## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2024-04-14 - DataFrame Row Iteration
**Learning:** Using `df.apply(func, axis=1)` in Pandas is severely slow (~3x slower) compared to a list comprehension that zips columns with `df.itertuples(index=False, name=None)` because it creates a new Series for each row.
**Action:** Always prefer `[func(dict(zip(cols, row))) for row in df.itertuples(index=False, name=None)]` over `df.apply(func, axis=1)` when row-level dict-like data processing is required.
