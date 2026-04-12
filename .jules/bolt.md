## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2026-02-17 - Iteration vs `.apply()` for DataFrame Row Operations
**Learning:** `parsed_df.apply(..., axis=1)` is notoriously slow in pandas because it converts each row to a Series. Using a list comprehension with `zip(columns, row)` over `parsed_df.itertuples(index=False, name=None)` avoids this overhead and yields ~2x performance gains for row-wise operations that process string or dict data.
**Action:** Replace `df.apply(func, axis=1)` with list comprehensions over `itertuples` when row-level function execution is necessary.
