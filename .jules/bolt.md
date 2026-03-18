## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).
## 2026-02-18 - Pandas iterrows/apply(axis=1) vs itertuples
**Learning:** Pandas `apply(axis=1)` creates a Series for every row, which adds massive overhead. `itertuples(index=False)` yields lightweight namedtuples and is ~6x faster than `apply(axis=1)` for applying a complex function to each row, provided the function uses direct attribute access (e.g., `row.column_name`).
**Action:** Replace `df.apply(func, axis=1)` with list comprehensions over `df.itertuples(index=False)` when applying a custom row-wise function to a DataFrame, and use direct attribute access in the function.
