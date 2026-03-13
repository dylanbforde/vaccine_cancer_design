## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2026-02-17 - pandas.DataFrame.apply(axis=1) vs itertuples(index=False)
**Learning:** `pandas.DataFrame.apply(axis=1)` is significantly slower than iterating over rows as NamedTuples via list comprehensions with `itertuples(index=False)`, mostly due to the overhead of creating `pd.Series` objects on each row.
**Action:** Use list comprehensions with `itertuples(index=False)` when executing row-level operations, taking care to extract columns from rows using `getattr(row, col_name, default)` to make the function compatible with `NamedTuples` while optionally falling back to `.get()` for backwards compatibility with `pd.Series`.
