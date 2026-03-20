## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2026-03-20 - Pandas apply/iterrows bottleneck
**Learning:** `df.apply(func, axis=1)` and `df.iterrows()` are massive performance bottlenecks when iterating through large pandas DataFrames because they convert rows into `pd.Series` objects under the hood, carrying large overheads. Tests show that using `df.itertuples(index=False)` inside a list comprehension reduces runtime by up to ~10x.
**Action:** Replace `apply` and `iterrows` with `itertuples(index=False)` combined with list comprehensions when processing DataFrame rows in Python. When parsing the NamedTuple from `itertuples`, use dot notation (e.g., `row.col_name` instead of `row['col_name']`). Avoid fallback logic like `.get` and instead use `getattr(row, "col_name", default)` for optional fields.
