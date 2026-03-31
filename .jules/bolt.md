## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2024-05-17 - Pandas Row Iteration Bottleneck
**Learning:** `parsed_df.apply(func, axis=1)` is a massive performance bottleneck because pandas constructs a `Series` object for every single row. Using `.to_dict('records')` is also slow. The fastest generic way to iterate rows as dictionaries is using `[func(dict(zip(cols, row))) for row in df.itertuples(index=False, name=None)]`.
**Action:** Always replace `.apply(axis=1)` with the `itertuples(index=False, name=None)` + `dict(zip(cols, row))` list comprehension pattern for row-wise operations in pandas DataFrames when vectorization is not possible. Ensure the function signature accepts `Any` or `dict` instead of `pd.Series`.
