## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2024-05-19 - Pandas Row Iteration Performance Bottleneck
**Learning:** Using `apply(axis=1)` in Pandas is severely bottlenecked due to creating a new `pd.Series` object for every row.
**Action:** Replace `df.apply(func, axis=1)` with list comprehensions over `df.itertuples(index=False)`, and update the target function's signature and access patterns (e.g., using `row.attribute` instead of `row['column']` and avoiding fallback `getattr()`). Ensure tests use `collections.namedtuple` instead of dictionaries or `pd.Series` to mock row inputs accurately.
