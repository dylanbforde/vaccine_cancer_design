## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2024-05-19 - Pandas DataFrame apply(axis=1) Iteration Performance Fix
**Learning:** `DataFrame.apply(axis=1)` is notoriously slow because it instantiates a `pd.Series` object for every single row. During testing, iterating over 10,000 rows with a moderate-sized dictionary payload using `apply(axis=1)` was roughly 4-5x slower than optimal iteration loops.
**Action:** Always prefer `itertuples(index=False)` with a list comprehension instead of `apply(axis=1)`. To maintain compatibility for downstream functions expecting standard dictionary indexing via row string keys, you can add `row = row._asdict() if hasattr(row, '_asdict') else row` at the start of your row processing function to cheaply convert the incoming `namedtuple` back into a dictionary. This provides significant performance gains while mitigating breaking logic changes!
