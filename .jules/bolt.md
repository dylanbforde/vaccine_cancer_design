## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2026-04-05 - Pandas DataFrame Row Iteration Optimization
**Learning:** For applying row-wise operations that process complex custom logic via a custom function (like `generate_peptides`), the native `pd.DataFrame.apply(axis=1)` method has substantial overhead. Using a list comprehension with `df.itertuples(index=False, name=None)` and converting rows to dicts using `dict(zip(df.columns, row))` is significantly faster.
**Action:** Replace `df.apply(func, axis=1)` with list comprehensions and `to_dict('records')` for compute-heavy row iterations where vectorization is not possible, but ensure `func` gracefully handles the dict input instead of expecting a pandas Series.
