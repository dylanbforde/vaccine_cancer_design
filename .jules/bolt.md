## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2026-03-15 - Vectorized Pandas Row Processing with itertuples
**Learning:** `pd.DataFrame.apply(func, axis=1)` is notoriously slow because it creates a new Pandas Series object for every single row. Using a list comprehension over `itertuples(index=False)` avoids this overhead. However, namedtuples from `itertuples` require attribute access (`row.col`) instead of dictionary access (`row["col"]`). Implementing `getattr` with a `try...except AttributeError` fallback allows a function to support both fast vectorization via `itertuples` and traditional `.apply()`/dict processing, yielding a ~4.5x speedup for complex row-wise generation like `generate_peptides`.
**Action:** Replace `df.apply(func, axis=1)` with `[func(row) for row in df.itertuples(index=False)]` when processing millions of rows and ensure the function logic uses dot-notation attribute access with a `try...except` fallback for `dict`/`pd.Series` compatibility.
