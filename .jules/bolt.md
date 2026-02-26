## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2026-02-17 - JSON Serialization and Pandas Isna Crash
**Learning:** JSON `load` deserializes tuples as lists, causing `isinstance(..., tuple)` checks to fail and caches to be ignored. Furthermore, `pd.isna(list_obj)` crashes with `ValueError` because it returns an array of booleans.
**Action:** Allow `(tuple, list)` in type checks for cached data, convert back to tuple at the boundary if needed, and always check type (e.g., `not isinstance(x, (tuple, list))`) *before* or *instead of* calling `pd.isna(x)` to avoid crashes.
