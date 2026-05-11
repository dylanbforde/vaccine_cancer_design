## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).
## 2026-02-18 - Pandas String Operations Performance
**Learning:** For Pandas `Series` containing strings, using vectorized string methods like `.str.contains()` and `.str.len()` is much faster (approx 3.65x speedup in tests) than using `.apply()` with a lambda function for validation logic.
**Action:** Always prefer `.str` vectorized methods over `.apply()` when processing text/string columns in Pandas DataFrames or Series.
