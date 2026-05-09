## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).
## 2026-05-09 - Pandas Vectorization for Data Validation
**Learning:** Using `.apply()` with lambda functions for string checks (like length and character matching) in pandas is significantly slower than using native vectorized string methods (`.str.contains()` and `.str.len()`). However, care must be taken to cast boolean `.sum()` results to `int()` to prevent serialization issues downstream.
**Action:** Always prefer native `.str` accessors for element-wise string operations on Series over `.apply()`. Handle missing values carefully with `.notna()` and cast counts to native python types when inserting into dictionaries.
