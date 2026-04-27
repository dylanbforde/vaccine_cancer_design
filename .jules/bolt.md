## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).
## 2025-02-20 - Vectorized Operations for Series

**Learning:** When dealing with string evaluations inside Pandas Series (like amino acid validations), native Python `lambda` functions inside `.apply()` create a severe bottleneck. The overhead of iterating through Python objects is massive compared to delegating the looping to the underlying C implementations using Pandas' `.str` methods.

**Action:** Whenever iterating over Pandas Series containing strings, search for `.apply(lambda ...)` and systematically replace them with `.str.contains()`, `.str.len()`, or other vectorized `.str` equivalents. Make sure to properly handle missing values (e.g., using `fillna(False)` or bitwise `notna() & ...`) to retain identical evaluation semantics.
