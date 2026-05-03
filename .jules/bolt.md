## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2026-02-18 - Vectorizing Pandas String Validations
**Learning:** Using `.apply(lambda x: ...)` for string validations (e.g., checking amino acids and sequence length) across large Pandas Series is a significant performance bottleneck. Replacing these with vectorized operations like `.str.contains` and `.str.len` yields massive speedups (e.g., from ~0.93s to ~0.25s for 650k rows).
**Action:** Always prefer Pandas vectorized string operations over `.apply()` lambda loops for row-wise evaluations on text data.
