## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-21 - Regex Compilation in Loops
**Learning:** Defining regex patterns inside a function called repeatedly (130k+ times) caused significant overhead due to re-compilation/cache lookups. Pre-compiling them as module constants yielded ~45% speedup on the function.
**Action:** Move regex patterns to module-level constants for high-frequency parsing functions.
