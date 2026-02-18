## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-03-01 - Regex Pre-compilation
**Learning:** Compiling regex patterns at module level instead of inside a frequently called function (100k+ calls) yielded a ~35% speedup (3.1µs -> 2.0µs per call).
**Action:** Always pre-compile regex patterns when they are used inside loops or frequently called functions.
