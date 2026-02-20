## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-18 - DNA Translation Loop Optimization
**Learning:** Replacing a loop with list comprehension for DNA translation was significantly SLOWER when stop codons are present. The loop breaks early, while list comprehension processes the entire string.
**Action:** When translating DNA sequences where premature stop codons are expected (e.g., frameshifts), use a loop with early exit.
