## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2026-02-17 - Iterrows Optimization applied
**Learning:** `iterrows` has been successfully replaced with `itertuples(index=False)` in `model/training.py` and `vaccine_design/pipeline.py` which led to faster iteration over DataFrames without series creation overhead.
**Action:** Always prefer `itertuples(index=False)` or vectorization over `iterrows` in the future for iterating pandas DataFrames.
