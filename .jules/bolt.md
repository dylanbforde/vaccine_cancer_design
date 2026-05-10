## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2026-02-18 - Pandas Element-wise Operations vs Vectorized Operations
**Learning:** For row-wise or element-wise string evaluations in Pandas Series, native `apply(lambda ...)` can create severe bottlenecks due to looping. Using vectorized string functions (e.g., `str.contains(..., na=False)` and `str.len()`) significantly boosts performance. Benchmarking showed an improvement of ~70% (~1.4s -> ~0.4s) on 1 million iterations.
**Action:** Replace `apply` and lambdas in string validation/processing over pandas Series with vectorized `.str` functions wherever possible. Use parameters like `na=False` and combinators like `& peptides.notna()` to gracefully handle null values instead of slow logic operators.
