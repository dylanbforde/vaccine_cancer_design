## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).
## 2026-02-18 - Pandas Apply vs Vectorized String Methods\n**Learning:** In pandas, replacing `.apply(lambda ...)` loops with vectorized string operations (like `.str.contains()` and `.str.len()`) for large series significantly improves processing speed. Benchmarking showed a 3x speedup on evaluating amino acid strings and lengths, resolving a bottleneck in data validation.\n**Action:** Use vectorized string operations over `.apply()` loops for Series string evaluations to optimize performance.
