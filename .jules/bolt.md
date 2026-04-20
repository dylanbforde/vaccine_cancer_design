## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).
## 2026-02-18 - Pandas String Operations Over apply(lambda)
**Learning:** Using `.apply(lambda ...)` on Pandas Series of strings is significantly slower than using native vectorizer operations (`.str.contains()`, `.str.len()`). Benchmarking showed ~3x speedup. Care must be taken to safely handle `np.nan` values which throw boolean indexing errors unless specifically mitigated via `.fillna(False)` or `& .notna()`.
**Action:** Default to using pandas `.str` accessor methods when dealing with column-level string manipulations and validations.
