## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2024-05-24 - Pandas Row Iteration vs NamedTuples with Zip
**Learning:** `iterrows` in Pandas creates immense performance overhead when iterating row-by-row to generate features, taking around 0.13 seconds for just 1000 rows. A significant boost is achieved using `itertuples(index=False, name=None)` mapped back via `dict(zip(df.columns, row))`, cutting execution time by over 50%. Graph generation via networkx/torch in sequence encoding was also identified as a CPU intensive hot-path that benefits immensely from LRU caching.
**Action:** When iterating rows directly (instead of using vectorization), avoid `iterrows()`. Pair `itertuples(index=False, name=None)` with column zipping. Cache predictable expensive struct/graph creations based on determinable constraints (like peptide length).
