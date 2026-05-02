## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2026-05-02 - Caching Edge Index Generation
**Learning:** In graph neural network processing for peptides, graph topologies based strictly on sequence length are highly repetitive. Generating `edge_index` dynamically for every sequence iteration introduces massive overhead.
**Action:** Always memoize graph topology generation when the output only depends on standard static parameters (like sequence length). Using `@staticmethod` with `@functools.lru_cache` perfectly eliminates this overhead without polluting instance caches.
