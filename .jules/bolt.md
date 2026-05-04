## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2026-02-18 - Caching GNN Graph Structures
**Learning:** In Graph Neural Networks (GNNs) handling fixed-size structures (like uniform 9-mer peptides), creating identical `edge_index` tensors repeatedly creates a severe CPU bottleneck. Logic-only benchmarking showed caching `edge_index` with `functools.lru_cache` gives an approximately 99% speedup for tensor creation. To prevent memory leaks by caching the `self` instance, the cached method MUST be converted to a `@staticmethod`.
**Action:** Always extract invariant tensor creation logic into `@staticmethod`s and apply `@functools.lru_cache` to drastically reduce CPU overhead during graph data generation.
