## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).
## 2026-02-18 - PyTorch Geometric Edge Index Caching
**Learning:** PyTorch Geometric graph data generation inside loops (like `pipeline.py` processing mutations) involves repeatedly recalculating static properties like `edge_index` for peptides of identical lengths. Because `edge_index` describes sequence proximity (e.g. residues within 3 positions), it is deterministic based purely on the length of the peptide. Repeatedly building nested Python loops and instantiating new PyTorch tensors per iteration is a significant performance bottleneck.
**Action:** Extract length-dependent deterministic tensor creations into `@staticmethod`s and cache them using `@functools.lru_cache()`. This provides near instantaneous $O(1)$ lookup instead of $O(N^2)$ recalculation. Always look for purely deterministic struct calculations inside data generation pipelines.
