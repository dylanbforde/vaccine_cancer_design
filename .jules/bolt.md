## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2026-02-17 - PyTorch Edge Index Tensor Caching
**Learning:** PyTorch Geometric graph pipeline creates massive numbers of identically sized tensors when processing standard inputs (like 9-mer peptides). Generating redundant `edge_index` tensors repeatedly wastes substantial CPU cycles.
**Action:** Extract standard graph generation methods to `@staticmethod` and use `@functools.lru_cache` to cache deterministic tensor topologies based on input length.
