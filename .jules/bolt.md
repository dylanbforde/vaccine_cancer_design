## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2026-02-18 - Caching edge_index generation in GNN pipeline
**Learning:** Generating the PyTorch Geometric `edge_index` dynamically via loop and tensor allocation is extremely slow for identical graph structures (e.g., thousands of 9-mer peptides). Caching this structure per sequence length provides >100x speedup per call since most peptides have identical lengths.
**Action:** When working with graphs that have static topologies based on size (like sequence data), always cache the adjacency matrix/`edge_index` generation step to avoid redundant structure creation and tensor allocations.
