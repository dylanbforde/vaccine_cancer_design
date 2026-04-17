## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2026-02-18 - Pandas Graph Construction Bottleneck
**Learning:** `iterrows()` is incredibly slow for constructing torch graph Data objects. Coupling this with redundant graph edge creation (via `create_edge_index`) for uniform sequence lengths (e.g., 9-mer peptides) wastes substantial CPU time. Combined, removing iterrows and caching the `edge_index` by peptide length provides an ~80% speedup.
**Action:** When converting Pandas rows to graph networks or PyTorch objects, always use `itertuples(index=False, name=None)` with `zip(columns, row)` and cache invariant properties like adjacency matrices (edge_index) based on static dimensions like length.
