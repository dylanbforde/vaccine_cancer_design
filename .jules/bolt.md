## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2024-05-01 - Empty Tensor Shape Mismatch with PyTorch Geometric
**Learning:** PyTorch Geometric strictly expects `edge_index` to be of shape `[2, num_edges]`. When creating a graph for a sequence with no edges (e.g., length 0 or 1), calling `.t()` on an empty tensor list via `torch.tensor([], dtype=torch.long).t()` results in a 1D tensor of shape `[0]`, which crashes the pipeline later.
**Action:** Always explicitly return `torch.empty((2, 0), dtype=torch.long)` when generating an empty `edge_index` for PyTorch Geometric, rather than relying on transposed empty lists.
