## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).
## 2025-02-28 - PyTorch Edge Index Caching Over Iterrows
**Learning:** In PyTorch Geometric data processing within tight loops (like `create_training_data` and `process_mutations`), `df.iterrows()` combined with repeated `create_edge_index()` calls for fixed-length peptide sequences causes severe bottlenecks. `df.itertuples(index=False, name=None)` with `dict(zip(cols, row))` reduces pandas overhead significantly, while caching the `edge_index` tensor by `len(peptide)` avoids redundant graph structure computation and tensor allocations.
**Action:** Always replace `iterrows()` with `itertuples()` for row iteration, and cache repetitive PyTorch tensor instantiations (like `edge_index` based on peptide length) across iterations.
