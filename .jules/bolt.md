## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2024-05-19 - Pandas Iteration Performance
**Learning:** Pandas `.iterrows()` creates a massive overhead compared to `.itertuples()` for iterating through rows, and pairing `.itertuples(index=False, name=None)` with `dict(zip(df.columns, row))` allows dictionary-like access to rows nearly 10x faster than iterrows without the namedtuple mangling issues.
**Action:** Use `.itertuples(index=False, name=None)` instead of `.iterrows()` when row iteration is strictly necessary in Pandas DataFrames, or use vectorization if possible.

## 2024-05-19 - Torch Graph Data Generation
**Learning:** Calling `create_edge_index` repeatedly to generate standard tensor data in GNN loops (`vaccine_design/pipeline.py` & `model/training.py`) wastes significant CPU cycles. Since the edge index depends strictly on the peptide length, caching it leads to a sizable performance improvement.
**Action:** Cache edge_index tensors by sequence length to avoid redundant, costly CPU operations when building graph data from sequences in standard formats.
