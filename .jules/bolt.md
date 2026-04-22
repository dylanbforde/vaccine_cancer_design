## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).
## 2024-05-24 - Pandas .iterrows() Overhead in ML Pipelines
**Learning:** PyTorch/Graph ML pipelines processing thousands of rows (e.g., peptide sequences) suffer massive overhead from Pandas `.iterrows()`. Benchmarking showed `.iterrows()` takes 9.6x longer than `.itertuples(index=False, name=None)` combined with `dict(zip(columns, row))` when reconstructing dictionaries for downstream parsing.
**Action:** Always replace `.iterrows()` in high-throughput data preparation loops with `itertuples(index=False, name=None)` to bypass Pandas Series boxing overhead.

## 2024-05-24 - PyTorch Graph Edge Calculation Redundancy
**Learning:** Generating static node topology (e.g., adjacent edges for a sequence of length N) repeatedly within an iteration loop causes significant CPU bottlenecking during data processing. For `PeptideEncoder.create_edge_index`, calculating identical edges for 100,000 sequences took 2.24s un-cached vs 0.008s cached.
**Action:** Use `@staticmethod` combined with `@functools.lru_cache(maxsize=32)` on topology generation functions that solely depend on deterministic inputs like sequence length.
