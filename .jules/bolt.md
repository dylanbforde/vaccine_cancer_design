## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2024-05-18 - Faster Pandas Row Iteration using List Comprehensions
**Learning:** `apply(axis=1)` is notoriously slow. Iterating through rows via `itertuples(index=False)` combined with list comprehensions instead yields a roughly 6x speedup. `itertuples` returns a NamedTuple, meaning column items must be fetched as attributes (`getattr(row, "field", None)`) instead of via bracket index (`row["field"]`).
**Action:** Replace `apply(axis=1)` with list comprehensions using `itertuples(index=False)` in Pandas operations for row-wise iterations, when vectorization isn't viable.
