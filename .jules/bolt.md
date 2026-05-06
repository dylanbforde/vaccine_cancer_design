## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2024-05-06 - Replacing iterrows with itertuples in Pandas
**Learning:** Using `mutations_df.iterrows()` in `VaccineDesignPipeline.process_mutations` resulted in significantly slow iteration over rows (~0.39s for 1000 rows). By replacing it with `mutations_df.itertuples(index=False, name=None)` and dictionary generation using `dict(zip(columns, row_tuple))`, we achieved an 8.7x performance improvement (~0.045s for 1000 rows).
**Action:** When performing row-by-row iteration over Pandas DataFrames where vectorization isn't easily possible (like creating complex nested object graphs like `Data(x=..., edge_index=...)`), strongly prefer `itertuples()` + `zip(columns, ...)` over `.iterrows()` for significantly reduced overhead while preserving `dict` mapping usability.
