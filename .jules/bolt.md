## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2026-02-17 - Itertuples Dict Zipping
**Learning:** `parsed_df.apply(func, axis=1)` is known to be slow, but `parsed_df.to_dict('records')` can also be sub-optimal. Using a list comprehension with `zip(df.columns, row)` and `df.itertuples(index=False, name=None)` gives an incredibly fast row dictionary creation (avoiding NamedTuple overhead), showing a 2-3x speedup in `generate_peptides` generation over `apply(axis=1)`.
**Action:** When a function requires dictionary-like access per row across a large DataFrame, prefer zipping columns with `itertuples(index=False, name=None)` inside a list comprehension.
