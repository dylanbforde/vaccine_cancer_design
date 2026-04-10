## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).
## 2026-04-10 - Pandas Iteration vs Vectorization
**Learning:** `.apply(axis=1)` in Pandas creates a new Series object for every row, causing severe overhead. Benchmarks showed it taking >1s for 100k rows compared to ~0.5s for list comprehensions using `dict(zip())` with `.itertuples(index=False, name=None)`.
**Action:** Always use vectorized operations if possible. If row-by-row iteration is necessary and returns non-scalar or complex data, use list comprehensions over `.itertuples()` paired with `dict(zip())` to simulate named access instead of `.apply(axis=1)`.

## 2026-04-10 - Pandas Vectorized String Operations
**Learning:** Using `.apply(lambda x: ...)` for simple string checks (like length or regex matching) in a Pandas Series is very slow. Native Pandas vectorized string methods like `.str.match()` and `.str.len()` execute at C-speed and are significantly faster (e.g., 0.029s vs 0.076s for length checks).
**Action:** Always use Pandas `.str` accessors (`.str.contains`, `.str.match`, `.str.len`, etc.) combined with bitwise operators (`&`, `|`, `~`) for boolean masks instead of `.apply()` when validating or manipulating strings in a DataFrame or Series.
