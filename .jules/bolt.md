## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2026-02-17 - Lru_cache on Instance Methods
**Learning:** Applying `@lru_cache` directly to an instance method causes `self` to be included in the cache key. This bypasses the cache across different instances of the class and can lead to memory leaks by holding references to those instances.
**Action:** When caching utility-like class methods (e.g., generating sequence graphs), convert them to `@staticmethod` first so the cache is shared globally and instances can be garbage collected.
