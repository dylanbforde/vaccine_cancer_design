## 2026-02-17 - Pandas Iteration vs Vectorization
**Learning:** `iterrows` in pandas is significantly slower than `apply` or vectorized operations, especially for large datasets (130k rows). Benchmarking showed ~20x speedup.
**Action:** Always prefer `apply` or vectorization when processing DataFrame columns. Avoid `iterrows` loop.

## 2026-02-17 - Pre-compiled Regex Overhead
**Learning:** Defining a list of regex patterns inside a frequently called function (`parse_protein_change`) causes significant overhead due to repeated list construction and regex compilation.
**Action:** Move constant regex patterns to module scope and pre-compile them using `re.compile`. This yielded a ~35% performance improvement (0.28s -> 0.18s for 90k calls).

## 2026-02-18 - JSON Deserialization and Cache Invalidations
**Learning:** `json.load()` converts tuples to lists. Checking `isinstance(..., tuple)` on JSON-loaded data caused a 100% cache miss rate in `get_sequences`, forcing repeated network calls.
**Action:** When validating cache entries loaded from JSON, always check `isinstance(val, (list, tuple))`.

## 2026-02-18 - Repeated File I/O in Loops
**Learning:** Calling `json.load()` inside a frequently called function (`get_sequences` called in a batch loop) caused significant I/O overhead (2.5s for 10 calls on 50k items).
**Action:** Use a module-level dictionary (`_MEMORY_CACHES`) to store loaded file content in memory, preventing repeated disk reads. This yielded an ~8.5x speedup.
