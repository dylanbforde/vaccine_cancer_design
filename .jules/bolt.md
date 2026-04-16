## 2024-06-12 - Pandas apply Performance Bottleneck
**Learning:** In data processing pipelines, Pandas `.apply(axis=1)` imposes severe overhead by constructing a new `Series` object for every row.
**Action:** Always replace `.apply(axis=1)` with list comprehensions combined with `zip(df.columns, row)` and `df.itertuples(index=False, name=None)` for significant performance gains when calling arbitrary functions row-by-row.
