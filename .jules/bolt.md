## 2025-04-25 - Pandas DataFrame Iteration Bottleneck
**Learning:** Iterating over Pandas DataFrames using `.iterrows()` introduces severe performance overhead because it instantiates a new `pd.Series` object for every row.
**Action:** Always replace `.iterrows()` with `itertuples(index=False, name=None)` paired with `dict(zip(cols, row_tuple))` for row-wise iteration when dictionary-like access is needed. This maintains identical syntax for the loop body while speeding up iteration by 10x-20x.
