# DP common Mantel-Haenszel odds ratio from one fixed capsule table

Purely post-processes one validated `ds.vertDPContingency` whose rows
are public strata and whose four columns are publicly mapped 2-by-2
cells. It computes the finite-snapshot common Mantel-Haenszel odds ratio
and a conservative simultaneous mechanism-noise region by interval
arithmetic. It never computes a classical CMH p-value or sampling
confidence interval. The capsule must use add/remove adjacency with the
sole `consistent_joint_cell_else_exclude_v1` rule: each admitted unit
contributes to exactly one global stratum-by-cell coordinate, giving
block L1 sensitivity one.

## Usage

``` r
ds.vertDPMantelHaenszel(x, cell_map = NULL, level = 0.95)
```

## Arguments

- x:

  One released `ds.vertDPContingency` with K rows and four columns.

- cell_map:

  A public canonical cell mapping as described for
  [`ds.vertMantelHaenszel()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMantelHaenszel.md);
  optional only for canonical column names.

- level:

  Simultaneous coverage for the DP mechanism-noise cell box.

## Value

A `ds.vertDPMantelHaenszel` object. The estimate is explicitly typed as
finite, zero, infinite, or non-estimable. The function performs zero
server calls and consumes zero additional privacy budget.
