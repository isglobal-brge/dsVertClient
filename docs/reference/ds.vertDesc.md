# Federated descriptive statistics with approximate quantiles

Compute a federated
[`summary()`](https://rdrr.io/r/base/summary.html)-style table of
numeric variables held across vertically partitioned servers. Exact
means, standard deviations, minimums, and maximums are obtained from
each server's local moments; quantiles are interpolated from histogram
bucket counts (user-controlled granularity). No observation-level
quantity is reconstructed at the client.

## Usage

``` r
ds.vertDesc(
  data_name,
  variables = NULL,
  probs = c(0.25, 0.5, 0.75),
  n_buckets = 100L,
  verbose = TRUE,
  datasources = NULL
)
```

## Arguments

- data_name:

  Character. Name of the aligned data frame held on each server.

- variables:

  Optional: a character vector of variables of interest, a named list
  mapping server -\> variables, or `NULL` (default) to auto-detect.

- probs:

  Numeric vector of quantile probabilities to report. Defaults to the
  usual tertiles (0.25, 0.5, 0.75).

- n_buckets:

  Integer. Number of uniform histogram buckets used for quantile
  interpolation. Higher values give tighter quantile resolution at the
  cost of one extra aggregate call per variable; 100 is a reasonable
  default for clinical cohort descriptives.

- verbose:

  Logical. Print per-variable progress when TRUE.

- datasources:

  DataSHIELD connections. If `NULL`, auto-detected.

## Value

A data frame with one row per variable containing columns:

- `server`: which server holds the variable

- `variable`: column name

- `n`: number of non-missing observations

- `n_na`: number of missing observations

- `mean`, `sd`, `min`, `max`: local moments

- one column per requested quantile (named `q25`, `q50`, `q75` etc. by
  default)

## Details

Per-variable flow:

1.  `dsvertLocalMomentsDS` returns the scalar moments in a single
    aggregate call.

2.  Edges for `dsvertHistogramDS` are derived from the moments as
    `seq(min, max, length.out = n_buckets + 1)`; a tiny padding is added
    to guarantee `max` itself falls in the last bucket.

3.  Client-side linear interpolation within the target bucket converts
    the cumulative count to the requested quantile.

The quantile estimate inherits the histogram's resolution; for
`n_buckets = 100` on a roughly normal cohort, the median is typically
accurate within 0.5\\ at linear extra cost.

## Examples

``` r
if (FALSE) { # \dontrun{
conns <- DSI::datashield.login(logindata)
ds.vertDesc("DA", variables = c("age", "bmi", "glu", "bp"),
            probs = c(0.25, 0.5, 0.75, 0.9))
} # }
```
