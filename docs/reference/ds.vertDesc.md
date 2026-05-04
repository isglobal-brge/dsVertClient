# Federated descriptive statistics with approximate quantiles

Compute a federated
[`summary()`](https://rdrr.io/r/base/summary.html)-style table of
numeric variables held across vertically partitioned servers. Exact
means and standard deviations are obtained from each server's local
moments; exact extrema are suppressed by default because they can
disclose outliers. Quantiles are interpolated from disclosure-checked
histogram bucket counts. No observation-level quantity is reconstructed
at the client.

## Usage

``` r
ds.vertDesc(
  data_name,
  variables = NULL,
  probs = c(0.25, 0.5, 0.75),
  n_buckets = 100L,
  range_sd = getOption("dsvert.desc_range_sd", 4),
  exact_extrema = getOption("dsvert.allow_exact_extrema", FALSE),
  open_ended = getOption("dsvert.desc_open_ended", TRUE),
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
  cost of one extra aggregate call per variable; when
  `datashield.privacyLevel > 0`, the effective bucket count is capped by
  the cohort size and privacy threshold.

- range_sd:

  Numeric. When exact extrema are not released, histogram edges are
  built over `mean +/- range_sd * sd`. Values outside that range are
  retained as under/overflow aggregate counts.

- exact_extrema:

  Logical. Return exact min/max and use them for the histogram range.
  Defaults to `getOption("dsvert.allow_exact_extrema", FALSE)`.

- open_ended:

  Logical. When exact extrema are suppressed, use open-ended first/last
  histogram buckets so outliers are absorbed into coarse tail aggregates
  rather than disclosed as separate under/overflow counts.

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

- `mean`, `sd`: local moments

- `min`, `max`: exact extrema only when `exact_extrema = TRUE`;
  otherwise `NA`

- `range_low`, `range_high`: histogram range used for approximate
  quantiles

- `quantile_status`: `"ran"` or the suppression reason

- one column per requested quantile (named `q25`, `q50`, `q75` etc. by
  default)

## Details

Per-variable flow:

1.  `dsvertLocalMomentsDS` returns scalar moments in a single aggregate
    call; exact extrema are omitted unless explicitly enabled.

2.  Edges for `dsvertHistogramDS` are derived from either exact extrema
    (opt-in) or a moment-bounded range `mean +/- range_sd * sd`.

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
