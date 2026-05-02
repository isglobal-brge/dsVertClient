# Pure helper: compute chi-square statistics from observed counts + margins. Exposed as an internal function so unit tests can exercise the math path without a DataSHIELD round trip.

Pure helper: compute chi-square statistics from observed counts +
margins. Exposed as an internal function so unit tests can exercise the
math path without a DataSHIELD round trip.

## Usage

``` r
.dsvert_chisq_compute(observed, row_margins, col_margins, n, correct = TRUE)
```

## Arguments

- observed:

  Integer matrix of cell counts.

- row_margins:

  Integer vector of row sums.

- col_margins:

  Integer vector of column sums.

- n:

  Total observation count.

- correct:

  Apply Yates continuity correction for 2x2 tables.

## Value

list with statistic, df, p_value, observed, expected, residuals, n,
correct.
