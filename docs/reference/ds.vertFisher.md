# Federated Fisher exact test (same-server case)

Exact conditional p-value for a 2x2 (or small 2xK) contingency table of
two variables held at the same server. Uses the
[`stats::fisher.test`](https://rdrr.io/r/stats/fisher.test.html)
implementation on the released count matrix (itself already subject to
datashield.privacyLevel small-cell suppression), so this is a safe
exact-inference complement to ds.vertChisq when any expected cell count
is below the chi-square validity floor.

## Usage

``` r
ds.vertFisher(
  data_name,
  var1,
  var2,
  server = NULL,
  alternative = c("two.sided", "greater", "less"),
  conf.int = TRUE,
  conf.level = 0.95,
  datasources = NULL
)
```

## Arguments

- data_name, var1, var2, server, datasources:

  Same semantics as in `ds.vertChisq`.

- alternative:

  One of "two.sided" (default), "greater", "less"; only consulted for
  2x2 tables.

- conf.int:

  Logical. Return an odds-ratio confidence interval for 2x2.

- conf.level:

  Confidence level for 2x2 OR interval.

## Value

An object of class `ds.vertFisher` (a `htest`-compatible list plus table
metadata).
