# Differentially private fixed-domain contingency table

Reshapes one purpose-bound, one-contribution-per-unit DP histogram. It
does not apply ordinary chi-square or Fisher inference, whose reference
distributions are not calibrated for DP-noised cells. A unit contributes
once only when all of its valid fixed-domain rows occupy the same cell;
conflicting valid cells contribute zero, while missing and out-of-domain
rows are ignored for consistency. The public custodian-owned
aggregation-policy identifier is returned with every release. Both
variables must currently be on one peer; this is not a cross-peer
vertical association mechanism.

## Usage

``` r
ds.vertDPContingency(
  data_name,
  row_var,
  col_var,
  server = NULL,
  datasources = NULL
)

# S3 method for class 'ds.vertDPContingency'
print(x, ...)
```

## Arguments

- data_name:

  Name of the registered protected data frame.

- row_var, col_var:

  Fixed-domain categorical variables. They may be held by one custodian
  or by two custodians declared in a signed `vertical_cross_specs`
  categorical-pair descriptor.

- server:

  Optional owner assertion. For a cross-owner table it may name either
  declared source owner and never triggers column discovery.

- datasources:

  Named DataSHIELD connections.

- x:

  A `ds.vertDPContingency` object.

- ...:

  Additional print arguments.

## Value

A validated non-negative DP table with marginal and simultaneous
mechanism-noise accuracy metadata.
