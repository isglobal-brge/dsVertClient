# Differentially private privacy-unit count

Under add/remove adjacency, the two custodian-designated pinned peers
create one sticky discrete-Laplace draw and one bounded clamp inside
exact GC; the client accepts only byte-identical releases signed by both
peers. Under fixed-cohort replace-one adjacency, the count is the
zero-cost public policy constant unanimously returned by all connected
peers and checked against the consensus policy capacity. Capsule reuse
has no request quota, remaining budget, or history-dependent accuracy
change. Reported accuracy covers mechanism noise, not population
sampling uncertainty.

## Usage

``` r
ds.vertDPCount(data_name, server = NULL, datasources = NULL)

# S3 method for class 'ds.vertDPCount'
print(x, ...)
```

## Arguments

- data_name:

  Name of the registered protected data frame.

- server:

  Optional returned-source label for the fixed-cohort public-policy
  path. Every connected peer must return the same value and it must
  equal the consensus policy capacity. On the joint add/remove path it
  is only a compatibility assertion.

- datasources:

  Named DataSHIELD connections.

- x:

  A `ds.vertDPCount` object.

- ...:

  Additional print arguments.

## Value

A validated DP count release with mechanism, accuracy, sticky-noise, and
composition metadata.
