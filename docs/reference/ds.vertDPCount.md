# Differentially private privacy-unit count

Every connected peer signs one canonical current aligned-snapshot
contract. Under add/remove adjacency, exactly two pinned authorities
combine persistent sticky randomness with the protected count inside
exact MPC; the client sees only one bounded release signed by the
finalizer. Under fixed-cohort replace-one adjacency, all K peers sign the
same public cohort size and PSI-run binding, so no MPC session or DP
noise is needed. Privacy metadata is per canonical signed artifact:
distinct artifacts compose and no finite global composition claim is
made. Reported accuracy covers mechanism noise, not population sampling
uncertainty.

## Usage

``` r
ds.vertDPCount(data_name, server = NULL, datasources = NULL)

# S3 method for class 'ds.vertDPCount'
print(x, ...)
```

## Arguments

- data_name:

  Name of the protected data frame.

- server:

  Optional connected-peer assertion. For add/remove Count the signed
  contract selects the finalizer. For fixed-cohort Count this is the
  label attached to the K-consensus public result.

- datasources:

  Named DataSHIELD connections.

- x:

  A `ds.vertDPCount` object.

- ...:

  Additional print arguments.

## Value

A signed Count release with mechanism, accuracy, and per-artifact privacy
metadata.
