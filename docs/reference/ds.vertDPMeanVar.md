# Differentially private bounded mean and variance

Clips finite records to public bounds, averages within privacy unit,
quantises normalised moments, and selects either sequential
granular-Laplace coordinates or one joint TV-accounted
approximate-Gaussian vector. Coordinate radii cover mechanism noise
only. The compact payload does not expose the noisy sums needed for a
client-verifiable propagated ratio interval and provides no sampling
confidence interval.

## Usage

``` r
ds.vertDPMeanVar(data_name, variable, server = NULL, datasources = NULL)

# S3 method for class 'ds.vertDPMeanVar'
print(x, ...)
```

## Arguments

- data_name:

  Name of the registered protected data frame.

- variable:

  Numeric variable with custodian-owned bounds.

- server:

  Optional datasource name.

- datasources:

  Named DataSHIELD connections.

- x:

  A `ds.vertDPMeanVar` object.

- ...:

  Additional print arguments.

## Value

A validated release containing a noisy effective count, bounded mean,
population central second moment, typed estimability state, and
mechanism metadata.
