# Differentially Private Fixed-Grid Descriptive Statistics

Public descriptive component of the package's unified disclosure-safe
interface. One sticky server release supplies every requested statistic;
quantile probabilities are client-side post-processing.

## Usage

``` r
ds.vertDPDescribe(data_name, analysis_id,
  probs = c(0.25, 0.5, 0.75), server = NULL,
  datasources = NULL)

# S3 method for class 'ds.vertDPDescribe'
print(x, ...)
```

## Arguments

- data_name:

  Name of the registered protected data frame.

- analysis_id:

  Custodian-owned describe specification id.

- probs:

  Public quantile probabilities strictly between zero and one. They only
  alter local post-processing and never create a different server query.

- server:

  Optional datasource name. If omitted, the lexicographically first
  datasource is selected deterministically.

- datasources:

  DataSHIELD connections.

- x:

  A `ds.vertDPDescribe` object.

- ...:

  Additional print arguments (currently unused).

## Value

A `ds.vertDPDescribe` object containing DP-noisy effective counts,
bounded mean, population variance and standard deviation, fixed
histograms, and histogram-CDF quantiles. The object also includes
conservative simultaneous mean/variance regions propagated from
mechanism noise and deterministic quantisation, the effective epsilon of
each mechanism family, simultaneous histogram-noise radii, public bounds
and grids, and typed estimability status.

## Details

The method never requests or reconstructs observed minima or maxima and
never uses data-adaptive bins. Mean and variance are projected onto the
geometry of the public bounds, so returned point estimates are
internally consistent. When the DP-noisy count is positive but its
simultaneous lower bound is zero, the point remains available with
status `dp_point_available_count_not_certified_positive` rather than
being silently presented as fully certified.

The release validates both add/remove and fixed-cohort replace-one
adjacency. In the latter case, histogram sensitivity accounts for moving
a unit between two bins.

Quantiles are upper endpoints of the fixed public histogram bins. Their
bands are obtained from a simultaneous box around every DP histogram
coordinate and include the corresponding public grid interval. They
cover DP mechanism noise under the release certificate. They are not
sampling confidence intervals, do not provide standard errors or
p-values, and should not be interpreted as covering sampling or model
uncertainty.

Mean and variance regions propagate the same simultaneous coordinate box
for count, quantised sum and quantised sum of squares, together with the
declared per-unit quantisation bound and bounded-moment geometry. They
are conditional on a positive effective count and likewise exclude
sampling uncertainty.

The client uses a fail-closed DSI aggregate contract: partial or null
results and transport callbacks are rejected, remote error text is not
propagated, and a failed release is not retried.
