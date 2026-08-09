# Differentially private one-way frequency distribution

Extracts one signed categorical-marginal block from the single reusable
joint-DP capsule. Each admitted patient contributes to at most one
public level after the custodian-defined repeated-record collapse. The
returned count and proportion regions cover mechanism noise only; no
exact marginal, missingness count, or patient-level value is requested.

## Usage

``` r
ds.vertDPFrequency(data_name, variable, server = NULL, datasources = NULL)
```

## Arguments

- data_name:

  Signed protected dataset name.

- variable:

  Signed fixed-domain categorical variable.

- server:

  Optional expected owner server.

- datasources:

  DataSHIELD connections.

## Value

A `ds.vertDPFrequency` object.
