# Disclosure-safe descriptive statistics compatibility adapter

Return the historical `ds.vertDesc` data-frame shape from one
custodian-owned, sticky `ds.vertDPDescribe` capsule artifact. Counts,
moments and quantiles are differentially private. Observed extrema and
data-adaptive histogram queries are never requested.

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
  datasources = NULL,
  analysis_id = NULL
)
```

## Arguments

- data_name:

  Character. Name of the aligned data frame held on each server.

- variables:

  Optional character vector, or a named list mapping the artifact-owning
  server to variables, used only to filter the variables already present
  in the signed artifact. `NULL` returns every signed variable. A
  variable absent from the artifact is an error.

- probs:

  Numeric vector of quantile probabilities to report. Defaults to the
  usual quartiles (0.25, 0.5, 0.75).

- n_buckets:

  Legacy compatibility argument. If explicitly supplied, it must equal
  the fixed signed grid length of every selected variable. It never
  changes the server workload.

- range_sd:

  Legacy data-adaptive range argument. Its omitted default is ignored.
  Explicit values are rejected because ranges are fixed by the
  custodian-owned artifact.

- exact_extrema:

  Logical legacy argument. `TRUE` is always rejected: exact observed
  minima and maxima are not disclosure-safe releases.

- open_ended:

  Legacy adaptive-tail argument. Its omitted default is ignored.
  Explicit values are rejected because the signed public bounds, grids
  and invalid bin define the histogram semantics.

- verbose:

  Logical. Print per-variable progress when TRUE.

- datasources:

  DataSHIELD connections. If `NULL`, auto-detected.

- analysis_id:

  Custodian-owned describe specification id. It must name an existing
  artifact in the signed capsule. The client never discovers or chooses
  an id remotely.

## Value

A data frame with one row per variable containing columns:

- `server`: which server holds the variable

- `variable`: column name

- `n`: DP-noisy effective count

- `n_na`: DP-noisy invalid/missing-unit count

- `mean`, `sd`: bounded DP mean and square root of the DP population
  central second moment (not the usual sample SD)

- `min`, `max`: always `NA`

- `range_low`, `range_high`: signed public bounds

- `quantile_status`: fixed-grid DP post-processing status

- one column per requested quantile (named `q25`, `q50`, `q75` etc. by
  default)

## Details

The function makes one call to `ds.vertDPDescribe`. Quantile
probabilities and variable selection are then pure client-side
post-processing of that single immutable release. The returned object
has class `ds.vertDesc` for compatibility and carries `dp_release`,
`dp_descriptives`, and `dp_quantile_bands` attributes with the formal
release and mechanism-uncertainty metadata. These regions exclude
sampling uncertainty. In particular, `n_na` is the noised fixed invalid
bin rather than an exact row-level NA count, and quantiles are
fixed-grid upper endpoints rather than interpolation from
analyst-selected bins.

## Examples

``` r
if (FALSE) { # \dontrun{
conns <- DSI::datashield.login(logindata)
ds.vertDesc("DA", variables = c("age", "bmi", "glu", "bp"),
            probs = c(0.25, 0.5, 0.75, 0.9),
            analysis_id = "baseline_describe_v1")
} # }
```
