# Epidemiological measures with simultaneous DP-noise regions

Post-processes one already released 2-by-2 DP histogram. The regions
cover the certified granular-Laplace or approximate-Gaussian mechanism
noise with at least `level` probability by a union bound over all cells.
They are not sampling confidence intervals and do not make an ordinary
chi-square or Fisher p-value valid for noisy counts. Returned
finite-snapshot measures include group and population risks, risk
difference, risk and odds ratios, attributable fractions among the
exposed and in the population, and a direction-typed number needed to
benefit or harm.

## Usage

``` r
ds.vertDPEpi2x2(x, exposed = 2L, event = 2L, level = 0.95)
```

## Arguments

- x:

  A released `ds.vertDPContingency` with a 2-by-2 table.

- exposed, event:

  Exposed row and event column by index or level name.

- level:

  Simultaneous coverage for DP mechanism noise.

## Value

A `ds.vertDPEpi2x2` object. No server call is made.
