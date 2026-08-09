# Epidemiological effect measures from a protected 2-by-2 table

Computes risks, risk difference, risk ratio, odds ratio, attributable
fractions among the exposed and in the population, and number needed to
benefit/harm from an existing aggregate 2-by-2 table. This function is
entirely client-side: it performs no additional DataSHIELD query and
therefore does not enlarge the disclosure surface of the table-producing
method.

## Usage

``` r
ds.vertEpi2x2(
  x,
  exposed = 2L,
  event = 2L,
  level = 0.95,
  zero_correction = c("none", "if_zero", "always")
)
```

## Arguments

- x:

  A 2-by-2 count matrix, or an object such as `ds.vertChisq` containing
  an `observed` matrix.

- exposed:

  Row identifying the exposed group, as a name or index. Defaults to the
  second row.

- event:

  Column identifying the event, as a name or index. Defaults to the
  second column.

- level:

  Confidence level for Wald intervals on an ordinary aggregate, or
  simultaneous DP-mechanism coverage when `x` is a DP-aware
  `ds.vertChisq` result.

- zero_correction:

  One of `"none"` (default), `"if_zero"`, or `"always"`. The latter two
  explicitly request a 0.5 Haldane–Anscombe correction for ratio
  estimates and report that choice in the result.

## Value

For an ordinary aggregate, a `ds.vertEpi2x2` object with standard
epidemiological estimates and Wald intervals. For a DP-aware
`ds.vertChisq` result, a `ds.vertDPEpi2x2` object with simultaneous
mechanism-noise regions and no ordinary sampling interval.

## Details

A DP-aware `ds.vertChisq` object is routed through its validated source
capsule. Ordinary Wald formulae and continuity corrections are not
applied to DP-noised cells, and no additional server call is made.
