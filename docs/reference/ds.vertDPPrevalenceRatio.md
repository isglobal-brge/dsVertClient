# Cross-sectional prevalence effects from one DP 2-by-2 release

Provides an explicitly cross-sectional naming view of
[`ds.vertDPEpi2x2()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPEpi2x2.md).
It performs no new statistical calculation: prevalence difference,
prevalence ratio, prevalence odds ratio, attributable prevalence
fractions, and the reciprocal prevalence-difference summary are aliases
of the corresponding risk-named values and regions in the same validated
result.

## Usage

``` r
ds.vertDPPrevalenceRatio(x, exposed, prevalent, level = 0.95)

ds.vertDPPrevalenceRatioInference(
  x,
  exposed,
  prevalent,
  level = 0.95,
  mechanism_alpha_share = 0.5
)
```

## Arguments

- x:

  A released `ds.vertDPContingency` with a 2-by-2 table.

- exposed:

  Required row name or index for the exposed group. There is no default
  because orientation is a scientific declaration.

- prevalent:

  Required column name or index for prevalent disease or status. There
  is no default because orientation is a scientific declaration.

- level:

  Simultaneous coverage for DP mechanism noise, passed unchanged to
  [`ds.vertDPEpi2x2()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPEpi2x2.md).

- mechanism_alpha_share:

  Fraction of total non-coverage allocated to the DP mechanism event,
  passed unchanged to
  [`ds.vertDPEpi2x2Inference()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPEpi2x2Inference.md).

## Value

A `ds.vertDPPrevalenceRatio` object that retains every field of the
underlying `ds.vertDPEpi2x2` result and adds prevalence-named aliases.
The aliases are numerically identical to the source fields. No server
call or additional privacy release is made.

`ds.vertDPPrevalenceRatioInference()` returns a
`ds.vertDPPrevalenceRatioInference` object retaining every field and
coverage statement of the underlying `ds.vertDPEpi2x2Inference` result,
plus prevalence-named aliases. It makes no server call and consumes no
additional privacy.

## Details

Calling this function declares a cross-sectional interpretation; the
design is not and cannot be inferred from the released table. The method
adds no causal, temporal, or population-transportability claim. Its
mechanism regions exclude population-sampling uncertainty; use
`ds.vertDPPrevalenceRatioInference()` when its iid cross-sectional
sampling model is scientifically justified.
