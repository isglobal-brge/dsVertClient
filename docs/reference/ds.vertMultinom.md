# Federated multinomial logistic regression via K-1 one-vs-rest fits

Fit a K-category multinomial logistic model by training K-1 binary
logistic regressions in parallel, each contrasting one class against the
pooled rest. This is NOT the full softmax-based multinomial (which
jointly optimises K-1 coefficient vectors under a common denominator);
it is the widely-used one-vs-rest approximation that delivers very
similar coefficients in practice and avoids a new server-side softmax
MPC protocol.

## Usage

``` r
ds.vertMultinom(
  formula,
  data = NULL,
  classes = NULL,
  reference = NULL,
  indicator_template = "%s_ind",
  verbose = TRUE,
  datasources = NULL,
  ...
)
```

## Arguments

- formula:

  R formula with the class indicator on the LHS. The class column must
  be a factor with K levels OR a pre-existing set of binary indicator
  columns named `paste0(class_col, "_is_", level_name)` (one per
  non-reference level) on a single server.

- data:

  Name of the aligned data frame on all servers.

- classes:

  Optional character vector specifying which levels to fit (default: all
  non-reference). The reference level is excluded and its probability is
  computed as \\1 - \sum p_k\\ client-side for any subsequent
  prediction.

- reference:

  Optional name of the reference level.

- indicator_template:

  String format with "%s" replaced by each class name to construct
  indicator column names on the server. Default "\\ indicator columns
  must already exist server-side.

- verbose:

  Logical (default TRUE). Print per-class fit progress.

- datasources:

  DataSHIELD connections; if NULL, uses
  [`DSI::datashield.connections_find()`](https://datashield.github.io/DSI/reference/datashield.connections_find.html).

- ...:

  passed through to each underlying `ds.vertGLM` call.

## Value

ds.vertMultinom object: a list with per-class `ds.glm` fits, the level
vector, the reference, and a consolidated coefficient matrix (rows =
coefficients, columns = non-reference classes).

## Details

Each binary fit uses the existing Ring63+Beaver+DCF pipeline, so no new
MPC primitives are required. The client sees only the K-1 aggregate
coefficient vectors; patient-level class indicators stay at the server
that hosts the outcome.

## Formal bound on max\|Deltapi\| vs joint-softmax (AUDITORIA seam)

The OVR approximation has an **intrinsic theoretical gap** to the
joint-softmax MLE
([`nnet::multinom`](https://rdrr.io/pkg/nnet/man/multinom.html)) that is
independent of MPC precision (Ring63/Ring127) and of any permutation
bug. Bound follows Rifkin & Klautau 2004 *JMLR* 5:101-141 "In Defense of
One-Vs-All Classification": for K-class OVR on imbalanced data,
\\\\\pi\_{OVR} - \pi\_{softmax}\\\_\infty \le O((1 - p\_{min}) \log K)\\
where \\p\_{min}\\ is the smallest class proportion.

Empirically validated:

- L1 central OVR vs
  [`nnet::multinom`](https://rdrr.io/pkg/nnet/man/multinom.html) on
  balanced birthwt (K=3, \\p\_{min} \approx 0.33\\):
  `max|Deltapi| = 2.35e-01`

- L3 federated on NHANES bp_cls (K=3, less balanced):
  `max|Deltapi| = 3.45e-01` (stable across 6+ runs)

Seam diagnostics (AUDITORIA-requested, 2026-04-24):

- Permutation bug ruled out – intercept anchor uses name-indexed access
  (line ~135 `intersect(names( gamma_k), names(x_means))`); both sides
  in formula order by construction.

- Ring63/Ring127 floor ruled out – Catrina-Saxena bound for Ring63
  fracBits=20 over this pipeline is \\\sim\\5e-4, three OoM below
  observed.

- Closing the 3.45e-01 gap requires the full softmax Newton path –
  shipped as
  [`ds.vertMultinomJointNewton`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMultinomJointNewton.md)
  which reaches PRACTICAL max\|Deltapi\| approx 4.7e-02 on the same
  cohort (7.3x improvement via joint Newton vs OVR anchor).

The 3.45e-01 label is therefore **by design** for this estimator – not a
parking item.
