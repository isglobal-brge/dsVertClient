\name{dp_grouped_grid_cross}
\alias{dp_grouped_grid_cross}
\alias{dp_lmm_grid}
\alias{dp_glmm_grid}
\alias{dp_gee_grid}
\title{Select a Signed Cross-owner Grouped Finite-grid Candidate}
\usage{
dp_lmm_grid(outcome, predictors, contract, policy, schema_manifest,
  datasources = NULL)
dp_glmm_grid(outcome, predictors, contract, policy, schema_manifest,
  datasources = NULL)
dp_gee_grid(outcome, predictors, contract, policy, schema_manifest,
  datasources = NULL)
}
\arguments{
\item{outcome}{One owner-qualified outcome, such as \code{site_a$y}.}
\item{predictors}{Owner-qualified covariates in signed canonical order.}
\item{contract}{The complete grouped-grid contract signed by both authorities.}
\item{policy}{Public custodian policy with the two pinned compute-and-noise authorities.}
\item{schema_manifest}{Doubly signed aligned dataset schema.}
\item{datasources}{Optional DataSHIELD connections for the integrated reader.}
}
\description{Validates a finite-grid grouped contract and requests its sticky
joint-DP artifact. The production reader currently fails closed pending fused
producer integration. Existing same-owner adapters are unchanged.}
\details{LMM accepts a signed fixed covariance or an explicit finite ML variance grid,
ordered by variance then coefficient candidate. The ML objective includes the
private log determinant and a candidate-independent count shift; it does not
provide REML.
GLMM uses fixed nonadaptive GH5 quadrature with signed variance. Its target is a
finite-quadrature surrogate, and numerical circuit error is separate from its
Gaussian-integral approximation error. No protected optimisation is performed.
GEE retains likelihood, bread, and componentwise clipped cluster-score meat;
selection uses likelihood. Independence, exchangeable, and AR(1) correlation
parameters are public and signed. The privacy unit is one admitted patient with
one analysis row. Sensitivity includes the change of the entire affected cluster.
Both custodians must sign the schema and the complete candidate contract. An
unavailable producer, invalid signature, stale evidence, or malformed descriptor
fails with the standard transcript-safe error. No plaintext fallback exists.}
\value{After integration, a list containing the selected candidate, transformed
fixed-effect coefficients, signed parameters and release provenance. No standard
errors or p-values.}
