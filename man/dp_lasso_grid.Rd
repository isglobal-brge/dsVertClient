\name{dp_lasso_grid}
\alias{dp_lasso_grid}
\title{Finite Signed Cross-Owner LASSO Grid}
\usage{
dp_lasso_grid(covariates, signed_contract, base_contract, policy,
  schema_manifest, datasources = NULL)
}
\arguments{
\item{covariates}{Owner-qualified predictor references signed by both custodians.}
\item{signed_contract}{Both custodians' signed LASSO extension contract.}
\item{base_contract}{Authenticated base cross-owner release contract.}
\item{policy}{Public peer pins, capacity, lattice and adjacency policy.}
\item{schema_manifest}{Signed public schema.}
\item{datasources}{DataSHIELD connections for the authenticated base release.}
}
\description{Adds a public slope-only L1 penalty to a signed finite candidate
release. Binomial and Poisson models reuse the certified piecewise-polynomial
loss vector evaluated inside MPC; Gaussian
models reuse the existing cross-owner sufficient statistics and use half RSS.
All objectives use public capacity rather than a protected complete-case count.
No optimization on protected data and no standard errors are provided.}
\details{This namespace-internal integration entry uses the existing authenticated
base release and remains unpromoted pending complete family validation. A contract state or signature alone does
not prove execution. No production plaintext fallback or caller-provided noisy
vector is accepted. The Gaussian adapter must use the existing manifest and
certificate validators. The nonlinear producer must supply its signed profile,
error certificate and approximation-aware caps; the historical Step-1 numeric
profile is not a production fallback. See INTEGRATION_LASSO.md in the companion
server source.}
\value{Selected public candidates, original-scale and normalized paths, the
signed extension, and a transport-anchored base release certificate.}
