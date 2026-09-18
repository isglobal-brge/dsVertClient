\name{dp_ordinal_grid}
\alias{dp_ordinal_grid}
\title{Select a signed cross-owner ordinal cumulative-logit grid candidate}
\usage{
dp_ordinal_grid(formula, data, analysis_id, signed_contract,
  policy, schema_manifest, datasources = NULL)
}
\arguments{
\item{formula}{An additive formula of qualified owner-column references,
for example \code{site_a$y ~ site_a$x + site_b$z}.}
\item{data}{Signed logical dataset name.}
\item{analysis_id}{Signed analysis identifier.}
\item{signed_contract}{Complete cross-owner grid contract signed by both
pinned custodians.}
\item{policy}{Pinned public policy with two computation/noise authorities,
adjacency, capacity and numeric output precision.}
\item{schema_manifest}{Public schema signed by both custodians.}
\item{datasources}{DataSHIELD connections; reserved for release integration.}
}
\value{After authenticated release integration, a finite-grid candidate with
coefficients and thresholds and no standard errors. Currently raises the fixed
transcript-safe unavailable error.}
\description{Validates both custodians' signatures and the complete public
numeric and source contract. The release boundary fails closed before any
DataSHIELD call until Step 2 supplies its fused producer and authenticated
joint-DP release integration.}
\details{The ordinal domain has 2--8 ordered classes, 1--16 predictors and at
most 256 signed candidates. Intercepts are fixed at zero; candidate slopes
have L1 norm at most 8, thresholds have magnitude at most 8 and adjacent
thresholds differ by at least 1/16. These restrictions give every outcome
probability the bound \eqn{p \geq \min(\sigma(-T),\,d e^{-T}/(1+e^{-T})^2)},
where \eqn{T=A+\max|c|\leq16} and \eqn{d\geq1/16} is the minimum threshold
gap. This is strictly greater than \eqn{2^{-28}}, keeping the log loss finite
(the two-class case uses only the endpoint bound). Normalized predictor
coefficients and thresholds are converted
to original predictor units only after DP selection. No optimizer accesses
protected data. Candidate ties use signed canonical order.

The signed log-sigmoid profile uses 64 quadratic pieces and 32-bit words
with 16 fractional bits. Its uniform loss error is at most 0.00012 per row,
plus \eqn{2^{-g-1}} output rounding if \eqn{g<16} (zero additional rounding
otherwise). Selection targets the certified profile. The exact loss total
differs by at most the admitted row count times the certified row error;
this approximation allowance is included conservatively in the signed caps
and sensitivities. Privacy parameters and noise calibration are unchanged.

Candidate losses are integer coordinates on the signed output lattice
\eqn{S=2^g}, with \eqn{8\leq g\leq18}. The two-authority DP vector release
must be authenticated before candidate post-processing. A signed descriptor
or either state string alone is never release evidence. Existing same-owner
\code{ds.vertOrdinal} routes retain their original behavior.}
\seealso{\code{\link{dp_multinomial_grid}}, \code{\link{ds.vertOrdinal}}}

\note{This integration entry point is namespace-internal until the shared public method registry and authenticated release path are wired together.}
\keyword{internal}
