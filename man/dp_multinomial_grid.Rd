\name{dp_multinomial_grid}
\alias{dp_multinomial_grid}
\title{Select a signed cross-owner multinomial grid candidate}
\usage{
dp_multinomial_grid(formula, data, analysis_id, signed_contract,
  policy, schema_manifest, datasources = NULL)
}
\arguments{
\item{formula}{An additive formula of qualified owner-column references,
for example \code{site_a$y ~ site_a$x + site_b$z}.}
\item{data}{Signed logical dataset name.}
\item{analysis_id}{Signed analysis identifier.}
\item{signed_contract}{Complete cross-owner grid contract signed by all
pinned custodians.}
\item{policy}{Pinned public policy with two computation/noise authorities,
adjacency, capacity and numeric output precision.}
\item{schema_manifest}{Public schema signed by all pinned custodians.}
\item{datasources}{Authenticated DataSHIELD connections.}
}
\value{A finite-grid candidate with coefficients and no standard errors,
after verifying the authenticated DP release certificate.}
\description{Validates all pinned custodians' signatures and the complete
public numeric and source contract, then uses the shared two-authority DP
release lifecycle. This internal adapter remains pending full family
promotion validation.}
\details{The multinomial domain has 2--8 classes, 1--16 predictors and at most
256 signed candidates. Each nonreference class has an intercept and slopes
with component magnitude at most 8 and L1 norm at most 16. Class predictors
are added across owners before the nonlinear loss. No optimizer accesses
protected data. Candidate ties use signed canonical order.

The signed nonlinearity uses 64 quadratic pieces and 32-bit words with 16
fractional bits. Its uniform loss error is at most 0.0011 per row, plus
\eqn{2^{-g-1}} output rounding if \eqn{g<16} (zero additional rounding
otherwise). Selection targets the certified profile. The exact loss total
differs by at most the admitted row count times the certified row error;
this approximation allowance is included conservatively in the signed caps
and sensitivities. Privacy parameters and noise calibration are unchanged.

Candidate losses are integer coordinates on the signed output lattice
\eqn{S=2^g}, with \eqn{8\leq g\leq18}. The two-authority DP vector release
must be authenticated before candidate post-processing. A signed descriptor
or either state string alone is never release evidence. Existing same-owner
\code{ds.vertMultinom} routes retain their original behavior.}
\seealso{\code{\link{dp_ordinal_grid}}, \code{\link{ds.vertMultinom}}}

\note{This entry point is namespace-internal pending full promotion evidence.}
\keyword{internal}
