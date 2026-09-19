\name{dp_nb_grid}
\alias{dp_nb_grid}
\title{Validate a signed cross-owner negative-binomial finite grid}
\usage{
dp_nb_grid(formula, data_name, analysis_id, signed_contract,
  policy, schema_manifest, datasources = NULL)
}
\arguments{
\item{formula}{An intercept-bearing formula with owner-qualified outcome and predictors, for example \code{peer_a$y ~ peer_a$x + peer_b$z}.}
\item{data_name}{The signed, prealigned logical dataset identifier.}
\item{analysis_id}{The signed analysis identifier.}
\item{signed_contract}{The complete NB2 cross-owner contract signed by all pinned custodians.}
\item{policy}{Public pinned peer, adjacency, capacity and precision policy.}
\item{schema_manifest}{The matching schema with every pinned custodian's signature.}
\item{datasources}{Reserved DataSHIELD connections for the authenticated release adapter.}
}
\description{
Admits the versioned \code{nb_grid_cross_v1} public contract for a finite paired
beta/theta grid whose covariates may span both custodians. The correct NB2
negative log-PMF is evaluated on the shared linear predictor by a separately
registered secure producer using a certified 32-bit piecewise-quadratic
softplus profile. The signed row caps include approximation error. Theta is signed and restricted to powers of two
from one eighth through 128. The complete contract, bounds, predictor order,
PSI binding, numeric certificate and both signatures are checked locally.
}
\details{
This isolated family extension does not yet connect the fused secure producer
and authenticated sticky joint DP release. The entry therefore fails closed
before any server request, including for a valid contract. It provides no
plaintext fallback or caller-supplied evaluation callback. Integration must
check bilateral producer/release evidence before selecting a candidate from
the DP vector. The same-owner negative-binomial routes remain independent.
No protected-data optimizer, standard errors, covariance or p-values are used.
}
\value{Until authenticated release wiring is installed, a fixed transcript-safe
\code{dsvert_dp_public_failure} error.}
