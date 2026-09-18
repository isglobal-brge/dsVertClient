\name{dp_cox_grid}
\alias{dp_cox_grid}
\title{Select a signed cross-owner Breslow Cox candidate}
\usage{dp_cox_grid(formula, data, analysis_id, datasources = NULL)}
\arguments{
\item{formula}{An additive \code{Surv(time, event)} formula with owner-qualified
predictors, for example \code{site_a$x + site_b$z}.}
\item{data}{A signed logical dataset name or vertical federation object.}
\item{analysis_id}{The signed finite-grid analysis identifier.}
\item{datasources}{Exactly two custodians, also the compute/noise authorities.}
}
\description{Selects the first canonical minimum of a sticky jointly noised
loss vector from a finite coefficient grid signed by both custodians.}
\details{This family entry is currently namespace-internal and exposed to the
integration session through its single registration function. Public export
requires the package method-inventory and maturity-registry registration.
The production release runner currently fails closed pending connection
to the authenticated fused producer and existing joint-DP vector lifecycle.
This function does not enable a plaintext evaluator or a general nonlinear RPC.
Time and event must have the same owner. Predictors may span both custodians.
Observed binary64 times are clipped to signed bounds and retain exact ordering
and ties; ties use Breslow denominators, including censored observations at the
event time. The execution contract uses one private packed-row permutation,
exact prefix and loss additions on shares, and certified exp/log profile kernels.
Log batches have a fixed public shape and mask event status privately.
Missing inputs exclude the patient privately. There is no intercept,
protected-data optimizer, baseline hazard, or standard error. Existing same-owner
\code{ds.vertCox} routes are unchanged.}
\value{After release integration, a finite-grid point estimate with raw-unit
coefficients and hazard ratios, without sampling inference.}
\keyword{internal}
