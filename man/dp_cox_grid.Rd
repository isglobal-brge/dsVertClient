\name{dp_cox_grid}
\alias{dp_cox_grid}
\title{Select a signed cross-owner Breslow Cox candidate}
\usage{dp_cox_grid(formula, data, analysis_id, datasources = NULL)}
\arguments{
\item{formula}{An additive \code{Surv(time, event)} formula with owner-qualified
predictors, for example \code{site_a$x + site_b$z}.}
\item{data}{A signed logical dataset name or vertical federation object.}
\item{analysis_id}{The signed finite-grid analysis identifier.}
\item{datasources}{Two to five signed source owners, including the two
compute/noise authorities.}
}
\description{Selects the first canonical minimum of a sticky jointly noised
loss vector from a finite coefficient grid signed by all source owners.}
\details{The release uses the authenticated staged Synopsis lifecycle and
requires signed evidence from both pinned computation authorities. Durable
cold reads verify the public compilation and publication without rereading sources.
This function does not enable a plaintext evaluator or a general nonlinear RPC.
Time and event must have the same owner. Predictors may span all source owners.
Observed binary64 times are clipped to signed bounds and retain exact ordering
and ties; ties use Breslow denominators, including censored observations at the
event time. The execution contract uses one private packed-row permutation,
exact prefix and loss additions on shares, and certified exp/log profile kernels.
Log batches have a fixed public shape and mask event status privately.
All source owners sign the public admission envelope. The integrated release
is scoped to at most 400 rows and 50 candidates.
The client and server reject requests exceeding that envelope before source
access; signatures cannot override the implementation limits.
Missing inputs exclude the patient privately. There is no intercept,
protected-data optimizer, baseline hazard, or standard error. Existing same-owner
\code{ds.vertCox} routes are unchanged.}
\value{A finite-grid point estimate with raw-unit
coefficients and hazard ratios, without sampling inference.}
