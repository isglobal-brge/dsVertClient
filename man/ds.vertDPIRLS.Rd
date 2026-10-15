\name{ds.vertDPIRLS}
\alias{ds.vertDPIRLS}
\title{Execute a signed differential-private IRLS path}
\usage{ds.vertDPIRLS(analysis_id, datasources = NULL, verbose = TRUE)}
\arguments{
\item{analysis_id}{Custodian analysis selector.}
\item{datasources}{Named DataSHIELD connections.}
\item{verbose}{Emit fixed public progress messages.}}
\description{Selects a custodian-provisioned analysis with every statistical
choice in its unanimously signed contract. Configure the independently verified
named peer pinset in the session option \code{dsvert.dp_irls.trusted_pins} before
use. No beta, predecessor, budget, bounds, family or seed is accepted.}
\value{Verified public chain, coefficients and signed-translation releases.
Sampling inference and standard errors are unavailable pending B4.}
