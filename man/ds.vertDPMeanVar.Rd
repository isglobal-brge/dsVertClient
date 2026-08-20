\name{ds.vertDPMeanVar}
\alias{ds.vertDPMeanVar}
\alias{print.ds.vertDPMeanVar}
\title{Differentially private bounded mean and variance}
\usage{
ds.vertDPMeanVar(data_name, variable, server = NULL, datasources = NULL)

\method{print}{ds.vertDPMeanVar}(x, ...)
}
\arguments{
\item{data_name}{Name of the registered protected data frame.}
\item{variable}{Numeric variable with custodian-owned bounds.}
\item{server}{Optional datasource name.}
\item{datasources}{Named DataSHIELD connections.}
\item{x}{A \code{ds.vertDPMeanVar} object.}
\item{...}{Additional print arguments.}
}
\value{
A validated release containing a noisy effective count, bounded mean,
population central second moment, typed estimability state, and mechanism
metadata bound to the signed canonical Synopsis artifact.
}
\description{
One canonical signed Synopsis artifact clips one contribution per protected
unit to immutable bounds and releases the effective count, normalised quantised
sum, and normalised quantised sum of squares in a common Ring128 vector. Replay
of the same durable artifact is unlimited post-processing and is never denied
by lifetime history. Distinct artifacts compose and no finite global
composition bound is claimed. The client applies only bounded-moment and exact
affine post-processing. Reported regions cover mechanism noise and deterministic
quantisation, not sampling uncertainty.
}
