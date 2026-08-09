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
metadata.
}
\description{
Clips finite records to public bounds, averages within privacy unit, quantises
normalised moments, and selects either sequential granular-Laplace coordinates
or one joint TV-accounted approximate-Gaussian vector. Coordinate radii cover
mechanism noise only. The compact payload does not expose the noisy sums needed
for a client-verifiable propagated ratio interval and provides no sampling
confidence interval.
}
