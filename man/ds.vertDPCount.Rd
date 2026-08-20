\name{ds.vertDPCount}
\alias{ds.vertDPCount}
\alias{print.ds.vertDPCount}
\title{Differentially private privacy-unit count}
\usage{
ds.vertDPCount(data_name, server = NULL, datasources = NULL)

\method{print}{ds.vertDPCount}(x, ...)
}
\arguments{
\item{data_name}{Name of the registered protected data frame.}
\item{server}{Optional connected-peer assertion. When supplied it must match
the signed Synopsis finalizer.}
\item{datasources}{Named DataSHIELD connections.}
\item{x}{A \code{ds.vertDPCount} object.}
\item{...}{Additional print arguments.}
}
\value{
A signed Count release with mechanism, accuracy, and per-artifact privacy
metadata.
}
\description{
The signed durable Synopsis contains one bounded \code{admitted_count}
coordinate for its canonical aligned snapshot. Every connected peer validates
that same artifact, while exactly two pinned authorities supply the sticky
mechanism randomness. This function reads and validates that coordinate only:
it neither draws additional noise nor opens a second MPC output. Replaying the
artifact is post-processing. The privacy metadata is per canonical signed
artifact: distinct artifacts compose and no finite global composition claim is
made. Reported accuracy covers mechanism noise, not population sampling
uncertainty.
}
