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
\item{server}{Optional connected-peer assertion. For add/remove Count the
signed contract selects the finalizer. For fixed-cohort Count this is the
label attached to the K-consensus public result.}
\item{datasources}{Named DataSHIELD connections.}
\item{x}{A \code{ds.vertDPCount} object.}
\item{...}{Additional print arguments.}
}
\value{
A signed Count release with mechanism, accuracy, and per-artifact privacy
metadata.
}
\description{
Every connected peer signs one canonical current aligned-snapshot contract.
Under add/remove adjacency, exactly two pinned authorities combine persistent
sticky randomness with the protected count inside exact MPC; the client sees
only one bounded release signed by the finalizer. Under fixed-cohort
replace-one adjacency, all K peers sign the same public cohort size and PSI-run
binding, so no MPC session or DP noise is needed. The privacy metadata is per
canonical signed artifact: distinct artifacts compose and no finite global
composition claim is made. Reported accuracy covers mechanism noise, not
population sampling uncertainty.
}
