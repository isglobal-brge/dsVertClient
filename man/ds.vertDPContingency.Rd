\name{ds.vertDPContingency}
\alias{ds.vertDPContingency}
\alias{print.ds.vertDPContingency}
\title{Differentially private fixed-domain contingency table}
\usage{
ds.vertDPContingency(
  data_name,
  row_var,
  col_var,
  server = NULL,
  datasources = NULL
)

\method{print}{ds.vertDPContingency}(x, ...)
}
\arguments{
\item{data_name}{Name of the registered protected data frame or a version-2
\code{ds.vertFederation}.}
\item{row_var, col_var}{Fixed-domain categorical variables authorized as one
unambiguous signed same- or cross-owner pair. With a federation, bare names
must have one public-schema owner; otherwise use exact
\code{server$column} references.}
\item{server}{Optional same-owner assertion; it never triggers column
discovery.}
\item{datasources}{Named DataSHIELD connections.}
\item{x}{A \code{ds.vertDPContingency} object.}
\item{...}{Additional print arguments.}
}
\value{
A validated non-negative DP table with marginal and simultaneous
mechanism-noise accuracy metadata.
}
\description{
Reshapes one purpose-bound, one-contribution-per-unit DP histogram from a
canonical reduced Synopsis artifact selected from signed same- or cross-owner catalog
metadata. The artifact contains only the authorized pair and its two
marginals, has sticky unlimited replay, and has no lifetime admission gate. It does not
apply ordinary chi-square or Fisher inference, whose reference distributions
are not calibrated for DP-noised cells. A unit contributes once only when all
of its valid fixed-domain rows occupy the same cell; conflicting valid cells
contribute zero, while missing and out-of-domain rows are ignored for
consistency. The public custodian-owned aggregation-policy identifier is
returned with every release. For a cross-owner pair, private exact-GC shares
of the aligned product remain unopened until the single noisy release.
}
