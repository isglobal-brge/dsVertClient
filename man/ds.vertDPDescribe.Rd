\name{ds.vertDPDescribe}
\alias{ds.vertDPDescribe}
\alias{print.ds.vertDPDescribe}
\title{Differentially Private Fixed-Grid Descriptive Statistics}
\usage{
ds.vertDPDescribe(data_name, analysis_id,
  probs = c(0.25, 0.5, 0.75), server = NULL,
  datasources = NULL, resume = NULL)

\method{print}{ds.vertDPDescribe}(x, ...)
}
\arguments{
  \item{data_name}{Name of the registered protected data frame.}
  \item{analysis_id}{Custodian-owned describe specification id.}
  \item{probs}{Public quantile probabilities strictly between zero and one.
  They only alter local post-processing and never create a different server
  query.}
  \item{server}{Optional datasource name. If omitted, the
  lexicographically first datasource is selected deterministically.}
  \item{datasources}{DataSHIELD connections.}
  \item{resume}{Optional portable resume token returned in the \code{resume}
  field of an earlier \code{ds.vertDPDescribe()} result, or that earlier result
  itself. The token is authenticated by the pinned peers' signed public
  bootstrap and manifest evidence and contains no connection handles.}
  \item{x}{A \code{ds.vertDPDescribe} object.}
  \item{...}{Additional print arguments (currently unused).}
}
\value{
A \code{ds.vertDPDescribe} object containing DP-noisy effective counts,
bounded mean, population variance and standard deviation, fixed histograms,
and histogram-CDF quantiles. The object also includes conservative simultaneous
mean/variance regions propagated from mechanism noise and deterministic
quantisation, the effective epsilon of each mechanism family, simultaneous
histogram-noise radii, public bounds and grids, and typed estimability status.
Synopsis results also include the portable \code{resume} token and the logical
\code{cleanup_pending} state.
}
\description{
Public descriptive component of the package's unified disclosure-safe
interface. One sticky server release supplies every requested statistic;
quantile probabilities are client-side post-processing.
}
\details{
The method never requests or reconstructs observed minima or maxima and never
uses data-adaptive bins. Mean and variance are projected onto the geometry of
the public bounds, so returned point estimates are internally consistent.
When the DP-noisy count is positive but its simultaneous lower bound is zero,
the point remains available with status
\code{dp_point_available_count_not_certified_positive} rather than being
silently presented as fully certified.

The release validates both add/remove and fixed-cohort replace-one adjacency.
In the latter case, histogram sensitivity accounts for moving a unit between
two bins.

Quantiles are upper endpoints of the fixed public histogram bins. Their bands
are obtained from a simultaneous box around every DP histogram coordinate and
include the corresponding public grid interval. They cover DP mechanism noise
under the release certificate. They are not sampling confidence intervals,
do not provide standard errors or p-values, and should not be interpreted as
covering sampling or model uncertainty.

Mean and variance regions propagate the same simultaneous coordinate box for
count, quantised sum and quantised sum of squares, together with the declared
per-unit quantisation bound and bounded-moment geometry. They are conditional
on a positive effective count and likewise exclude sampling uncertainty.

The Synopsis route has no request, rate, catalogue-count or lifetime budget
that can deny a future analysis. Its privacy claim is per canonical signed
artifact. Distinct artifacts compose, and no finite global composition bound is
claimed. A resume call rebinds the saved signed bootstrap to the live peer
identities and replays the durable publication without rerunning Claim,
Compile, source transfer or session execution.

This release currently supports same-owner Describe artifacts. A cold exact-GC
Synopsis execution fails before Claim or Compile; an already published exact
artifact may still be replayed through the publication fast path.

The client uses a fail-closed DSI aggregate contract: partial or null results
and transport callbacks are rejected, remote error text is not propagated, and
a failed release is not retried.
}
