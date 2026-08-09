\name{ds.vertDPCumulativeIncidence}
\alias{ds.vertDPCumulativeIncidence}
\title{Competing-Risks Cumulative Incidence from One DP Release}
\usage{
ds.vertDPCumulativeIncidence(x, cause = NULL)
}
\arguments{
  \item{x}{A released \code{ds.vertDPSurvival} object.}
  \item{cause}{Optional public event-cause label. \code{NULL} returns every
  released cause.}
}
\value{
A long data frame with public grid time, cause, and the Aalen--Johansen
cumulative-incidence estimate derived from the DP histogram.
}
\description{
Pure post-processing of the same histogram used for Kaplan--Meier and
Nelson--Aalen. It makes no server request and consumes no additional privacy
allocation.
}
