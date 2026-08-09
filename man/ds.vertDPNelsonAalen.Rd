\name{ds.vertDPNelsonAalen}
\alias{ds.vertDPNelsonAalen}
\title{Nelson--Aalen Curve from One DP Survival Release}
\usage{
ds.vertDPNelsonAalen(x)
}
\arguments{
  \item{x}{A released \code{ds.vertDPSurvival} object.}
}
\value{
A data frame with public grid time and the DP-histogram-derived cumulative
hazard. No sampling confidence interval is claimed.
}
\description{
Pure post-processing; this function makes no server request and consumes no
additional privacy allocation.
}
