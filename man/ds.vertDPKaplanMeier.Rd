\name{ds.vertDPKaplanMeier}
\alias{ds.vertDPKaplanMeier}
\title{Kaplan--Meier Curve from One DP Survival Release}
\usage{
ds.vertDPKaplanMeier(x)
}
\arguments{
  \item{x}{A released \code{ds.vertDPSurvival} object.}
}
\value{
A data frame with public grid time and the DP-histogram-derived survival
estimate. Its \code{uncertainty_scope} attribute states that only DP mechanism
noise is covered.
}
\description{
Pure post-processing; this function makes no server request and consumes no
additional privacy allocation.
}
