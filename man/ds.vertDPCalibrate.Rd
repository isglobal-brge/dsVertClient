\name{ds.vertDPCalibrate}
\alias{ds.vertDPCalibrate}
\title{Calibrate deployed DP mechanism candidates}
\usage{
ds.vertDPCalibrate(
  capsule_epsilon = c(1, 3),
  peer_count = 2L,
  sensitivity = 1,
  confidence = 0.95,
  capsule_delta = .DSVERT_DP_DEFAULT_CAPSULE_DELTA,
  coordinate_count = 1L,
  gaussian_l2_sensitivity = sensitivity,
  objective = "auto"
)
}
\arguments{
\item{capsule_epsilon}{Positive fixed consortium epsilon values assigned
once to immutable capsules. Reusing a capsule does not change epsilon.}
\item{peer_count}{Positive number of pinned vertical peers. This is planning
metadata only and never divides capsule epsilon.}
\item{sensitivity}{Positive integer per-coordinate L1 sensitivity.}
\item{confidence}{Central probability for the separately reported planning
radius. Candidate selection always uses its declared 95-percent objective.}
\item{capsule_delta}{Fixed consortium delta assigned to each immutable
capsule. Zero disables Gaussian.}
\item{coordinate_count}{Public vector dimension.}
\item{gaussian_l2_sensitivity}{Public joint L2 sensitivity.}
\item{objective}{One of \code{"auto"}, \code{"marginal_95_abs"}, or
\code{"simultaneous_95_abs"}.}
}
\value{
A data frame with fixed per-capsule parameters, granular-Laplace and
TV-accounted approximate-Gaussian utility values, the utility winner, the
actually deployed formal backend, support flags, and explicit
no-operation-quota indicators. Gaussian is calibration-only until its joint
formal plan is selected by the signed server manifest.
}
\description{
Reproduces the public deterministic utility comparison without accessing
protected data. Selection remains fail-closed to deployed formal backends. It
quantifies DP mechanism noise, not sampling uncertainty and not a universally
optimal privacy setting.
}
