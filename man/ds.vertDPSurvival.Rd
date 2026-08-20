\name{ds.vertDPSurvival}
\alias{ds.vertDPSurvival}
\alias{print.ds.vertDPSurvival}
\title{Differentially Private Non-Parametric Survival Curves}
\usage{
ds.vertDPSurvival(data_name, analysis_id, server = NULL,
  datasources = NULL)

\method{print}{ds.vertDPSurvival}(x, ...)
}
\arguments{
  \item{data_name}{Name of the registered protected data frame.}
  \item{analysis_id}{Custodian-owned survival specification id.}
  \item{server}{Optional datasource name. If omitted, the
  lexicographically first connected datasource is used deterministically.}
  \item{datasources}{DataSHIELD connections.}
  \item{x}{A \code{ds.vertDPSurvival} object.}
  \item{...}{Additional print arguments (currently unused).}
}
\value{
A \code{ds.vertDPSurvival} object containing the single released DP
histogram, DP-derived risk/event/censor values, Kaplan--Meier survival,
Nelson--Aalen cumulative hazard, cause-specific cumulative incidence, signed
Synopsis provenance, and per-artifact privacy metadata.
}
\description{
Reads one canonical signed Synopsis release and derives all three
non-parametric estimators by post-processing that same sticky fixed-grid DP
histogram. Exact replay is unlimited and does not use lifetime admission.
}
\details{
Time is discretised on a public, fixed grid selected by the custodian. With
multiple records per patient, the estimand uses the first event or, when no
event exists, the latest censoring record. Delayed entry is supported by a
second fixed histogram margin.

Reported accuracy radii cover DP mechanism noise in histogram coordinates.
They are not sampling confidence intervals, do not include public-grid
discretisation error, and do not justify ordinary p-values. Partial, null, or
errored DSI results are rejected without exposing remote error text or retrying
the release.
}
