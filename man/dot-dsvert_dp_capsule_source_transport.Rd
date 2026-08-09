\name{dot-dsvert_dp_capsule_source_transport}
\alias{.dsvert_dp_capsule_source_transport}
\title{Relay an immutable biomedical capsule source to its pinned peers}
\usage{
.dsvert_dp_capsule_source_transport(
  manifest_json,
  datasources,
  status = NULL,
  allocation_openings = NULL,
  .aggregate = DSI::datashield.aggregate
)
}
\arguments{
\item{manifest_json}{Canonical immutable biomedical capsule manifest.}

\item{datasources}{Complete named DataSHIELD connection set.}

\item{status}{Optional already validated full capsule-status handshake.}

\item{allocation_openings}{Canonical cross-signed allocation openings,
named exactly by the two designated peers.}

\item{.aggregate}{Injectable DSI aggregate implementation for tests.}
}
\value{
A redacted internal sampler-handoff receipt.
}
\description{
Internal orchestration only. It performs a complete reusable-capsule status
handshake, relays paired authenticated Ring128 ciphertexts in canonical
owner/chunk order, and returns public completion metadata. It never returns a
coordinate, share, mask, seed, plaintext, or patient-derived digest.
}
\keyword{internal}
