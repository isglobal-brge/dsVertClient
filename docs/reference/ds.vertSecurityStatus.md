# Verify the consortium disclosure profile

Checks that every connected server enforces dsVert's single
disclosure-safe profile, a matching custodian-owned attestation of the
dedicated logical dsVert surface and, by default, the coherent
joint-capsule policy, pinset, privacy epoch and authenticated
lifetime-budget control-plane handshake. The availability of a new
capsule reservation is reported separately in the joint-DP status
telemetry and does not redefine consortium readiness. The surface
attestation is an administrative assertion, not live introspection: it
must be renewed after an administrator changes the effective callable
surface. Opal and Armadillo/Rock use the same logical contract and
token; connector-specific admin tooling provisions it in the server
profile or server/container environment only. This client sends no
attestation option, token or claim in the DSI expression. Under
security-profile schema v4, the returned route map reports biomedical
joint-DP profile eligibility and authenticated control-plane readiness
separately. It explicitly does not evaluate route-specific dataset
admission, manifest construction, numeric runtime capabilities, or a
live release. The top-level `ready` value and server compatibility alias
`formal_dp_claim_eligible` apply only to that biomedical route; they
never promote formal GLM or formal Cox, whose route-specific `ready`
fields remain false while sealed.

## Usage

``` r
ds.vertSecurityStatus(datasources = NULL, require_ready = TRUE)
```

## Arguments

- datasources:

  DataSHIELD connections; active connections by default.

- require_ready:

  Fail if the custodian-attested remote surface or the biomedical
  joint-DP consortium policy/authenticated control plane is not ready.
  This is not a route-specific execution preflight and never promotes
  formal GLM or Cox.

## Value

A `ds.vertSecurityStatus` object with explicit route readiness.
