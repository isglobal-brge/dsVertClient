# Connector-neutral remote-surface attestation (2026-08-08)

## Contract boundary

The runtime attestation identifies one logical dsVert disclosure-safe surface;
it does not identify a DSI connector or an Opal profile. The token is derived
from the exact aggregate/assign name, type and `dsVert::function` mapping in the
installed dsVert `DESCRIPTION`:

```
dsvert-custodian-surface-attestation-v1:dsvert-disclosure-safe-v1:<sha256-of-canonical-method-map>
```

The token is public metadata, not a secret and not live introspection. It is a
custodian-owned operational assertion that the *effective callable inventory*
was independently checked and exactly matches that contract: no aliases, no
extra exact/generic endpoints and no missing purpose-bound endpoint. Installing
dsVert, reading its `DESCRIPTION`, or starting an R process never creates this
assertion automatically.

## Opal administration

`provision_opal.R` remains deliberately Opal-specific. It disables and
restricts the dedicated non-default `dsvert` profile, removes its complete
prior aggregate/assign inventory, installs and registers the staged tarball's
exact contract, verifies the resulting inventory, and only then persists the
connector-neutral token as the server-side profile option
`dsvert.remote_surface_attestation`. A final inventory check occurs before the
profile is enabled. The generic token does not weaken these Opal controls.

## Armadillo/Rock administration

For every Rock service/image behind Armadillo, the server administrator must:

1. deploy an isolated dsVert DataSHIELD surface and install the intended dsVert
   server tarball;
2. inspect the effective Rock aggregate/assign registration using the
   deployment's administrative interface and compare it exactly with the
   installed dsVert `AggregateMethods` and `AssignMethods` contract;
3. only after that comparison succeeds, calculate the deterministic expected
   token inside the same server image with
   `Rscript --vanilla -e 'cat(dsVert:::.dsvert_remote_surface_attestation_expected())'`;
4. set `DSVERT_REMOTE_SURFACE_ATTESTATION` in the Rock server/container
   environment, not in the analyst/client R session, and restart the service;
5. verify `ds.vertSecurityStatus(require_ready = FALSE)` reports
   `verified_custodian_attestation` for that server.

The helper in step 3 computes the expected contract token only. It does not
inspect Rock registration and must never be used as evidence that step 2
succeeded. If Rock exposes any additional callable method, the administrator
must correct that surface rather than set the token.

## Fail-closed lifecycle

- A missing, malformed, stale or conflicting attestation keeps readiness false.
- If both the Opal-style server option and the Rock-style server environment
  variable are present, they must be byte-identical; disagreement fails closed.
- The public client call is always the zero-argument
  `dsvertSecurityProfileDS()` expression. dsVertClient never forwards a client
  option, environment value or caller-supplied attestation.
- Before changing packages or registrations, remove the deployment token or
  disable the surface. Reconcile the complete inventory, derive the new token,
  provision it server-side, and only then re-enable service.
- For a DSV1 wire upgrade, keep release traffic paused and complete that
  deploy/reconcile/attest sequence on every server before deploying the paired
  client. The intermediate mixed generation is expected to fail at the first
  framed, data-free manifest phase; it is not a compatibility window.
- `verified_custodian_attestation` establishes local surface eligibility only.
  Consortium policy/pinset readiness and route-specific dataset/numeric
  preflight remain separate checks; formal GLM and Cox remain not ready.
