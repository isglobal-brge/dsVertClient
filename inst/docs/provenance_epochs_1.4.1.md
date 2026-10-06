# Provenance epochs in dsVertClient 1.4.1

Count, Frequency and Synopsis identify each release using public source and
analysis semantics plus authenticated authorization epochs. Every new explicit
custodian authorization allocates an independent 256-bit random epoch, including
an authorization of identical bytes. Every participating PSI peer contributes
an epoch; multi-dataset Synopsis binds the full authenticated dependency closure,
including sources which affect membership without contributing a selected
column. Equal stable identities therefore describe the same admitted frozen
material, rather than testing equality of private aligned content.

## Coordinated upgrade and custodian workflow

Install dsVert and dsVertClient 1.4.1 together. The client performs a public
capability check across all peers before PSI initialization, Count compilation,
Frequency Claim or Synopsis Bootstrap can materialize private sources. Older or
mixed protocols fail closed. Analyst arguments are unchanged. The public
capability check is data-free. Synopsis Bootstrap and capsule planning may read
and durably capture an authorized source on first registration, and advertise
`data_access = TRUE`; they do not release source statistics.

A raw padded-PSI source can register automatically once when it matches its
custodian-approved raw digest. Explicit descriptors with authenticated
pre-alignment provenance may also register automatically. The server records a
custodian notice. Re-running PSI, reading a descriptor and restarting a process
are not new authorization events.

**Legacy generic aligned imports, and descriptors without reconstructible
authenticated provenance, require one explicit migration.** Every analysis of
an unregistered import fails with the same migration error regardless of its
content. In the trusted server R process, approve the descriptor for the current
frame using the existing provisioning route and run:

```r
dsVert::dsvertAuthorizeSource(
  data = aligned, descriptor = getOption("dsvert.dp.datasets")$aligned,
  event_id = "migration-1.4.1", kind = "dataset", patient_column = "id")
```

Use the configured privacy-unit column and repeat for every imported dependency.
The approved private descriptor must match the frame. Migration durably captures
the frame as a materialized root; it does not reconstruct unknown historical
upstream PSI authorizations. The root's explicit authorization defines the new
boundary. Keep descriptor digests and authorization records private.

For a later publication, approve the new descriptor and call
`dsVert::dsvertReauthorizeSource()` with the same arguments and a new `event_id`.
For raw PSI use `kind = "psi"` and the approved raw-source descriptor. A retry of
the same event must reuse its original `event_id`; changing that ID deliberately
allocates a new epoch, even for identical data. Custodians can inspect retained
metadata with `dsVert::dsvertListSourceAuthorizations()`. These three functions
are local administrative APIs and are not on the DataSHIELD analyst allowlist.

## Staleness, replay and recovery

Registered source material is frozen durably. Editing a source in place does
not select a fresh identity or expose a data-dependent freshness flag. Existing
public identities continue to use their captured material and committed answers;
a request that cannot use intact authorized material fails closed. Reauthorization
is the publication event that allows fresh data. Answers can remain stale
indefinitely if custodians never reauthorize, and ordinary DP uncertainty
intervals do not include staleness error.

Retain the authorization store together with the node identity, noise roots and
existing artifact/replay state across restarts and upgrades. Missing, corrupted
or inconsistent retained state requires custodian recovery. After provisioning,
retain `dsVert::dsvertListSourceAuthorizations()$uuid` outside the state volume
and set the server option `dsvert.dp.provenance_store_uuid` to that value across
restarts. Without this external pin, loss of both the store and its marker is
indistinguishable from a first installation. A MAC cannot detect rollback of an
entire valid state snapshot: trusted storage and consistent backups remain
required. Deleting state is not a refresh procedure. Store capacity is a public
admission limit and does not evict old identities or grant a reroll.

Count, Frequency and Synopsis results include `stable_receipt`, containing the
semantic artifact key and canonical public semantics. It is reconstructed only
after release/compilation authentication and must remain identical on exact
replay. Session IDs, PSI attestation IDs and run signatures remain separate
operational evidence and can change across independent PSI runs. Stable receipt
equality is not a test for equality of private content across authorizations.

## Guarantee and limits

With intact authenticated state, a fixed complete epoch vector and fixed public
analysis semantics admit one effective input and retain one sticky release.
Every new authorization changes all affected artifact/noise contexts and stable
receipts independently of whether bytes changed. Numeric outputs can still
coincide by chance. Content digests remain private admission/integrity checks.
The existing samplers, sensitivity bounds, epsilon/delta allocations and backend
selection are unchanged. Legacy capsule and formal GLM/Cox paths retain their
existing contracts.

This closes the specified identity/receipt equality channel, not the full
DataSHIELD transcript. Authorization timing, raw-source first-load mismatch
refusals, sealed aligned-object tampering, infrastructure/integrity failures,
PSI leakage and longitudinal
composition retain their existing limits. Publication scheduling based on
private data needs its own analysis. A new authorization is a new release for
composition; replay is not a new privacy release. See isglobal-brge/dsVert#25.
