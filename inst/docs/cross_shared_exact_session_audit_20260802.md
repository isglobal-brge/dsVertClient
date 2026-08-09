# Shared cross-family exact-MPC session audit (2026-08-02)

## Decision

One immutable capsule that contains both cross-owner categorical and Gaussian
artifacts uses one authenticated exact-MPC peer session and one private
alignment gate. This reuse is safe under the package threat model because the
alignment predicate is bound to the capsule source contract, private layout,
ordered source set, designated pinned computation pair and complete coordinate
geometry. It is not an analysis-family predicate.

Categorical and Gaussian computations do not share operation authority. Every
subsequent operation retains its existing distinct producer, purpose prefix,
analysis identifier, stage identifier, artifact digest, transcript digest and
numeric certificate. The server remains authoritative: a client-side shared
context cannot bypass the requirement for one completed alignment batch in the
same authenticated session and for the same capsule and source contract.

## Performance effect

For a dual-family capsule, exact transport setup calls decrease from two to
one and full alignment-gate executions decrease from two to one. Source
transport was already common. Family-specific exact multiplications and signed
result receipts are unchanged, so numerical output and DP accounting are
unchanged. The total speedup depends on padded coordinate count, source count,
network latency and the relative cost of the family computations; no fixed
wall-clock factor is claimed.

Single-family capsules retain the previous one-session, one-gate path.

## Retry and tamper behavior

An identical shared context is replay-exact. A changed manifest, capsule,
source contract, source purpose, private layout, ordered source set, pinned
peer set, coordinate geometry, chunk geometry or family-domain contract is
rejected before a family starts. A client-visible batch identifier is not an
authority: even a syntactically valid changed identifier cannot select or
authorize server data because categorical and Gaussian bind endpoints locate
the unique completed capsule/contract batch inside the pinned session.

If either family fails, the common session is cleaned once through the
peer-specific signed capability minted for the fixed cross-owner cleanup
purpose. The generic cleanup endpoint remains blocked, and the capability can
name only the exact session and pinned-peer binding that issued it. The generic
endpoint is also unregistered and is not used as a promoted-route fallback. A
successful first family remains represented by its signed, idempotent server
receipt; retrying the capsule cannot turn a changed second-family request into
the first family's domain.

## Validation evidence

- K=2, K=3, K=4 and K=5 tests assert exactly one setup, one alignment gate, one
  cleanup and the same bound session at both family orchestrators.
- Replay tests accept an identical shared context and reject contract, pin,
  source-order, manifest, geometry, extra-field and domain tampering.
- Real categorical orchestration for K=2/3/4/5 and Gaussian orchestration for
  K=3/4/5 consume a supplied shared session without reopening transport,
  rerunning alignment or cleaning a session they do not own.
- The checked-multiplication client runner is executed for K=2/3/4/5 while
  every emitted claim/start/validity/receipt/commit name traverses the actual
  server production guard. Generic vecmul bind, GLM softplus and generic MPC
  cleanup remain rejected by that guard.
- The complete alignment, categorical, Gaussian and joint-vector capsule test
  bundle remains green. Invalid cross-signed allocation is rejected before
  source access or the shared exact gate.

The security boundary remains two pinned, semi-honest, non-colluding
computation peers with an untrusted analyst relay. Malicious computation
peers, peer collusion, host compromise and availability are not claimed.
