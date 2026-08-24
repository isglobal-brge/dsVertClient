# Disclosure audit — formal Cox public result and legacy computation

`ds.vertCox()`, `ds.vert.cox()` and `ds.vert.coxph()` accept a
custodian-configured `formal_analysis_id` only to read one already completed,
two-authority-signed sticky formal Cox opening. The read-only response contains
the certified coefficient lattice values and hazard-ratio point/range values;
it cannot start a Cox computation or reveal source records, shares, paths,
keys, retries or intermediate state.

Security-profile schema v5 continues to report
`route_claims$formal_cox_ready = FALSE` with state
`sealed_no_recipient_encrypted_r_dsi_lifecycle_or_end_to_end_numeric_certificate`:
new formal Cox execution is not ready. It separately reports the completed
read-only result boundary as
`route_claims$formal_cox_public_result_ready = TRUE`. The top-level client
`ready` value and server compatibility alias `formal_dp_claim_eligible` do not
authorize a new Cox analysis.

Without `formal_analysis_id`, the compatibility frontdoors fail before any DSI
call. `ds.vertCoxProfileNonDisclosive()` accepts the same read-only certificate
selector for compatibility. `ds.vertCoxDiscreteNonDisclosive()` and
`ds.vert.cox(method = "discrete")` accept a different, custodian-configured
selector only to read a completed two-authority binomial certificate bound to
the supplied `Surv()` formula and one fixed pooled-logistic time grid. That
reader has no caller-controlled grid, expansion, optimisation or privacy
parameter; it cannot start a new release and returns no covariance or sampling
inference. It is never an interchangeable Cox partial-likelihood result. The
legacy calculation endpoints remain absent from `AggregateMethods`.

## Why new Cox computation remains unavailable

The retained research implementation used exact score, information and
likelihood aggregates and, in older variants, row-order risk-set metadata.
Encryption and pinned peer identities protect transport and peer
authentication, but they do not turn deliberate exact analyst openings into a
contribution-bounded DP release. Repeated or adaptive exact model aggregates
also have no non-reconstruction guarantee.

The discrete-time public reader targets a fixed-grid pooled-logistic hazard
estimand, not a Cox partial-likelihood estimand, and therefore cannot be
presented as an interchangeable fallback. A fresh discrete-time computation
remains unavailable.

## Remaining promotion gate

A releasable Cox route needs one purpose-bound, signed and sticky joint-DP
capsule with explicit clipping and contribution bounds; private risk-set
evaluation; certified ties, strata and delayed-entry semantics; covariance and
identifiability certificates; and independent multi-process DSI validation for
K=2 and K>=3. Internal prototype components are not release evidence until all
of those gates are closed together.

This note records the boundary: it is not authorization to call the legacy
implementation or register any retired server endpoint.
