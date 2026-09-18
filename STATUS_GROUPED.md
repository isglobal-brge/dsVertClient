# Grouped client status

Resumed 2026-09-18 after quota interruption. Base cb26ecd.
Recovered signed contract mirror, dp_lmm_grid/dp_glmm_grid/dp_gee_grid,
documentation and contract tests. Production reader deliberately fails closed.
Authoritative lane milestones and gates: ../dsVert/STATUS_GROUPED.md.

Final resumed implementation commits: 83512b6, 872bb6d. Client focused suite:
5 tests / 151 assertions, no failures/errors/warnings. All production entry
points remain fail-closed. Grouped production/full-envelope delivery is blocked
by the arithmetic gate documented in ../dsVert/BLOCKED_GROUPED.md.
Pod R CMD check remains active at /workspace/dsvert/grouped/final/client-check.log;
initial unavailable dependencies and current check state are recorded in the
server STATUS_GROUPED.md. No completed clean package-check claim.


## Resumed 2026-09-18 — revised arithmetic gate

Clean at resume. Commit 016e400 binds the new certified range-reduced profile
and its manifest hash, matching server 4f2d314. Focused client contracts:
151 assertions, no failures/errors/warnings. Server-side paired-source tests
also pass the 15 synthetic DSLite/lmer/glmer/geeglm epsilon comparisons.
The old <=2000-AND blocker is void; the new <=5000 scalar gate passes.
Whole-release traffic is independently excessive for the current kernels;
see the server's revised status/benchmark evidence. Client reader stays closed.
R CMD check results are being collected from the prior snapshot; those do not
constitute a final clean check of this hash update.
