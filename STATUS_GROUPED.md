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

## Resumed 2026-09-18 — clause 5 and package-check recovery

The old composed-GC traffic blocker is withdrawn; no measured release capacity
is yet admitted. Server components and unresolved exact-arithmetic/fusion
interface are documented in ../dsVert/STATUS_GROUPED.md and
../dsVert/CLAUSE5_ARITHMETIC_GROUPED.md. Reader remains fail closed.

Inherited client R CMD check finished with 1 ERROR, 3 WARNINGs, 1 NOTE. Tests:
23868 pass, 1 failure, 51 skipped. The error was OUR missing method inventory
for the new exports, not a pre-existing failure. Fixed by 429b2a7: three
quarantine entries via additive shared-registry lines; new grouped helper and
tests. Existing inventory/maturity expectations now explicitly cover the
quarantined methods without weakening assertions for existing methods.
Focused registry/inventory/maturity/contract tests: 1237 assertions pass,
zero failures, errors or warnings. Existing warnings remain documented in
server status; no clean full-check claim.

New-snapshot R CMD build passed; full post-fix R CMD check is running at
/workspace/dsvert/grouped/clause5/client-check.log (exit: client-check.exit).
Read these results before relaunching. No production release enabled.

## Addendum 3 — exact arithmetic components

Client 91ea81b mirrors server's new `grouped-lmm-stats-f264-q64-v2` profile and
exact f100/no-intermediate-rounding semantics. Focused contracts: 151 assertions
pass. Paired server synthetic DSLite comparison now calls the matching tagged
LMM statistic and passes 202 assertions across all five families/three epsilons.
The reader remains fail closed. Backend ownership is resolved; see server
STATUS_GROUPED.md for measured products, remaining GEE whitening/certification,
full-release matrix and authenticated Step-2 fusion. No admitted capacity or
completed current-snapshot R CMD check is claimed. Historical traffic and
ownership blocker statements above are superseded.
