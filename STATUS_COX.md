# Cox client handoff

Resumed 2026-09-18 from the quota-interrupted worktree. Preserved and committed
the recovered client in 19b0e02, added the synthetic DSLite/oracle/coxph comparison
and server/client canonical parity in 6b993e2.

The implementation is in R/dp_cox_grid_cross.R; one registration function,
.dsvert_dp_cox_grid_cross_client_register(), returns its `entry` (dp_cox_grid),
validator and disabled release function. The entry is namespace-internal while
fused release and public inventory/maturity registration are pending. Both
real custodian signatures are required; all production release attempts fail
with the fixed public failure. No protected optimizer, standard errors or
baseline-hazard path is enabled.

The initial full R check detected our premature NAMESPACE export through the
existing exact-surface inventory test. Removed that export and preserved the
existing inventory/status code and tests. Cox contract/parity checks now pass
39 expectations; the full inventory suite also passes. The two-peer PUBLIC
synthetic DSLite test passed 23 expectations at epsilon 1,4,8, comparing a
test-tagged Go integer evaluator with pure R and survival::coxph (Breslow).
Its reference noise is not evidence of production joint-DP/sticky execution.

The authoritative numeric proof, cost model, full validation outcomes, commits
and remaining gates are in the sibling server checkout: STATUS_COX.md,
DECISIONS_COX.md, NUMERIC_CERTIFICATE_COX.md and INTEGRATION_COX.md. The current
full-size resource model misses the 30 GB envelope; promotion remains blocked.

Final functional commit: c35322e. The initial full client check reached tests
and caught the export/inventory mismatch; the focused repair passed. Static
checks show three pre-existing warning categories (MI non-ASCII, codoc mismatch,
ordinal duplicate argument). The long full checks were stopped at the blocked
handoff, not reported as passes. The server BLOCKED_COX.md and machine-readable
package_check_status.json distinguish completed Cox tests from partial broad
coverage. No modified pre-existing files remain relative to cb26ecd.

## 2026-09-18 resumed — revised share-based lane

The earlier 30 GB gate/handoff is superseded by GATE_REVISION clauses 4 and 5.
Client commits 7ffacd0, a1966dc and df6b7eb mirror the single packed permutation,
share prefix/loss sums, streaming checked OT and private-key receipt mode. No
public export or release lifecycle is enabled; registration remains internal.

Full client R CMD check completed: **25425 passing expectations, 0 failures,
0 test warnings, 48 skips**. Overall **3 pre-existing WARNING categories** remain:
MI non-ASCII source, codoc mismatches, and duplicate ordinal Rd argument.
The explicit tagged two-custodian DSLite/pooled-coxph test separately passed all
23 expectations without skips. Server/client signature/canonical checks remain
part of that completed suite. The server Go family gate passed 21 top-level
and 7 nested tests after correcting its fixture working directory.

The new family kernel's first completed pod envelope, 2000 rows x 16 candidates,
is 8.473 GB and 21.64 minutes including cold compilation, with exact oracle
equality. The remaining matrix and final signed admission are still pending;
this is not a production source-to-joint-DP claim. The sibling server STATUS,
INTEGRATION and NUMERIC_CERTIFICATE documents contain the authoritative evidence.

Resumed admission validation: client targeted suite now passes 50 expectations;
server targeted suite passes its original 626 plus 11 admission expectations.
Client Rd validation passes. Paired installed baseline/current synopsis-artifact
logs reproduce the same missing test-helper error (source-loaded tests pass).
Related installed synopsis failures are being checked against that baseline.
