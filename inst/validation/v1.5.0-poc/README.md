# DP-IRLS proof of concept (isglobal-brge/dsVert#26): end-to-end run, 15 Oct 2026

This is one complete fixed-path DP Newton/IRLS release through DataSHIELD. It is evidence that the mechanism works
end to end, not a qualification of utility or inference.

**Setup.**
- Two custodians in separate R/DSLite processes, with real PSI alignment.
- Real two-party exact garbled circuits with IKNP/KOS OT.
- v4 two-full-draw discrete-Laplace noise on both lanes.
- Dual-signed stage publications, exact rational successor updates and an offline client chain verification.
- Model: binomial/logit `y ~ x + z`, n = 64 synthetic persons, ε = 8 over three releases (two optimizer stages and
  the terminal audit, each with score and curvature lanes).

| Stage | Wall seconds | Protected compute | Noise | Publication |
|---|---:|---:|---:|---:|
| 0 (β = 0) | 268.5 | 188.6 | 6.1 | 31.6 |
| 1 | 453.3 | 367.4 | 6.4 | 30.5 |
| 2 (terminal audit) | 492.9 | 398.1 | 8.0 | 33.9 |

Path 1,219 s. Setup and PSI 53 s. Repeated client calls returned byte-identical chains.

| Term | DP estimate | glm() on the pooled synthetic data |
|---|---:|---:|
| (Intercept) | −0.363 | −0.455 |
| x | 0.757 | 1.414 |
| z | −0.317 | −1.222 |

**Reading.**
- The signs agree.
- At n = 64 the DP noise and the two-step schedule (stage-0 trust truncation, stage-1 floor shift) leave the
  estimates far from the MLE.
- The refusal policy correctly marked sampling inference **unavailable**.
- Utility at realistic n is the subject of the #26 simulations.
- Runtime is dominated by DSI relay round trips, not by the circuit: the same stage runs in seconds in-process.

The files are `poc-report.json`, `v3-poc8-evidence.json`, `v3-poc8-smoke-ledger.jsonl` and `poc-n64-run8.log`.

To reproduce: `DSVERT_B3_POC=1 tools/validate_dp_irls_b3.R`, then `tools/poc_dp_irls_report.R`.
