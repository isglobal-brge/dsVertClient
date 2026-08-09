# GLM centralized-reference validation — 2026-07-31

This is a historical numerical record from the removed
`test-glm-central-reference-integration.R` compatibility harness. The harness
depended on exact releases and an unattested numeric escape hatch, so it was
deleted when dsVert moved to its single disclosure-safe profile. These figures
are useful only as a legacy regression baseline; they are not evidence for the
current formal route, a numeric certificate, a disclosure proof, a DP claim or
a production attestation.

## Scope and profile

- Two logical DSLite custodians with distinct, deterministic Ed25519
  identities, name-bound reciprocal peer pins, trusted-peer enforcement, and
  legacy peer fallback disabled.
- Vertically split synthetic data with a common, already aligned row order.
  PSI alignment itself is outside this oracle; it has separate protocol tests.
- The historical run used the now-removed exact-release selector and an
  unattested numeric escape hatch. Neither can enable a server operation in
  the current package.
- The comparison oracle is base R `glm()`. Coefficient error is reported both
  in model units and divided by the corresponding reference standard error;
  standard-error error is also reported relative to the reference standard
  error.
- No server method is evaluated directly by the fixture: all calls traverse
  the registered DSI/DataSHIELD surface.

## Observed Gaussian accuracy

The figures below are repeated-run ranges on the arm64 validation host. The
tests use fixed tolerances with additional headroom; they do not assert elapsed
time because that is hardware dependent.

| Case | Maximum coefficient error | Coefficient error / reference SE | Maximum SE error | Relative SE error | Deviance error | Test limits |
|---|---:|---:|---:|---:|---:|---|
| Baseline, n=80, p=2 | 1.24–1.40e-8 | 3.59–4.02e-7 | 3.07–3.45e-9 | 8.96–10.0e-8 | 1.17–1.32e-6 | coefficient <1e-6; SE <1e-7; deviance <2e-5 |
| Extreme signed/scaled predictors, n=120, p=3 | 0.00332–0.0211 | 1.86–2.38e-7 | scale dependent | 1.41–2.11e-7 | 1.13–1.70e-6 | coefficient/SE <1e-4; relative SE <1e-4; deviance <1e-3 |
| Small n/p, n=24, p=1 | scale dependent | 7.98–8.20e-7 | scale dependent | 4.2–5.42e-7 | 1.9–2.48e-7 | coefficient/SE, relative SE, deviance <1e-4 |
| Weighted, n=48, p=2 | 2.18–2.32e-6 | 2.49–2.65e-5 | not requested | not requested | 9.72–9.97e-6 | coefficient/SE and deviance <1e-4 |

Repeating the baseline on identical protected data with fresh protocol
randomness changed no coefficient by more than `2e-6` and no standard error by
more than `2e-7`.

An exactly collinear design (`x2 = 2*x1`) returned the typed condition
`non_identifiable` with reason `rank_deficient_design_matrix`. The exercised
unpenalized route did not silently substitute ridge, Firth, or a pseudoinverse.

## Weighted-route communication evidence

The weighted fit used 15 fixed iterations, took 38.3 seconds, and made 1,028
client-level DSI aggregate invocations: 68.53 per iteration. Each aggregate
invocation normally addresses the two custodians; this table counts client
protocol invocations, not individual HTTP requests.

| Aggregate method | Calls |
|---|---:|
| `dsvertColNamesDS` | 3 |
| `getObsCountDS` | 1 |
| `psiAlignmentManifestDS` | 1 |
| `dsvertNumericPolicyDS` | 1 |
| `glmRing63TransportInitDS` | 1 |
| `mpcStoreTransportKeysDS` | 1 |
| unclassified expression | 1 |
| `k2ShareWeightsDS` | 1 |
| `mpcStoreBlobDS` | 355 |
| `k2ReceiveWeightSharesDS` | 1 |
| `k2ShareInputDS` | 2 |
| `k2ReceiveShareDS` | 2 |
| `k2ComputeEtaShareDS` | 32 |
| `k2IdentityLinkDS` | 32 |
| `k2PrepareWeightedResidualShareDS` | 32 |
| `dsvertBeaverPolicyDS` | 2 |
| `k2OtBeaverSampleDS` | 64 |
| `k2IknpBaseReceiverSetupDS` | 2 |
| `k2IknpBaseSenderChoicesDS` | 2 |
| `k2IknpBaseReceiverEncryptDS` | 2 |
| `k2IknpBaseSenderFinalizeDS` | 2 |
| `k2IknpReceiverExtendDS` | 64 |
| `k2IknpSenderEncryptDS` | 64 |
| `k2IknpReceiverDecryptDS` | 64 |
| `k2OtBeaverFinalizeDS` | 64 |
| `k2BeaverVecmulConsumeTripleDS` | 32 |
| `k2BeaverVecmulR1DS` | 32 |
| `k2BeaverVecmulR2DS` | 32 |
| `k2FinalizeWeightedResidualShareDS` | 32 |
| `k2StoreGradTripleDS` | 32 |
| `k2GradientR1DS` | 32 |
| `k2GradientR2DS` | 32 |
| `mpcGcDS` | 4 |
| `glmRing63PrepDevianceDS` | 2 |
| `mpcCleanupDS` | 2 |

The regression test records a target of at most 20 aggregate invocations per
iteration as known unmet debt. Persistent DSI connections avoid connection
setup overhead, but do not eliminate these protocol turns; reducing them needs
a protocol-level consolidation, not request throttling.

## Open nonlinear blocker

Binomial, Poisson, offsets, weights combined with nonlinear links, separation,
and their standard errors do not yet have defensible full-stack tolerances on
the current exact-F2 route. In the attempted binomial fixture (`n=180`, `p=3`),
the first checked multiply/truncate was still running after more than two
minutes. A complete sigmoid evaluation needs many such operations, so allowing
that run to continue would not provide a practical regression test.

The test suite records this as an explicit skip rather than using infinite or
invented tolerances. These families must not be described as full-stack
reference-validated until a practical batched/vectorized exact nonlinear path
exists and the oracle measurements are rerun. Exact arithmetic primitives
alone are insufficient evidence for end-to-end accuracy or performance.
