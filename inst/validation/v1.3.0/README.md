# 1.3.0 RC1 run-at-pin evidence

Executed against the immutable annotated `v1.3.0-rc1` source pair:

- dsVert: `f9d4a19980d9e78b6ab49fc64e6f8ebbd8e7c927`
- dsVertClient: `5e4c3ea7cd72b71bfac5ef91c21ee77094071888`

The package code remains at this pair (both DESCRIPTION versions are 1.3.0).
The source-only harness and worked-example scaffold include the state-path fix
described below; their SHA256s are recorded in the campaign manifest. Historical heavy cross-owner
campaigns remain evidence at their own recorded ancestor pairs; this directory
does not claim to rerun or complete those campaigns. No release tag is created.

## Execution

One fresh RunPod CPU pod, `lyzk3uwe3wlds7`, uses `cpu3c-32-64`, Ubuntu 22.04,
a 50 GB container disk and no volume. The committed Linux AMD64 MPC runtime
is checked against the committed `inst/bin/SHA256SUMS`; it is not rebuilt.
R and system dependencies follow the existing CPU provisioning recipe.
The supplied R dependency archives are supplemented with Ubuntu's
`r-cran-devtools` because the battery requires it.

Use sibling source trees with `dsVert` at the server ref above and
`dsVertClient` at this evidence commit, which supplies the corrected harnesses
while keeping the client package implementation at rc1. Set
`DSVERT_VALIDATION_ROOT` to their parent directory, install the packages and
validation dependencies, and run:

```sh
export DSVERT_SERVER_SOURCE="$DSVERT_VALIDATION_ROOT/dsVert"
export R_LIBS_USER="$DSVERT_VALIDATION_ROOT/R-library"
Rscript inst/validation/v1.3.0/run_promoted_battery.R
timeout 900 Rscript tools/validate_formal_dp_e2e.R --k=2,3,5 --output=inst/validation/v1.3.0/formal_dp_e2e_v130_2026-09-20.json
Rscript inst/validation/v1.3.0/export_registry.R
DSVERT_EXAMPLE_N=50 timeout 2400 Rscript -e 'source("inst/validation/v1.3.0/worked_example_driver.R", echo=TRUE, max.deparse.length=Inf)'
DSVERT_EXAMPLE_N=500 timeout 2400 Rscript -e 'source("inst/validation/v1.3.0/worked_example_driver.R", echo=TRUE, max.deparse.length=Inf)'
```

The evidence-only battery driver evaluates the committed historical runner's
`devtools::test` call with `stop_on_failure=FALSE`. It uses the v1.2.1 manifest's
split: full Synopsis at K=2; Count and Gaussian LASSO focus at K=2,3,5;
describe/Gaussian/GLM-grid/survival at K=3,5; cross-owner tamper at K=3.
The battery has a six-hour wall-time cap. Its CSV checkpoints completed gates;
the technical log also preserves completed assertions within an interrupted gate.
A driver exit code alone does not establish a passing test suite: consult the
structured failed/error/skipped columns and the manifest.

## Worked examples

The scripts adapt `../v1.2.0/worked_example_{driver,custodian}.R`.
Server A retains y, x1 and x2; Server B retains x3 and x4. The model is
`y ~ x1 + x3`, written using rc1's required owner-qualified references as
`site_a$y ~ site_a$x1 + site_b$x3`. Both custodians sign the grid contract,
using the procedure in dsVert's `inst/cross-grid-v2/validate_dslite.R`.
The signed normalized coefficient candidates, ordered as intercept/x1/x3, are
`(-0.5, 0.8, -0.6)`, `(0, 0, 0)`, and `(0.5, -0.8, 0.6)`.
The reported source-scale coefficients apply the signed [-4,4] predictor
bounds. The analysis id is `glm_primary`, epsilon is 1 and delta is 2^-100.

The original seed 20260829 and data generator are retained, including the
outcome generator's dependence on x1 and x2. The pooled reference fit uses
the requested x1/x3 model; the first grid candidate is therefore not labelled
as the true coefficient vector for that model. Each sample size starts fresh
custodial identities and obtains its own sticky draw. PSI output reports
alignment status without disclosing intersection cardinality. The stated n
comes from the synthetic generator, not a PSI release.

The original capacity of 1024 is retained; PSI capacity is made explicit and
the grid uses rc1's supported 16-bit release quantization and certified
piecewise-v2 profile. Custodian signing is distinct from the exported analyst
calls. No private keys or protected runtime state are included here.

## Multi-process bootstrap correction

The original attempts failed at `ds.getIdentityPks(conns)` with the client's
aggregate transport rejection. Calling `dsvertIdentityPkDS()` directly through
one isolated peer's DSLite server exposed the underlying error verbatim:

```
The DP noise root must be outside temporary and installed-library trees
```

The verified allowlist already contained all 97 required endpoints. Identity
initialization follows `dsvertIdentityPkDS()` → `.get_identity_keypair()` →
`.get_identity_seed()` → `.dsvert_dp_noise_root_for_identity_recovery()` and
reaches the server's DP noise-root path validation; both old
scaffolds placed `dsvert.state_dir` under `tempdir()`. The rc1 service-state
contract rejects that location. The fix creates owner-only, unique state
directories beside the sibling source trees, outside temporary and installed
library trees. Run from a persistent, writable validation root. The formal
harness retains its guarded cleanup and explicit `--keep-state` behavior; the
worked example removes its state on completion and registers exit cleanup.
No server policy, calibration, admission limit, package implementation, release manifest,
or fleet configuration was changed.

The worked-example driver asserts successful fit and Count releases and exact
certificate and whole-release digest equality on replay. Synthetic reference
fits are for orientation; they are not private-data releases. Fresh custodial
secrets produce independent DP draws for each sample size.

The first corrected n=50 attempt passed identity pinning, PSI and grid signing,
but its 900-second command timeout interrupted the cross-owner fit. The signed
1024-slot workload uses 32-row batches; native completion markers showed forward
progress. Both final captures use a 2400-second command timeout, with the same
capacity, grid, numerical profile and privacy policy. The interrupted attempt's
private state and orphan workers were removed before starting fresh captures.

## Final results

The retained promoted-historical battery has 48 executed test blocks and 1,898
passing expectations, with zero failures or errors. The retained registry has
114 rows: 108 promoted, 2 provisional and 4 quarantine, covering 97 endpoints.
The repaired formal DSLite harness passes K=2,3,5, including identity persistence,
status-only PSI, sticky Count replay, no additional sampler after client restart,
and signed-release revalidation after service restart. Both worked examples pass
with byte-identical certificate and whole-release digest replay and a DP Count
release. All five campaign exit statuses are zero. The existing source-harness
cleanup tests also pass (30 expectations).

The campaign manifest records UTC execution dates, package refs,
runtime checksums and evidence hashes. These results retain the original battery
and registry byte-for-byte and replace the three failed multi-process captures.
