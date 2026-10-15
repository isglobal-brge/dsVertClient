# Cross-owner correlation and PCA in 1.4.2

Correlation and PCA are post-processing of the existing signed Gaussian
Synopsis. No regression fit is required for offline input. The registered set
must contain an intercept and at least two bounded numeric variables. All pairs
use the **joint complete cases of the full registered set**, even if the analyst
requests a subset. This differs from the default pairwise-complete `ds.vertDPCor`.

## Custodian configuration

These are administrator-side examples, not analyst privacy controls. Provision
separate persistent identities/state directories, pin every peer out of band,
and authorize local raw inputs for `ds.psiAlign` using the existing
`dsvertPSISourceDescriptor` workflow. `DA` below is the resulting private padded
alignment. The same dataset identity/version and alignment group refer to the
logical cohort at all peers. Use scientifically approved bounds and capacities;
the illustrative bounds and epsilon do not promise useful correlation precision.

The following common configuration is applied independently in each custodian
process, with its own `local_columns` and `bounds`. It is frozen when registered.

```r
configure_correlation <- function(peer, peers, local_columns, bounds, spec) {
  options(
    dsvert.peer_name = peer,
    dsvert.dp.total_epsilon = 3,
    dsvert.dp.total_delta = 2^-100,
    dsvert.dp.domain = "correlation-study-v1",
    dsvert.dp.cohort_id = "cohort-v1",
    dsvert.dp.adjacency = "add_remove_patient",
    dsvert.dp.patient_column = "patient_id",
    dsvert.dp.unit_capacity = 1024L,
    dsvert.dp.max_records_per_unit = 1L,
    dsvert.dp.numeric_grid_bits = 16L,
    dsvert.dp.overflow_policy = "reject_snapshot",
    dsvert.dp.designated_noise_peers = peers[1:2],
    dsvert.dp.datasets = list(DA = list(id = "correlation-cohort", version = "v1")),
    dsvert.dp.capsule_dataset_mapping = list(DA = local_columns),
    dsvert.dp.numeric_bounds = bounds,
    dsvert.dp.workload_scope = list(
      mode = "catalog_v1", numeric_moments = character(),
      categorical_marginals = character(), categorical_pairs = list(),
      correlations = list()),
    dsvert.dp.gaussian_specs = spec)
}
```

Two-owner correlation-only registration (run each call in the named process):

```r
peers <- c("site_a", "site_b")
# The smallest owner-qualified reference, site_a$age, is the fixed anchor.
configure_correlation("site_a", peers, "age", list(age = c(18, 90)),
  list(cor_primary = list(version = "v2", dataset = "DA",
    outcome = "age", predictors = "height", intercept = TRUE)))
configure_correlation("site_b", peers, "height", list(height = c(120, 220)), list())
```

Three-owner registration with both same-owner and cross-owner pairs (a separate
study configuration; do not change the first registered artifact in place):

```r
peers <- c("site_a", "site_b", "site_c")
configure_correlation("site_a", peers, c("age", "bmi"),
  list(age = c(18, 90), bmi = c(15, 50)),
  list(cor_mixed = list(version = "v2", dataset = "DA",
    outcome = "age",
    predictors = c("bmi", "height", "marker"),
    intercept = TRUE)))
configure_correlation("site_b", peers, "height", list(height = c(120, 220)), list())
configure_correlation("site_c", peers, "marker", list(marker = c(0, 100)), list())
```

Registration fragments must use the existing server's **canonical schema labels**:
bare names for globally unique columns, `owner$column` for physical names shared
by different owners. For example, two `x` columns use `site_a$x` and `site_b$x`.
The signed schema still binds every unique bare name to its actual owner/dataset.
The public anchor identity (`site_a$age` here) and analyst references remain
owner-qualified. A raw registration fragment using qualified aliases for unique
columns is rejected by the existing exact local-fragment check; 1.4.2 does not
relax that check or introduce a new registration protocol.

Only the anchor/outcome owner registers the global analysis fragment. All peers
approve and sign the global schema, workload, source and compile receipts; every
source owner supplies its protected inputs. Two designated authorities produce
and attest the joint noisy release. The other owners are not omitted from
provenance dependencies. With a scientific outcome, retain that registered role;
without one, choose the smallest reference using R's `sort(..., method="radix")`
at registration. Changing the anchor or labels is a **new artifact/release** and
may compose with the old one. Do not optimize the anchor using private utility.

Source data stay frozen under the 1.4.1 authorization epoch workflow. A new
publication needs explicit custodian reauthorization; retries retain the same
event/epoch and sticky answer. Generic imported alignments require the existing
explicit migration. No new authorization bypass or public content hash is used.

## Analyst use

```r
# The current server capabilities and signed workload must validate.
r <- ds.vertCor("DA", c("site_a$age", "site_b$height"),
                analysis_id = "cor_primary", datasources = conns)
r_explicit <- ds.vertDPCor("DA", analysis_id = "cor_primary",
                artifact_family = "gaussian_models", datasources = conns)
p <- ds.vertPCA(cor_result = r)

# An already acquired current Gaussian result OR its certificate is enough.
# trusted_pins is obtained independently, never copied from an untrusted file.
r_offline <- ds.vertCor(gaussian_result = saved_current_certificate,
                       trusted_pinset = trusted_pins)
p_offline <- ds.vert.pca(gaussian_result = saved_current_certificate,
                        trusted_pinset = trusted_pins)
```

Offline inputs must be caller-anchored or anchored by a current authenticated
transport session. Every release, replay receipt, provenance dependency,
coordinate and cross authority attestation is revalidated. Supplied matrices,
coefficients and validity flags are not evidence. Conversion makes no DSI call,
requires no regression fit and costs no extra epsilon/delta. Online acquisition
can call the existing release/replay workflow. PCA returns eigenvalues, loadings
and variance, never individual scores; unstable eigendirections remain marked.

Clipping, complete-case exclusions, protected product quantization, the full
catalog's joint DP noise, clamping and PSD projection affect utility. The
Gaussian block contains all pair products, but is asymmetric in its storage
roles and also includes marginal coordinates. It is not an optimally calibrated
symmetric correlation-only release. Report the signed bounds, roles, registered
cohort and mechanism regions; these do not cover population sampling uncertainty.
