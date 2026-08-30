# =========================================================================
# dsVert / dsVertClient v1.2.0 -- chapter-4 worked example (binomial GLM)
# Vertically partitioned synthetic data over DSLite, exported analyst API.
# Custodian provisioning uses process-isolated DSLite peers exactly as in
# the released validation harness (tools/validate_formal_dp_e2e.R); see
# out/worked_example_custodian.R for the sourced scaffold.
# =========================================================================

WR <- "/private/tmp/claude-501/-Users-david-Documents-GitHub-thesis/a3f58816-b8a7-4aef-b778-fbf9c9639dd7/scratchpad/dsvert_campaign"
SERVER_DIR <- "/private/tmp/claude-501/-Users-david-Documents-GitHub-thesis/a3f58816-b8a7-4aef-b778-fbf9c9639dd7/scratchpad/snapshots-current/dsVert-120"
MPC_BINARY <- file.path(SERVER_DIR, "inst", "bin", "darwin-arm64", "dsvert-mpc")

source(file.path(WR, "out", "worked_example_custodian.R"))
library(dsVertClient)
packageVersion("dsVert")
packageVersion("dsVertClient")

## ---- Synthetic vertically partitioned data (CUSTODIAN side) -------------
## n = 500 common patients; site_a additionally holds 3 local-only records
## and site_b 4, so PSI has a real intersection to compute.  Outcome and two
## covariates live on site_a; two further covariates live on site_b.
set.seed(20260829)
n <- 50L
patient_ids <- sprintf("patient-%04d", seq_len(n))
x1 <- runif(n, -3, 3)
x2 <- runif(n, -3, 3)
x3 <- runif(n, -3, 3)
x4 <- runif(n, -3, 3)
y  <- rbinom(n, 1L, plogis(-0.5 + 0.8 * x1 - 0.6 * x2))

raw_site_a <- data.frame(
  patient_id = c(patient_ids, sprintf("site-a-extra-%02d", 1:3)),
  y  = c(y,  rbinom(3L, 1L, 0.5)),
  x1 = c(x1, runif(3, -3, 3)),
  x2 = c(x2, runif(3, -3, 3)),
  stringsAsFactors = FALSE)
raw_site_a <- raw_site_a[sample(nrow(raw_site_a)), ]

raw_site_b <- data.frame(
  patient_id = c(patient_ids, sprintf("site-b-extra-%02d", 1:4)),
  x3 = c(x3, runif(4, -3, 3)),
  x4 = c(x4, runif(4, -3, 3)),
  stringsAsFactors = FALSE)
raw_site_b <- raw_site_b[sample(nrow(raw_site_b)), ]

## Local (non-federated) reference fit on the pooled data, for orientation
## only -- the federated route below never sees these rows together.
coef(glm(y ~ x1 + x2, family = binomial()))

## ---- Custodian DP Synopsis policy (per peer, server-owned options) ------
dataset_id <- "worked-example-binomial-v120"
state_root <- file.path(tempdir(), "we-state")
dir.create(state_root, recursive = TRUE, mode = "0700")

make_spec <- function(name) {
  state_dir <- file.path(state_root, name)
  dir.create(state_dir, recursive = TRUE, mode = "0700")
  list(name = name, state_dir = state_dir,
       dslite_home = file.path(state_dir, "dslite"),
       mpc_binary = MPC_BINARY,
       dataset_id = dataset_id, dataset_version = "v1")
}
specs <- list(site_a = make_spec("site_a"), site_b = make_spec("site_b"))

base_policy <- function(spec) list(
  total_epsilon = 1,
  total_delta = 2^-100,
  domain = "worked-example-binomial",
  cohort_id = "worked-example-cohort-v1",
  adjacency = "add_remove_patient",
  patient_column = "patient_id",
  unit_capacity = 1024L,
  max_records_per_unit = 1L,
  fixed_cohort_size = NULL,
  overflow_policy = "reject_snapshot",
  datasets = list(DA = list(id = spec$dataset_id,
                            version = spec$dataset_version)),
  workload_scope = list(
    mode = "catalog_v1", numeric_moments = character(),
    categorical_marginals = character(), categorical_pairs = list(),
    correlations = list()),
  designated_noise_peers = c("site_a", "site_b"),
  synopsis_state_path = file.path(spec$state_dir, "dp-synopsis"))

## site_a additionally signs a finite binomial likelihood grid: the
## custodian-owned candidate coefficient vectors (Intercept, x1, x2) that
## the promoted ds.vertGLM binomial route is allowed to select among.
policy_a <- c(base_policy(specs$site_a), list(
  numeric_bounds = list(y = c(0, 1), x1 = c(-4, 4), x2 = c(-4, 4)),
  capsule_dataset_mapping = list(DA = c("y", "x1", "x2")),
  gaussian_specs = list(glm_primary = list(
    version = "binomial_grid_v1", dataset = "DA",
    outcome = "y", predictors = c("x1", "x2"), intercept = TRUE,
    beta_grid = list(c(-0.5, 0.8, -0.6),
                     c(0, 0, 0),
                     c(0.5, -0.8, 0.6))))))
policy_b <- c(base_policy(specs$site_b), list(
  numeric_bounds = list(x3 = c(-4, 4), x4 = c(-4, 4)),
  capsule_dataset_mapping = list(DA = c("x3", "x4")),
  gaussian_specs = list()))

## ---- Custodian provisioning: identity bootstrap, then pinned restart ----
peers <- list(
  site_a = we_boot_peer(specs$site_a, SERVER_DIR, raw = raw_site_a,
                        policy = policy_a),
  site_b = we_boot_peer(specs$site_b, SERVER_DIR, raw = raw_site_b,
                        policy = policy_b))
conns <- we_connections(peers)

## [ANALYST API] persistent peer identity discovery
identities <- ds.getIdentityPks(conns)
str(identities)

## Custodians pin each other's identity keys and restart their services.
pins <- list(
  site_a = c(site_b = unname(unlist(identities$site_b))),
  site_b = c(site_a = unname(unlist(identities$site_a))))
we_close_peers(peers)
peers <- list(
  site_a = we_boot_peer(specs$site_a, SERVER_DIR, pins = pins$site_a,
                        raw = raw_site_a, policy = policy_a),
  site_b = we_boot_peer(specs$site_b, SERVER_DIR, pins = pins$site_b,
                        raw = raw_site_b, policy = policy_b))
conns <- we_connections(peers)
identical(ds.getIdentityPks(conns), identities)

## ---- [ANALYST API] 1. record alignment via pinned fixed-capacity PSI ----
t_align <- system.time(
  ds.psiAlign("D", "patient_id", "DA", datasources = conns, verbose = FALSE))
t_align

aligned <- ds.isPsiAligned("DA", datasources = conns)
str(aligned)

## ---- [ANALYST API] 2. promoted binomial GLM (signed finite grid) --------
t_fit <- system.time(
  fit <- tryCatch(
    ds.vertGLM(y ~ x1 + x2, data = "DA", family = "binomial",
               analysis_id = "glm_primary", datasources = conns,
               verbose = FALSE),
    error = function(e) e))
t_fit
if (inherits(fit, "error")) {
  cat("ds.vertGLM (finite-grid binomial) FAILED with:\n")
  print(conditionMessage(fit))
} else {
  print(class(fit))
  print(fit$coefficients)
  print(fit$selected_candidate)
  str(fit, max.level = 1)
}

## ---- [ANALYST API] 3. replay: identical call, sticky signed release -----
if (!inherits(fit, "error")) {
  replay <- ds.vertGLM(y ~ x1 + x2, data = "DA", family = "binomial",
                       analysis_id = "glm_primary", datasources = conns,
                       verbose = FALSE)
  cat("coefficients identical:            ",
      identical(fit$coefficients, replay$coefficients), "\n")
  cat("certificate_sha256 identical:      ",
      identical(fit$certificate_sha256, replay$certificate_sha256), "\n")
  cat("provenance_certificate identical:  ",
      identical(fit$provenance_certificate, replay$provenance_certificate), "\n")
  cat("digest(fit):    ", digest::digest(fit), "\n")
  cat("digest(replay): ", digest::digest(replay), "\n")
  cat("whole-release digest identical:    ",
      identical(digest::digest(fit), digest::digest(replay)), "\n")
}

## ---- [ANALYST API] 4. cohort size as a DP release (same aligned object) --
count <- tryCatch(ds.vertDPCount("DA", datasources = conns),
                  error = function(e) e)
if (inherits(count, "error")) {
  cat("ds.vertDPCount FAILED with:\n"); print(conditionMessage(count))
} else {
  print(c(value = count$value, accuracy_95_abs = count$accuracy_95_abs))
  cat("count released:", count$released, " sticky_replay:", count$sticky_replay, "\n")
}

we_close_peers(peers)
sessionInfo()$otherPkgs$dsVertClient$Version
