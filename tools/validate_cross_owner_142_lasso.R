# Custodian-side companion to validate_cross_owner_142.R.
# Uses the released DSLite process scaffold, real PSI, custodian signatures and
# the production Synopsis lifecycle. No oracle, source/noise capture or mocks.
# Sourcing this file does not start a federation or release any artifact.

vertcor142_lasso_register <- function(peers, conns, identities, family,
    domain, cohort_id, capacity, bits, outcome, predictors, maximum) {
  stopifnot(family %in% c("binomial", "poisson"), length(peers) == 2L)
  cf <- function(name) get(name, asNamespace("dsVertClient"), inherits = FALSE)
  bootstrap <- cf(".dsvert_dp_synopsis_bootstrap_build_v1")(conns)
  schema <- jsonlite::fromJSON(bootstrap$manifest_bundle$schema_json,
    simplifyVector = FALSE)
  all_pins <- stats::setNames(unname(unlist(identities)), names(peers))
  policy <- list(peer_pinset = all_pins,
    peer_pinset_sha256 = cf(".dsvert_dp_capsule_source_hash")(as.list(all_pins)),
    designated_noise_peers = names(peers), unit_capacity = capacity,
    numeric_grid_bits = bits, adjacency = "add_remove_patient")
  beta <- list(c(-0.25, 0.25), c(0.25, 0.25))
  beta <- beta[order(vapply(beta, function(value)
    cf(".dsvert_joint_dp_client_json")(as.list(value)), character(1L)), method = "radix")]
  make_contract <- function(schema) {
    authenticated <- cf(".dsvert_dp_glm_grid_cross_schema_validate")(
      policy, schema$logical_snapshot, schema, cf(".dsvert_dp_glm_grid_cross_verify"))
    group <- schema$datasets$DA$alignment_group
    alignment <- list(version = "existing_prealigned_logical_dataset_v1",
      method = "pinned_psi_ordered_manifest_v1", alignment_group = group,
      public_alignment_contract_sha256 = cf(".dsvert_dp_capsule_source_hash")(list(
        logical_snapshot = schema$logical_snapshot, alignment_group = group,
        method = "pinned_psi_ordered_manifest_v1")), public_patient_dependent_hash = FALSE)
    raw <- list(version = paste0(family, "_grid_cross_v1"), analysis_id = "grid",
      dataset = "DA", outcome = outcome, predictor_order = predictors,
      beta_grid = beta, max_outcome = maximum, alignment = alignment)
    spec <- cf(".dsvert_dp_glm_grid_cross_spec")(raw, policy, authenticated, "piecewise_v2")
    artifact <- cf(".dsvert_dp_glm_grid_cross_artifact")(spec)
    list(version = "dsvert-cross-owner-grid-signed-contract-v1", spec = spec,
      artifact = artifact,
      source_contract = cf(".dsvert_dp_glm_grid_cross_source_contract")(spec, artifact))
  }
  contract <- make_contract(schema)
  workload <- list(version = "dsvert-biomedical-capsule-workload-contract-v2",
    describe = list(), survival = list(), vertical_cross = list(),
    gaussian = list(grid = list(owner_peer = names(peers)[[1L]], spec = list(
      version = contract$spec$version, dataset = "DA", contract = contract))))
  # The snapshot commits to the embedded workload. Sign its fixed point before
  # enabling that workload, using each custodian's private key in its own process.
  schema$signatures <- NULL
  snapshot_context <- list(servers = names(peers)[[1L]],
    status = stats::setNames(list(list(policy = list(domain = domain,
      cohort_id = cohort_id, peer_pinset_sha256 = policy$peer_pinset_sha256))),
      names(peers)[[1L]]))
  schema$logical_snapshot <- cf(".dsvert_dp_capsule_manifest_expected_snapshot")(
    snapshot_context, schema$datasets, schema$logical_snapshot$alignment_protocol_version,
    workload)
  schema$signatures <- lapply(peers, function(peer) peer$worker$run(function(schema) {
    identity <- dsVert:::.get_identity_keypair()
    dsVert:::.dsvert_relay_sign_message(dsVert:::.dsvert_dp_capsule_schema_message(schema),
      identity$identity_sk)
  }, args = list(schema)))
  contract <- make_contract(schema)
  contract$signatures <- lapply(peers, function(peer) peer$worker$run(function(contract) {
    identity <- dsVert:::.get_identity_keypair()
    dsVert:::.dsvert_relay_sign_message(dsVert:::.dsvert_dp_glm_grid_cross_message(contract),
      identity$identity_sk)
  }, args = list(contract)))
  cf(".dsvert_dp_glm_grid_profile_admit")(contract, policy, schema)
  workload$gaussian$grid$spec$contract <- contract
  stopifnot(identical(cf(".dsvert_dp_capsule_manifest_expected_snapshot")(
    snapshot_context, schema$datasets, schema$logical_snapshot$alignment_protocol_version,
    workload), schema$logical_snapshot))
  peers[[1L]]$worker$run(function(contract) {
    options(dsvert.dp.gaussian_specs = list(grid = list(version = contract$spec$version,
      dataset = "DA", contract = dsVert:::.dsvert_dp_canonical_json(
        dsVert:::.dsvert_dp_canonical_query_value(contract)))))
    TRUE
  }, args = list(contract))
  base <- cf(".dsvert_dp_lasso_cross_base")(contract, policy)
  candidates <- unlist(lapply(c(0.5, 0), function(lambda) lapply(beta,
    function(value) list(lambda = lambda, beta = value))), recursive = FALSE)
  lasso <- list(version = "dsvert-cross-owner-lasso-signed-contract-v1",
    spec = cf(".dsvert_dp_lasso_cross_spec")(list(version = "lasso_grid_cross_v1",
      analysis_id = "lasso", candidate_grid = candidates), base))
  lasso$signatures <- lapply(peers, function(peer) peer$worker$run(function(lasso) {
    identity <- dsVert:::.get_identity_keypair()
    dsVert:::.dsvert_relay_sign_message(dsVert:::.dsvert_dp_lasso_cross_message(lasso),
      identity$identity_sk)
  }, args = list(lasso)))
  # Fresh public admission checks the registered workload and all epoch bindings.
  admitted <- cf(".dsvert_dp_grid_cross_bootstrap")(conns)
  stopifnot(isTRUE(cf(".dsvert_dp_glm_grid_cross_equal")(schema, admitted$schema)))
  list(base = contract, extension = lasso, schema = admitted$schema,
    policy = admitted$policy, predictors = predictors)
}

vertcor142_lasso_smoke <- function(family, root, run_analysis,
    peer_hook = function(peers) invisible(NULL)) {
  stopifnot(length(family) >= 1L, length(family) <= 2L, !anyDuplicated(family),
    all(family %in% c("binomial", "poisson")),
    is.function(run_analysis), is.function(peer_hook))
  families <- family
  combined <- length(families) == 2L
  case_label <- paste(families, collapse = "+")
  started <- Sys.time()
  server <- file.path(root, "dsVert")
  out <- file.path(root, ".tmp", "vertcor")
  state <- tempfile(paste0("smoke-lasso-", case_label, "-"), tmpdir = out)
  dir.create(state, mode = "0700", recursive = TRUE)
  peers <- list()
  on.exit({ we_close_peers(peers); unlink(state, recursive = TRUE) }, add = TRUE)
  owners <- c("site_a", "site_b")
  # Tiny protocol smoke only: exercises a genuine signed finite-grid release,
  # with no claim about statistical accuracy or useful model selection.
  n <- 1L
  bits <- 16L
  domain <- paste0("lasso142-", if (combined) "combined" else families[[1L]])
  cohort_id <- "synthetic-v1"
  outcomes <- if (combined) data.frame(y_binomial = 1, y_poisson = 2) else
    data.frame(y = if (identical(families[[1L]], "binomial")) 1 else 2)
  raw <- list(
    site_a = cbind(data.frame(patient_id = sprintf("u%03d", seq_len(n))), outcomes),
    site_b = data.frame(patient_id = sprintf("u%03d", seq_len(n)),
      x2 = 0.75))
  specs <- stats::setNames(lapply(owners, function(owner) {
    path <- file.path(state, owner)
    dir.create(path, mode = "0700")
    list(name = owner, state_dir = path, dslite_home = file.path(path, "dslite"),
      mpc_binary = Sys.getenv("DS_VERT_MPC"), dataset_id = domain, dataset_version = "v1")
  }), owners)
  policies <- stats::setNames(lapply(owners, function(owner) {
    columns <- setdiff(names(raw[[owner]]), "patient_id")
    bounds <- stats::setNames(lapply(columns, function(column)
      c(0, if (identical(column, "y_poisson") ||
        (identical(column, "y") && identical(families[[1L]], "poisson"))) 4 else 1)), columns)
    # The two planned releases compose to epsilon eight. Retain the compatible
    # lifetime policy bound too; never reset a ledger or claim a global budget.
    list(total_epsilon = if (combined) 4 else 8, total_delta = 2^-100,
      lifetime_max_distinct_capsules = length(families), domain = domain, cohort_id = cohort_id,
      adjacency = "add_remove_patient", patient_column = "patient_id", unit_capacity = n,
      numeric_grid_bits = bits, max_records_per_unit = 1L, overflow_policy = "reject_snapshot",
      designated_noise_peers = owners,
      datasets = list(DA = list(id = specs[[owner]]$dataset_id, version = "v1")),
      workload_scope = list(mode = "catalog_v1", numeric_moments = character(),
        categorical_marginals = character(), categorical_pairs = list(), correlations = list()),
      numeric_bounds = bounds, capsule_dataset_mapping = list(DA = columns),
      gaussian_specs = list(), synopsis_state_path = file.path(specs[[owner]]$state_dir, "synopsis"))
  }), owners)
  boot <- function(pins = NULL) stats::setNames(lapply(owners, function(owner)
    we_boot_peer(specs[[owner]], server,
      if (is.null(pins)) NULL else pins[setdiff(owners, owner)],
      raw[[owner]], policies[[owner]])), owners)
  peers <- boot()
  identities <- ds.getIdentityPks(we_connections(peers))
  pins <- stats::setNames(unname(unlist(identities)), owners)
  we_close_peers(peers)
  peers <- boot(pins)
  peer_hook(peers)
  if (exists("report_peer_errors", mode = "function", inherits = TRUE)) {
    on.exit(report_peer_errors(peers), add = TRUE, after = FALSE)
  }
  for (peer in peers) peer$worker$run(function() {
    options(dsvert.psi.max_input_ids = 64L)
    TRUE
  })
  conns <- we_connections(peers)
  cat("LASSO", case_label, "PEER_PIDS",
    vapply(peers, function(peer) peer$worker$get_pid(), integer(1L)), "\n")
  ds.psiAlign("D", "patient_id", "DA", datasources = conns, verbose = FALSE)
  if (exists("record_event", mode = "function", inherits = TRUE)) {
    record_event(paste0("K2 ", case_label, " LASSO setup and PSI"), "passed",
      as.numeric(difftime(Sys.time(), started, units = "secs")))
  }
  results <- list()
  for (family in families) {
    registration_started <- Sys.time()
    # Custodian administration changes only the public registered workload. The
    # aligned snapshot, identities, authorization epochs and release ledger stay.
    for (peer in peers) peer$worker$run(function() {
      options(dsvert.dp.gaussian_specs = list())
      TRUE
    })
    maximum <- if (identical(family, "binomial")) 1 else 4
    outcome <- paste0("site_a$", if (combined) paste0("y_", family) else "y")
    registered <- vertcor142_lasso_register(peers, conns, identities, family,
      domain, cohort_id, n, bits, outcome, "site_b$x2", maximum)
    if (exists("record_event", mode = "function", inherits = TRUE)) {
      record_event(paste0("K2 ", family, " LASSO authenticated registration"), "passed",
        as.numeric(difftime(Sys.time(), registration_started, units = "secs")))
    }
    cat("LASSO", family, "AUTHENTICATED_REGISTRATION_COMPLETE", "\n")
    flush.console()
    label <- if (identical(family, "binomial")) "K2 binomial LASSO grid" else
      "K2 Poisson LASSO grid"
    # A passed analysis event includes certificate authentication and result checks;
    # process startup, PSI and custodian registration have their own timed events.
    fit <- run_analysis(label, {
      result <- if (identical(family, "binomial")) {
        ds.vertLASSOIter(stats::as.formula(paste(outcome, "~ site_b$x2")), "DA", family = family,
          analysis_id = "lasso", cross_contract = registered$extension,
          base_contract = registered$base, datasources = conns, verbose = FALSE)
      } else {
        dp_lasso_grid(registered$predictors, registered$extension, registered$base,
          registered$policy, registered$schema, datasources = conns)
      }
      checked <- ds.validateDPGaussianCertificate(result$provenance_certificate,
        trusted_pinset = pins)
      stopifnot(isTRUE(checked$integrity_valid),
        checked$authenticity %in% c("caller_anchored", "session_transport_anchored"),
        identical(result$implementation_state, "cross_owner_exact_gc_materialized"),
        length(result$paths) == 2L, is.null(result$standard_errors))
      saveRDS(result, file.path(out, paste0("smoke-lasso-", family, ".rds")))
      result
    })
    elapsed <- as.numeric(difftime(Sys.time(), registration_started, units = "secs"))
    cat("COMPLETE LASSO", family, "elapsed", elapsed, "\n")
    results[[family]] <- list(family = family, elapsed = elapsed,
      certificate_sha256 = fit$certificate_sha256, paths = length(fit$paths))
  }
  if (combined) stopifnot(!identical(results[[1L]]$certificate_sha256,
    results[[2L]]$certificate_sha256))
  invisible(if (combined) results else results[[1L]])
}
