# Signed-grid projection into the existing capsule coordinate/source machinery.
.dsvert_dp_glm_grid_cross_artifacts <- function(manifest) {
  artifacts <- manifest$workload$families$gaussian_models$artifacts
  if (!is.list(artifacts)) return(list())
  artifacts[vapply(artifacts, function(artifact) is.list(artifact) &&
    artifact$version %in% unname(.DSVERT_CLIENT_DP_GLM_GRID_CROSS_ARTIFACT_VERSIONS),
    logical(1L))]
}

.dsvert_dp_glm_grid_cross_workload_artifact <- function(contract) {
  spec <- contract$spec
  artifact <- contract$artifact
  artifact$transcript$producer <- artifact$transcript$operation
  artifact$transcript$operation <- NULL
  descriptors <- lapply(spec$predictors, function(x)
    x[c("column", "dataset", "owner_peer", "lower", "upper")])
  extra <- list(
    dataset = spec$dataset, family = spec$family,
    outcome = spec$outcome[c("column", "dataset", "owner_peer", "lower", "upper")],
    predictors = descriptors,
    predictor_order = unlist(spec$predictor_order, use.names = FALSE),
    input_variable_order = unlist(spec$input_variable_order, use.names = FALSE),
    design_terms = unlist(spec$design_terms, use.names = FALSE),
    beta_grid = spec$beta_grid, intercept = TRUE,
    observation_capacity = spec$observation_capacity,
    max_outcome = spec$max_outcome,
    statistic_maximum = spec$sensitivity$maximum_coordinates,
    source_raw_l1_sensitivity = spec$sensitivity$raw_l1_sensitivity,
    source_raw_l2_sensitivity = spec$sensitivity$raw_l2_sensitivity,
    natural_l1_sensitivity = spec$sensitivity$natural_l1_sensitivity,
    natural_l2_sensitivity = spec$sensitivity$natural_l2_sensitivity,
    candidate_loss_bounds = lapply(spec$sensitivity$candidate_bounds,
                                   `[[`, "loss_bound"),
    numeric_certificate = spec$numeric_contract,
    adjacency = spec$adjacency,
    estimation_scope = "certified_bounded_cross_owner_signed_grid_only_v2",
    signed_contract = .dsvert_joint_dp_client_json(.dsvert_joint_dp_client_canonical(contract)))
  c(artifact, extra)
}

.dsvert_dp_glm_grid_cross_source_blocks <- function(artifact, cursor) {
  spec <- .dsvert_dp_glm_grid_cross_embedded_contract(artifact)$spec
  blocks <- list()
  for (variable in artifact$input_variable_order) {
    outcome <- identical(variable, spec$outcome$reference)
    descriptor <- if (outcome) artifact$outcome else artifact$predictors[[variable]]
    for (kind in c("value", "validity")) {
      size <- spec$observation_capacity
      end <- cursor + size - 1
      if (end > .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_MAX_TRANSPORT_COORDINATES) {
        .dsvert_dp_glm_grid_cross_fail()
      }
      key <- paste(artifact$analysis_id, variable, kind, sep = "::")
      blocks[[key]] <- c(list(input_family = "glm_grid",
        analysis_id = artifact$analysis_id, variable = descriptor$column,
        kind = kind, outcome = outcome), descriptor[c(
          "dataset", "owner_peer", "lower", "upper")], list(
        start = as.integer(cursor), end = as.integer(end), length = as.integer(size),
        fraction_bits = if (identical(kind, "validity") || outcome) 0L else 50L,
        maximum = if (identical(kind, "validity")) 1 else
          if (outcome) spec$max_outcome else 2^50))
      cursor <- end + 1
    }
  }
  list(blocks = blocks, cursor = cursor)
}

.dsvert_dp_glm_grid_cross_embedded_contract <- function(artifact) {
  if (!is.character(artifact$signed_contract) || length(artifact$signed_contract) != 1L ||
      is.na(artifact$signed_contract)) .dsvert_dp_glm_grid_cross_fail()
  .dsvert_joint_dp_client_decode(artifact$signed_contract,
    "signed cross-grid contract", .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_MANIFEST_BYTES)
}

# The signed schema and signatures are checked by preflight and certificate
# reconstruction. Rebuild all public arithmetic fields again at use time.
.dsvert_dp_glm_grid_cross_client_artifact <- function(artifact, data_name,
    analysis_id, owner_peer, adjacency, scale, capacity, family) {
  contract <- .dsvert_dp_glm_grid_cross_embedded_contract(artifact)
  spec <- contract$spec
  .dsvert_dp_glm_grid_cross_equal(artifact,
    .dsvert_dp_glm_grid_cross_workload_artifact(contract))
  if (!identical(spec$family, family) || !identical(spec$dataset, data_name) ||
      !identical(spec$analysis_id, analysis_id) ||
      (!is.null(owner_peer) && !identical(owner_peer, spec$owner_peer)) ||
      !identical(spec$adjacency, adjacency) ||
      !isTRUE(all.equal(spec$observation_capacity, capacity)) ||
      !isTRUE(all.equal(2^spec$numeric_grid_bits, scale)) ||
      !identical(spec$numeric_contract$profile_identity, .DSVERT_DP_GLM_GRID_PROFILE_V2)) {
    .dsvert_dp_glm_grid_cross_fail()
  }
  .dsvert_dp_glm_grid_cross_equal(spec, .dsvert_dp_glm_grid_profile_spec(spec))
  artifact$beta_grid <- lapply(spec$beta_grid, function(x) unlist(x, use.names = FALSE))
  artifact$statistic_maximum <- unlist(artifact$statistic_maximum, use.names = FALSE)
  artifact$candidate_loss_bounds <- unlist(artifact$candidate_loss_bounds, use.names = FALSE)
  artifact
}

.dsvert_dp_glm_grid_formula_reference <- function(value) {
  if (is.symbol(value)) return(as.character(value))
  if (is.call(value) && length(value) == 3L && identical(value[[1L]], as.name("$")) &&
      is.symbol(value[[2L]]) && is.symbol(value[[3L]])) {
    return(paste(as.character(value[[2L]]), as.character(value[[3L]]), sep = "$"))
  }
  NULL
}

# Snapshot identity binds the public plan, excluding only fields derived from
# that identity and signatures. The complete contract remains bound separately
# by workload/admission hashes and unanimous signatures.
.dsvert_dp_glm_grid_cross_snapshot_workload <- function(workload) {
  for (id in names(workload$gaussian)) {
    raw <- workload$gaussian[[id]]$spec
    if (!raw$version %in% unname(.DSVERT_CLIENT_DP_GLM_GRID_CROSS_SPEC_VERSIONS)) next
    contract <- .dsvert_dp_glm_grid_cross_raw_contract(raw)
    plan <- contract$spec
    plan$schema_sha256 <- NULL
    plan$logical_snapshot <- NULL
    plan$alignment$public_alignment_contract_sha256 <- NULL
    workload$gaussian[[id]]$spec <- list(version = raw$version,
      dataset = raw$dataset, public_candidate_plan = plan)
  }
  workload
}

.dsvert_dp_glm_grid_cross_raw_contract <- function(raw) {
  .dsvert_dp_glm_grid_cross_fields(raw, c("version", "dataset", "contract"))
  value <- raw$contract
  if (is.character(value) && length(value) == 1L && !is.na(value)) {
    value <- .dsvert_joint_dp_client_decode(value, "signed cross-grid contract",
      .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_MANIFEST_BYTES)
  }
  if (!is.list(value) || !identical(value$spec$version, raw$version) ||
      !identical(value$spec$dataset, raw$dataset)) .dsvert_dp_glm_grid_cross_fail()
  value
}

.dsvert_dp_glm_grid_formula_check <- function(formula, artifact) {
  outcome <- .dsvert_dp_glm_grid_formula_reference(formula[[2L]])
  predictors <- attr(stats::terms(formula), "term.labels")
  expected_outcome <- if (artifact$version %in%
      unname(.DSVERT_CLIENT_DP_GLM_GRID_CROSS_ARTIFACT_VERSIONS)) {
    .dsvert_dp_glm_grid_cross_embedded_contract(artifact)$spec$outcome$reference
  } else artifact$outcome$column
  if (!identical(expected_outcome, outcome) ||
      !setequal(artifact$predictor_order, predictors)) {
    stop("formula must match the signed finite GLM grid artifact", call. = FALSE)
  }
  invisible(TRUE)
}

.dsvert_dp_glm_grid_cross_noise_policy <- function(manifest) {
  if (length(.dsvert_dp_glm_grid_cross_artifacts(manifest))) {
    "dsvert-cross-grid-exact-gc-cost-policy-v2"
  } else .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_COST_POLICY_VERSION
}
