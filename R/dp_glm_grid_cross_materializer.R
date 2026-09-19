# Signed-grid projection into the existing capsule coordinate/source machinery.
.dsvert_dp_glm_grid_cross_artifacts <- function(manifest) {
  artifacts <- manifest$workload$families$gaussian_models$artifacts
  if (!is.list(artifacts)) return(list())
  artifacts[vapply(artifacts, function(artifact) is.list(artifact) &&
    artifact$version %in% c(unname(.DSVERT_CLIENT_DP_GLM_GRID_CROSS_ARTIFACT_VERSIONS),
      "bounded-lmm-cross-grid-v1"),
    logical(1L))]
}

.dsvert_dp_glm_grid_cross_workload_artifact <- function(contract) {
  spec <- contract$spec
  artifact <- contract$artifact
  artifact$transcript$producer <- artifact$transcript$operation
  artifact$transcript$operation <- NULL
  descriptors <- lapply(spec$predictors, function(x)
    x[c("column", "dataset", "owner_peer", "lower", "upper")])
  categorical <- spec$family %in% c("multinomial", "ordinal")
  extra <- list(
    dataset = spec$dataset, family = spec$family,
    outcome = spec$outcome[c("column", "dataset", "owner_peer",
      if (categorical) "levels" else c("lower", "upper"))],
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
  if (identical(spec$family, "nb")) extra$theta_grid <- spec$theta_grid
  if (categorical) {
    extra$class_count <- spec$class_count
    extra$class_order <- spec$class_order
  }
  c(artifact, extra)
}

.dsvert_dp_glm_grid_cross_source_blocks <- function(artifact, cursor) {
  if (identical(artifact$version, "bounded-lmm-cross-grid-v1")) {
    return(.dsvert_dp_grouped_cross_source_blocks(artifact, cursor))
  }
  spec <- .dsvert_dp_glm_grid_cross_embedded_contract(artifact)$spec
  blocks <- list()
  for (variable in artifact$input_variable_order) {
    outcome <- identical(variable, spec$outcome$reference)
    descriptor <- if (outcome) artifact$outcome else artifact$predictors[[variable]]
    categorical_outcome <- outcome && spec$family %in% c("multinomial", "ordinal")
    if (categorical_outcome) {
      descriptor$lower <- 0
      descriptor$upper <- spec$max_outcome
    }
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
      if (categorical_outcome) blocks[[key]]$levels <- unlist(spec$class_order, use.names = FALSE)
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
  if (identical(artifact$version, "bounded-lmm-cross-grid-v1")) {
    return(.dsvert_dp_grouped_cross_client_artifact(artifact, data_name,
      analysis_id, owner_peer, adjacency, scale, capacity, family))
  }
  contract <- .dsvert_dp_glm_grid_cross_embedded_contract(artifact)
  spec <- contract$spec
  .dsvert_dp_glm_grid_cross_equal(artifact,
    .dsvert_dp_glm_grid_cross_workload_artifact(contract))
  if (!identical(spec$family, family) || !identical(spec$dataset, data_name) ||
      !identical(spec$analysis_id, analysis_id) ||
      (!is.null(owner_peer) && !identical(owner_peer, spec$owner_peer)) ||
      !identical(spec$adjacency, adjacency) ||
      !isTRUE(all.equal(spec$observation_capacity, capacity)) ||
      !isTRUE(all.equal(2^spec$numeric_grid_bits, scale))) {
    .dsvert_dp_glm_grid_cross_fail()
  }
  if (identical(family, "nb")) {
    .dsvert_dp_glm_grid_cross_equal(spec$numeric_contract, .dsvert_dp_nb_grid_cross_numeric())
    .dsvert_dp_glm_grid_cross_equal(spec$sensitivity,
      .dsvert_dp_nb_grid_cross_sensitivity(spec$beta_grid,
        unlist(spec$theta_grid, use.names = FALSE), spec$max_outcome,
        spec$numeric_grid_bits, capacity, adjacency))
  } else if (family %in% c("multinomial", "ordinal")) {
    .dsvert_dp_glm_grid_cross_equal(spec$numeric_contract,
      .dsvert_dp_categorical_grid_cross_numeric(family))
    .dsvert_dp_glm_grid_cross_equal(spec$sensitivity,
      .dsvert_dp_categorical_grid_cross_sensitivity(
        if (identical(family, "ordinal")) spec$candidate_grid else spec$beta_grid,
        family, spec$class_count, length(spec$predictor_order),
        spec$numeric_grid_bits, capacity, adjacency))
  } else {
    if (!identical(spec$numeric_contract$profile_identity, .DSVERT_DP_GLM_GRID_PROFILE_V2)) {
      .dsvert_dp_glm_grid_cross_fail()
    }
    .dsvert_dp_glm_grid_cross_equal(spec, .dsvert_dp_glm_grid_profile_spec(spec))
  }
  artifact$beta_grid <- lapply(spec$beta_grid, function(x) unlist(x, use.names = FALSE))
  artifact$statistic_maximum <- unlist(artifact$statistic_maximum, use.names = FALSE)
  artifact$candidate_loss_bounds <- unlist(artifact$candidate_loss_bounds, use.names = FALSE)
  artifact
}

.dsvert_dp_grouped_cross_workload_artifact <- function(contract) {
  artifact <- .dsvert_dp_glm_grid_cross_workload_artifact(contract)
  spec <- contract$spec
  artifact$transcript$producer <- paste0("dp.", spec$family, "-grid-cross.v1")
  artifact$outcome_encoding <- spec$outcome_encoding
  artifact$predictor_encoding <- spec$predictor_encoding
  artifact$routing_inputs <- spec$routing_inputs
  artifact$grouping <- spec$grouping
  artifact$parameters <- spec$parameters
  artifact$candidate_loss_bounds <- lapply(spec$sensitivity$candidate_bounds,
                                           `[[`, "per_cluster_caps")
  if (identical(spec$family, "lmm")) {
    artifact$source_coordinate_scaling <-
      "all_coordinates_already_on_common_numeric_lattice_v1"
    # Match the manifest reader's scalar/array representation while preserving
    # the complete signed contract as its original opaque JSON string.
    artifact$candidate_loss_bounds <- unlist(artifact$candidate_loss_bounds, use.names = FALSE)
    artifact <- .dsvert_joint_dp_client_canonical(jsonlite::fromJSON(
      .dsvert_joint_dp_client_json(artifact), simplifyVector = TRUE,
      simplifyDataFrame = FALSE, simplifyMatrix = FALSE))
    artifact$participating_peers <- as.list(artifact$participating_peers)
    artifact$computation_peers <- as.list(artifact$computation_peers)
  }
  artifact
}

.dsvert_dp_grouped_cross_source_blocks <- function(artifact, cursor) {
  spec <- .dsvert_dp_glm_grid_cross_embedded_contract(artifact)$spec
  layout <- .dsvert_dp_grouped_grid_cross_layout(spec)
  blocks <- list()
  for (signed_block in layout$blocks) {
    variable <- signed_block$reference
    outcome <- identical(variable, spec$outcome$reference)
    routing <- isTRUE(signed_block$private_routing_input)
    descriptor <- if (routing) spec$routing_inputs[[variable]] else {
      if (outcome) spec$outcome else spec$predictors[[variable]]
    }
    for (kind in c("value", "validity")) {
      size <- signed_block$length
      end <- cursor + size - 1
      if (end > .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_MAX_TRANSPORT_COORDINATES) {
        .dsvert_dp_grouped_grid_cross_fail()
      }
      key <- paste(artifact$analysis_id, variable, kind, sep = "::")
      block <- list(input_family = "grouped_grid", analysis_id = artifact$analysis_id,
        variable = descriptor$column, kind = kind, outcome = outcome,
        dataset = descriptor$dataset, owner_peer = descriptor$owner_peer,
        start = as.integer(cursor), end = as.integer(end), length = as.integer(size),
        fraction_bits = if (kind == "validity") 0L else signed_block$value_fraction_bits,
        maximum = if (kind == "validity") 1 else signed_block$value_maximum,
        private_routing_input = routing)
      if (routing) block$levels <- unlist(descriptor$levels, use.names = FALSE) else {
        block$lower <- descriptor$lower
        block$upper <- descriptor$upper
        if (outcome) block$outcome_encoding <- spec$outcome_encoding
      }
      blocks[[key]] <- block
      cursor <- end + 1
    }
  }
  list(blocks = blocks, cursor = cursor)
}

.dsvert_dp_grouped_cross_client_artifact <- function(artifact, data_name,
    analysis_id, owner_peer, adjacency, scale, capacity, family = "lmm") {
  contract <- .dsvert_dp_glm_grid_cross_embedded_contract(artifact)
  spec <- contract$spec
  if (!identical(family, "lmm") || !identical(spec$family, "lmm") ||
      !identical(spec$parameters$objective, "ml") ||
      !identical(spec$dataset, data_name) || !identical(spec$analysis_id, analysis_id) ||
      (!is.null(owner_peer) && !identical(spec$owner_peer, owner_peer)) ||
      !identical(spec$adjacency, adjacency) ||
      !isTRUE(all.equal(spec$observation_capacity, capacity)) ||
      !isTRUE(all.equal(2^spec$numeric_grid_bits, scale))) {
    .dsvert_dp_grouped_grid_cross_fail()
  }
  parameters <- .dsvert_dp_grouped_grid_cross_ml_parameters(spec$parameters)
  .dsvert_dp_glm_grid_cross_equal(spec$numeric_contract,
    .dsvert_dp_grouped_grid_cross_ml_numeric(spec$grouping, parameters))
  .dsvert_dp_glm_grid_cross_equal(spec$sensitivity,
    .dsvert_dp_grouped_grid_cross_ml_sensitivity(spec$beta_grid,
      spec$numeric_grid_bits, spec$grouping, parameters, adjacency, spec$numeric_contract))
  candidates <- unlist(lapply(seq_along(parameters$variance_grid), function(v) {
    lapply(seq_along(spec$beta_grid), function(b) list(variance_index = v, beta_index = b))
  }), recursive = FALSE)
  .dsvert_dp_glm_grid_cross_equal(spec$candidate_grid, candidates)
  .dsvert_dp_glm_grid_cross_equal(spec$candidate_order, lapply(candidates, function(candidate) {
    .dsvert_dp_capsule_source_hash(list(objective = "ml",
      variance = parameters$variance_grid[[candidate$variance_index]],
      beta = spec$beta_grid[[candidate$beta_index]]))
  }))
  .dsvert_dp_grouped_grid_cross_artifact_validate(contract$artifact, spec)
  .dsvert_dp_glm_grid_cross_equal(contract$source_contract,
    .dsvert_dp_grouped_grid_cross_source_contract(spec, contract$artifact))
  .dsvert_dp_glm_grid_cross_equal(artifact,
    .dsvert_dp_grouped_cross_workload_artifact(contract))
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
    # Snapshot identity must not recursively include its own signatures/hashes.
    # Recognizing a public plan here does not authorize its producer or reader.
    versions <- c(unname(.DSVERT_CLIENT_DP_GLM_GRID_CROSS_SPEC_VERSIONS),
      paste0(.DSVERT_CLIENT_DP_GROUPED_GRID_CROSS_FAMILIES, "_grid_cross_v1"), "cox_grid_cross_v1")
    if (!raw$version %in% versions) next
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
    if (!identical(manifest$workload$capsule_mechanism$mechanism, "discrete-laplace")) {
      .dsvert_dp_glm_grid_cross_fail()
    }
    "dsvert-cross-grid-exact-gc-cost-policy-v2"
  } else .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_COST_POLICY_VERSION
}
