.DSVERT_CLIENT_DP_CAPSULE_PLAN_VERSION <-
  "dsvert-biomedical-capsule-plan-v1"

.dsvert_dp_capsule_plan_invalid <- function(what) {
  stop("The server-authoritative DP capsule plan has an invalid ", what,
       call. = FALSE)
}

.dsvert_dp_capsule_plan_mechanism <- function(manifest, context, scope,
                                               layout) {
  workload <- manifest$workload
  mechanism <- workload$capsule_mechanism
  selection <- workload$mechanism_selection
  lattice <- workload$release_lattice
  mechanism_fields <- c(
    "release_scope", "capability_id", "producer", "purpose",
    "source_context_hash", "mechanism", "mechanism_version", "sampler",
    "sensitivity_norm", "sensitivity", "coordinate_count", "uses_delta",
    "clipping_hash", "ring_bits", "frac_bits")
  lattice_fields <- c(
    "version", "transform_rule", "output_lattice_bits",
    "output_lattice_scale", "natural_l1_sensitivity",
    "integer_l1_sensitivity_steps", "natural_l2_sensitivity",
    "integer_l2_sensitivity_steps")
  selection_fields <- c(
    "version", "selector", "objective", "tie_break", "coordinate_count",
    "epsilon", "allocated_delta", "gaussian_eligible",
    "positive_delta_reserved", "deployed_backends",
    "gaussian_backend_available", "gaussian_unavailable_reason",
    "gaussian_calibration_request", "gaussian_calibration_rounding",
    "gaussian_plan", "gaussian_plan_sha256",
    "laplace_simultaneous_95_abs", "gaussian_simultaneous_95_abs",
    "deployment_rule", "winner", "utility_winner", "decision",
    "canonical_selector_invoked", "selector_certificate_sha256")
  invalid <- function() .dsvert_dp_capsule_plan_invalid(
    "mechanism or calibration contract")
  scalar <- function(value, positive = FALSE) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value >= 0 && (!isTRUE(positive) || value > 0)
  }
  if (!.dsvert_dp_has_exact_names(mechanism, mechanism_fields) ||
      !.dsvert_dp_has_exact_names(lattice, lattice_fields) ||
      !.dsvert_dp_has_exact_names(selection, selection_fields)) {
    invalid()
  }
  reference <- tryCatch(
    context$status[[context$servers[[1L]]]]$policy,
    error = function(error) NULL)
  l1 <- lattice$integer_l1_sensitivity_steps
  l2 <- lattice$integer_l2_sensitivity_steps
  scale <- lattice$output_lattice_scale
  bits <- lattice$output_lattice_bits
  expected_winner <- if (identical(
      mechanism$mechanism,
      .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM)) "gaussian" else "laplace"
  expected_sensitivity <- if (identical(
      mechanism$sensitivity_norm, "l2")) l2 else l1
  profile <- tryCatch(
    .dsvert_vector_profile(mechanism, selection),
    error = function(error) NULL)
  backends <- tryCatch(
    .dsvert_dp_capsule_manifest_string_array(
      selection$deployed_backends, "deployed capsule backends"),
    error = function(error) character())
  rounding_fields <- c(
    "declared_policy_values", "epsilon", "delta", "l2_sensitivity")
  rounding <- selection$gaussian_calibration_rounding
  rounding_valid <- .dsvert_dp_has_exact_names(rounding, rounding_fields) &&
    all(vapply(rounding, .dsvert_dp_is_string, logical(1L)))
  plan_valid <- if (is.null(selection$gaussian_plan)) {
    is.null(selection$gaussian_plan_sha256)
  } else {
    is.list(selection$gaussian_plan) &&
      .dsvert_dp_capsule_source_hex(selection$gaussian_plan_sha256) &&
      identical(selection$gaussian_plan_sha256,
                .dsvert_dp_capsule_source_hash(selection$gaussian_plan))
  }
  gaussian_available <- is.logical(selection$gaussian_backend_available) &&
    length(selection$gaussian_backend_available) == 1L &&
    !is.na(selection$gaussian_backend_available) &&
    isTRUE(selection$gaussian_backend_available)
  expected_backends <- if (gaussian_available) c(
    "joint-discrete-laplace-v3", "dyadic-discrete-gaussian-tv-bounded-v2")
  else "joint-discrete-laplace-v3"
  gaussian_availability_valid <- if (gaussian_available) {
    is.null(selection$gaussian_unavailable_reason) &&
      is.list(selection$gaussian_calibration_request) &&
      is.list(selection$gaussian_plan) &&
      scalar(selection$gaussian_simultaneous_95_abs)
  } else {
    .dsvert_dp_is_string(selection$gaussian_unavailable_reason) &&
      is.null(selection$gaussian_simultaneous_95_abs)
  }
  utility_valid <- scalar(selection$laplace_simultaneous_95_abs) &&
    if (gaussian_available) {
      if (identical(expected_winner, "gaussian")) {
        selection$gaussian_simultaneous_95_abs <
          selection$laplace_simultaneous_95_abs
      } else {
        selection$gaussian_simultaneous_95_abs >=
          selection$laplace_simultaneous_95_abs
      }
    } else {
      identical(expected_winner, "laplace")
    }
  valid <- is.list(reference) && is.list(profile) &&
    identical(mechanism$release_scope, "joint_mpc_single_opening") &&
    identical(mechanism$capability_id, "joint_mpc_single_opening_v1") &&
    .dsvert_dp_capsule_source_hex(mechanism$source_context_hash) &&
    .dsvert_dp_capsule_source_hex(mechanism$clipping_hash) &&
    mechanism$sensitivity_norm %in% c("l1", "l2") &&
    scalar(mechanism$sensitivity, TRUE) &&
    identical(as.numeric(mechanism$coordinate_count),
              as.numeric(layout$coordinate_count)) &&
    is.logical(mechanism$uses_delta) &&
    length(mechanism$uses_delta) == 1L && !is.na(mechanism$uses_delta) &&
    .dsvert_dp_is_integer(mechanism$ring_bits, 128, 128) &&
    .dsvert_dp_is_integer(mechanism$frac_bits, 0, 0) &&
    identical(lattice$version, "biomedical-capsule-common-lattice-v1") &&
    identical(
      lattice$transform_rule,
      "raw_coordinate_left_shift_to_common_numeric_grid_v1") &&
    .dsvert_dp_is_integer(bits, 8, 18) && scalar(scale, TRUE) &&
    scalar(lattice$natural_l1_sensitivity, TRUE) && scalar(l1, TRUE) &&
    scalar(lattice$natural_l2_sensitivity, TRUE) && scalar(l2, TRUE) &&
    .dsvert_dp_num_equal(scale, 2^as.numeric(bits), 128) &&
    .dsvert_dp_num_equal(
      l1, lattice$natural_l1_sensitivity * scale, 128) &&
    .dsvert_dp_num_equal(
      l2, lattice$natural_l2_sensitivity * scale, 128) &&
    .dsvert_dp_num_equal(
      l1, scope$projected_cost$projected_integer_l1_sensitivity, 128) &&
    .dsvert_dp_num_equal(
      l2, scope$projected_cost$projected_integer_l2_sensitivity, 128) &&
    .dsvert_dp_num_equal(mechanism$sensitivity,
                         expected_sensitivity, 128) &&
    identical(selection$version,
              "biomedical-capsule-noise-selection-v4") &&
    identical(
      selection$selector,
      "formal_fixed_work_backend_minimum_simultaneous_95_radius_v2") &&
    identical(selection$objective, "simultaneous_95_abs") &&
    identical(
      selection$tie_break,
      "laplace_unless_fixed_work_gaussian_strictly_improves") &&
    identical(as.numeric(selection$coordinate_count),
              as.numeric(layout$coordinate_count)) &&
    scalar(selection$epsilon, TRUE) &&
    scalar(selection$allocated_delta) && selection$allocated_delta < 1 &&
    .dsvert_dp_num_equal(
      selection$epsilon, reference$capsule_epsilon, 128) &&
    .dsvert_dp_num_equal(
      selection$allocated_delta, reference$capsule_delta, 128) &&
    selection$winner %in% c("laplace", "gaussian") &&
    identical(selection$winner, selection$utility_winner) &&
    identical(selection$winner, expected_winner) &&
    identical(
      mechanism$sensitivity_norm,
      if (identical(expected_winner, "gaussian")) "l2" else "l1") &&
    identical(mechanism$uses_delta, selection$allocated_delta > 0) &&
    is.logical(selection$gaussian_eligible) &&
    identical(selection$gaussian_eligible,
              selection$allocated_delta > 0) &&
    is.logical(selection$positive_delta_reserved) &&
    identical(selection$positive_delta_reserved,
              selection$allocated_delta > 0) &&
    identical(backends, expected_backends) && rounding_valid && plan_valid &&
    gaussian_availability_valid && utility_valid &&
    (if (selection$allocated_delta > 0) {
      is.list(selection$gaussian_calibration_request)
    } else {
      is.null(selection$gaussian_calibration_request)
    }) &&
    .dsvert_dp_is_string(selection$decision) &&
    identical(selection$deployment_rule,
              "formal_backend_or_explicit_laplace_fallback") &&
    is.logical(selection$canonical_selector_invoked) &&
    identical(selection$canonical_selector_invoked, TRUE) &&
    .dsvert_dp_capsule_source_hex(selection$selector_certificate_sha256)
  if (!isTRUE(valid)) invalid()

  calibration_fields <- c(
    "version", "selector", "objective", "tie_break", "winner",
    "decision", "coordinate_count", "epsilon", "allocated_delta",
    "gaussian_eligible", "positive_delta_reserved", "deployed_backends",
    "gaussian_backend_available", "gaussian_unavailable_reason",
    "laplace_simultaneous_95_abs", "gaussian_simultaneous_95_abs",
    "gaussian_plan_sha256", "deployment_rule",
    "canonical_selector_invoked", "selector_certificate_sha256")
  calibration <- stats::setNames(lapply(
    calibration_fields, function(field) selection[[field]]),
    calibration_fields)
  list(
    mechanism = list(
      mechanism = mechanism$mechanism,
      mechanism_version = mechanism$mechanism_version,
      sampler = mechanism$sampler,
      sensitivity_norm = mechanism$sensitivity_norm,
      sensitivity = as.numeric(mechanism$sensitivity),
      uses_delta = mechanism$uses_delta,
      ring_bits = as.integer(mechanism$ring_bits),
      frac_bits = as.integer(mechanism$frac_bits),
      release_lattice = lattice),
    calibration = calibration)
}

.dsvert_dp_capsule_plan_families <- function(manifest, layout) {
  family_names <- c(
    "admitted_count", "numeric_moments", "numeric_pair_moments",
    "gaussian_models", "fixed_numeric_histograms",
    "categorical_marginals", "categorical_pairs",
    "correlation_artifacts", "describe_artifacts", "survival_artifacts")
  blocks <- layout$blocks
  block_family <- if (length(blocks)) vapply(
    blocks, `[[`, character(1L), "family") else character()
  block_length <- if (length(blocks)) vapply(
    blocks, `[[`, numeric(1L), "length") else numeric()
  metadata_counts <- c(
    correlation_artifacts = length(
      manifest$workload$families$correlation_artifacts),
    describe_artifacts = length(
      manifest$workload$families$describe_artifacts))
  stats::setNames(lapply(family_names, function(family) {
    direct <- which(block_family == family)
    count <- if (family %in% names(metadata_counts)) {
      unname(metadata_counts[[family]])
    } else {
      length(direct)
    }
    list(
      artifact_count = as.integer(count),
      coordinate_count = as.integer(sum(block_length[direct])))
  }), family_names)
}

.dsvert_dp_capsule_plan_cost_validate <- function(manifest, scope, layout) {
  cost <- scope$projected_cost
  families <- manifest$workload$families
  sensitivity <- manifest$workload$sensitivity
  blocks <- layout$blocks
  coordinate_count <- function(family) sum(vapply(
    blocks[vapply(blocks, function(block) identical(block$family, family),
                  logical(1L))],
    `[[`, numeric(1L), "length"))
  expected <- c(
    numeric_moment_coordinate_count = coordinate_count("numeric_moments"),
    numeric_pair_coordinate_count =
      coordinate_count("numeric_pair_moments"),
    categorical_marginal_coordinate_count =
      coordinate_count("categorical_marginals"),
    categorical_pair_coordinate_count =
      coordinate_count("categorical_pairs"),
    gaussian_model_coordinate_count = coordinate_count("gaussian_models"),
    projected_coordinate_count = as.numeric(layout$coordinate_count))
  valid <- is.list(sensitivity) &&
    all(vapply(names(expected), function(field) {
      identical(as.numeric(cost[[field]]), as.numeric(expected[[field]]))
    }, logical(1L))) &&
    identical(
      as.numeric(cost$included_cross_categorical_pair_count),
      as.numeric(length(families$categorical_pairs$cross_artifacts))) &&
    .dsvert_dp_num_equal(
      cost$projected_integer_l1_sensitivity, sensitivity$l1, 128) &&
    .dsvert_dp_num_equal(
      cost$projected_integer_l2_sensitivity, sensitivity$l2, 128)
  if (!isTRUE(valid)) {
    .dsvert_dp_capsule_plan_invalid("projected workload cost")
  }
  invisible(TRUE)
}

#' Dry-run one server-authoritative DP capsule
#'
#' Performs only the reusable-capsule status handshake and the three signed
#' manifest phases (draft, global-schema signature, and byte-identical build).
#' The selection is owned by custodian policy and signed workload
#' specifications: this function has no analyst-controlled expansion input.
#' It validates the current manifest, primitive-scope, coordinate-layout, and
#' mechanism contracts, but never resolves a protected snapshot, materialises
#' coordinates, invokes a producer, or creates a DP release.
#'
#' @param datasources Complete named DataSHIELD connection set. If `NULL`, use
#'   the active connections.
#' @param status Optional result from `ds.vertDPStatus()`. It is validated
#'   again against `datasources`; if omitted, the handshake is performed.
#' @return A compact `ds.vertDPCapsulePlan` object containing the signed
#'   capsule identity, immutable primitive selection, projected coordinate and
#'   sensitivity cost, family inventory, mechanism calibration, pinset, and
#'   explicit zero-access/zero-release guarantees for this dry-run.
#' @export
ds.vertDPCapsulePlan <- function(datasources = NULL, status = NULL) {
  .dsvert_dp_capsule_plan_impl(
    datasources, status = status, .aggregate = DSI::datashield.aggregate)
}

.dsvert_dp_capsule_plan_impl <- function(
    datasources = NULL, status = NULL,
    .aggregate = DSI::datashield.aggregate) {
  datasources <- .dsvert_dp_datasources(datasources)
  bundle <- .dsvert_dp_capsule_manifest_build(
    datasources, status = status, .aggregate = .aggregate)
  if (!is.list(bundle) || !is.list(bundle$context) ||
      !.dsvert_dp_is_string(bundle$manifest_json)) {
    .dsvert_dp_capsule_plan_invalid("manifest bundle")
  }
  validated <- .dsvert_dp_capsule_source_manifest(
    bundle$manifest_json, bundle$context)
  manifest <- if (is.list(validated)) validated$value else NULL
  if (!is.list(manifest)) {
    .dsvert_dp_capsule_plan_invalid("manifest")
  }
  scope <- .dsvert_dp_capsule_primitive_scope_validate(manifest)
  layout <- .dsvert_dp_capsule_vector_layout(manifest)
  context <- bundle$context
  servers <- context$servers
  reference <- tryCatch(
    context$status[[servers[[1L]]]]$policy,
    error = function(error) NULL)
  bundle_valid <- is.character(servers) && length(servers) >= 2L &&
    !anyNA(servers) && !anyDuplicated(servers) &&
    identical(servers, sort(servers, method = "radix")) &&
    is.character(context$pinset) &&
    identical(names(context$pinset), servers) &&
    is.character(context$designated) && length(context$designated) == 2L &&
    identical(context$designated,
              sort(context$designated, method = "radix")) &&
    all(context$designated %in% servers) && is.list(reference) &&
    .dsvert_dp_capsule_source_hex(reference$peer_pinset_sha256) &&
    identical(.dsvert_dp_pinset_hash(context$pinset),
              reference$peer_pinset_sha256) &&
    .dsvert_dp_capsule_source_hex(bundle$manifest_sha256) &&
    identical(bundle$manifest_sha256, digest::digest(
      bundle$manifest_json, "sha256", serialize = FALSE)) &&
    .dsvert_dp_capsule_source_hex(bundle$capsule_id) &&
    identical(validated$capsule_id, bundle$capsule_id) &&
    identical(.dsvert_joint_dp_client_canonical(bundle$logical_snapshot),
              .dsvert_joint_dp_client_canonical(manifest$logical_snapshot)) &&
    .dsvert_dp_is_integer(bundle$artifact_commitment_count, 0,
                          .DSVERT_DP_MAX_COORDINATES) &&
    .dsvert_dp_capsule_source_hex(bundle$artifact_commitments_root) &&
    identical(bundle$deterministic_replay, TRUE) &&
    identical(bundle$operation_or_request_limit, FALSE) &&
    identical(bundle$history_can_deny_operation, FALSE)
  if (!isTRUE(bundle_valid)) {
    .dsvert_dp_capsule_plan_invalid("manifest bundle")
  }
  .dsvert_dp_capsule_plan_cost_validate(manifest, scope, layout)
  mechanism <- .dsvert_dp_capsule_plan_mechanism(
    manifest, context, scope, layout)
  families <- .dsvert_dp_capsule_plan_families(manifest, layout)
  primitive_fields <- c(
    "version", "mode", "authority", "analyst_expandable",
    "client_query_can_add_coordinates", "consensus", "mismatch_behavior",
    "compatibility_default", "recommended_deployment_mode")
  primitive_scope <- scope[primitive_fields]
  explicit <- scope$selection$explicit_catalog
  references <- scope$selection$referenced_by_signed_specs
  included <- scope$selection$included
  selection <- list(
    sha256 = scope$selection_sha256,
    mode = scope$selection$mode,
    explicit_catalog_counts = list(
      numeric_moments = as.integer(length(explicit$numeric_moments)),
      categorical_marginals =
        as.integer(length(explicit$categorical_marginals)),
      categorical_pairs = as.integer(length(explicit$categorical_pairs)),
      correlations = as.integer(length(explicit$correlations))),
    signed_spec_reference_counts = list(
      numeric = as.integer(length(references$numeric)),
      categorical = as.integer(length(references$categorical)),
      describe = as.integer(length(references$describe)),
      survival = as.integer(length(references$survival)),
      gaussian = as.integer(length(references$gaussian)),
      vertical_cross = as.integer(length(references$vertical_cross))),
    included_counts = list(
      numeric_moments = as.integer(length(included$numeric_moments)),
      categorical_marginals =
        as.integer(length(included$categorical_marginals)),
      same_owner_categorical_pairs =
        as.integer(length(included$same_owner_categorical_pairs)),
      same_owner_correlations =
        as.integer(length(included$same_owner_correlations))))
  result <- list(
    version = .DSVERT_CLIENT_DP_CAPSULE_PLAN_VERSION,
    capsule = list(
      capsule_id = bundle$capsule_id,
      manifest_sha256 = bundle$manifest_sha256,
      logical_snapshot = manifest$logical_snapshot),
    primitive_scope = primitive_scope,
    selection = selection,
    projected_cost = scope$projected_cost,
    artifacts = list(
      commitment_count = as.integer(bundle$artifact_commitment_count),
      commitments_root = bundle$artifact_commitments_root,
      coordinate_layout = list(
        version = layout$version,
        sha256 = layout$sha256,
        coordinate_count = as.integer(layout$coordinate_count))),
    families = families,
    mechanism = mechanism$mechanism,
    calibration = mechanism$calibration,
    consortium = list(
      K = as.integer(length(servers)), peers = servers,
      designated_noise_peers = context$designated,
      peer_pinset_sha256 = reference$peer_pinset_sha256),
    guarantees = list(
      data_access = FALSE, release_created = FALSE,
      privacy_cost = c(epsilon = 0, delta = 0),
      operation_limit = FALSE, request_limit = FALSE,
      history_can_deny_operation = FALSE))
  class(result) <- c("ds.vertDPCapsulePlan", "list")
  result
}

#' @export
print.ds.vertDPCapsulePlan <- function(x, ...) {
  cat("dsVert DP capsule plan (dry-run; no data access; no release)\n")
  cat("capsule:", x$capsule$capsule_id, "\n")
  cat("manifest:", x$capsule$manifest_sha256, "\n")
  cat("scope:", x$primitive_scope$mode, "| coordinates:",
      x$artifacts$coordinate_layout$coordinate_count, "\n")
  cat("mechanism:", x$mechanism$mechanism, "| epsilon:",
      format(x$calibration$epsilon), "| delta:",
      format(x$calibration$allocated_delta), "\n")
  cat("peers:", x$consortium$K, "| pinset:",
      x$consortium$peer_pinset_sha256, "\n")
  cat("planning privacy cost: epsilon=0, delta=0; limits/history gates: false\n")
  invisible(x)
}
