# Internal consumer for the formal GLM joint-DP inferential artifact.
#
# The artifact is evaluated at an already released DP coefficient vector.  It
# does not reconstruct protected fitted values, scores or patient rows.  The
# complete stacked vector is opened once, after hidden validity has been
# consumed and independent full peer noise has been applied.  This file is not
# exported and deliberately contains no DSI request path.

.DSVERT_FORMAL_GLM_INFERENCE_VERSION <-
  "dsvert-formal-glm-joint-inference-release-v1"
.DSVERT_FORMAL_GLM_INFERENCE_DOMAIN <-
  "dsVert/formal-glm/joint-inference-release/v1|"
.DSVERT_FORMAL_GLM_INFERENCE_COORDINATES <-
  "dsvert-formal-glm-joint-inference-coordinate-contract-v1"
.DSVERT_FORMAL_GLM_INFERENCE_BASE_BINDING <-
  "dsvert-formal-glm-joint-inference-base-binding-v1"
.DSVERT_FORMAL_GLM_INFERENCE_PROTECTED_BINDING <-
  "dsvert-formal-glm-joint-inference-protected-binding-v1"
.DSVERT_FORMAL_GLM_INFERENCE_MAX_RELEASE_BYTES <- 8L * 1024L^2
.DSVERT_FORMAL_GLM_INFERENCE_RELEASED <- "released"
.DSVERT_FORMAL_GLM_INFERENCE_NOT_PROMOTED <-
  "not_promoted_sensitivity_not_certified"

.dsvert_formal_glm_inference_unavailable <- function(
    statistic, reason, release_sha256 = NULL) {
  structure(list(
    statistic = statistic,
    status = "joint_dp_inference_capability_unavailable",
    reason = reason,
    additional_server_calls = 0L,
    fallback_used = FALSE,
    public_release_sha256 = release_sha256),
    class = c("dsvert_formal_glm_inference_unavailable", "list"))
}

.dsvert_formal_glm_inference_hash <- function(value) {
  digest::digest(.dsvert_joint_dp_client_json(value), algo = "sha256",
                 serialize = FALSE)
}

.dsvert_formal_glm_inference_triangle <- function(coefficients, block) {
  result <- character()
  p <- length(coefficients)
  for (row in seq_len(p)) {
    for (column in seq.int(row, p)) {
      result <- c(result, paste0(
        block, "[", row, ",", column, "]:",
        coefficients[[row]], "|", coefficients[[column]]))
    }
  }
  result
}

.dsvert_formal_glm_inference_capability_fields <- function() {
  c("information_upper_triangle", "score_meat_upper_triangle",
    "score_vector", "canonical_bounded_log_likelihood_at_dp_beta",
    "integrated_pwl_surrogate_loss_at_dp_beta", "admitted_n")
}

.dsvert_formal_glm_inference_validate_capabilities <- function(value) {
  fields <- .dsvert_formal_glm_inference_capability_fields()
  .dsvert_formal_glm_release_exact_names(
    value, fields, "inference block capabilities")
  if (!identical(value$information_upper_triangle,
                 .DSVERT_FORMAL_GLM_INFERENCE_RELEASED)) {
    .dsvert_formal_glm_release_stop(
      "The formal inference artifact must contain its information matrix",
      "inference_information_missing")
  }
  allowed <- c(.DSVERT_FORMAL_GLM_INFERENCE_RELEASED,
               .DSVERT_FORMAL_GLM_INFERENCE_NOT_PROMOTED)
  valid <- vapply(value, function(item) {
    is.character(item) && length(item) == 1L && !is.na(item) &&
      item %in% allowed
  }, logical(1L))
  if (any(!valid)) {
    .dsvert_formal_glm_release_stop(
      "The formal inference artifact promoted an unknown capability",
      "inference_capability_invalid")
  }
  value
}

.dsvert_formal_glm_inference_coordinate_layout <- function(
    coefficients, capabilities) {
  capabilities <- .dsvert_formal_glm_inference_validate_capabilities(
    capabilities)
  blocks <- list(
    information_upper_triangle =
      .dsvert_formal_glm_inference_triangle(coefficients, "H"),
    score_meat_upper_triangle =
      .dsvert_formal_glm_inference_triangle(coefficients, "J"),
    score_vector = paste0("score[", seq_along(coefficients), "]:",
                          coefficients),
    canonical_bounded_log_likelihood_at_dp_beta =
      "canonical_bounded_log_likelihood_at_dp_beta",
    integrated_pwl_surrogate_loss_at_dp_beta =
      "integrated_pwl_surrogate_loss_at_dp_beta",
    admitted_n = "admitted_n")
  unname(unlist(blocks[vapply(names(blocks), function(field) {
    identical(capabilities[[field]], .DSVERT_FORMAL_GLM_INFERENCE_RELEASED)
  }, logical(1L))], use.names = FALSE))
}

.dsvert_formal_glm_inference_integer_bound <- function(value, scale,
                                                       direction) {
  value <- .dsvert_glm_rat_mul(value, scale)
  quotient <- value$numerator %/% value$denominator
  remainder <- value$numerator %% value$denominator
  increment <- remainder != .dsvert_glm_bn_zero()
  if (identical(direction, "floor")) {
    magnitude <- if (value$sign < 0L && increment) {
      quotient + .dsvert_glm_bn_one()
    } else quotient
  } else {
    magnitude <- if (value$sign > 0L && increment) {
      quotient + .dsvert_glm_bn_one()
    } else quotient
  }
  if (value$sign < 0L && magnitude != .dsvert_glm_bn_zero()) {
    paste0("-", as.character(magnitude))
  } else {
    as.character(magnitude)
  }
}

.dsvert_formal_glm_inference_rational_sum <- function(values) {
  Reduce(.dsvert_glm_rat_add, values, init = .dsvert_glm_rat("0"))
}

.dsvert_formal_glm_inference_server_plan <- function(
    point_result, capabilities, output_lattice_scale) {
  .dsvert_formal_glm_inference_validate_point(point_result)
  capabilities <- .dsvert_formal_glm_inference_validate_capabilities(
    capabilities)
  scale_text <- .dsvert_formal_glm_release_integer_text(
    output_lattice_scale, "inference plan lattice scale", positive = TRUE)
  scale <- .dsvert_glm_rat(scale_text)
  schema <- point_result$schema
  coefficient_names <- unlist(
    schema$estimand$coefficient_order, use.names = FALSE)
  coefficients <- schema$estimand$coefficients
  p <- length(coefficients)
  capacity <- .dsvert_glm_rat(
    as.character(schema$estimand$capacity_normalizer))
  gamma <- .dsvert_glm_rat(
    as.character(schema$adjacency$triangle_factor_gamma))
  beta_scale <- .dsvert_glm_rat(point_result$output_lattice_scale)
  beta <- lapply(point_result$coefficient_lattice_steps, function(value) {
    .dsvert_glm_rat_div(value, beta_scale)
  })

  x_lower <- x_upper <- x_abs <- vector("list", p)
  for (index in seq_len(p)) {
    term <- coefficients[[index]]$term
    if (identical(term$kind, "intercept")) {
      x_lower[[index]] <- x_upper[[index]] <- x_abs[[index]] <-
        .dsvert_glm_rat("1")
    } else if (identical(term$kind, "numeric")) {
      x_lower[[index]] <- .dsvert_glm_rat(term$clipping_lower)
      x_upper[[index]] <- .dsvert_glm_rat(term$clipping_upper)
      x_abs[[index]] <- .dsvert_glm_rat_max(
        .dsvert_glm_rat_abs(x_lower[[index]]),
        .dsvert_glm_rat_abs(x_upper[[index]]))
    } else if (identical(term$kind, "categorical_indicator")) {
      x_lower[[index]] <- .dsvert_glm_rat("0")
      x_upper[[index]] <- x_abs[[index]] <- .dsvert_glm_rat("1")
    } else {
      .dsvert_formal_glm_release_stop(
        "The inference plan found an unsupported signed design term",
        "inference_sensitivity_derivation_invalid")
    }
  }
  offset <- schema$estimand$offset
  eta_lower <- .dsvert_glm_rat(offset$lower)
  eta_upper <- .dsvert_glm_rat(offset$upper)
  for (index in seq_len(p)) {
    endpoints <- list(
      .dsvert_glm_rat_mul(beta[[index]], x_lower[[index]]),
      .dsvert_glm_rat_mul(beta[[index]], x_upper[[index]]))
    eta_lower <- .dsvert_glm_rat_add(
      eta_lower, Reduce(.dsvert_glm_rat_min, endpoints))
    eta_upper <- .dsvert_glm_rat_add(
      eta_upper, Reduce(.dsvert_glm_rat_max, endpoints))
  }
  eta_abs <- .dsvert_glm_rat_max(
    .dsvert_glm_rat_abs(eta_lower), .dsvert_glm_rat_abs(eta_upper))
  mean_abs <- Reduce(.dsvert_glm_rat_max, lapply(
    schema$link_surrogate$values, function(value) {
      .dsvert_glm_rat_abs(value)
    }))
  slope_abs <- Reduce(.dsvert_glm_rat_max, lapply(
    schema$link_surrogate$slopes, function(value) {
      .dsvert_glm_rat_abs(value)
    }))
  response_upper <- if (identical(schema$estimand$family, "binomial")) {
    .dsvert_glm_rat("1")
  } else {
    .dsvert_glm_rat(schema$estimand$column_registry[[
      schema$estimand$response]]$upper)
  }
  residual_abs <- .dsvert_glm_rat_add(response_upper, mean_abs)
  weight_abs <- .dsvert_glm_rat(
    schema$estimand$weights$maximum_patient_weight)
  weight_squared <- .dsvert_glm_rat_pow(weight_abs, 2L)
  residual_squared <- .dsvert_glm_rat_pow(residual_abs, 2L)

  h_patient <- j_patient <- list()
  for (row in seq_len(p)) {
    for (column in seq.int(row, p)) {
      x_product <- .dsvert_glm_rat_mul(
        x_abs[[row]], x_abs[[column]])
      h_patient[[length(h_patient) + 1L]] <- .dsvert_glm_rat_mul(
        .dsvert_glm_rat_mul(weight_abs, slope_abs), x_product)
      j_patient[[length(j_patient) + 1L]] <- .dsvert_glm_rat_mul(
        .dsvert_glm_rat_mul(weight_squared, residual_squared), x_product)
    }
  }
  score_patient <- lapply(x_abs, function(value) {
    .dsvert_glm_rat_mul(
      .dsvert_glm_rat_mul(weight_abs, residual_abs), value)
  })
  loglik_patient <- if (identical(schema$estimand$family, "binomial")) {
    log_two_upper <- .dsvert_glm_rat_log_interval(
      "2", schema$numeric$reference_precision_bits)$upper
    .dsvert_glm_rat_mul(weight_abs, .dsvert_glm_rat_add(
      .dsvert_glm_rat_mul("2", eta_abs), log_two_upper))
  } else {
    exponential_upper <- .dsvert_glm_rat_exp_interval(
      eta_abs, schema$numeric$reference_precision_bits)$upper
    log_y_upper <- if (.dsvert_glm_rat_cmp(response_upper, "1") > 0L) {
      .dsvert_glm_rat_log_interval(
        response_upper, schema$numeric$reference_precision_bits)$upper
    } else .dsvert_glm_rat("0")
    .dsvert_glm_rat_mul(weight_abs,
      .dsvert_formal_glm_inference_rational_sum(list(
        .dsvert_glm_rat_mul(response_upper, eta_abs),
        exponential_upper,
        .dsvert_glm_rat_mul(response_upper, log_y_upper))))
  }
  loss_patient <- .dsvert_glm_rat_mul(
    weight_abs, .dsvert_glm_rat_mul(
      .dsvert_glm_rat_add(mean_abs, response_upper), eta_abs))

  patient_blocks <- list(
    information_upper_triangle = h_patient,
    score_meat_upper_triangle = j_patient,
    score_vector = score_patient,
    canonical_bounded_log_likelihood_at_dp_beta = list(loglik_patient),
    integrated_pwl_surrogate_loss_at_dp_beta = list(loss_patient),
    admitted_n = list(.dsvert_glm_rat("1")))
  patient_abs <- unname(unlist(patient_blocks[vapply(
    names(patient_blocks), function(field) {
      identical(capabilities[[field]], .DSVERT_FORMAL_GLM_INFERENCE_RELEASED)
    }, logical(1L))], recursive = FALSE, use.names = FALSE))

  ridge <- lapply(coefficients, function(coefficient) {
    .dsvert_glm_rat(coefficient$ridge)
  })
  beta_squared <- lapply(beta, .dsvert_glm_rat_pow, exponent = 2L)
  ridge_loss <- .dsvert_glm_rat_mul(
    .dsvert_glm_rat_div(capacity, "2"),
    .dsvert_formal_glm_inference_rational_sum(Map(
      .dsvert_glm_rat_mul, ridge, beta_squared)))
  aggregate_lower <- aggregate_upper <- list()
  triangle_position <- 1L
  for (row in seq_len(p)) {
    for (column in seq.int(row, p)) {
      spread <- .dsvert_glm_rat_mul(
        capacity, h_patient[[triangle_position]])
      deterministic <- if (row == column) {
        .dsvert_glm_rat_mul(capacity, ridge[[row]])
      } else .dsvert_glm_rat("0")
      aggregate_lower[[length(aggregate_lower) + 1L]] <-
        .dsvert_glm_rat_sub(deterministic, spread)
      aggregate_upper[[length(aggregate_upper) + 1L]] <-
        .dsvert_glm_rat_add(deterministic, spread)
      triangle_position <- triangle_position + 1L
    }
  }
  append_symmetric_bounds <- function(bounds) {
    for (bound in bounds) {
      aggregate_lower[[length(aggregate_lower) + 1L]] <<-
        .dsvert_glm_rat_neg(.dsvert_glm_rat_mul(capacity, bound))
      aggregate_upper[[length(aggregate_upper) + 1L]] <<-
        .dsvert_glm_rat_mul(capacity, bound)
    }
  }
  if (identical(capabilities$score_meat_upper_triangle,
                .DSVERT_FORMAL_GLM_INFERENCE_RELEASED)) {
    append_symmetric_bounds(j_patient)
  }
  if (identical(capabilities$score_vector,
                .DSVERT_FORMAL_GLM_INFERENCE_RELEASED)) {
    append_symmetric_bounds(score_patient)
  }
  if (identical(
      capabilities$canonical_bounded_log_likelihood_at_dp_beta,
      .DSVERT_FORMAL_GLM_INFERENCE_RELEASED)) {
    append_symmetric_bounds(list(loglik_patient))
  }
  if (identical(capabilities$integrated_pwl_surrogate_loss_at_dp_beta,
                .DSVERT_FORMAL_GLM_INFERENCE_RELEASED)) {
    spread <- .dsvert_glm_rat_mul(capacity, loss_patient)
    aggregate_lower[[length(aggregate_lower) + 1L]] <-
      .dsvert_glm_rat_sub(ridge_loss, spread)
    aggregate_upper[[length(aggregate_upper) + 1L]] <-
      .dsvert_glm_rat_add(ridge_loss, spread)
  }
  if (identical(capabilities$admitted_n,
                .DSVERT_FORMAL_GLM_INFERENCE_RELEASED)) {
    aggregate_lower[[length(aggregate_lower) + 1L]] <- .dsvert_glm_rat("0")
    aggregate_upper[[length(aggregate_upper) + 1L]] <- capacity
  }
  coordinate_order <- .dsvert_formal_glm_inference_coordinate_layout(
    coefficient_names, capabilities)
  if (length(patient_abs) != length(coordinate_order) ||
      length(aggregate_lower) != length(coordinate_order) ||
      length(aggregate_upper) != length(coordinate_order)) {
    .dsvert_formal_glm_release_stop(
      "The inference sensitivity derivation has an inconsistent shape",
      "inference_sensitivity_derivation_invalid")
  }
  lower_steps <- vapply(aggregate_lower, function(value) {
    .dsvert_formal_glm_inference_integer_bound(value, scale, "floor")
  }, character(1L))
  absolute_upper_steps <- vapply(aggregate_upper, function(value) {
    .dsvert_formal_glm_inference_integer_bound(value, scale, "ceiling")
  }, character(1L))
  shifted_upper_steps <- vapply(seq_along(lower_steps), function(index) {
    .dsvert_glm_rat_json(.dsvert_glm_rat_sub(
      absolute_upper_steps[[index]], lower_steps[[index]]))$numerator
  }, character(1L))
  l1_sensitivity <- .dsvert_glm_rat_mul(
    gamma, .dsvert_formal_glm_inference_rational_sum(patient_abs))
  sensitivity_steps <- .dsvert_formal_glm_inference_integer_bound(
    l1_sensitivity, scale, "ceiling")
  fields <- list(
    version = "dsvert-formal-glm-joint-inference-server-plan-v1",
    schema_sha256 = point_result$schema_sha256,
    beta_public_release_sha256 = point_result$public_release_sha256,
    adjacency = schema$adjacency$definition,
    adjacency_triangle_factor = schema$adjacency$triangle_factor_gamma,
    coordinate_order_sha256 = .dsvert_formal_glm_inference_hash(
      as.list(coordinate_order)),
    output_lattice_scale = scale_text,
    sensitivity_norm_bound =
      "exact_rational_stacked_l1_upper_bounds_joint_l2",
    per_patient_coordinate_abs_bounds = lapply(
      patient_abs, .dsvert_glm_rat_json),
    aggregate_coordinate_lower_bounds = lapply(
      aggregate_lower, .dsvert_glm_rat_json),
    aggregate_coordinate_upper_bounds = lapply(
      aggregate_upper, .dsvert_glm_rat_json),
    stacked_l1_sensitivity_upper = .dsvert_glm_rat_json(l1_sensitivity),
    l2_sensitivity_steps = sensitivity_steps,
    coordinate_lower_bounds_steps = as.list(lower_steps),
    shifted_upper_bounds_steps = as.list(shifted_upper_steps),
    beta_is_fixed_public_input = TRUE,
    ridge_terms_are_public_deterministic_zero_sensitivity = TRUE,
    invalid_rows_are_zero_contribution_before_noise = TRUE)
  c(fields, list(certificate_sha256 =
    .dsvert_formal_glm_inference_hash(fields)))
}

.dsvert_formal_glm_inference_base_expected <- function(point_result) {
  certificate <- point_result$release_certificate
  phase <- certificate$phase19_binding
  identity <- certificate$release_identity
  schema <- point_result$schema
  list(
    version = .DSVERT_FORMAL_GLM_INFERENCE_BASE_BINDING,
    public_release_sha256 = point_result$public_release_sha256,
    release_instance_id = point_result$release_instance_id,
    release_contract_sha256 = point_result$release_contract_sha256,
    privacy_epoch_sha256 = point_result$privacy_epoch_sha256,
    materialization_root_sha256 = identity$materialization_root_sha256,
    post_execution_root_sha256 = phase$post_execution_root_sha256,
    final_vector_root_sha256 = identity$final_vector_root_sha256,
    coordinate_order_sha256 = identity$coordinate_order_sha256,
    logical_snapshot_sha256 = .dsvert_formal_glm_inference_hash(
      schema$authority$logical_snapshot),
    output_lattice_scale = point_result$output_lattice_scale,
    shifted_coefficient_lattice_steps =
      as.list(unname(point_result$shifted_coefficient_lattice_steps)),
    shifted_upper_bounds_steps =
      as.list(unname(point_result$shifted_upper_bounds_steps)),
    beta_is_public_fixed_adaptive_input = TRUE,
    same_logical_snapshot = TRUE)
}

.dsvert_formal_glm_inference_validate_point <- function(point_result) {
  schema_hash <- if (is.list(point_result$schema)) tryCatch(
    digest::digest(.dsvert_joint_dp_client_json(point_result$schema),
                   algo = "sha256", serialize = FALSE),
    error = function(error) NULL) else NULL
  if (!inherits(point_result, "dsvert_formal_dp_glm") ||
      !is.list(point_result$release_certificate) ||
      !is.list(point_result$schema) ||
      !is.character(point_result$public_release_sha256) ||
      length(point_result$public_release_sha256) != 1L ||
      !identical(point_result$schema_sha256,
                 point_result$release_certificate$schema_sha256) ||
      !identical(schema_hash, point_result$schema_sha256) ||
      !identical(names(point_result$coefficients),
                 unlist(point_result$schema$estimand$coefficient_order,
                        use.names = FALSE)) ||
      !identical(point_result$protected_path$intermediate_openings, 0L) ||
      !identical(point_result$protected_path$final_openings, 1L)) {
    .dsvert_formal_glm_release_stop(
      "Invalid formal DP GLM point release for joint inference",
      "inference_base_release_invalid")
  }
  .dsvert_formal_glm_release_hash(
    point_result$public_release_sha256, "base public release hash")
  scale <- tryCatch(.dsvert_glm_rat(point_result$output_lattice_scale),
                    error = function(error) NULL)
  signed_steps <- tryCatch(lapply(
    point_result$coefficient_lattice_steps, .dsvert_glm_rat),
    error = function(error) NULL)
  decoded <- if (is.null(scale) || is.null(signed_steps)) numeric() else
    vapply(signed_steps, function(value) {
      .dsvert_glm_rat_double(.dsvert_glm_rat_div(value, scale))
    }, numeric(1L))
  names(decoded) <- names(point_result$coefficients)
  tolerance <- if (is.null(scale)) Inf else
    .dsvert_glm_rat_double(.dsvert_glm_rat_div("1", scale)) /
      2 + .Machine$double.eps
  if (!identical(names(decoded), names(point_result$coefficients)) ||
      length(decoded) != length(point_result$coefficients) ||
      any(!is.finite(decoded)) ||
      any(abs(decoded - point_result$coefficients) > tolerance)) {
    .dsvert_formal_glm_release_stop(
      "The formal DP beta object was modified after signed decoding",
      "inference_base_release_invalid")
  }
  invisible(TRUE)
}

.dsvert_formal_glm_inference_verify_signatures <- function(
    release, unsigned, designated, pinset) {
  signatures <- release$signatures
  if (!is.list(signatures) || !is.null(names(signatures)) ||
      length(signatures) != 2L) {
    .dsvert_formal_glm_release_stop(
      "The formal inference artifact lacks both compute-peer signatures",
      "inference_signature_set_invalid")
  }
  signers <- vapply(signatures, function(item) {
    valid <- is.list(item) &&
      identical(sort(names(item), method = "radix"),
                c("signature", "signer")) &&
      is.character(item$signer) && length(item$signer) == 1L &&
      !is.na(item$signer) &&
      is.character(item$signature) && length(item$signature) == 1L &&
      !is.na(item$signature) &&
      nchar(item$signature, type = "bytes") == 86L &&
      grepl("^[A-Za-z0-9_-]{86}$", item$signature)
    if (valid) item$signer else NA_character_
  }, character(1L))
  if (anyNA(signers) || anyDuplicated(signers) ||
      !identical(signers, designated)) {
    .dsvert_formal_glm_release_stop(
      "The formal inference signature set is not canonical",
      "inference_signature_set_invalid")
  }
  message <- charToRaw(paste0(
    .DSVERT_FORMAL_GLM_INFERENCE_DOMAIN,
    .dsvert_joint_dp_client_json(unsigned)))
  for (index in seq_along(signatures)) {
    peer <- signers[[index]]
    public <- tryCatch(.dsvert_joint_dp_client_b64url(
      unname(pinset[[peer]]), 32L, "formal inference identity key"),
      error = function(error) raw())
    signature <- tryCatch(.dsvert_joint_dp_client_b64url(
      signatures[[index]]$signature, 64L, "formal inference signature"),
      error = function(error) raw())
    verified <- length(public) == 32L && length(signature) == 64L &&
      tryCatch(openssl::ed25519_verify(
        message, signature, openssl::read_ed25519_pubkey(public)),
        error = function(error) FALSE)
    if (!isTRUE(verified)) {
      .dsvert_formal_glm_release_stop(
        paste0("Invalid formal inference signature from pinned peer ", peer),
        "inference_signature_invalid")
    }
  }
  invisible(TRUE)
}

.dsvert_formal_glm_inference_scientific_contract_expected <- function() {
  list(
    beta_evaluation = "already_released_dp_beta_fixed_public_input",
    information_definition = paste0(
      "sum_i_weight_i_pwl_slope_eta_i_x_i_x_i_transpose_plus_",
      "capacity_times_ridge_diagonal"),
    information_scale = "patient_sum_scale",
    score_definition = "sum_i_weight_i_x_i_times_y_i_minus_pwl_mean_eta_i",
    score_meat_definition = paste0(
      "sum_i_patient_score_i_patient_score_i_transpose_uncentered_hc0"),
    pwl_derivative_tie_rule =
      "left_segment_at_internal_knots_endpoints_use_adjacent_segment",
    surrogate_loss_definition = paste0(
      "sum_i_weight_i_integrated_pwl_loss_plus_capacity_half_",
      "ridge_quadratic"),
    canonical_log_likelihood_definition = paste0(
      "bounded_canonical_family_log_likelihood_at_dp_beta_",
      "not_optimized_pwl_loglik"),
    ridge_estimand = "fully_ridge_regularized_target_not_unpenalized_mle",
    box_estimand = "box_constrained_target_not_unbounded_mle",
    regularity_scope = paste0(
      "working_covariance_conditional_on_no_positive_asymptotic_mass_",
      "at_pwl_knots"),
    penalized_curvature_inverse = paste0(
      "information_inverse_working_laplace_scale_not_sampling_",
      "covariance_or_standard_errors"),
    model_based_covariance = paste0(
      "unavailable_requires_unpenalized_likelihood_bread_and_",
      "model_variance_meat"),
    robust_covariance =
      "information_inverse_score_meat_information_inverse_hc0",
    score_meat_psd_postprocessing =
      "symmetric_eigen_projection_negative_eigenvalues_to_zero",
    beta_mechanism_region_requires_strict_box_interior = TRUE,
    dp_uncertainty_separate_from_sampling_uncertainty = TRUE,
    wald = "unavailable_no_noise_aware_null_calibration",
    likelihood_ratio =
      "unavailable_requires_joint_nested_same_snapshot_contrast_artifact",
    aic = "unavailable_not_classical_for_fixed_iteration_ridge_box_target",
    deviance = "unavailable_not_reconstructed_from_patient_residuals")
}

.dsvert_formal_glm_inference_scientific_contract <- function(value) {
  expected <- .dsvert_formal_glm_inference_scientific_contract_expected()
  .dsvert_formal_glm_release_exact_names(
    value, names(expected), "inference scientific contract")
  valid <- vapply(names(expected), function(field) {
    identical(value[[field]], expected[[field]])
  }, logical(1L))
  if (any(!valid)) {
    .dsvert_formal_glm_release_stop(
      "The formal inference artifact changed its scientific estimand",
      "inference_scientific_contract_invalid")
  }
  value
}

.dsvert_formal_glm_inference_triangle_matrix <- function(
    values, coefficients) {
  p <- length(coefficients)
  result <- matrix(0, p, p, dimnames = list(coefficients, coefficients))
  position <- 1L
  for (row in seq_len(p)) {
    for (column in seq.int(row, p)) {
      result[row, column] <- result[column, row] <- values[[position]]
      position <- position + 1L
    }
  }
  result
}

.dsvert_formal_glm_inference_block_values <- function(
    values, coefficients, capabilities) {
  p <- length(coefficients)
  triangle <- p * (p + 1L) / 2L
  position <- 1L
  take <- function(count) {
    selected <- values[seq.int(position, length.out = count)]
    position <<- position + count
    selected
  }
  result <- list(information_upper_triangle = take(triangle))
  for (field in .dsvert_formal_glm_inference_capability_fields()[-1L]) {
    if (identical(capabilities[[field]],
                  .DSVERT_FORMAL_GLM_INFERENCE_RELEASED)) {
      count <- switch(field,
        score_meat_upper_triangle = triangle,
        score_vector = p,
        canonical_bounded_log_likelihood_at_dp_beta = 1L,
        integrated_pwl_surrogate_loss_at_dp_beta = 1L,
        admitted_n = 1L)
      result[[field]] <- take(count)
    } else {
      result[[field]] <- NULL
    }
  }
  if (position != length(values) + 1L) {
    .dsvert_formal_glm_release_stop(
      "The formal inference coordinate blocks do not exhaust the vector",
      "inference_coordinate_shape_invalid")
  }
  result
}

.dsvert_formal_glm_inference_psd <- function(matrix) {
  symmetric <- (matrix + t(matrix)) / 2
  decomposition <- eigen(symmetric, symmetric = TRUE)
  clipped <- pmax(decomposition$values, 0)
  projected <- decomposition$vectors %*% (clipped *
    t(decomposition$vectors))
  projected <- (projected + t(projected)) / 2
  dimnames(projected) <- dimnames(matrix)
  list(
    matrix = projected,
    raw_eigenvalues = decomposition$values,
    projected_eigenvalues = clipped,
    projection_applied = any(decomposition$values < 0),
    negative_eigenvalue_mass = sum(abs(pmin(decomposition$values, 0))))
}

.dsvert_formal_glm_joint_inference_result <- function(
    point_result, release, trusted_pinset) {
  .dsvert_formal_glm_inference_validate_point(point_result)
  pinset <- .dsvert_formal_glm_release_validate_pinset(
    structure(list(unsigned_schema = point_result$schema),
              class = "dsvert_formal_glm_compilation"), trusted_pinset)
  top_fields <- c(
    "version", "schema_sha256", "family", "coefficient_order",
    "base_release_binding", "coordinate_contract",
    "shifted_coordinate_lattice_steps", "coordinate_lower_bounds_steps",
    "shifted_upper_bounds_steps", "output_lattice_scale",
    "server_derivation",
    "protected_binding", "authority", "release_identity", "mechanism",
    "numeric_certificate", "scientific_contract", "signatures")
  .dsvert_formal_glm_release_exact_names(
    release, top_fields, "joint inference release")
  if (!identical(release$version, .DSVERT_FORMAL_GLM_INFERENCE_VERSION)) {
    .dsvert_formal_glm_release_stop(
      "Unknown formal joint inference release version",
      "inference_version_invalid")
  }
  unsigned <- release[setdiff(names(release), "signatures")]
  canonical <- tryCatch(.dsvert_joint_dp_client_json(unsigned),
                        error = function(error) NULL)
  if (is.null(canonical) || nchar(canonical, type = "bytes") >
      .DSVERT_FORMAL_GLM_INFERENCE_MAX_RELEASE_BYTES) {
    .dsvert_formal_glm_release_stop(
      "The formal inference release is non-canonical or oversized",
      "inference_schema_invalid")
  }

  coefficients <- unlist(
    point_result$schema$estimand$coefficient_order, use.names = FALSE)
  released_order <- .dsvert_formal_glm_release_array(
    release$coefficient_order, "inference coefficient order")
  if (!identical(release$schema_sha256, point_result$schema_sha256) ||
      !identical(release$family, point_result$family) ||
      !release$family %in% c("binomial", "poisson") ||
      !identical(released_order, coefficients)) {
    .dsvert_formal_glm_release_stop(
      "The formal inference release changed its signed model",
      "inference_schema_binding_invalid")
  }

  expected_base <- .dsvert_formal_glm_inference_base_expected(point_result)
  if (!identical(.dsvert_joint_dp_client_canonical(
    release$base_release_binding),
    .dsvert_joint_dp_client_canonical(expected_base))) {
    .dsvert_formal_glm_release_stop(
      "The formal inference release is not bound to this DP beta release",
      "inference_base_binding_invalid")
  }

  coordinate_contract <- release$coordinate_contract
  .dsvert_formal_glm_release_exact_names(coordinate_contract, c(
    "version", "coordinate_order", "coordinate_order_sha256",
    "upper_triangle_order", "block_capabilities"),
    "inference coordinate contract")
  capabilities <- .dsvert_formal_glm_inference_validate_capabilities(
    coordinate_contract$block_capabilities)
  expected_coordinates <- .dsvert_formal_glm_inference_coordinate_layout(
    coefficients, capabilities)
  coordinates <- .dsvert_formal_glm_release_array(
    coordinate_contract$coordinate_order, "inference coordinate order")
  if (!identical(coordinate_contract$version,
                 .DSVERT_FORMAL_GLM_INFERENCE_COORDINATES) ||
      !identical(coordinate_contract$upper_triangle_order,
                 "row_major_i_then_j_for_i_less_than_or_equal_to_j") ||
      !identical(coordinates, expected_coordinates) ||
      !identical(coordinate_contract$coordinate_order_sha256,
                 .dsvert_formal_glm_inference_hash(as.list(coordinates)))) {
    .dsvert_formal_glm_release_stop(
      "The formal inference coordinate layout is not canonical",
      "inference_coordinate_contract_invalid")
  }

  shifted <- .dsvert_formal_glm_release_array(
    release$shifted_coordinate_lattice_steps,
    "shifted inference coordinate lattice")
  lower <- .dsvert_formal_glm_release_array(
    release$coordinate_lower_bounds_steps,
    "inference coordinate lower bounds")
  upper <- .dsvert_formal_glm_release_array(
    release$shifted_upper_bounds_steps,
    "inference coordinate shifted upper bounds")
  if (length(shifted) != length(coordinates) ||
      length(lower) != length(coordinates) ||
      length(upper) != length(coordinates)) {
    .dsvert_formal_glm_release_stop(
      "The formal inference vector has the wrong shape",
      "inference_coordinate_shape_invalid")
  }
  shifted <- vapply(seq_along(shifted), function(index) {
    .dsvert_formal_glm_release_integer_text(
      shifted[[index]], paste0("shifted inference coordinate ", index))
  }, character(1L))
  lower <- vapply(seq_along(lower), function(index) {
    .dsvert_formal_glm_release_integer_text(
      lower[[index]], paste0("inference lower bound ", index))
  }, character(1L))
  upper <- vapply(seq_along(upper), function(index) {
    .dsvert_formal_glm_release_integer_text(
      upper[[index]], paste0("inference shifted upper bound ", index))
  }, character(1L))
  scale_text <- .dsvert_formal_glm_release_integer_text(
    release$output_lattice_scale, "inference output lattice scale",
    positive = TRUE)
  scale <- .dsvert_glm_rat(scale_text)
  expected_server_plan <- .dsvert_formal_glm_inference_server_plan(
    point_result, capabilities, scale_text)
  if (!identical(.dsvert_joint_dp_client_canonical(
      release$server_derivation),
      .dsvert_joint_dp_client_canonical(expected_server_plan)) ||
      !identical(as.list(lower),
                 expected_server_plan$coordinate_lower_bounds_steps) ||
      !identical(as.list(upper),
                 expected_server_plan$shifted_upper_bounds_steps)) {
    .dsvert_formal_glm_release_stop(
      paste0("The inference bounds or joint sensitivity were not derived ",
             "from the signed schema and fixed public beta"),
      "inference_sensitivity_derivation_invalid")
  }
  shifted_rat <- lapply(shifted, .dsvert_glm_rat)
  upper_rat <- lapply(upper, .dsvert_glm_rat)
  in_bounds <- vapply(seq_along(shifted_rat), function(index) {
    .dsvert_glm_rat_cmp(shifted_rat[[index]], "0") >= 0L &&
      .dsvert_glm_rat_cmp(upper_rat[[index]], "0") >= 0L &&
      .dsvert_glm_rat_cmp(shifted_rat[[index]], upper_rat[[index]]) <= 0L
  }, logical(1L))
  if (any(!in_bounds)) {
    .dsvert_formal_glm_release_stop(
      "A formal inference coordinate escaped its server-owned bound",
      "inference_projection_invalid")
  }
  decoded_rat <- lapply(seq_along(shifted), function(index) {
    .dsvert_glm_rat_div(
      .dsvert_glm_rat_add(lower[[index]], shifted[[index]]), scale)
  })
  decoded <- vapply(decoded_rat, .dsvert_glm_rat_double, numeric(1L))
  names(decoded) <- coordinates
  if (any(!is.finite(decoded))) {
    .dsvert_formal_glm_release_stop(
      "A formal inference coordinate is not representable",
      "inference_numeric_unrepresentable")
  }

  protected <- release$protected_binding
  .dsvert_formal_glm_release_exact_names(protected, c(
    "version", "protected_path", "source_materialization_root_sha256",
    "base_post_execution_root_sha256",
    "inference_post_execution_root_sha256",
    "execution_receipt_pair_sha256", "same_snapshot_as_base",
    "protected_data_e2e_verified", "hidden_execution_validity_consumed",
    "invalid_execution_maps_to_neutral_before_noise",
    "intermediate_openings", "final_openings"),
    "inference protected binding")
  if (!identical(protected$version,
                 .DSVERT_FORMAL_GLM_INFERENCE_PROTECTED_BINDING) ||
      !identical(protected$protected_path,
                 "phase19_verified_beta_to_joint_inference_v1") ||
      !identical(protected$source_materialization_root_sha256,
                 expected_base$materialization_root_sha256) ||
      !identical(protected$base_post_execution_root_sha256,
                 expected_base$post_execution_root_sha256)) {
    .dsvert_formal_glm_release_stop(
      "The formal inference release bypassed its protected base snapshot",
      "inference_protected_path_invalid")
  }
  for (field in c("inference_post_execution_root_sha256",
                  "execution_receipt_pair_sha256")) {
    .dsvert_formal_glm_release_hash(protected[[field]], field)
  }
  for (field in c("same_snapshot_as_base", "protected_data_e2e_verified",
                  "hidden_execution_validity_consumed",
                  "invalid_execution_maps_to_neutral_before_noise")) {
    .dsvert_formal_glm_release_bool(
      protected[[field]], TRUE, gsub("_", " ", field))
  }
  intermediate_openings <- .dsvert_formal_glm_release_integer(
    protected$intermediate_openings, 0L, 0L,
    "inference intermediate openings")
  final_openings <- .dsvert_formal_glm_release_integer(
    protected$final_openings, 1L, 1L, "inference final openings")

  authority <- release$authority
  .dsvert_formal_glm_release_exact_names(authority, c(
    "pinset_sha256", "custodian_peers", "designated_compute_peers",
    "all_k_manifest_admission_verified",
    "all_k_manifest_signature_set_sha256", "server_derived_bounds",
    "server_derived_sensitivity", "server_derived_contribution_caps",
    "server_derived_block_availability", "operation_limit",
    "request_limit", "history_can_deny_operation"),
    "inference authority")
  custodians <- .dsvert_formal_glm_release_array(
    authority$custodian_peers, "inference custodians")
  designated <- .dsvert_formal_glm_release_array(
    authority$designated_compute_peers, "inference compute peers")
  schema_authority <- point_result$schema$authority
  if (!identical(custodians, unlist(
      schema_authority$custodian_peers, use.names = FALSE)) ||
      !identical(designated, unlist(
        schema_authority$designated_peers, use.names = FALSE)) ||
      length(custodians) < 2L || length(designated) != 2L ||
      !all(designated %in% custodians) || anyDuplicated(custodians) ||
      !identical(authority$pinset_sha256,
                 schema_authority$pinset_sha256)) {
    .dsvert_formal_glm_release_stop(
      "The formal inference release changed the pinned consortium",
      "inference_authority_invalid")
  }
  .dsvert_formal_glm_release_bool(
    authority$all_k_manifest_admission_verified, TRUE,
    "inference all-K manifest admission")
  .dsvert_formal_glm_release_hash(
    authority$all_k_manifest_signature_set_sha256,
    "inference all-K manifest signature set")
  for (field in c("server_derived_bounds", "server_derived_sensitivity",
                  "server_derived_contribution_caps",
                  "server_derived_block_availability")) {
    .dsvert_formal_glm_release_bool(
      authority[[field]], TRUE, gsub("_", " ", field))
  }
  for (field in c("operation_limit", "request_limit",
                  "history_can_deny_operation")) {
    .dsvert_formal_glm_release_bool(
      authority[[field]], FALSE, gsub("_", " ", field))
  }

  identity <- release$release_identity
  .dsvert_formal_glm_release_exact_names(identity, c(
    "release_instance_id", "release_contract_sha256",
    "privacy_epoch_sha256", "materialization_root_sha256",
    "base_public_release_sha256", "final_vector_root_sha256",
    "coordinate_order_sha256", "sticky_absolute_coordinates",
    "retry_restart_rechunk_same_release",
    "root_rotation_starts_visible_composed_release"),
    "inference release identity")
  for (field in c("release_instance_id", "release_contract_sha256",
                  "privacy_epoch_sha256", "materialization_root_sha256",
                  "base_public_release_sha256", "final_vector_root_sha256",
                  "coordinate_order_sha256")) {
    .dsvert_formal_glm_release_hash(identity[[field]], field)
  }
  if (!identical(identity$materialization_root_sha256,
                 expected_base$materialization_root_sha256) ||
      !identical(identity$base_public_release_sha256,
                 point_result$public_release_sha256) ||
      !identical(identity$coordinate_order_sha256,
                 coordinate_contract$coordinate_order_sha256)) {
    .dsvert_formal_glm_release_stop(
      "The formal inference release identity is inconsistent",
      "inference_identity_invalid")
  }
  for (field in c("sticky_absolute_coordinates",
                  "retry_restart_rechunk_same_release",
                  "root_rotation_starts_visible_composed_release")) {
    .dsvert_formal_glm_release_bool(
      identity[[field]], TRUE, gsub("_", " ", field))
  }

  mechanism <- release$mechanism
  .dsvert_formal_glm_release_exact_names(mechanism, c(
    "backend", "mechanism", "sampler", "epsilon", "delta",
    "l2_sensitivity_steps", "sigma_squared_steps",
    "mechanism_variance_upper_steps", "simultaneous_95_abs_steps",
    "nominal_variance_multiplier", "release_delta_aggregation",
    "complete_epsilon_delta_per_peer", "epsilon_divided_by_peer_count",
    "finite_support_transfer_charged", "one_global_vector_release",
    "automatic_fallback_used", "composition_rule", "base_epsilon",
    "base_delta", "composed_epsilon", "composed_delta",
    "beta_release_is_public_adaptive_input", "ledger_accounting_only",
    "gaussian_core_delta", "finite_support_transfer_delta",
    "calibration_rule", "minimum_sigma_squared_steps",
    "simultaneous_radius_rule"),
    "inference DP mechanism")
  if (!identical(mechanism$backend, .DSVERT_FORMAL_GLM_DP_BACKEND) ||
      !identical(mechanism$mechanism, .DSVERT_FORMAL_GLM_DP_MECHANISM) ||
      !identical(mechanism$sampler, .DSVERT_FORMAL_GLM_DP_SAMPLER) ||
      !identical(mechanism$release_delta_aggregation,
                 "max_per_peer_not_sum") ||
      !identical(mechanism$composition_rule,
                 "basic_sequential_adaptive_composition_exact_sum") ||
      !identical(mechanism$calibration_rule, paste0(
        "sigma_squared_gte_sensitivity_squared_times_2_log_",
        "5_over_4_delta_core_over_epsilon_squared")) ||
      !identical(mechanism$simultaneous_radius_rule,
                 "union_chebyshev_95_percent_all_coordinates")) {
    .dsvert_formal_glm_release_stop(
      "The formal inference release substituted its DP mechanism",
      "inference_mechanism_invalid")
  }
  epsilon <- .dsvert_formal_glm_release_rat(
    mechanism$epsilon, "inference epsilon", positive = TRUE)
  delta <- .dsvert_formal_glm_release_rat(
    mechanism$delta, "inference delta", positive = TRUE)
  base_epsilon <- .dsvert_formal_glm_release_rat(
    mechanism$base_epsilon, "base epsilon", positive = TRUE)
  base_delta <- .dsvert_formal_glm_release_rat(
    mechanism$base_delta, "base delta", positive = TRUE)
  composed_epsilon <- .dsvert_formal_glm_release_rat(
    mechanism$composed_epsilon, "composed epsilon", positive = TRUE)
  composed_delta <- .dsvert_formal_glm_release_rat(
    mechanism$composed_delta, "composed delta", positive = TRUE)
  core_delta <- .dsvert_formal_glm_release_rat(
    mechanism$gaussian_core_delta, "Gaussian core delta", positive = TRUE)
  transfer_delta <- .dsvert_formal_glm_release_rat(
    mechanism$finite_support_transfer_delta,
    "finite-support transfer delta", positive = TRUE)
  base_mechanism <- point_result$release_certificate$mechanism
  if (!identical(mechanism$base_epsilon, base_mechanism$epsilon) ||
      !identical(mechanism$base_delta, base_mechanism$delta) ||
      .dsvert_glm_rat_cmp(
        composed_epsilon,
        .dsvert_glm_rat_add(base_epsilon, epsilon)) != 0L ||
      .dsvert_glm_rat_cmp(
        composed_delta, .dsvert_glm_rat_add(base_delta, delta)) != 0L ||
      .dsvert_glm_rat_cmp(
        delta, .dsvert_glm_rat_add(core_delta, transfer_delta)) != 0L ||
      .dsvert_glm_rat_cmp(epsilon, "1") > 0L ||
      .dsvert_glm_rat_cmp(delta, "1") >= 0L ||
      .dsvert_glm_rat_cmp(composed_delta, "1") >= 0L) {
    .dsvert_formal_glm_release_stop(
      "The formal inference release has an invalid privacy composition",
      "inference_composition_invalid")
  }
  sensitivity <- .dsvert_formal_glm_release_integer_text(
    mechanism$l2_sensitivity_steps, "inference L2 sensitivity",
    positive = TRUE)
  if (!identical(sensitivity,
                 expected_server_plan$l2_sensitivity_steps)) {
    .dsvert_formal_glm_release_stop(
      "The inference mechanism changed the derived joint sensitivity",
      "inference_sensitivity_derivation_invalid")
  }
  minimum_sigma <- .dsvert_formal_glm_release_rat(
    mechanism$minimum_sigma_squared_steps,
    "minimum inference sigma squared", positive = TRUE)
  log_upper <- .dsvert_glm_rat_log_interval(
    .dsvert_glm_rat_div("5", .dsvert_glm_rat_mul("4", core_delta)),
    point_result$schema$numeric$reference_precision_bits)$upper
  expected_minimum_sigma <- .dsvert_glm_rat_div(
    .dsvert_glm_rat_mul(
      .dsvert_glm_rat_pow(sensitivity, 2L),
      .dsvert_glm_rat_mul("2", log_upper)),
    .dsvert_glm_rat_pow(epsilon, 2L))
  if (.dsvert_glm_rat_cmp(minimum_sigma, expected_minimum_sigma) != 0L) {
    .dsvert_formal_glm_release_stop(
      "The inference Gaussian calibration certificate is inconsistent",
      "inference_mechanism_invalid")
  }
  sigma_squared <- .dsvert_formal_glm_release_rat(
    mechanism$sigma_squared_steps, "inference sigma squared",
    positive = TRUE)
  variance_upper <- .dsvert_formal_glm_release_rat(
    mechanism$mechanism_variance_upper_steps,
    "inference variance upper bound", positive = TRUE)
  if (.dsvert_glm_rat_cmp(sigma_squared, minimum_sigma) < 0L ||
      .dsvert_glm_rat_cmp(
      variance_upper, .dsvert_glm_rat_mul(sigma_squared, "2")) < 0L) {
    .dsvert_formal_glm_release_stop(
      "The inference variance bound omitted the second full peer draw",
      "inference_mechanism_invalid")
  }
  radius_steps <- .dsvert_formal_glm_release_integer_text(
    mechanism$simultaneous_95_abs_steps,
    "inference simultaneous radius", positive = TRUE)
  simultaneous_threshold <- .dsvert_glm_rat_mul(
    as.character(20L * length(coordinates)), variance_upper)
  if (.dsvert_glm_rat_cmp(
      .dsvert_glm_rat_pow(radius_steps, 2L),
      simultaneous_threshold) < 0L) {
    .dsvert_formal_glm_release_stop(
      "The inference simultaneous radius lacks 95% joint coverage",
      "inference_mechanism_invalid")
  }
  .dsvert_formal_glm_release_integer(
    mechanism$nominal_variance_multiplier, 2L, 2L,
    "inference nominal variance multiplier")
  for (field in c("complete_epsilon_delta_per_peer",
                  "finite_support_transfer_charged",
                  "one_global_vector_release",
                  "beta_release_is_public_adaptive_input",
                  "ledger_accounting_only")) {
    .dsvert_formal_glm_release_bool(
      mechanism[[field]], TRUE, gsub("_", " ", field))
  }
  for (field in c("epsilon_divided_by_peer_count",
                  "automatic_fallback_used")) {
    .dsvert_formal_glm_release_bool(
      mechanism[[field]], FALSE, gsub("_", " ", field))
  }

  numeric <- release$numeric_certificate
  .dsvert_formal_glm_release_exact_names(numeric, c(
    "ring_bits", "no_wrap_certified", "final_projection",
    "coordinate_bounds_sha256", "stacked_sensitivity_sha256",
    "arithmetic_l2_error_upper", "output_quantization_l2_error_upper",
    "total_deterministic_l2_error_upper", "certificate_sha256"),
    "inference numeric certificate")
  ring_bits <- .dsvert_formal_glm_release_integer(
    numeric$ring_bits, 128L, 128L, "inference ring width")
  .dsvert_formal_glm_release_bool(
    numeric$no_wrap_certified, TRUE, "inference no-wrap certificate")
  if (!identical(numeric$final_projection,
                 "nonnegative_shifted_public_coordinate_box_v1") ||
      !identical(numeric$coordinate_bounds_sha256,
                 .dsvert_formal_glm_inference_hash(list(
                   lower = as.list(lower), shifted_upper = as.list(upper)))) ||
      !identical(numeric$stacked_sensitivity_sha256,
                 .dsvert_formal_glm_inference_hash(list(
                   coordinate_order_sha256 =
                     coordinate_contract$coordinate_order_sha256,
                   l2_sensitivity_steps = sensitivity,
                   output_lattice_scale = scale_text)))) {
    .dsvert_formal_glm_release_stop(
      "The inference numeric certificate changed its bounds or sensitivity",
      "inference_numeric_certificate_invalid")
  }
  arithmetic_error <- .dsvert_formal_glm_release_rat(
    numeric$arithmetic_l2_error_upper,
    "inference arithmetic error")
  quantization_error <- .dsvert_formal_glm_release_rat(
    numeric$output_quantization_l2_error_upper,
    "inference quantization error")
  total_error <- .dsvert_formal_glm_release_rat(
    numeric$total_deterministic_l2_error_upper,
    "inference total deterministic error")
  if (.dsvert_glm_rat_cmp(
      total_error,
      .dsvert_glm_rat_add(arithmetic_error, quantization_error)) < 0L) {
    .dsvert_formal_glm_release_stop(
      "The inference numeric error certificate omits a component",
      "inference_numeric_certificate_invalid")
  }
  numeric_hash <- .dsvert_formal_glm_release_hash(
    numeric$certificate_sha256, "inference numeric certificate hash")
  scientific <- .dsvert_formal_glm_inference_scientific_contract(
    release$scientific_contract)

  .dsvert_formal_glm_inference_verify_signatures(
    release, unsigned, designated, pinset)

  blocks <- .dsvert_formal_glm_inference_block_values(
    decoded, coefficients, capabilities)
  information <- .dsvert_formal_glm_inference_triangle_matrix(
    blocks$information_upper_triangle, coefficients)
  meat <- if (is.null(blocks$score_meat_upper_triangle)) NULL else
    .dsvert_formal_glm_inference_triangle_matrix(
      blocks$score_meat_upper_triangle, coefficients)
  score <- if (is.null(blocks$score_vector)) NULL else
    stats::setNames(blocks$score_vector, coefficients)

  scale_double <- .dsvert_glm_rat_double(scale)
  radius <- .dsvert_glm_rat_double(
    .dsvert_glm_rat_div(radius_steps, scale))
  variance <- .dsvert_glm_rat_double(variance_upper) / scale_double^2
  deterministic_error <- .dsvert_glm_rat_double(total_error)
  if (!is.finite(scale_double) || scale_double <= 0 ||
      !is.finite(radius) || radius <= 0 ||
      !is.finite(variance) || variance <= 0 ||
      !is.finite(deterministic_error) || deterministic_error < 0) {
    .dsvert_formal_glm_release_stop(
      "The inference accuracy certificate is not representable",
      "inference_numeric_unrepresentable")
  }
  p <- length(coefficients)
  information_eigen <- eigen(information, symmetric = TRUE,
                             only.values = TRUE)$values
  lambda_min <- min(information_eigen)
  information_operator_radius <- p * radius + deterministic_error
  spectral_margin <- lambda_min - information_operator_radius
  base_scale <- .dsvert_glm_rat(
    point_result$release_certificate$output_lattice_scale)
  box <- vapply(
    point_result$release_certificate$shifted_upper_bounds_steps,
    function(value) .dsvert_glm_rat_double(.dsvert_glm_rat_div(
      value, .dsvert_glm_rat_mul("2", base_scale))), numeric(1L))
  names(box) <- coefficients
  base_radius <- .dsvert_glm_rat_double(.dsvert_glm_rat_div(
    point_result$release_certificate$mechanism$simultaneous_95_abs_steps,
    base_scale))
  beta <- point_result$coefficients[coefficients]
  region_lower <- pmax(beta - base_radius, -box)
  region_upper <- pmin(beta + base_radius, box)
  box_interior <- all(region_lower > -box & region_upper < box)
  stable <- is.finite(spectral_margin) && spectral_margin >
    sqrt(.Machine$double.eps) * max(1, max(abs(information)))
  release_hash <- digest::digest(
    .dsvert_joint_dp_client_json(release), algo = "sha256",
    serialize = FALSE)
  common_reason <- NULL
  if (!box_interior) {
    common_reason <- paste0(
      "The simultaneous DP mechanism region for beta touches the signed ",
      "coefficient box; boundary-constrained working covariance is not ",
      "silently treated as an interior covariance")
  } else if (!stable) {
    common_reason <- paste0(
      "The released information matrix is not stably positive definite ",
      "under its simultaneous DP and deterministic error region")
  }
  information_inverse <- NULL
  if (is.null(common_reason)) {
    information_inverse <- tryCatch(solve(information),
                                    error = function(error) NULL)
    if (is.null(information_inverse) ||
        any(!is.finite(information_inverse))) {
      common_reason <- "The information matrix could not be stably inverted"
      information_inverse <- NULL
    } else {
      information_inverse <- (information_inverse +
                                t(information_inverse)) / 2
      dimnames(information_inverse) <- list(coefficients, coefficients)
    }
  }

  penalized_curvature <- if (!is.null(common_reason)) {
    .dsvert_formal_glm_inference_unavailable(
      "penalized_curvature_inverse", common_reason, release_hash)
  } else {
    structure(list(
      status = "available_penalized_curvature_inverse",
      inverse = information_inverse,
      formula = "H_total_inverse",
      scale = "patient_sum_scale",
      estimand = paste0(
        "fully_ridge_regularized_box_constrained_fixed_iteration_",
        "pwl_surrogate_target"),
      sampling_covariance = FALSE,
      sampling_standard_errors = FALSE,
      interpretation = "working_laplace_or_penalized_curvature_scale_only"),
      class = c("dsvert_formal_glm_penalized_curvature", "list"))
  }
  model_based <- .dsvert_formal_glm_inference_unavailable(
    "model_based_covariance",
    paste0("H_total includes the signed ridge penalty, so H_total^-1 is ",
           "only a penalized curvature inverse. Model-based sampling ",
           "covariance requires a separately released likelihood bread and ",
           "model-variance meat."), release_hash)

  meat_psd <- if (is.null(meat)) NULL else
    .dsvert_formal_glm_inference_psd(meat)
  robust <- if (is.null(meat)) {
    .dsvert_formal_glm_inference_unavailable(
      "robust_covariance",
      "The server did not promote a sensitivity-certified score-meat block",
      release_hash)
  } else if (!is.null(common_reason)) {
    .dsvert_formal_glm_inference_unavailable(
      "robust_covariance", common_reason, release_hash)
  } else {
    covariance <- information_inverse %*% meat_psd$matrix %*%
      information_inverse
    covariance <- (covariance + t(covariance)) / 2
    dimnames(covariance) <- list(coefficients, coefficients)
    standard_errors <- sqrt(pmax(diag(covariance), 0))
    names(standard_errors) <- coefficients
    structure(list(
      status = "available_conditional_robust_hc0_covariance",
      covariance = covariance,
      standard_errors = standard_errors,
      formula = "H_inverse_J_psd_H_inverse",
      score_meat_psd_projection = meat_psd,
      finite_sample_correction = "none_hc0",
      regularity_scope = scientific$regularity_scope,
      dp_mechanism_uncertainty_included = FALSE),
      class = c("dsvert_formal_glm_robust_covariance", "list"))
  }

  scalar_or_unavailable <- function(field, label) {
    value <- blocks[[field]]
    if (is.null(value)) {
      return(.dsvert_formal_glm_inference_unavailable(
        field,
        paste0("The server did not promote a sensitivity-certified ", label,
               " block"), release_hash))
    }
    unname(value[[1L]])
  }
  score_output <- if (is.null(score)) {
    .dsvert_formal_glm_inference_unavailable(
      "score_vector",
      "The server did not promote a sensitivity-certified score block",
      release_hash)
  } else score

  composed_privacy <- list(
    base_epsilon = .dsvert_glm_rat_double(base_epsilon),
    base_delta = .dsvert_glm_rat_double(base_delta),
    inference_epsilon = .dsvert_glm_rat_double(epsilon),
    inference_delta = .dsvert_glm_rat_double(delta),
    gaussian_core_delta = .dsvert_glm_rat_double(core_delta),
    finite_support_transfer_delta = .dsvert_glm_rat_double(transfer_delta),
    l2_sensitivity_steps = sensitivity,
    calibration_rule = mechanism$calibration_rule,
    composed_epsilon = .dsvert_glm_rat_double(composed_epsilon),
    composed_delta = .dsvert_glm_rat_double(composed_delta),
    composition_rule = mechanism$composition_rule,
    operation_limit = FALSE,
    request_limit = FALSE,
    history_can_deny_operation = FALSE,
    history_is_accounting_only = TRUE)
  mechanism_uncertainty <- list(
    scope = paste0(
      "simultaneous_95_percent_dp_mechanism_region_separate_from_",
      "sampling_working_covariance"),
    coordinate_abs_radius = radius,
    coordinate_variance_upper = variance,
    information_operator_abs_radius_upper = information_operator_radius,
    information_lambda_min_released = lambda_min,
    information_spectral_margin = spectral_margin,
    information_inverse_operator_perturbation_upper = if (stable) {
      information_operator_radius /
        (lambda_min * (lambda_min - information_operator_radius))
    } else Inf,
    deterministic_l2_error_upper = deterministic_error,
    covariance_not_combined_with_sampling = TRUE)
  statistic_status <- list(
    penalized_curvature_inverse = penalized_curvature$status,
    model_based_covariance = model_based$status,
    robust_covariance = robust$status,
    score_vector = if (inherits(
      score_output, "dsvert_formal_glm_inference_unavailable")) {
      score_output$status
    } else "released_dp_noisy_sum_scale_score",
    canonical_bounded_log_likelihood_at_dp_beta = if (is.null(
      blocks$canonical_bounded_log_likelihood_at_dp_beta)) {
      "joint_dp_inference_capability_unavailable"
    } else "released_dp_noisy_not_optimized_loglik",
    integrated_pwl_surrogate_loss_at_dp_beta = if (is.null(
      blocks$integrated_pwl_surrogate_loss_at_dp_beta)) {
      "joint_dp_inference_capability_unavailable"
    } else "released_dp_noisy_surrogate_loss",
    admitted_n = if (is.null(blocks$admitted_n)) {
      "joint_dp_inference_capability_unavailable"
    } else "released_dp_noisy_active_patient_count",
    wald = scientific$wald,
    likelihood_ratio = scientific$likelihood_ratio,
    aic = scientific$aic,
    deviance = scientific$deviance)

  result <- list(
    status = "dp_joint_inference_artifact",
    family = point_result$family,
    coefficients = point_result$coefficients,
    information = information,
    score_meat = if (is.null(meat_psd)) NULL else meat_psd$matrix,
    score_meat_raw_dp_release = meat,
    score = score_output,
    penalized_curvature = penalized_curvature,
    model_based = model_based,
    robust = robust,
    canonical_bounded_log_likelihood_at_dp_beta = scalar_or_unavailable(
      "canonical_bounded_log_likelihood_at_dp_beta",
      "canonical bounded log-likelihood"),
    integrated_pwl_surrogate_loss_at_dp_beta = scalar_or_unavailable(
      "integrated_pwl_surrogate_loss_at_dp_beta",
      "integrated PWL surrogate loss"),
    admitted_n = scalar_or_unavailable("admitted_n", "admitted-n"),
    statistic_status = statistic_status,
    mechanism_uncertainty = mechanism_uncertainty,
    privacy = composed_privacy,
    scientific_contract = scientific,
    capabilities = capabilities,
    beta_public_release_sha256 = point_result$public_release_sha256,
    public_release_sha256 = release_hash,
    release_instance_id = identity$release_instance_id,
    privacy_epoch_sha256 = identity$privacy_epoch_sha256,
    schema_sha256 = point_result$schema_sha256,
    logical_snapshot_sha256 = expected_base$logical_snapshot_sha256,
    materialization_root_sha256 =
      expected_base$materialization_root_sha256,
    protected_path = list(
      path = protected$protected_path,
      protected_data_e2e_verified = TRUE,
      hidden_execution_validity_consumed = TRUE,
      invalid_execution_maps_to_neutral_before_noise = TRUE,
      intermediate_openings = intermediate_openings,
      final_openings = final_openings,
      custodian_count = length(custodians),
      designated_compute_peers = designated,
      all_k_manifest_admission_verified = TRUE),
    numeric_certificate = list(
      ring_bits = ring_bits,
      no_wrap_certified = TRUE,
      certificate_sha256 = numeric_hash),
    root_rotation = "new_visible_composed_release_never_blocks",
    release_certificate = unsigned,
    production_ready = FALSE)
  class(result) <- c("dsvert_formal_glm_joint_inference", "list")
  result
}

.dsvert_formal_glm_joint_inference_statistic <- function(result, statistic) {
  if (!inherits(result, "dsvert_formal_glm_joint_inference")) {
    stop("result must be an internal formal joint-inference object",
         call. = FALSE)
  }
  allowed <- c("penalized_curvature_inverse",
               "model_based_covariance", "model_based_standard_errors",
               "robust_covariance", "robust_standard_errors", "wald",
               "likelihood_ratio", "aic", "deviance")
  if (!is.character(statistic) || length(statistic) != 1L ||
      is.na(statistic) || !statistic %in% allowed) {
    stop("Unknown formal joint-inference statistic", call. = FALSE)
  }
  if (identical(statistic, "penalized_curvature_inverse")) {
    if (inherits(result$penalized_curvature,
                 "dsvert_formal_glm_inference_unavailable")) {
      return(result$penalized_curvature)
    }
    return(result$penalized_curvature$inverse)
  }
  if (statistic %in% c("model_based_covariance",
                       "model_based_standard_errors")) {
    return(result$model_based)
  }
  if (statistic %in% c("robust_covariance", "robust_standard_errors")) {
    if (inherits(result$robust,
                 "dsvert_formal_glm_inference_unavailable")) {
      return(result$robust)
    }
    return(if (grepl("standard_errors$", statistic)) {
      result$robust$standard_errors
    } else result$robust$covariance)
  }
  .dsvert_formal_glm_inference_unavailable(
    statistic,
    switch(statistic,
      likelihood_ratio = paste0(
        "LR requires a separately calibrated, jointly signed nested-model ",
        "contrast on the same protected snapshot; two marginal releases are ",
        "never subtracted"),
      wald = "Wald tests require a separately validated noise-aware null calibration",
      aic = "Classical AIC is not defined for this fixed-iteration ridge-box target",
      deviance = "Protected residual deviance was not reconstructed"),
    result$public_release_sha256)
}

.dsvert_formal_glm_joint_nested_compatibility <- function(full, reduced) {
  if (!inherits(full, "dsvert_formal_glm_joint_inference") ||
      !inherits(reduced, "dsvert_formal_glm_joint_inference")) {
    stop("full and reduced must be formal joint-inference objects",
         call. = FALSE)
  }
  full_names <- names(full$coefficients)
  reduced_names <- names(reduced$coefficients)
  reasons <- character()
  if (!identical(full$family, reduced$family)) {
    reasons <- c(reasons, "different_family")
  }
  if (!identical(full$logical_snapshot_sha256,
                 reduced$logical_snapshot_sha256) ||
      !identical(full$materialization_root_sha256,
                 reduced$materialization_root_sha256)) {
    reasons <- c(reasons, "different_protected_snapshot")
  }
  if (!all(reduced_names %in% full_names) ||
      length(reduced_names) >= length(full_names)) {
    reasons <- c(reasons, "not_strictly_nested_coefficient_space")
  }
  if (length(reasons)) {
    return(structure(list(
      status = "incompatible_for_likelihood_ratio",
      reasons = unique(reasons), likelihood_ratio = NULL,
      additional_server_calls = 0L),
      class = c("dsvert_formal_glm_nested_compatibility", "list")))
  }
  structure(list(
    status = "requires_joint_nested_contrast_artifact",
    reasons = "marginal_dp_loglik_releases_must_not_be_subtracted",
    likelihood_ratio = NULL,
    same_snapshot = TRUE,
    structurally_nested = TRUE,
    additional_server_calls = 0L),
    class = c("dsvert_formal_glm_nested_compatibility", "list"))
}
