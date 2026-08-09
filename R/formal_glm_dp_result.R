# Internal consumer for the formal binomial/Poisson coefficient release.
#
# This file deliberately registers no public method.  It accepts only a
# two-peer-signed release bound to the complete trusted K-peer pinset, the
# Phase-1.8 -> Phase-1.9 protected path and the server-derived joint-DP
# certificate.  Sampling inference is not reconstructed from a noisy point
# estimate: the only uncertainty matrix exposed here is an upper bound for the
# mechanism noise itself.

.DSVERT_FORMAL_GLM_DP_RELEASE_VERSION <-
  "dsvert-formal-glm-dp-public-release-v1"
.DSVERT_FORMAL_GLM_DP_RELEASE_DOMAIN <-
  "dsVert/formal-glm/dp-public-release/v1|"
.DSVERT_FORMAL_GLM_DP_PHASE19_BINDING_VERSION <-
  "dsvert-formal-glm-phase19-public-binding-v1"
.DSVERT_FORMAL_GLM_DP_MAX_RELEASE_BYTES <- 8L * 1024L^2
.DSVERT_FORMAL_GLM_DP_BACKEND <- paste0(
  "independent_full_global_dyadic_discrete_gaussian_",
  "tv_bounded_ring128_v2")
.DSVERT_FORMAL_GLM_DP_MECHANISM <-
  "dyadic_discrete_gaussian_truncated_tv_bounded"
.DSVERT_FORMAL_GLM_DP_SAMPLER <- paste0(
  "cks-target-outward-rational-dyadic-cdf-hkdf-sha256-",
  "chacha20-coordinate-domain-v2")
.DSVERT_FORMAL_GLM_DP_INFERENCE_REQUIRED <-
  "requires_joint_dp_inference_artifact"

.dsvert_formal_glm_release_condition <- function(message, reason,
                                                 certificate = NULL) {
  structure(
    list(message = message, call = NULL, reason = reason,
         certificate = certificate),
    class = c("dsvert_formal_glm_release_error",
              "dsvert_numeric_condition", "error", "condition"))
}

.dsvert_formal_glm_release_stop <- function(message, reason,
                                            certificate = NULL) {
  stop(.dsvert_formal_glm_release_condition(
    message, reason = reason, certificate = certificate))
}

.dsvert_formal_glm_release_exact_names <- function(value, expected, what) {
  if (!is.list(value) || is.null(names(value)) || anyNA(names(value)) ||
      any(!nzchar(names(value))) || anyDuplicated(names(value)) ||
      !identical(sort(names(value), method = "radix"),
                 sort(expected, method = "radix"))) {
    .dsvert_formal_glm_release_stop(
      paste0("Invalid formal GLM ", what, " fields"),
      "release_schema_invalid")
  }
  invisible(TRUE)
}

.dsvert_formal_glm_release_string <- function(value, what,
                                              maximum = 4096L) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value) || nchar(value, type = "bytes") > maximum ||
      grepl("\\x00", value)) {
    .dsvert_formal_glm_release_stop(
      paste0("Invalid formal GLM ", what), "release_schema_invalid")
  }
  enc2utf8(value)
}

.dsvert_formal_glm_release_hash <- function(value, what) {
  value <- .dsvert_formal_glm_release_string(value, what, 64L)
  if (!grepl("^[0-9a-f]{64}$", value)) {
    .dsvert_formal_glm_release_stop(
      paste0("Invalid formal GLM ", what), "release_binding_invalid")
  }
  value
}

.dsvert_formal_glm_release_bool <- function(value, expected, what) {
  if (!is.logical(value) || length(value) != 1L || is.na(value) ||
      !identical(value, expected)) {
    .dsvert_formal_glm_release_stop(
      paste0("Invalid formal GLM ", what), "release_contract_invalid")
  }
  value
}

.dsvert_formal_glm_release_integer <- function(value, lower, upper, what) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value != floor(value) || value < lower || value > upper) {
    .dsvert_formal_glm_release_stop(
      paste0("Invalid formal GLM ", what), "release_schema_invalid")
  }
  as.integer(value)
}

.dsvert_formal_glm_release_integer_text <- function(value, what,
                                                    positive = FALSE) {
  value <- .dsvert_formal_glm_release_string(value, what, 4096L)
  pattern <- if (isTRUE(positive)) {
    "^[1-9][0-9]*$"
  } else {
    "^(0|-[1-9][0-9]*|[1-9][0-9]*)$"
  }
  if (!grepl(pattern, value)) {
    .dsvert_formal_glm_release_stop(
      paste0("Invalid canonical integer for ", what),
      "release_numeric_certificate_invalid")
  }
  value
}

.dsvert_formal_glm_release_array <- function(value, what) {
  if (is.character(value) && is.null(names(value))) {
    result <- unname(value)
  } else if (is.list(value) && is.null(names(value)) &&
             all(vapply(value, function(item) {
               is.character(item) && length(item) == 1L && !is.na(item)
             }, logical(1L)))) {
    result <- unname(unlist(value, use.names = FALSE))
  } else {
    result <- character()
  }
  if (!length(result) || anyNA(result) ||
      any(!nzchar(result)) || any(nchar(result, type = "bytes") > 4096L)) {
    .dsvert_formal_glm_release_stop(
      paste0("Invalid formal GLM ", what), "release_schema_invalid")
  }
  enc2utf8(result)
}

.dsvert_formal_glm_release_rat <- function(value, what,
                                           positive = FALSE,
                                           allow_zero = TRUE) {
  parsed <- tryCatch(.dsvert_glm_rat(value), error = function(error) NULL)
  canonical <- if (is.null(parsed)) NULL else .dsvert_glm_rat_json(parsed)
  valid <- !is.null(parsed) && identical(value, canonical)
  if (isTRUE(valid)) {
    comparison <- .dsvert_glm_rat_cmp(parsed, "0")
    valid <- if (isTRUE(positive)) comparison > 0L else
      comparison > 0L || (isTRUE(allow_zero) && comparison == 0L)
  }
  if (!isTRUE(valid)) {
    .dsvert_formal_glm_release_stop(
      paste0("Invalid exact rational for ", what),
      "release_numeric_certificate_invalid")
  }
  parsed
}

.dsvert_formal_glm_release_validate_pinset <- function(compilation, pinset) {
  schema <- compilation$unsigned_schema
  custodians <- unlist(schema$authority$custodian_peers, use.names = FALSE)
  if (!is.character(pinset) || is.null(names(pinset))) {
    .dsvert_formal_glm_release_stop(
      "The formal GLM release requires the trusted complete pinset",
      "release_pinset_invalid")
  }
  pinset <- pinset[order(names(pinset), method = "radix")]
  if (!isTRUE(.dsvert_dp_validate_pinset(pinset)) ||
      !identical(names(pinset), custodians) ||
      !identical(.dsvert_dp_pinset_hash(pinset),
                 schema$authority$pinset_sha256)) {
    .dsvert_formal_glm_release_stop(
      "The formal GLM release pinset does not match the signed schema",
      "release_pinset_invalid")
  }
  pinset
}

.dsvert_formal_glm_release_verify_signatures <- function(
    release, unsigned, designated, pinset) {
  signatures <- release$signatures
  if (!is.list(signatures) || !is.null(names(signatures)) ||
      length(signatures) != 2L ||
      any(vapply(signatures, is.null, logical(1L)))) {
    .dsvert_formal_glm_release_stop(
      "The formal GLM release lacks both compute-peer signatures",
      "release_signature_set_invalid")
  }
  signers <- vapply(signatures, function(item) {
    if (!is.list(item) || !identical(sort(names(item), method = "radix"),
                                     c("signature", "signer")) ||
        !is.character(item$signer) || length(item$signer) != 1L ||
        is.na(item$signer) || !is.character(item$signature) ||
        length(item$signature) != 1L || is.na(item$signature) ||
        nchar(item$signature, type = "bytes") != 86L ||
        !grepl("^[A-Za-z0-9_-]{86}$", item$signature)) {
      return(NA_character_)
    }
    item$signer
  }, character(1L))
  if (anyNA(signers) || anyDuplicated(signers) ||
      !identical(signers, designated)) {
    .dsvert_formal_glm_release_stop(
      "The formal GLM release signature set is not canonical",
      "release_signature_set_invalid")
  }
  message <- charToRaw(paste0(
    .DSVERT_FORMAL_GLM_DP_RELEASE_DOMAIN,
    .dsvert_joint_dp_client_json(unsigned)))
  for (index in seq_along(signatures)) {
    peer <- signers[[index]]
    public <- tryCatch(.dsvert_joint_dp_client_b64url(
      unname(pinset[[peer]]), 32L, "formal GLM release identity key"),
      error = function(error) raw())
    signature <- tryCatch(.dsvert_joint_dp_client_b64url(
      signatures[[index]]$signature, 64L, "formal GLM release signature"),
      error = function(error) raw())
    verified <- length(public) == 32L && length(signature) == 64L &&
      tryCatch(openssl::ed25519_verify(
        message, signature, openssl::read_ed25519_pubkey(public)),
        error = function(error) FALSE)
    if (!isTRUE(verified)) {
      .dsvert_formal_glm_release_stop(
        paste0("Invalid formal GLM release signature from pinned peer ",
               peer), "release_signature_invalid")
    }
  }
  invisible(TRUE)
}

.dsvert_formal_glm_release_postprocessing_contract <- function(value) {
  fields <- c(
    "coefficient_point_estimate", "mechanism_covariance",
    "sampling_covariance", "log_likelihood", "deviance", "aic",
    "wald", "likelihood_ratio", "prediction",
    "protected_residual_diagnostics")
  .dsvert_formal_glm_release_exact_names(
    value, fields, "post-processing contract")
  expected <- list(
    coefficient_point_estimate = "available",
    mechanism_covariance = "upper_bound_available",
    sampling_covariance = .DSVERT_FORMAL_GLM_DP_INFERENCE_REQUIRED,
    log_likelihood = .DSVERT_FORMAL_GLM_DP_INFERENCE_REQUIRED,
    deviance = .DSVERT_FORMAL_GLM_DP_INFERENCE_REQUIRED,
    aic = .DSVERT_FORMAL_GLM_DP_INFERENCE_REQUIRED,
    wald = .DSVERT_FORMAL_GLM_DP_INFERENCE_REQUIRED,
    likelihood_ratio = .DSVERT_FORMAL_GLM_DP_INFERENCE_REQUIRED,
    prediction = "public_covariates_only",
    protected_residual_diagnostics = "unavailable")
  matches <- vapply(names(expected), function(field) {
    identical(value[[field]], expected[[field]])
  }, logical(1L))
  if (any(!matches)) {
    .dsvert_formal_glm_release_stop(
      "The formal GLM release attempted to promote unsupported inference",
      "release_inference_contract_invalid")
  }
  value
}

.dsvert_formal_glm_dp_result <- function(compilation, release,
                                         trusted_pinset) {
  .dsvert_formal_glm_validate_compilation(compilation)
  pinset <- .dsvert_formal_glm_release_validate_pinset(
    compilation, trusted_pinset)
  top_fields <- c(
    "version", "schema_sha256", "family", "coefficient_order",
    "shifted_coefficient_lattice_steps", "shifted_upper_bounds_steps",
    "output_lattice_scale",
    "phase19_binding", "authority", "release_identity", "mechanism",
    "numeric_certificate", "postprocessing_contract", "signatures")
  .dsvert_formal_glm_release_exact_names(release, top_fields, "release")
  if (!identical(release$version,
                 .DSVERT_FORMAL_GLM_DP_RELEASE_VERSION)) {
    .dsvert_formal_glm_release_stop(
      "Unknown formal GLM DP release version", "release_version_invalid")
  }
  unsigned <- release[setdiff(names(release), "signatures")]
  canonical <- tryCatch(.dsvert_joint_dp_client_json(unsigned),
                        error = function(error) NULL)
  if (is.null(canonical) ||
      nchar(canonical, type = "bytes") >
        .DSVERT_FORMAL_GLM_DP_MAX_RELEASE_BYTES) {
    .dsvert_formal_glm_release_stop(
      "The formal GLM DP release is non-canonical or oversized",
      "release_schema_invalid")
  }

  schema <- compilation$unsigned_schema
  order <- .dsvert_formal_glm_release_array(
    release$coefficient_order, "coefficient order")
  expected_order <- unlist(
    schema$estimand$coefficient_order, use.names = FALSE)
  if (!identical(release$schema_sha256, compilation$sha256) ||
      !identical(release$family, schema$estimand$family) ||
      !release$family %in% c("binomial", "poisson") ||
      !identical(order, expected_order)) {
    .dsvert_formal_glm_release_stop(
      "The formal GLM release does not match its compiled estimand",
      "release_schema_binding_invalid")
  }
  shifted_steps <- .dsvert_formal_glm_release_array(
    release$shifted_coefficient_lattice_steps,
    "shifted coefficient lattice")
  shifted_upper <- .dsvert_formal_glm_release_array(
    release$shifted_upper_bounds_steps, "shifted coefficient bounds")
  if (length(shifted_steps) != length(order) ||
      length(shifted_upper) != length(order)) {
    .dsvert_formal_glm_release_stop(
      "The formal GLM coefficient vector has the wrong shape",
      "release_coordinate_shape_invalid")
  }
  shifted_steps <- vapply(seq_along(shifted_steps), function(index) {
    .dsvert_formal_glm_release_integer_text(
      shifted_steps[[index]], paste0("shifted coefficient step ", index))
  }, character(1L))
  shifted_upper <- vapply(seq_along(shifted_upper), function(index) {
    .dsvert_formal_glm_release_integer_text(
      shifted_upper[[index]], paste0("shifted upper bound ", index),
      positive = TRUE)
  }, character(1L))
  scale_text <- .dsvert_formal_glm_release_integer_text(
    release$output_lattice_scale, "output lattice scale", positive = TRUE)
  scale <- .dsvert_glm_rat(scale_text)

  phase <- release$phase19_binding
  .dsvert_formal_glm_release_exact_names(phase, c(
    "version", "protected_path", "post_execution_root_sha256",
    "execution_receipt_pair_sha256",
    "global_materialization_root_sha256", "protected_data_e2e_verified",
    "hidden_execution_validity_consumed",
    "invalid_execution_maps_to_neutral_before_noise",
    "intermediate_openings", "final_openings"), "Phase-1.9 binding")
  if (!identical(phase$version,
                 .DSVERT_FORMAL_GLM_DP_PHASE19_BINDING_VERSION) ||
      !identical(phase$protected_path, "phase18_v2_to_phase19_v2")) {
    .dsvert_formal_glm_release_stop(
      "The formal GLM release bypassed the protected materializer",
      "release_protected_path_invalid")
  }
  invisible(lapply(c(
    "post_execution_root_sha256", "execution_receipt_pair_sha256",
    "global_materialization_root_sha256"), function(field) {
      .dsvert_formal_glm_release_hash(phase[[field]], field)
    }))
  .dsvert_formal_glm_release_bool(
    phase$protected_data_e2e_verified, TRUE,
    "protected-data E2E attestation")
  .dsvert_formal_glm_release_bool(
    phase$hidden_execution_validity_consumed, TRUE,
    "hidden execution-validity consumption")
  .dsvert_formal_glm_release_bool(
    phase$invalid_execution_maps_to_neutral_before_noise, TRUE,
    "invalid-execution neutralization")
  intermediate_openings <- .dsvert_formal_glm_release_integer(
    phase$intermediate_openings, 0L, 0L, "intermediate opening count")
  final_openings <- .dsvert_formal_glm_release_integer(
    phase$final_openings, 1L, 1L, "final opening count")

  authority <- release$authority
  .dsvert_formal_glm_release_exact_names(authority, c(
    "pinset_sha256", "custodian_peers", "designated_compute_peers",
    "all_k_manifest_admission_verified",
    "all_k_manifest_signature_set_sha256", "server_derived_bounds",
    "server_derived_sensitivity", "server_derived_contribution_caps",
    "operation_limit", "request_limit", "history_can_deny_operation"),
    "release authority")
  custodians <- .dsvert_formal_glm_release_array(
    authority$custodian_peers, "custodian peers")
  designated <- .dsvert_formal_glm_release_array(
    authority$designated_compute_peers, "designated compute peers")
  expected_custodians <- unlist(
    schema$authority$custodian_peers, use.names = FALSE)
  expected_designated <- unlist(
    schema$authority$designated_peers, use.names = FALSE)
  if (!identical(custodians, expected_custodians) ||
      !identical(designated, expected_designated) ||
      length(custodians) < 2L || length(designated) != 2L ||
      anyDuplicated(custodians) || anyDuplicated(designated) ||
      !all(designated %in% custodians) ||
      !identical(authority$pinset_sha256,
                 schema$authority$pinset_sha256)) {
    .dsvert_formal_glm_release_stop(
      "The formal GLM release changed the pinned consortium",
      "release_authority_invalid")
  }
  .dsvert_formal_glm_release_bool(
    authority$all_k_manifest_admission_verified, TRUE,
    "all-K manifest admission")
  .dsvert_formal_glm_release_hash(
    authority$all_k_manifest_signature_set_sha256,
    "all-K manifest signature-set hash")
  for (field in c("server_derived_bounds", "server_derived_sensitivity",
                  "server_derived_contribution_caps")) {
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
    "final_vector_root_sha256", "coordinate_order_sha256",
    "sticky_absolute_coordinates", "retry_restart_rechunk_same_release",
    "root_rotation_starts_visible_composed_release"), "release identity")
  hashes <- c(
    "release_instance_id", "release_contract_sha256",
    "privacy_epoch_sha256", "materialization_root_sha256",
    "final_vector_root_sha256", "coordinate_order_sha256")
  invisible(lapply(hashes, function(field) {
    .dsvert_formal_glm_release_hash(identity[[field]], field)
  }))
  if (!identical(identity$materialization_root_sha256,
                 phase$global_materialization_root_sha256)) {
    .dsvert_formal_glm_release_stop(
      "The formal GLM release changed its materialization root",
      "release_identity_invalid")
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
    "automatic_fallback_used"), "DP mechanism")
  if (!identical(mechanism$backend, .DSVERT_FORMAL_GLM_DP_BACKEND) ||
      !identical(mechanism$mechanism,
                 .DSVERT_FORMAL_GLM_DP_MECHANISM) ||
      !identical(mechanism$sampler, .DSVERT_FORMAL_GLM_DP_SAMPLER) ||
      !identical(mechanism$release_delta_aggregation,
                 "max_per_peer_not_sum")) {
    .dsvert_formal_glm_release_stop(
      "The formal GLM release substituted its DP backend",
      "release_mechanism_invalid")
  }
  epsilon <- .dsvert_formal_glm_release_rat(
    mechanism$epsilon, "epsilon", positive = TRUE)
  delta <- .dsvert_formal_glm_release_rat(
    mechanism$delta, "delta", positive = TRUE)
  if (.dsvert_glm_rat_cmp(delta, "1") >= 0L) {
    .dsvert_formal_glm_release_stop(
      "The formal GLM delta must be below one",
      "release_mechanism_invalid")
  }
  sensitivity_steps <- .dsvert_formal_glm_release_integer_text(
    mechanism$l2_sensitivity_steps, "L2 sensitivity", positive = TRUE)
  sigma_squared <- .dsvert_formal_glm_release_rat(
    mechanism$sigma_squared_steps, "sigma squared", positive = TRUE)
  variance_upper <- .dsvert_formal_glm_release_rat(
    mechanism$mechanism_variance_upper_steps,
    "mechanism variance upper bound", positive = TRUE)
  minimum_variance <- .dsvert_glm_rat_mul(sigma_squared, "2")
  if (.dsvert_glm_rat_cmp(variance_upper, minimum_variance) < 0L) {
    .dsvert_formal_glm_release_stop(
      "The mechanism variance bound omitted the second full draw",
      "release_mechanism_invalid")
  }
  simultaneous_steps <- .dsvert_formal_glm_release_integer_text(
    mechanism$simultaneous_95_abs_steps,
    "simultaneous mechanism radius", positive = TRUE)
  .dsvert_formal_glm_release_integer(
    mechanism$nominal_variance_multiplier, 2L, 2L,
    "nominal variance multiplier")
  for (field in c("complete_epsilon_delta_per_peer",
                  "finite_support_transfer_charged",
                  "one_global_vector_release")) {
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
    "arithmetic_l2_error_upper", "optimizer_l2_error_upper",
    "link_l2_error_upper", "output_quantization_l2_error_upper",
    "total_deterministic_l2_error_upper", "certificate_sha256"),
    "numeric certificate")
  ring_bits <- .dsvert_formal_glm_release_integer(
    numeric$ring_bits, 128L, 128L, "ring width")
  .dsvert_formal_glm_release_bool(
    numeric$no_wrap_certified, TRUE, "no-wrap certificate")
  if (!identical(numeric$final_projection,
                 "nonnegative_shifted_public_coefficient_box_v1")) {
    .dsvert_formal_glm_release_stop(
      "The formal GLM release omitted its final coefficient projection",
      "release_numeric_certificate_invalid")
  }
  error_fields <- c(
    "arithmetic_l2_error_upper", "optimizer_l2_error_upper",
    "link_l2_error_upper", "output_quantization_l2_error_upper")
  errors <- lapply(error_fields, function(field) {
    .dsvert_formal_glm_release_rat(numeric[[field]], field)
  })
  total_error <- .dsvert_formal_glm_release_rat(
    numeric$total_deterministic_l2_error_upper,
    "total deterministic error")
  component_sum <- Reduce(.dsvert_glm_rat_add, errors,
                          init = .dsvert_glm_rat("0"))
  if (.dsvert_glm_rat_cmp(total_error, component_sum) < 0L) {
    .dsvert_formal_glm_release_stop(
      "The formal GLM total numeric error excludes a component",
      "release_numeric_certificate_invalid")
  }
  numeric_hash <- .dsvert_formal_glm_release_hash(
    numeric$certificate_sha256, "numeric certificate hash")
  postprocessing <- .dsvert_formal_glm_release_postprocessing_contract(
    release$postprocessing_contract)

  .dsvert_formal_glm_release_verify_signatures(
    release, unsigned, designated, pinset)

  source_boxes <- lapply(schema$estimand$coefficients, function(coefficient) {
    .dsvert_glm_rat(coefficient$box_abs)
  })
  two <- .dsvert_glm_bn("2")
  upper_rationals <- lapply(shifted_upper, .dsvert_glm_rat)
  even_bounds <- vapply(upper_rationals, function(value) {
    value$sign > 0L && value$denominator == .dsvert_glm_bn_one() &&
      value$numerator %% two == .dsvert_glm_bn_zero()
  }, logical(1L))
  if (any(!even_bounds)) {
    .dsvert_formal_glm_release_stop(
      "A shifted formal GLM coefficient box is not a positive even lattice",
      "release_projection_invalid")
  }
  box_steps <- lapply(upper_rationals, function(value) {
    .dsvert_glm_rat_div(value, "2")
  })
  shifted_rationals <- lapply(shifted_steps, .dsvert_glm_rat)
  signed_steps <- lapply(seq_along(shifted_rationals), function(index) {
    .dsvert_glm_rat_sub(shifted_rationals[[index]], box_steps[[index]])
  })
  coefficient_rationals <- lapply(signed_steps, function(value) {
    .dsvert_glm_rat_div(value, scale)
  })
  lattice_boxes <- lapply(box_steps, function(value) {
    .dsvert_glm_rat_div(value, scale)
  })
  one_step <- .dsvert_glm_rat_div("1", scale)
  valid_projection <- vapply(seq_along(coefficient_rationals), function(index) {
    shifted_in_range <-
      .dsvert_glm_rat_cmp(shifted_rationals[[index]], "0") >= 0L &&
      .dsvert_glm_rat_cmp(shifted_rationals[[index]],
                          upper_rationals[[index]]) <= 0L
    outward_box <-
      .dsvert_glm_rat_cmp(lattice_boxes[[index]],
                          source_boxes[[index]]) >= 0L &&
      .dsvert_glm_rat_cmp(
        .dsvert_glm_rat_sub(lattice_boxes[[index]],
                            source_boxes[[index]]), one_step) < 0L
    shifted_in_range && outward_box
  }, logical(1L))
  if (any(!valid_projection)) {
    .dsvert_formal_glm_release_stop(
      "A formal GLM coefficient escaped its shifted public box",
      "release_projection_invalid")
  }
  coefficients <- vapply(
    coefficient_rationals, .dsvert_glm_rat_double, numeric(1L))
  names(coefficients) <- order
  if (any(!is.finite(coefficients))) {
    .dsvert_formal_glm_release_stop(
      "A formal GLM coefficient is not representable as a finite double",
      "release_numeric_unrepresentable")
  }
  signed_step_text <- vapply(signed_steps, function(value) {
    .dsvert_glm_rat_json(value)$numerator
  }, character(1L))
  boxes_double <- vapply(
    lattice_boxes, .dsvert_glm_rat_double, numeric(1L))
  names(boxes_double) <- order
  scale_double <- .dsvert_glm_rat_double(scale)
  variance_double <- .dsvert_glm_rat_double(variance_upper) /
    (scale_double^2)
  radius <- .dsvert_glm_rat_double(
    .dsvert_glm_rat_div(simultaneous_steps, scale))
  if (!is.finite(scale_double) || scale_double <= 0 ||
      !is.finite(variance_double) || variance_double <= 0 ||
      !is.finite(radius) || radius <= 0) {
    .dsvert_formal_glm_release_stop(
      "The formal GLM mechanism accuracy is not representable",
      "release_numeric_unrepresentable")
  }
  mechanism_covariance <- diag(variance_double, length(order))
  dimnames(mechanism_covariance) <- list(order, order)
  mechanism_se <- stats::setNames(
    rep(sqrt(variance_double), length(order)), order)
  mechanism_lower <- pmax(coefficients - radius, -boxes_double)
  mechanism_upper <- pmin(coefficients + radius, boxes_double)

  signed_canonical <- .dsvert_joint_dp_client_json(release)
  public_release_sha256 <- digest::digest(
    signed_canonical, algo = "sha256", serialize = FALSE)
  accuracy <- list(
    arithmetic_l2_error_upper = .dsvert_glm_rat_double(errors[[1L]]),
    optimizer_l2_error_upper = .dsvert_glm_rat_double(errors[[2L]]),
    link_l2_error_upper = .dsvert_glm_rat_double(errors[[3L]]),
    output_quantization_l2_error_upper =
      .dsvert_glm_rat_double(errors[[4L]]),
    total_deterministic_l2_error_upper =
      .dsvert_glm_rat_double(total_error),
    mechanism_simultaneous_95_abs = radius,
    numeric_certificate_sha256 = numeric_hash,
    no_wrap_certified = TRUE,
    ring_bits = ring_bits)
  privacy <- list(
    mechanism = mechanism$mechanism,
    backend = mechanism$backend,
    sampler = mechanism$sampler,
    epsilon = .dsvert_glm_rat_double(epsilon),
    delta = .dsvert_glm_rat_double(delta),
    l2_sensitivity_steps = sensitivity_steps,
    nominal_variance_multiplier = 2L,
    operation_limit = FALSE,
    request_limit = FALSE,
    history_can_deny_operation = FALSE,
    history_is_accounting_only = TRUE)
  result <- list(
    status = "dp_regularized_point_estimate",
    family = release$family,
    link = schema$estimand$link,
    estimand = schema$estimand$target,
    estimability = "unique_regularized_box_constrained_target",
    coefficients = coefficients,
    coefficient_lattice_steps = signed_step_text,
    shifted_coefficient_lattice_steps = shifted_steps,
    shifted_upper_bounds_steps = shifted_upper,
    output_lattice_scale = scale_text,
    coefficient_box = boxes_double,
    mechanism_covariance_upper = mechanism_covariance,
    mechanism_std_error_upper = mechanism_se,
    mechanism_region = data.frame(
      estimate = unname(coefficients),
      lower = unname(mechanism_lower),
      upper = unname(mechanism_upper),
      row.names = order, check.names = FALSE),
    covariance = NULL,
    std_errors = NULL,
    deviance = NA_real_,
    log_likelihood = NA_real_,
    aic = NA_real_,
    fitted.values = NULL,
    residuals = NULL,
    statistic_status = postprocessing,
    accuracy = accuracy,
    privacy = privacy,
    privacy_epoch_sha256 = identity$privacy_epoch_sha256,
    release_instance_id = identity$release_instance_id,
    release_contract_sha256 = identity$release_contract_sha256,
    public_release_sha256 = public_release_sha256,
    protected_path = list(
      phase18_to_phase19 = phase$protected_path,
      protected_data_e2e_verified = TRUE,
      hidden_execution_validity_consumed = TRUE,
      invalid_execution_maps_to_neutral_before_noise = TRUE,
      intermediate_openings = intermediate_openings,
      final_openings = final_openings,
      custodian_count = length(custodians),
      designated_compute_peers = designated,
      all_k_manifest_admission_verified = TRUE),
    root_rotation = "new_visible_composed_release_never_blocks",
    schema_sha256 = compilation$sha256,
    schema = schema,
    release_certificate = unsigned,
    internal_review_ready = TRUE,
    production_ready = FALSE)
  class(result) <- c("dsvert_formal_dp_glm", "list")
  result
}

.dsvert_formal_glm_statistic <- function(result, statistic) {
  if (!inherits(result, "dsvert_formal_dp_glm")) {
    stop("result must be an internal formal DP GLM object", call. = FALSE)
  }
  supported <- c(
    "sampling_covariance", "standard_errors", "log_likelihood",
    "deviance", "aic", "wald", "likelihood_ratio",
    "protected_residuals")
  if (!is.character(statistic) || length(statistic) != 1L ||
      is.na(statistic) || !statistic %in% supported) {
    stop("Unknown formal GLM statistic", call. = FALSE)
  }
  reason <- if (identical(statistic, "protected_residuals")) {
    paste0("Protected fitted values and residuals are never reconstructed; ",
           "a separately registered bounded joint-DP diagnostic artifact ",
           "would be required")
  } else {
    paste0("The coefficient-only DP release cannot support ", statistic,
           "; a separately calibrated joint-DP information/likelihood ",
           "artifact is required")
  }
  structure(list(
    statistic = statistic,
    status = .DSVERT_FORMAL_GLM_DP_INFERENCE_REQUIRED,
    reason = reason,
    additional_server_calls = 0L,
    fallback_used = FALSE,
    public_release_sha256 = result$public_release_sha256),
    class = c("dsvert_formal_glm_unavailable", "list"))
}

.dsvert_formal_glm_effects <- function(result) {
  if (!inherits(result, "dsvert_formal_dp_glm")) {
    stop("result must be an internal formal DP GLM object", call. = FALSE)
  }
  region <- result$mechanism_region
  scale <- if (identical(result$family, "binomial")) {
    "odds_ratio"
  } else {
    "incidence_rate_ratio"
  }
  structure(list(
    scale = scale,
    estimate = stats::setNames(exp(region$estimate), rownames(region)),
    mechanism_lower = stats::setNames(exp(region$lower), rownames(region)),
    mechanism_upper = stats::setNames(exp(region$upper), rownames(region)),
    uncertainty = "simultaneous_dp_mechanism_only_not_sampling"),
    class = c("dsvert_formal_glm_effects", "list"))
}

.dsvert_formal_glm_public_design <- function(result, newdata,
                                             offset = NULL) {
  if (!is.data.frame(newdata) || nrow(newdata) < 1L) {
    stop("newdata must be a non-empty public data frame", call. = FALSE)
  }
  schema <- result$schema
  terms <- schema$estimand$coefficients
  design <- matrix(0, nrow = nrow(newdata), ncol = length(terms),
                   dimnames = list(NULL, names(result$coefficients)))
  clipped <- matrix(FALSE, nrow = nrow(newdata), ncol = length(terms),
                    dimnames = dimnames(design))
  for (index in seq_along(terms)) {
    term <- terms[[index]]$term
    if (identical(term$kind, "intercept")) {
      design[, index] <- 1
    } else if (identical(term$kind, "numeric")) {
      name <- term$source_column
      if (!name %in% names(newdata) || !is.numeric(newdata[[name]]) ||
          length(newdata[[name]]) != nrow(newdata) ||
          anyNA(newdata[[name]]) || any(!is.finite(newdata[[name]]))) {
        stop("Public prediction column '", name,
             "' must be finite numeric", call. = FALSE)
      }
      lower <- .dsvert_glm_rat_double(term$clipping_lower)
      upper <- .dsvert_glm_rat_double(term$clipping_upper)
      raw <- as.numeric(newdata[[name]])
      clipped[, index] <- raw < lower | raw > upper
      design[, index] <- pmax(lower, pmin(upper, raw))
    } else if (identical(term$kind, "categorical_indicator")) {
      name <- term$source_column
      if (!name %in% names(newdata) || length(newdata[[name]]) != nrow(newdata) ||
          anyNA(newdata[[name]])) {
        stop("Public prediction column '", name,
             "' must contain registered levels", call. = FALSE)
      }
      values <- as.character(newdata[[name]])
      registered <- unlist(
        schema$estimand$column_registry[[name]]$levels,
        use.names = FALSE)
      if (any(!values %in% registered)) {
        stop("Public prediction column '", name,
             "' contains an unregistered level", call. = FALSE)
      }
      design[, index] <- as.numeric(values == term$source_level)
    } else {
      stop("Unsupported formal GLM public prediction term", call. = FALSE)
    }
  }
  offset_policy <- schema$estimand$offset
  if (is.null(offset)) {
    if (!identical(offset_policy$mode, "none")) {
      stop("A public bounded offset is required for prediction", call. = FALSE)
    }
    offset <- rep(0, nrow(newdata))
    offset_clipped <- rep(FALSE, nrow(newdata))
  } else {
    if (identical(offset_policy$mode, "none")) {
      stop("This formal GLM has no registered prediction offset",
           call. = FALSE)
    }
    if (!is.numeric(offset) || !length(offset) %in% c(1L, nrow(newdata)) ||
        anyNA(offset) || any(!is.finite(offset))) {
      stop("offset must be one finite public value per prediction row",
           call. = FALSE)
    }
    offset <- rep(as.numeric(offset), length.out = nrow(newdata))
    lower <- .dsvert_glm_rat_double(offset_policy$lower)
    upper <- .dsvert_glm_rat_double(offset_policy$upper)
    offset_clipped <- offset < lower | offset > upper
    offset <- pmax(lower, pmin(upper, offset))
  }
  list(design = design, offset = offset, clipped = clipped,
       offset_clipped = offset_clipped)
}

.dsvert_formal_glm_pwl_public <- function(eta, table) {
  knots <- vapply(table$knots, .dsvert_glm_rat_double, numeric(1L))
  values <- vapply(table$values, .dsvert_glm_rat_double, numeric(1L))
  slopes <- vapply(table$slopes, .dsvert_glm_rat_double, numeric(1L))
  tolerance <- 16 * .Machine$double.eps *
    pmax(1, max(abs(knots)), abs(eta))
  if (any(!is.finite(eta)) ||
      any(eta < knots[[1L]] - tolerance) ||
      any(eta > knots[[length(knots)]] + tolerance)) {
    stop("Public prediction escaped the certified link-surrogate domain",
         call. = FALSE)
  }
  eta <- pmax(knots[[1L]], pmin(knots[[length(knots)]], eta))
  segment <- findInterval(
    eta, knots, rightmost.closed = TRUE, all.inside = TRUE)
  values[segment] + slopes[segment] * (eta - knots[segment])
}

.dsvert_formal_glm_predict_public <- function(
    result, newdata, type = c("response", "link"),
    link = c("surrogate", "canonical"), offset = NULL) {
  if (!inherits(result, "dsvert_formal_dp_glm")) {
    stop("result must be an internal formal DP GLM object", call. = FALSE)
  }
  type <- match.arg(type)
  link <- match.arg(link)
  public <- .dsvert_formal_glm_public_design(result, newdata, offset)
  eta <- as.vector(public$offset + public$design %*% result$coefficients)
  coefficient_lower <- result$mechanism_region$lower
  coefficient_upper <- result$mechanism_region$upper
  positive <- public$design >= 0
  eta_lower <- public$offset + rowSums(
    ifelse(positive,
           sweep(public$design, 2L, coefficient_lower, `*`),
           sweep(public$design, 2L, coefficient_upper, `*`)))
  eta_upper <- public$offset + rowSums(
    ifelse(positive,
           sweep(public$design, 2L, coefficient_upper, `*`),
           sweep(public$design, 2L, coefficient_lower, `*`)))
  transform <- function(value) {
    if (identical(type, "link")) return(value)
    if (identical(link, "surrogate")) {
      return(.dsvert_formal_glm_pwl_public(
        value, result$schema$link_surrogate))
    }
    if (identical(result$family, "binomial")) stats::plogis(value) else
      exp(value)
  }
  fit <- transform(eta)
  mechanism_lower <- transform(eta_lower)
  mechanism_upper <- transform(eta_upper)
  if (any(!is.finite(fit)) || any(!is.finite(mechanism_lower)) ||
      any(!is.finite(mechanism_upper))) {
    stop("Public prediction is not representable as a finite double",
         call. = FALSE)
  }
  structure(data.frame(
    fit = fit,
    mechanism_lower = mechanism_lower,
    mechanism_upper = mechanism_upper,
    public_input_clipped = rowSums(public$clipped) > 0L |
      public$offset_clipped,
    check.names = FALSE),
    type = type,
    link = link,
    uncertainty = "simultaneous_dp_mechanism_only_not_sampling",
    protected_data_used = FALSE,
    additional_server_calls = 0L)
}

# Machine-readable hand-off to the separately reviewed inferential phase.  It
# derives its shape solely from the signed compilation; it is not an analyst
# request builder and it does not assert that any listed coordinate is already
# materialized or released.
.dsvert_formal_glm_joint_inference_requirements <- function(compilation) {
  .dsvert_formal_glm_validate_compilation(compilation)
  schema <- compilation$unsigned_schema
  if (!schema$estimand$family %in% c("binomial", "poisson")) {
    .dsvert_formal_glm_release_stop(
      "The formal inference artifact supports binomial or Poisson only",
      "inference_artifact_family_invalid")
  }
  coefficients <- unlist(
    schema$estimand$coefficient_order, use.names = FALSE)
  p <- length(coefficients)
  triangle <- p * (p + 1) / 2
  blocks <- list(
    information_upper_triangle = list(
      coordinates = as.numeric(triangle),
      semantics = paste0(
        "patient_sum_penalized_surrogate_total_curvature_",
        "upper_triangle")),
    score_meat_upper_triangle = list(
      coordinates = as.numeric(triangle),
      semantics = paste0(
        "patient_sum_bounded_patient_score_outer_products_hc0_",
        "upper_triangle")),
    score_vector = list(
      coordinates = as.numeric(p),
      semantics = "patient_sum_bounded_score_at_fixed_public_dp_beta"),
    canonical_bounded_log_likelihood_at_dp_beta = list(
      coordinates = 1,
      semantics = paste0(
        "bounded_canonical_loglik_at_fixed_public_dp_beta_not_",
        "optimized_pwl_loglik")),
    fitted_surrogate_loss = list(
      coordinates = 1,
      semantics = "patient_sum_integrated_pwl_penalized_surrogate_loss"),
    admitted_n = list(
      coordinates = 1,
      semantics = "dp_noisy_active_patient_count"))
  coordinate_count <- sum(vapply(
    blocks, `[[`, numeric(1L), "coordinates"))
  .dsvert_joint_dp_client_canonical(list(
    version = "dsvert-formal-glm-joint-inference-requirements-v1",
    status = "requires_separate_protected_joint_dp_implementation",
    schema_sha256 = compilation$sha256,
    family = schema$estimand$family,
    coefficient_order = as.list(coefficients),
    coordinate_count = as.numeric(coordinate_count),
    blocks = blocks,
    derivation_authority = list(
      bounds = "server_from_cross_signed_schema",
      contribution_caps = "server_from_patient_level_materializer",
      sensitivity = paste0(
        "server_exact_rational_stacked_l1_upper_bound_for_joint_l2_",
        "at_fixed_public_dp_beta_with_outward_lattice_rounding"),
      analyst_sensitivity_fields = FALSE,
      analyst_bound_fields = FALSE,
      analyst_contribution_cap_fields = FALSE),
    protected_execution = list(
      required_input = "phase18_v2_to_phase19_v2_verified_blocks_only",
      hidden_execution_validity_consumed_before_noise = TRUE,
      invalid_execution_neutralized_before_noise = TRUE,
      intermediate_openings = 0,
      one_stacked_dp_vector_opening = TRUE,
      k_supported = "all_K_at_least_2_with_exactly_two_pinned_compute_peers"),
    scientific_labels = list(
      penalized_curvature_inverse = paste0(
        "H_total_inverse_working_laplace_scale_not_sampling_",
        "covariance_or_standard_errors"),
      sampling_covariance = paste0(
        "robust_h_total_inverse_J_h_total_inverse_hc0_conditional_",
        "research_estimate"),
      model_based_covariance = paste0(
        "unavailable_without_likelihood_bread_and_model_variance_",
        "meat"),
      loss = "integrated_pwl_surrogate_loss_not_ordinary_logLik",
      aic = "not_classical_aic_for_ridge_box_fixed_iteration_target",
      wald = "unavailable_until_noise_aware_coverage_validation",
      likelihood_ratio = paste0(
        "unavailable_for_penalized_surrogate_without_separate_",
        "mechanism_aware_null_calibration")),
    operation_limit = FALSE,
    request_limit = FALSE,
    history_can_deny_operation = FALSE,
    production_ready = FALSE))
}
