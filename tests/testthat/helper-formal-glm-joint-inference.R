formal_glm_joint_test_hash <- function(value) {
  digest::digest(value, algo = "sha256", serialize = FALSE)
}

formal_glm_joint_test_rat <- function(numerator, denominator = "1") {
  list(numerator = numerator, denominator = denominator)
}

formal_glm_joint_test_keys <- function(k) {
  peers <- paste0("site_", letters[seq_len(k)])
  private <- stats::setNames(lapply(peers, function(peer) {
    openssl::ed25519_keygen()
  }), peers)
  pinset <- vapply(private, function(key) {
    dsVertClient:::.dsvert_exact_gc_b64url_encode(
      tail(as.raw(as.list(key)$pubkey), 32L))
  }, character(1L))
  list(private = private,
       pinset = pinset[order(names(pinset), method = "radix")])
}

formal_glm_joint_test_model <- function(family = "binomial",
                                        formula = "y ~ x + group",
                                        capacity = 6L,
                                        weights = "unit",
                                        coefficient_box = "0.01",
                                        ridge = "500000") {
  columns <- list(
    y = if (identical(family, "binomial")) {
      list(kind = "binary", owner = "site_a")
    } else {
      list(kind = "count", owner = "site_a", upper = "8")
    },
    x = list(kind = "numeric", owner = "site_b",
             lower = "-1", upper = "1"),
    group = list(kind = "categorical", owner = "site_b",
                 levels = c("control", "treated"),
                 reference = "control", contrast = "treatment"))
  weight_policy <- list(mode = "unit", column = NULL)
  if (identical(weights, "bounded")) {
    columns$w <- list(kind = "weight", owner = "site_a", upper = "2")
    weight_policy <- list(mode = "bounded_column", column = "w")
  }
  list(
    formula = formula, family = family, capacity = capacity,
    columns = columns, adjacency = "add_remove",
    patient_collapse = list(
      unit = "aligned_patient", repeated_records = "reject_duplicates",
      row_order_invariant = TRUE, max_records_per_unit = 1L,
      conflict_policy = "zero_weight"),
    missingness = "complete_tuple_zero_weight",
    clipping = list(
      numeric = "clamp_then_quantize",
      categorical = "registered_level_or_zero_weight",
      binomial = "binary_or_zero_weight",
      poisson_count = "integer_then_patient_cap",
      offset = "clamp_then_quantize", weight = "clamp_then_quantize"),
    weights = weight_policy,
    offset = list(mode = "none", column = NULL),
    coefficient_box = coefficient_box, ridge = ridge,
    optimizer = list(
      alpha = "0.000001",
      iterations = 8L, start = "zero"),
    numeric = list(
      x_fraction_bits = 12L, offset_fraction_bits = 12L,
      beta_fraction_bits = 16L, link_fraction_bits = 20L,
      working_fraction_bits = 32L, reference_precision_bits = 64L),
    link = list(segments = 4L,
                construction = "uniform_monotone_dyadic_linear"),
    privacy = list(
      mechanism = "joint_discrete_gaussian_one_global_draw",
      epsilon = "1", delta = "1e-6",
      allocation = "one_stacked_capsule_vector"))
}

formal_glm_joint_test_authority <- function(k, pinset,
                                            snapshot = "snapshot-test") {
  peers <- paste0("site_", letters[seq_len(k)])
  list(
    consortium_id = "consortium-test", capsule_id = paste0("capsule-k", k),
    logical_snapshot = snapshot,
    policy_contract_sha256 = formal_glm_joint_test_hash(
      paste0("policy-", k)),
    custodian_peers = peers, designated_peers = peers[1:2],
    pinset_sha256 = dsVertClient:::.dsvert_dp_pinset_hash(pinset))
}

formal_glm_joint_test_data <- function(family = "binomial",
                                       weights = "unit") {
  result <- if (identical(family, "binomial")) {
    data.frame(
      y = c(0, 0, 1, 1, 0, 1), x = c(-1, -0.5, 0, 0.4, 0.8, 1),
      group = c("control", "treated", "control", "treated",
                "control", "treated"))
  } else {
    data.frame(
      y = c(0, 1, 1, 2, 3, 5), x = c(-1, -0.5, 0, 0.4, 0.8, 1),
      group = c("control", "treated", "control", "treated",
                "control", "treated"))
  }
  if (identical(weights, "bounded")) {
    result$w <- c(0.5, 1, 2, 1.5, 0.75, 1.25)
  }
  result
}

formal_glm_joint_test_point_release <- function(
    compilation, keys, coefficients, label = "point", epoch = "point-epoch") {
  schema <- compilation$unsigned_schema
  peers <- unlist(schema$authority$custodian_peers, use.names = FALSE)
  designated <- unlist(schema$authority$designated_peers, use.names = FALSE)
  coefficient_order <- unlist(
    schema$estimand$coefficient_order, use.names = FALSE)
  scale <- 2^16
  coefficients <- coefficients[coefficient_order]
  signed_steps <- round(coefficients * scale)
  box_steps <- vapply(schema$estimand$coefficients, function(coefficient) {
    ceiling(dsVertClient:::.dsvert_glm_rat_double(coefficient$box_abs) * scale)
  }, numeric(1L))
  shifted <- sprintf("%.0f", signed_steps + box_steps)
  shifted_upper <- sprintf("%.0f", 2 * box_steps)
  materialization <- formal_glm_joint_test_hash(
    paste0(label, "/materialization"))
  unsigned <- list(
    version = "dsvert-formal-glm-dp-public-release-v1",
    schema_sha256 = compilation$sha256,
    family = schema$estimand$family,
    coefficient_order = as.list(coefficient_order),
    shifted_coefficient_lattice_steps = as.list(unname(shifted)),
    shifted_upper_bounds_steps = as.list(unname(shifted_upper)),
    output_lattice_scale = sprintf("%.0f", scale),
    phase19_binding = list(
      version = "dsvert-formal-glm-phase19-public-binding-v1",
      protected_path = "phase18_v2_to_phase19_v2",
      post_execution_root_sha256 = formal_glm_joint_test_hash(
        paste0(label, "/post")),
      execution_receipt_pair_sha256 = formal_glm_joint_test_hash(
        paste0(label, "/receipts")),
      global_materialization_root_sha256 = materialization,
      protected_data_e2e_verified = TRUE,
      hidden_execution_validity_consumed = TRUE,
      invalid_execution_maps_to_neutral_before_noise = TRUE,
      intermediate_openings = 0, final_openings = 1),
    authority = list(
      pinset_sha256 = schema$authority$pinset_sha256,
      custodian_peers = as.list(peers),
      designated_compute_peers = as.list(designated),
      all_k_manifest_admission_verified = TRUE,
      all_k_manifest_signature_set_sha256 = formal_glm_joint_test_hash(
        paste0(label, "/all-k")),
      server_derived_bounds = TRUE, server_derived_sensitivity = TRUE,
      server_derived_contribution_caps = TRUE,
      operation_limit = FALSE, request_limit = FALSE,
      history_can_deny_operation = FALSE),
    release_identity = list(
      release_instance_id = formal_glm_joint_test_hash(
        paste0(label, "/", epoch)),
      release_contract_sha256 = formal_glm_joint_test_hash(
        paste0(label, "/contract/", epoch)),
      privacy_epoch_sha256 = formal_glm_joint_test_hash(epoch),
      materialization_root_sha256 = materialization,
      final_vector_root_sha256 = formal_glm_joint_test_hash(
        paste0(label, "/vector/", epoch)),
      coordinate_order_sha256 = formal_glm_joint_test_hash(
        paste(coefficient_order, collapse = "\n")),
      sticky_absolute_coordinates = TRUE,
      retry_restart_rechunk_same_release = TRUE,
      root_rotation_starts_visible_composed_release = TRUE),
    mechanism = list(
      backend = paste0(
        "independent_full_global_dyadic_discrete_gaussian_",
        "tv_bounded_ring128_v2"),
      mechanism = "dyadic_discrete_gaussian_truncated_tv_bounded",
      sampler = paste0(
        "cks-target-outward-rational-dyadic-cdf-hkdf-sha256-",
        "chacha20-coordinate-domain-v2"),
      epsilon = formal_glm_joint_test_rat("1"),
      delta = formal_glm_joint_test_rat("1", "1000000"),
      l2_sensitivity_steps = "128",
      sigma_squared_steps = formal_glm_joint_test_rat("16"),
      mechanism_variance_upper_steps = formal_glm_joint_test_rat("32"),
      simultaneous_95_abs_steps = "64", nominal_variance_multiplier = 2,
      release_delta_aggregation = "max_per_peer_not_sum",
      complete_epsilon_delta_per_peer = TRUE,
      epsilon_divided_by_peer_count = FALSE,
      finite_support_transfer_charged = TRUE,
      one_global_vector_release = TRUE, automatic_fallback_used = FALSE),
    numeric_certificate = list(
      ring_bits = 128, no_wrap_certified = TRUE,
      final_projection = "nonnegative_shifted_public_coefficient_box_v1",
      arithmetic_l2_error_upper = formal_glm_joint_test_rat("1", "65536"),
      optimizer_l2_error_upper = formal_glm_joint_test_rat("1", "1024"),
      link_l2_error_upper = formal_glm_joint_test_rat("1", "512"),
      output_quantization_l2_error_upper = formal_glm_joint_test_rat(
        "1", "32768"),
      total_deterministic_l2_error_upper = formal_glm_joint_test_rat(
        "1", "128"),
      certificate_sha256 = formal_glm_joint_test_hash(
        paste0(label, "/numeric"))),
    postprocessing_contract = list(
      coefficient_point_estimate = "available",
      mechanism_covariance = "upper_bound_available",
      sampling_covariance = "requires_joint_dp_inference_artifact",
      log_likelihood = "requires_joint_dp_inference_artifact",
      deviance = "requires_joint_dp_inference_artifact",
      aic = "requires_joint_dp_inference_artifact",
      wald = "requires_joint_dp_inference_artifact",
      likelihood_ratio = "requires_joint_dp_inference_artifact",
      prediction = "public_covariates_only",
      protected_residual_diagnostics = "unavailable"))
  message <- charToRaw(paste0(
    "dsVert/formal-glm/dp-public-release/v1|",
    dsVertClient:::.dsvert_joint_dp_client_json(unsigned)))
  signatures <- lapply(designated, function(peer) {
    list(signer = peer,
         signature = dsVertClient:::.dsvert_exact_gc_b64url_encode(
           openssl::ed25519_sign(message, keys$private[[peer]])))
  })
  c(unsigned, list(signatures = signatures))
}

formal_glm_joint_test_triangle <- function(matrix) {
  result <- numeric()
  for (row in seq_len(nrow(matrix))) {
    for (column in seq.int(row, ncol(matrix))) {
      result <- c(result, matrix[row, column])
    }
  }
  result
}

formal_glm_joint_test_capabilities <- function(
    j = TRUE, score = TRUE, loglik = TRUE, loss = TRUE,
    admitted_n = TRUE) {
  released <- "released"
  absent <- "not_promoted_sensitivity_not_certified"
  list(
    information_upper_triangle = released,
    score_meat_upper_triangle = if (j) released else absent,
    score_vector = if (score) released else absent,
    canonical_bounded_log_likelihood_at_dp_beta =
      if (loglik) released else absent,
    integrated_pwl_surrogate_loss_at_dp_beta =
      if (loss) released else absent,
    admitted_n = if (admitted_n) released else absent)
}

formal_glm_joint_test_release <- function(
    point, compilation, keys, oracle,
    capabilities = formal_glm_joint_test_capabilities(),
    label = "inference", epoch = "inference-epoch") {
  schema <- compilation$unsigned_schema
  peers <- unlist(schema$authority$custodian_peers, use.names = FALSE)
  designated <- unlist(schema$authority$designated_peers, use.names = FALSE)
  coefficients <- unlist(
    schema$estimand$coefficient_order, use.names = FALSE)
  coordinates <- dsVertClient:::.dsvert_formal_glm_inference_coordinate_layout(
    coefficients, capabilities)
  blocks <- list(
    information_upper_triangle = formal_glm_joint_test_triangle(
      oracle$information),
    score_meat_upper_triangle = formal_glm_joint_test_triangle(
      oracle$score_meat),
    score_vector = unname(oracle$score),
    canonical_bounded_log_likelihood_at_dp_beta =
      oracle$canonical_bounded_log_likelihood_at_dp_beta,
    integrated_pwl_surrogate_loss_at_dp_beta =
      oracle$integrated_pwl_surrogate_loss_at_dp_beta,
    admitted_n = oracle$admitted_n)
  values <- unname(unlist(blocks[vapply(names(blocks), function(field) {
    identical(capabilities[[field]], "released")
  }, logical(1L))], use.names = FALSE))
  scale <- 2^20
  scale_text <- sprintf("%.0f", scale)
  server_plan <- dsVertClient:::.dsvert_formal_glm_inference_server_plan(
    point, capabilities, scale_text)
  raw_steps <- round(values * scale)
  lower_text <- unlist(
    server_plan$coordinate_lower_bounds_steps, use.names = FALSE)
  upper_text <- unlist(
    server_plan$shifted_upper_bounds_steps, use.names = FALSE)
  shifted_text <- vapply(seq_along(raw_steps), function(index) {
    dsVertClient:::.dsvert_glm_rat_json(
      dsVertClient:::.dsvert_glm_rat_sub(
        sprintf("%.0f", raw_steps[[index]]), lower_text[[index]]))$numerator
  }, character(1L))
  coordinate_order_sha256 <- dsVertClient:::.dsvert_formal_glm_inference_hash(
    as.list(coordinates))
  materialization <- point$release_certificate$release_identity$
    materialization_root_sha256
  base <- dsVertClient:::.dsvert_formal_glm_inference_base_expected(point)
  sensitivity <- server_plan$l2_sensitivity_steps
  coordinate_bounds_sha256 <- dsVertClient:::.dsvert_formal_glm_inference_hash(
    list(lower = as.list(lower_text),
         shifted_upper = as.list(upper_text)))
  sensitivity_sha256 <- dsVertClient:::.dsvert_formal_glm_inference_hash(list(
    coordinate_order_sha256 = coordinate_order_sha256,
    l2_sensitivity_steps = sensitivity,
    output_lattice_scale = scale_text))
  coordinate_count <- length(coordinates)
  total_error <- dsVertClient:::.dsvert_glm_rat_json(
    dsVertClient:::.dsvert_glm_rat_div(
      as.character(coordinate_count), scale_text))
  epsilon <- dsVertClient:::.dsvert_glm_rat_div("1", "2")
  core_delta <- dsVertClient:::.dsvert_glm_rat_div("1", "2000000")
  calibration_log_upper <- dsVertClient:::.dsvert_glm_rat_log_interval(
    dsVertClient:::.dsvert_glm_rat_div(
      "5", dsVertClient:::.dsvert_glm_rat_mul("4", core_delta)),
    schema$numeric$reference_precision_bits)$upper
  minimum_sigma <- dsVertClient:::.dsvert_glm_rat_json(
    dsVertClient:::.dsvert_glm_rat_div(
      dsVertClient:::.dsvert_glm_rat_mul(
        dsVertClient:::.dsvert_glm_rat_pow(sensitivity, 2L),
        dsVertClient:::.dsvert_glm_rat_mul("2", calibration_log_upper)),
      dsVertClient:::.dsvert_glm_rat_pow(epsilon, 2L)))
  variance_upper <- dsVertClient:::.dsvert_glm_rat_mul(minimum_sigma, "2")
  radius_threshold <- dsVertClient:::.dsvert_glm_rat_mul(
    as.character(20L * coordinate_count), variance_upper)
  radius_steps_number <- ceiling(sqrt(
    dsVertClient:::.dsvert_glm_rat_double(radius_threshold)))
  while (dsVertClient:::.dsvert_glm_rat_cmp(
      dsVertClient:::.dsvert_glm_rat_pow(
        sprintf("%.0f", radius_steps_number), 2L),
      radius_threshold) < 0L) {
    radius_steps_number <- radius_steps_number + 1
  }
  radius_steps <- sprintf("%.0f", radius_steps_number)
  unsigned <- list(
    version = "dsvert-formal-glm-joint-inference-release-v1",
    schema_sha256 = compilation$sha256,
    family = schema$estimand$family,
    coefficient_order = as.list(coefficients),
    base_release_binding = base,
    coordinate_contract = list(
      version = "dsvert-formal-glm-joint-inference-coordinate-contract-v1",
      coordinate_order = as.list(coordinates),
      coordinate_order_sha256 = coordinate_order_sha256,
      upper_triangle_order =
        "row_major_i_then_j_for_i_less_than_or_equal_to_j",
      block_capabilities = capabilities),
    shifted_coordinate_lattice_steps = as.list(shifted_text),
    coordinate_lower_bounds_steps = as.list(lower_text),
    shifted_upper_bounds_steps = as.list(upper_text),
    output_lattice_scale = scale_text,
    server_derivation = server_plan,
    protected_binding = list(
      version = "dsvert-formal-glm-joint-inference-protected-binding-v1",
      protected_path = "phase19_verified_beta_to_joint_inference_v1",
      source_materialization_root_sha256 = materialization,
      base_post_execution_root_sha256 =
        point$release_certificate$phase19_binding$post_execution_root_sha256,
      inference_post_execution_root_sha256 = formal_glm_joint_test_hash(
        paste0(label, "/post/", epoch)),
      execution_receipt_pair_sha256 = formal_glm_joint_test_hash(
        paste0(label, "/receipts/", epoch)),
      same_snapshot_as_base = TRUE,
      protected_data_e2e_verified = TRUE,
      hidden_execution_validity_consumed = TRUE,
      invalid_execution_maps_to_neutral_before_noise = TRUE,
      intermediate_openings = 0, final_openings = 1),
    authority = list(
      pinset_sha256 = schema$authority$pinset_sha256,
      custodian_peers = as.list(peers),
      designated_compute_peers = as.list(designated),
      all_k_manifest_admission_verified = TRUE,
      all_k_manifest_signature_set_sha256 = formal_glm_joint_test_hash(
        paste0(label, "/all-k/", epoch)),
      server_derived_bounds = TRUE, server_derived_sensitivity = TRUE,
      server_derived_contribution_caps = TRUE,
      server_derived_block_availability = TRUE,
      operation_limit = FALSE, request_limit = FALSE,
      history_can_deny_operation = FALSE),
    release_identity = list(
      release_instance_id = formal_glm_joint_test_hash(
        paste0(label, "/", epoch)),
      release_contract_sha256 = formal_glm_joint_test_hash(
        paste0(label, "/contract/", epoch)),
      privacy_epoch_sha256 = formal_glm_joint_test_hash(epoch),
      materialization_root_sha256 = materialization,
      base_public_release_sha256 = point$public_release_sha256,
      final_vector_root_sha256 = formal_glm_joint_test_hash(
        paste0(label, "/vector/", epoch)),
      coordinate_order_sha256 = coordinate_order_sha256,
      sticky_absolute_coordinates = TRUE,
      retry_restart_rechunk_same_release = TRUE,
      root_rotation_starts_visible_composed_release = TRUE),
    mechanism = list(
      backend = paste0(
        "independent_full_global_dyadic_discrete_gaussian_",
        "tv_bounded_ring128_v2"),
      mechanism = "dyadic_discrete_gaussian_truncated_tv_bounded",
      sampler = paste0(
        "cks-target-outward-rational-dyadic-cdf-hkdf-sha256-",
        "chacha20-coordinate-domain-v2"),
      epsilon = formal_glm_joint_test_rat("1", "2"),
      delta = formal_glm_joint_test_rat("1", "1000000"),
      l2_sensitivity_steps = sensitivity,
      sigma_squared_steps = minimum_sigma,
      mechanism_variance_upper_steps =
        dsVertClient:::.dsvert_glm_rat_json(variance_upper),
      simultaneous_95_abs_steps = radius_steps,
      nominal_variance_multiplier = 2,
      release_delta_aggregation = "max_per_peer_not_sum",
      complete_epsilon_delta_per_peer = TRUE,
      epsilon_divided_by_peer_count = FALSE,
      finite_support_transfer_charged = TRUE,
      one_global_vector_release = TRUE, automatic_fallback_used = FALSE,
      composition_rule = "basic_sequential_adaptive_composition_exact_sum",
      base_epsilon = point$release_certificate$mechanism$epsilon,
      base_delta = point$release_certificate$mechanism$delta,
      composed_epsilon = formal_glm_joint_test_rat("3", "2"),
      composed_delta = formal_glm_joint_test_rat("1", "500000"),
      beta_release_is_public_adaptive_input = TRUE,
      ledger_accounting_only = TRUE,
      gaussian_core_delta = formal_glm_joint_test_rat("1", "2000000"),
      finite_support_transfer_delta =
        formal_glm_joint_test_rat("1", "2000000"),
      calibration_rule = paste0(
        "sigma_squared_gte_sensitivity_squared_times_2_log_",
        "5_over_4_delta_core_over_epsilon_squared"),
      minimum_sigma_squared_steps = minimum_sigma,
      simultaneous_radius_rule =
        "union_chebyshev_95_percent_all_coordinates"),
    numeric_certificate = list(
      ring_bits = 128, no_wrap_certified = TRUE,
      final_projection = "nonnegative_shifted_public_coordinate_box_v1",
      coordinate_bounds_sha256 = coordinate_bounds_sha256,
      stacked_sensitivity_sha256 = sensitivity_sha256,
      arithmetic_l2_error_upper = formal_glm_joint_test_rat("0"),
      output_quantization_l2_error_upper = total_error,
      total_deterministic_l2_error_upper = total_error,
      certificate_sha256 = formal_glm_joint_test_hash(
        paste0(label, "/numeric/", epoch))),
    scientific_contract =
      dsVertClient:::.dsvert_formal_glm_inference_scientific_contract_expected())
  message <- charToRaw(paste0(
    "dsVert/formal-glm/joint-inference-release/v1|",
    dsVertClient:::.dsvert_joint_dp_client_json(unsigned)))
  signatures <- lapply(designated, function(peer) {
    list(signer = peer,
         signature = dsVertClient:::.dsvert_exact_gc_b64url_encode(
           openssl::ed25519_sign(message, keys$private[[peer]])))
  })
  c(unsigned, list(signatures = signatures))
}

formal_glm_joint_test_resign <- function(release, compilation, keys) {
  unsigned <- release[setdiff(names(release), "signatures")]
  designated <- unlist(
    compilation$unsigned_schema$authority$designated_peers,
    use.names = FALSE)
  message <- charToRaw(paste0(
    "dsVert/formal-glm/joint-inference-release/v1|",
    dsVertClient:::.dsvert_joint_dp_client_json(unsigned)))
  release$signatures <- lapply(designated, function(peer) {
    list(signer = peer,
         signature = dsVertClient:::.dsvert_exact_gc_b64url_encode(
           openssl::ed25519_sign(message, keys$private[[peer]])))
  })
  release
}

formal_glm_joint_test_fixture <- function(
    k = 3L, family = "binomial", weights = "unit", data = NULL,
    capabilities = formal_glm_joint_test_capabilities(),
    point_coefficients = NULL, label = NULL,
    coefficient_box = "0.01", ridge = "500000") {
  keys <- formal_glm_joint_test_keys(k)
  authority <- formal_glm_joint_test_authority(k, keys$pinset)
  model <- formal_glm_joint_test_model(
    family = family, weights = weights,
    coefficient_box = coefficient_box, ridge = ridge)
  compilation <- dsVertClient:::.dsvert_formal_glm_compile(model, authority)
  if (is.null(data)) data <- formal_glm_joint_test_data(family, weights)
  point_oracle <- dsVertClient:::.dsvert_formal_glm_rational_oracle(
    compilation, data)
  if (is.null(point_coefficients)) {
    point_coefficients <- point_oracle$coefficients
  }
  if (is.null(label)) label <- paste0("k", k, "/", family, "/", weights)
  point_release <- formal_glm_joint_test_point_release(
    compilation, keys, point_coefficients,
    label = paste0(label, "/point"))
  point <- dsVertClient:::.dsvert_formal_glm_dp_result(
    compilation, point_release, keys$pinset)
  oracle <- dsVertClient:::.dsvert_formal_glm_joint_inference_oracle(
    point, compilation, data)
  release <- formal_glm_joint_test_release(
    point, compilation, keys, oracle, capabilities,
    label = paste0(label, "/inference"))
  list(keys = keys, compilation = compilation, data = data,
       point_oracle = point_oracle, point_release = point_release,
       point = point, oracle = oracle, release = release,
       capabilities = capabilities)
}
