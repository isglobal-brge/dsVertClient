formal_glm_dp_test_hash <- function(value) {
  digest::digest(value, algo = "sha256", serialize = FALSE)
}

formal_glm_dp_test_model <- function(family = "binomial", capacity = 6L,
                                     formula = "y ~ x + group") {
  list(
    formula = formula,
    family = family,
    capacity = capacity,
    columns = list(
      y = if (identical(family, "binomial")) {
        list(kind = "binary", owner = "site_a")
      } else {
        list(kind = "count", owner = "site_a", upper = "8")
      },
      x = list(kind = "numeric", owner = "site_b",
               lower = "-1", upper = "1"),
      group = list(kind = "categorical", owner = "site_b",
                   levels = c("control", "treated"),
                   reference = "control", contrast = "treatment")),
    adjacency = "add_remove",
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
      offset = "clamp_then_quantize",
      weight = "clamp_then_quantize"),
    weights = list(mode = "unit", column = NULL),
    offset = list(mode = "none", column = NULL),
    coefficient_box = "2",
    ridge = "1",
    optimizer = list(
      alpha = if (identical(family, "poisson")) "0.0001" else "0.1",
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

formal_glm_dp_test_authority <- function(k, pinset) {
  peers <- paste0("site_", letters[seq_len(k)])
  list(
    consortium_id = "consortium-test",
    capsule_id = paste0("capsule-k", k),
    logical_snapshot = "snapshot-test",
    policy_contract_sha256 = formal_glm_dp_test_hash(paste0("policy-", k)),
    custodian_peers = peers,
    designated_peers = peers[1:2],
    pinset_sha256 = dsVertClient:::.dsvert_dp_pinset_hash(pinset))
}

formal_glm_dp_test_keys <- function(k) {
  peers <- paste0("site_", letters[seq_len(k)])
  private <- stats::setNames(lapply(peers, function(peer) {
    openssl::ed25519_keygen()
  }), peers)
  pinset <- vapply(private, function(key) {
    dsVertClient:::.dsvert_exact_gc_b64url_encode(
      tail(as.raw(as.list(key)$pubkey), 32L))
  }, character(1L))
  pinset <- pinset[order(names(pinset), method = "radix")]
  list(private = private, pinset = pinset)
}

formal_glm_dp_test_rat <- function(numerator, denominator = "1") {
  list(numerator = numerator, denominator = denominator)
}

formal_glm_dp_test_release <- function(compilation, keys, coefficients,
                                       label = "release", epoch = "epoch-1") {
  schema <- compilation$unsigned_schema
  peers <- unlist(schema$authority$custodian_peers, use.names = FALSE)
  designated <- unlist(schema$authority$designated_peers, use.names = FALSE)
  order <- unlist(schema$estimand$coefficient_order, use.names = FALSE)
  coefficients <- coefficients[order]
  scale <- 2^16
  signed_steps <- round(coefficients * scale)
  box_steps <- vapply(schema$estimand$coefficients, function(coefficient) {
    ceiling(dsVertClient:::.dsvert_glm_rat_double(coefficient$box_abs) * scale)
  }, numeric(1L))
  shifted_steps <- sprintf("%.0f", signed_steps + box_steps)
  shifted_upper <- sprintf("%.0f", 2 * box_steps)
  unsigned <- list(
    version = "dsvert-formal-glm-dp-public-release-v1",
    schema_sha256 = compilation$sha256,
    family = schema$estimand$family,
    coefficient_order = as.list(order),
    shifted_coefficient_lattice_steps = as.list(unname(shifted_steps)),
    shifted_upper_bounds_steps = as.list(unname(shifted_upper)),
    output_lattice_scale = sprintf("%.0f", scale),
    phase19_binding = list(
      version = "dsvert-formal-glm-phase19-public-binding-v1",
      protected_path = "phase18_v2_to_phase19_v2",
      post_execution_root_sha256 = formal_glm_dp_test_hash(
        paste0(label, "/post")),
      execution_receipt_pair_sha256 = formal_glm_dp_test_hash(
        paste0(label, "/receipts")),
      global_materialization_root_sha256 = formal_glm_dp_test_hash(
        paste0(label, "/materialization")),
      protected_data_e2e_verified = TRUE,
      hidden_execution_validity_consumed = TRUE,
      invalid_execution_maps_to_neutral_before_noise = TRUE,
      intermediate_openings = 0,
      final_openings = 1),
    authority = list(
      pinset_sha256 = schema$authority$pinset_sha256,
      custodian_peers = as.list(peers),
      designated_compute_peers = as.list(designated),
      all_k_manifest_admission_verified = TRUE,
      all_k_manifest_signature_set_sha256 = formal_glm_dp_test_hash(
        paste0(label, "/all-k-manifest-signatures")),
      server_derived_bounds = TRUE,
      server_derived_sensitivity = TRUE,
      server_derived_contribution_caps = TRUE,
      operation_limit = FALSE,
      request_limit = FALSE,
      history_can_deny_operation = FALSE),
    release_identity = list(
      release_instance_id = formal_glm_dp_test_hash(paste0(label, "/", epoch)),
      release_contract_sha256 = formal_glm_dp_test_hash(
        paste0(label, "/contract/", epoch)),
      privacy_epoch_sha256 = formal_glm_dp_test_hash(epoch),
      materialization_root_sha256 = formal_glm_dp_test_hash(
        paste0(label, "/materialization")),
      final_vector_root_sha256 = formal_glm_dp_test_hash(
        paste0(label, "/vector/", epoch)),
      coordinate_order_sha256 = formal_glm_dp_test_hash(
        paste(order, collapse = "\n")),
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
      epsilon = formal_glm_dp_test_rat("1"),
      delta = formal_glm_dp_test_rat("1", "1000000"),
      l2_sensitivity_steps = "128",
      sigma_squared_steps = formal_glm_dp_test_rat("16"),
      mechanism_variance_upper_steps = formal_glm_dp_test_rat("32"),
      simultaneous_95_abs_steps = "64",
      nominal_variance_multiplier = 2,
      release_delta_aggregation = "max_per_peer_not_sum",
      complete_epsilon_delta_per_peer = TRUE,
      epsilon_divided_by_peer_count = FALSE,
      finite_support_transfer_charged = TRUE,
      one_global_vector_release = TRUE,
      automatic_fallback_used = FALSE),
    numeric_certificate = list(
      ring_bits = 128,
      no_wrap_certified = TRUE,
      final_projection = "nonnegative_shifted_public_coefficient_box_v1",
      arithmetic_l2_error_upper = formal_glm_dp_test_rat("1", "65536"),
      optimizer_l2_error_upper = formal_glm_dp_test_rat("1", "1024"),
      link_l2_error_upper = formal_glm_dp_test_rat("1", "512"),
      output_quantization_l2_error_upper = formal_glm_dp_test_rat("1", "32768"),
      total_deterministic_l2_error_upper = formal_glm_dp_test_rat("1", "128"),
      certificate_sha256 = formal_glm_dp_test_hash(
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

formal_glm_dp_test_resign <- function(release, compilation, keys) {
  unsigned <- release[setdiff(names(release), "signatures")]
  designated <- unlist(
    compilation$unsigned_schema$authority$designated_peers,
    use.names = FALSE)
  message <- charToRaw(paste0(
    "dsVert/formal-glm/dp-public-release/v1|",
    dsVertClient:::.dsvert_joint_dp_client_json(unsigned)))
  release$signatures <- lapply(designated, function(peer) {
    list(signer = peer,
         signature = dsVertClient:::.dsvert_exact_gc_b64url_encode(
           openssl::ed25519_sign(message, keys$private[[peer]])))
  })
  release
}

formal_glm_dp_test_fixture <- function(k = 3L, family = "binomial",
                                       data = NULL) {
  keys <- formal_glm_dp_test_keys(k)
  authority <- formal_glm_dp_test_authority(k, keys$pinset)
  model <- formal_glm_dp_test_model(family = family)
  # Only site_a/site_b own columns; additional custodians still participate in
  # the all-K protected fan-in and pinset.
  compilation <- dsVertClient:::.dsvert_formal_glm_compile(model, authority)
  if (is.null(data)) {
    data <- if (identical(family, "binomial")) {
      data.frame(y = c(0, 0, 1, 1, 0, 1),
                 x = c(-1, -0.5, 0, 0.4, 0.8, 1),
                 group = c("control", "treated", "control", "treated",
                           "control", "treated"))
    } else {
      data.frame(y = c(0, 1, 1, 2, 3, 5),
                 x = c(-1, -0.5, 0, 0.4, 0.8, 1),
                 group = c("control", "treated", "control", "treated",
                           "control", "treated"))
    }
  }
  oracle <- dsVertClient:::.dsvert_formal_glm_rational_oracle(
    compilation, data)
  release <- formal_glm_dp_test_release(
    compilation, keys, oracle$coefficients,
    label = paste0("k", k, "/", family))
  list(keys = keys, compilation = compilation, oracle = oracle,
       release = release, data = data)
}

test_that("formal DP GLM result validates K=2, K=3 and K=5 protected releases", {
  for (k in c(2L, 3L, 5L)) {
    fixture <- formal_glm_dp_test_fixture(k = k)
    result <- dsVertClient:::.dsvert_formal_glm_dp_result(
      fixture$compilation, fixture$release, fixture$keys$pinset)

    expect_s3_class(result, "dsvert_formal_dp_glm")
    expect_identical(result$family, "binomial")
    expect_identical(names(result$coefficients),
                     unlist(fixture$release$coefficient_order,
                            use.names = FALSE))
    expect_true(result$protected_path$protected_data_e2e_verified)
    expect_true(result$protected_path$hidden_execution_validity_consumed)
    expect_identical(result$protected_path$intermediate_openings, 0L)
    expect_identical(result$privacy$operation_limit, FALSE)
    expect_identical(result$privacy$request_limit, FALSE)
    expect_identical(result$privacy$history_can_deny_operation, FALSE)
    expect_null(result$std_errors)
    expect_null(result$covariance)
    expect_true(all(is.finite(result$mechanism_std_error_upper)))
    expect_true(all(diag(result$mechanism_covariance_upper) > 0))
    expect_identical(result$production_ready, FALSE)
  }
})

test_that("formal DP GLM decoding agrees with the central rational oracle", {
  for (family in c("binomial", "poisson")) {
    fixture <- formal_glm_dp_test_fixture(family = family)
    result <- dsVertClient:::.dsvert_formal_glm_dp_result(
      fixture$compilation, fixture$release, fixture$keys$pinset)
    quantization <- 1 / as.numeric(fixture$release$output_lattice_scale)
    expect_lte(max(abs(result$coefficients - fixture$oracle$coefficients)),
               quantization)
    oracle_lattice_error <- sqrt(sum(
      (result$coefficients - fixture$oracle$coefficients)^2))
    expect_lte(
      oracle_lattice_error,
      result$accuracy$output_quantization_l2_error_upper)

    prediction <- dsVertClient:::.dsvert_formal_glm_predict_public(
      result,
      data.frame(x = c(-0.75, 0.25),
                 group = c("control", "treated")),
      type = "response", link = "canonical")
    design <- cbind(1, c(0, 1), c(-0.75, 0.25))
    eta <- as.vector(design %*% result$coefficients)
    expected <- if (identical(family, "binomial")) {
      stats::plogis(eta)
    } else {
      exp(eta)
    }
    expect_equal(prediction$fit, expected, tolerance = 1e-12)
    expect_true(all(prediction$mechanism_lower <= prediction$fit))
    expect_true(all(prediction$mechanism_upper >= prediction$fit))
  }
})

test_that("formal DP GLM release is sticky across retry/restart/rechunk and rotates visibly", {
  fixture <- formal_glm_dp_test_fixture()
  first <- dsVertClient:::.dsvert_formal_glm_dp_result(
    fixture$compilation, fixture$release, fixture$keys$pinset)
  retry <- dsVertClient:::.dsvert_formal_glm_dp_result(
    fixture$compilation, fixture$release, fixture$keys$pinset)
  expect_identical(first$public_release_sha256, retry$public_release_sha256)
  expect_identical(first$coefficients, retry$coefficients)

  rotated <- formal_glm_dp_test_release(
    fixture$compilation, fixture$keys,
    fixture$oracle$coefficients + 1 / 2^16,
    label = "k3/binomial", epoch = "epoch-2")
  second <- dsVertClient:::.dsvert_formal_glm_dp_result(
    fixture$compilation, rotated, fixture$keys$pinset)
  expect_false(identical(first$public_release_sha256,
                         second$public_release_sha256))
  expect_false(identical(first$privacy_epoch_sha256,
                         second$privacy_epoch_sha256))
  expect_identical(second$privacy$history_can_deny_operation, FALSE)
  expect_identical(second$root_rotation,
                   "new_visible_composed_release_never_blocks")
})

test_that("formal DP GLM rejects tamper, forged pins, replay and quota fields", {
  fixture <- formal_glm_dp_test_fixture()
  validate <- function(release = fixture$release,
                       compilation = fixture$compilation,
                       pinset = fixture$keys$pinset) {
    dsVertClient:::.dsvert_formal_glm_dp_result(
      compilation, release, pinset)
  }
  mutations <- list(
    coefficient = function(value) {
      value$shifted_coefficient_lattice_steps[[1L]] <- "999"
      value
    },
    oversized_signature = function(value) {
      value$signatures[[1L]]$signature <- strrep("A", 1000L)
      value
    },
    phase = function(value) {
      value$phase19_binding$protected_path <- "legacy"
      value
    },
    validity = function(value) {
      value$phase19_binding$hidden_execution_validity_consumed <- FALSE
      value
    },
    quota = function(value) {
      value$authority$history_can_deny_operation <- TRUE
      value
    },
    statistic = function(value) {
      value$postprocessing_contract$wald <- "available"
      value
    })
  for (mutation in mutations) {
    expect_error(validate(mutation(fixture$release)),
                 class = "dsvert_formal_glm_release_error")
  }

  forged <- fixture$keys$pinset
  replacement <- formal_glm_dp_test_keys(2L)$pinset[[1L]]
  forged[[1L]] <- replacement
  expect_error(validate(pinset = forged),
               class = "dsvert_formal_glm_release_error")

  other <- formal_glm_dp_test_fixture(family = "poisson")$compilation
  expect_error(validate(compilation = other),
               class = "dsvert_formal_glm_release_error")
})

test_that("formal DP GLM rejects re-signed unsafe server contracts", {
  fixture <- formal_glm_dp_test_fixture()
  validate <- function(release) {
    release <- formal_glm_dp_test_resign(
      release, fixture$compilation, fixture$keys)
    dsVertClient:::.dsvert_formal_glm_dp_result(
      fixture$compilation, release, fixture$keys$pinset)
  }
  mutations <- list(
    outside_box = function(value) {
      value$shifted_coefficient_lattice_steps[[1L]] <- "999999999"
      value
    },
    variance_omits_peer = function(value) {
      value$mechanism$mechanism_variance_upper_steps <-
        formal_glm_dp_test_rat("31")
      value
    },
    numeric_total_too_low = function(value) {
      value$numeric_certificate$total_deterministic_l2_error_upper <-
        formal_glm_dp_test_rat("0")
      value
    },
    backend_substitution = function(value) {
      value$mechanism$backend <- "substituted"
      value
    },
    ring_width_substitution = function(value) {
      value$numeric_certificate$ring_bits <- 129
      value
    },
    hidden_validity = function(value) {
      value$phase19_binding$hidden_execution_validity_consumed <- FALSE
      value
    },
    opening_count = function(value) {
      value$phase19_binding$final_openings <- 2
      value
    },
    all_k = function(value) {
      value$authority$all_k_manifest_admission_verified <- FALSE
      value
    },
    history_gate = function(value) {
      value$authority$history_can_deny_operation <- TRUE
      value
    },
    false_wald = function(value) {
      value$postprocessing_contract$wald <- "available"
      value
    })
  for (name in names(mutations)) {
    expect_error(
      validate(mutations[[name]](fixture$release)),
      class = "dsvert_formal_glm_release_error",
      info = name)
  }
})

test_that("public prediction is bounded, explicit about clipping and offset-safe", {
  fixture <- formal_glm_dp_test_fixture()
  result <- dsVertClient:::.dsvert_formal_glm_dp_result(
    fixture$compilation, fixture$release, fixture$keys$pinset)

  clipped <- dsVertClient:::.dsvert_formal_glm_predict_public(
    result,
    data.frame(x = c(-99, 0.25), group = c("control", "treated")),
    type = "response", link = "surrogate")
  expect_identical(clipped$public_input_clipped, c(TRUE, FALSE))
  expect_true(all(is.finite(clipped$fit)))
  expect_error(
    dsVertClient:::.dsvert_formal_glm_predict_public(
      result, data.frame(x = 0, group = "unknown")),
    "unregistered level")
  expect_error(
    dsVertClient:::.dsvert_formal_glm_predict_public(
      result, data.frame(x = 0, group = "control"), offset = 0),
    "no registered prediction offset")
  expect_error(
    dsVertClient:::.dsvert_formal_glm_pwl_public(
      7, result$schema$link_surrogate),
    "escaped the certified")

  keys <- formal_glm_dp_test_keys(3L)
  model <- formal_glm_dp_test_model()
  model$columns$off <- list(
    kind = "offset", owner = "site_a", lower = "-1", upper = "1")
  model$offset <- list(mode = "bounded_offset", column = "off")
  compilation <- dsVertClient:::.dsvert_formal_glm_compile(
    model, formal_glm_dp_test_authority(3L, keys$pinset))
  data <- fixture$data
  data$off <- c(-0.5, 0, 0.25, 0.5, -0.25, 0.75)
  oracle <- dsVertClient:::.dsvert_formal_glm_rational_oracle(
    compilation, data)
  release <- formal_glm_dp_test_release(
    compilation, keys, oracle$coefficients, label = "bounded-offset")
  offset_result <- dsVertClient:::.dsvert_formal_glm_dp_result(
    compilation, release, keys$pinset)
  newdata <- data.frame(
    x = c(-0.5, 0.25, 0.75),
    group = c("control", "treated", "control"))
  prediction <- dsVertClient:::.dsvert_formal_glm_predict_public(
    offset_result, newdata, type = "link", offset = c(-2, 0.25, 2))
  design <- cbind(1, c(0, 1, 0), newdata$x)
  expected <- as.vector(
    design %*% offset_result$coefficients + c(-1, 0.25, 1))
  expect_equal(prediction$fit, expected, tolerance = 1e-12)
  expect_identical(prediction$public_input_clipped, c(TRUE, FALSE, TRUE))
})

test_that("sampling inference and protected diagnostics fail closed without a joint artifact", {
  fixture <- formal_glm_dp_test_fixture()
  result <- dsVertClient:::.dsvert_formal_glm_dp_result(
    fixture$compilation, fixture$release, fixture$keys$pinset)

  for (statistic in c("sampling_covariance", "standard_errors",
                      "log_likelihood", "deviance", "aic", "wald",
                      "likelihood_ratio", "protected_residuals")) {
    state <- dsVertClient:::.dsvert_formal_glm_statistic(result, statistic)
    expect_s3_class(state, "dsvert_formal_glm_unavailable")
    expect_identical(state$status,
                     "requires_joint_dp_inference_artifact")
    expect_identical(state$additional_server_calls, 0L)
  }
  expect_true(is.na(result$deviance))
  expect_true(is.na(result$log_likelihood))
  expect_true(is.na(result$aic))
  expect_null(result$fitted.values)
  expect_null(result$residuals)

  effects <- dsVertClient:::.dsvert_formal_glm_effects(result)
  expect_identical(effects$scale, "odds_ratio")
  expect_equal(effects$estimate, exp(result$coefficients))
  expect_true(all(effects$mechanism_lower <= effects$estimate))
  expect_true(all(effects$mechanism_upper >= effects$estimate))

  prediction_body <- paste(deparse(body(
    dsVertClient:::.dsvert_formal_glm_predict_public)), collapse = "\n")
  expect_false(grepl("datashield|aggregate|DSI", prediction_body,
                     ignore.case = TRUE))
})

test_that("regularized formal target remains explicit at ordinary degeneracies", {
  separated <- data.frame(
    y = c(0, 0, 0, 1, 1, 1), x = c(-1, -0.8, -0.4, 0.4, 0.8, 1),
    group = rep(c("control", "treated"), 3L))
  binomial <- formal_glm_dp_test_fixture(family = "binomial",
                                         data = separated)
  expect_error(
    dsVertClient:::.dsvert_formal_glm_ordinary_reference(
      binomial$compilation, separated),
    class = "non_identifiable")
  binomial_result <- dsVertClient:::.dsvert_formal_glm_dp_result(
    binomial$compilation, binomial$release, binomial$keys$pinset)
  expect_identical(binomial_result$estimability,
                   "unique_regularized_box_constrained_target")

  all_zero <- data.frame(
    y = rep(0, 6L), x = seq(-1, 1, length.out = 6L),
    group = rep(c("control", "treated"), 3L))
  poisson <- formal_glm_dp_test_fixture(family = "poisson", data = all_zero)
  expect_error(
    dsVertClient:::.dsvert_formal_glm_ordinary_reference(
      poisson$compilation, all_zero),
    class = "non_identifiable")
  poisson_result <- dsVertClient:::.dsvert_formal_glm_dp_result(
    poisson$compilation, poisson$release, poisson$keys$pinset)
  expect_identical(poisson_result$estimability,
                   "unique_regularized_box_constrained_target")
})

test_that("future joint inference contract is server-derived and scientifically typed", {
  for (family in c("binomial", "poisson")) {
    fixture <- formal_glm_dp_test_fixture(family = family)
    contract <- dsVertClient:::.dsvert_formal_glm_joint_inference_requirements(
      fixture$compilation)
    expect_identical(contract$family, family)
    expect_identical(contract$coordinate_count, 18)
    expect_identical(
      contract$blocks$information_upper_triangle$coordinates, 6)
    expect_identical(
      contract$blocks$score_meat_upper_triangle$coordinates, 6)
    expect_identical(contract$blocks$score_vector$coordinates, 3)
    expect_identical(
      contract$derivation_authority$analyst_sensitivity_fields, FALSE)
    expect_identical(contract$protected_execution$intermediate_openings, 0)
    expect_identical(contract$operation_limit, FALSE)
    expect_identical(contract$request_limit, FALSE)
    expect_identical(contract$history_can_deny_operation, FALSE)
    expect_match(contract$scientific_labels$loss,
                 "not_ordinary_logLik", fixed = TRUE)
    expect_match(contract$scientific_labels$penalized_curvature_inverse,
                 "not_sampling", fixed = TRUE)
    expect_match(contract$scientific_labels$aic,
                 "not_classical_aic", fixed = TRUE)
    expect_identical(contract$production_ready, FALSE)
  }
  expect_identical(
    names(formals(dsVertClient:::
      .dsvert_formal_glm_joint_inference_requirements)),
    "compilation")
})
