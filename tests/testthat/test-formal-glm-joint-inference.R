test_that("formal joint inference matches its central oracle for K=2,3,5", {
  for (family in c("binomial", "poisson")) {
    for (k in c(2L, 3L, 5L)) {
      fixture <- formal_glm_joint_test_fixture(k = k, family = family)
      result <- dsVertClient:::.dsvert_formal_glm_joint_inference_result(
        fixture$point, fixture$release, fixture$keys$pinset)
      scale <- as.numeric(fixture$release$output_lattice_scale)

      expect_s3_class(result, "dsvert_formal_glm_joint_inference")
      expect_identical(result$family, family)
      expect_lte(max(abs(result$information - fixture$oracle$information)),
                 1 / scale)
      expect_lte(max(abs(result$score_meat_raw_dp_release -
                           fixture$oracle$score_meat)), 1 / scale)
      expect_lte(max(abs(result$score - fixture$oracle$score)), 1 / scale)
      expect_equal(
        result$canonical_bounded_log_likelihood_at_dp_beta,
        fixture$oracle$canonical_bounded_log_likelihood_at_dp_beta,
        tolerance = 1 / scale)
      expect_equal(
        result$integrated_pwl_surrogate_loss_at_dp_beta,
        fixture$oracle$integrated_pwl_surrogate_loss_at_dp_beta,
        tolerance = 1 / scale)
      expect_equal(result$admitted_n, fixture$oracle$admitted_n,
                   tolerance = 1 / scale)
      expect_identical(result$penalized_curvature$status,
                       "available_penalized_curvature_inverse")
      expect_s3_class(result$model_based,
                      "dsvert_formal_glm_inference_unavailable")
      expect_identical(result$robust$status,
                       "available_conditional_robust_hc0_covariance")
      expect_equal(result$penalized_curvature$inverse,
                   solve(result$information), tolerance = 1e-10)
      expect_equal(
        result$robust$covariance,
        solve(result$information) %*% result$score_meat %*%
          solve(result$information), tolerance = 1e-10)
      expect_true(all(is.finite(result$robust$standard_errors)))
      expect_false(result$penalized_curvature$sampling_covariance)
      expect_false(result$penalized_curvature$sampling_standard_errors)
      expect_false(result$robust$dp_mechanism_uncertainty_included)
      expect_true(result$mechanism_uncertainty$
                    covariance_not_combined_with_sampling)
      expect_identical(result$protected_path$custodian_count, k)
      expect_identical(result$protected_path$intermediate_openings, 0L)
      expect_identical(result$protected_path$final_openings, 1L)
      expect_true(result$protected_path$hidden_execution_validity_consumed)
      expect_identical(result$privacy$operation_limit, FALSE)
      expect_identical(result$privacy$request_limit, FALSE)
      expect_identical(result$privacy$history_can_deny_operation, FALSE)
      expect_equal(result$privacy$composed_epsilon, 1.5)
      expect_equal(result$privacy$composed_delta, 2e-6)
      expect_identical(
        fixture$release$server_derivation$l2_sensitivity_steps,
        fixture$release$mechanism$l2_sensitivity_steps)
      expect_identical(
        fixture$release$server_derivation$sensitivity_norm_bound,
        "exact_rational_stacked_l1_upper_bounds_joint_l2")
      expect_identical(result$production_ready, FALSE)
    }
  }
})

test_that("formal inference capabilities fail closed without fabricating statistics", {
  capabilities <- formal_glm_joint_test_capabilities(
    j = FALSE, score = FALSE, loglik = FALSE, loss = FALSE,
    admitted_n = FALSE)
  fixture <- formal_glm_joint_test_fixture(capabilities = capabilities)
  result <- dsVertClient:::.dsvert_formal_glm_joint_inference_result(
    fixture$point, fixture$release, fixture$keys$pinset)

  expect_identical(result$penalized_curvature$status,
                   "available_penalized_curvature_inverse")
  expect_s3_class(result$model_based,
                  "dsvert_formal_glm_inference_unavailable")
  expect_equal(
    dsVertClient:::.dsvert_formal_glm_joint_inference_statistic(
      result, "penalized_curvature_inverse"),
    result$penalized_curvature$inverse)
  expect_s3_class(
    dsVertClient:::.dsvert_formal_glm_joint_inference_statistic(
      result, "model_based_standard_errors"),
    "dsvert_formal_glm_inference_unavailable")
  for (value in list(
      result$robust, result$score,
      result$canonical_bounded_log_likelihood_at_dp_beta,
      result$integrated_pwl_surrogate_loss_at_dp_beta,
      result$admitted_n)) {
    expect_s3_class(value, "dsvert_formal_glm_inference_unavailable")
    expect_identical(value$status,
                     "joint_dp_inference_capability_unavailable")
    expect_identical(value$additional_server_calls, 0L)
  }
  for (statistic in c("wald", "likelihood_ratio", "aic", "deviance")) {
    value <- dsVertClient:::.dsvert_formal_glm_joint_inference_statistic(
      result, statistic)
    expect_s3_class(value, "dsvert_formal_glm_inference_unavailable")
    expect_identical(value$fallback_used, FALSE)
  }
  expect_false("logLik" %in% names(result))
  expect_false("p.value" %in% names(result))
})

test_that("formal inference rejects tamper, replay, forged pins and gates", {
  fixture <- formal_glm_joint_test_fixture()
  validate <- function(release = fixture$release,
                       point = fixture$point,
                       pinset = fixture$keys$pinset) {
    dsVertClient:::.dsvert_formal_glm_joint_inference_result(
      point, release, pinset)
  }
  mutations <- list(
    coordinate = function(value) {
      value$shifted_coordinate_lattice_steps[[1L]] <- "1"
      value
    },
    signature = function(value) {
      value$signatures[[1L]]$signature <- strrep("A", 1000L)
      value
    },
    base = function(value) {
      value$base_release_binding$public_release_sha256 <-
        formal_glm_joint_test_hash("different-base")
      value
    },
    hidden_validity = function(value) {
      value$protected_binding$hidden_execution_validity_consumed <- FALSE
      value
    },
    quota = function(value) {
      value$authority$history_can_deny_operation <- TRUE
      value
    })
  for (name in names(mutations)) {
    expect_error(validate(mutations[[name]](fixture$release)),
                 class = "dsvert_formal_glm_release_error", info = name)
  }

  forged <- fixture$keys$pinset
  forged[[1L]] <- formal_glm_joint_test_keys(2L)$pinset[[1L]]
  expect_error(validate(pinset = forged),
               class = "dsvert_formal_glm_release_error")

  changed_point <- fixture$point
  changed_point$coefficients[[1L]] <-
    changed_point$coefficients[[1L]] + 0.25
  expect_error(validate(point = changed_point),
               class = "dsvert_formal_glm_release_error")
  changed_schema <- fixture$point
  changed_schema$schema$estimand$target <- "analyst-substituted-target"
  expect_error(validate(point = changed_schema),
               class = "dsvert_formal_glm_release_error")

  rotated_point_release <- formal_glm_joint_test_point_release(
    fixture$compilation, fixture$keys,
    fixture$point_oracle$coefficients,
    label = "k3/binomial/unit/point", epoch = "rotated-point")
  rotated_point <- dsVertClient:::.dsvert_formal_glm_dp_result(
    fixture$compilation, rotated_point_release, fixture$keys$pinset)
  expect_error(validate(point = rotated_point),
               class = "dsvert_formal_glm_release_error")
})

test_that("re-signed unsafe inferential contracts are still rejected", {
  fixture <- formal_glm_joint_test_fixture()
  validate <- function(release) {
    release <- formal_glm_joint_test_resign(
      release, fixture$compilation, fixture$keys)
    dsVertClient:::.dsvert_formal_glm_joint_inference_result(
      fixture$point, release, fixture$keys$pinset)
  }
  mutations <- list(
    quota = function(value) {
      value$authority$request_limit <- TRUE
      value
    },
    validity = function(value) {
      value$protected_binding$hidden_execution_validity_consumed <- FALSE
      value
    },
    composition = function(value) {
      value$mechanism$composed_epsilon <- formal_glm_joint_test_rat("1")
      value
    },
    peer_epsilon = function(value) {
      value$mechanism$epsilon_divided_by_peer_count <- TRUE
      value
    },
    variance = function(value) {
      value$mechanism$mechanism_variance_upper_steps <-
        formal_glm_joint_test_rat("1")
      value
    },
    ring = function(value) {
      value$numeric_certificate$ring_bits <- 127
      value
    },
    scientific = function(value) {
      value$scientific_contract$wald <- "available"
      value
    },
    bounds = function(value) {
      value$numeric_certificate$coordinate_bounds_sha256 <-
        formal_glm_joint_test_hash("forged-bounds")
      value
    },
    sensitivity = function(value) {
      value$mechanism$l2_sensitivity_steps <- "1"
      value
    },
    server_plan = function(value) {
      value$server_derivation$adjacency_triangle_factor <- 99
      value
    },
    simultaneous_radius = function(value) {
      value$mechanism$simultaneous_95_abs_steps <- "1"
      value
    },
    opening = function(value) {
      value$protected_binding$intermediate_openings <- 1
      value
    })
  for (name in names(mutations)) {
    expect_error(validate(mutations[[name]](fixture$release)),
                 class = "dsvert_formal_glm_release_error", info = name)
  }
})

test_that("joint inference is sticky and root rotation is visible but non-blocking", {
  fixture <- formal_glm_joint_test_fixture()
  first <- dsVertClient:::.dsvert_formal_glm_joint_inference_result(
    fixture$point, fixture$release, fixture$keys$pinset)
  retry <- dsVertClient:::.dsvert_formal_glm_joint_inference_result(
    fixture$point, fixture$release, fixture$keys$pinset)
  expect_identical(first$public_release_sha256,
                   retry$public_release_sha256)
  expect_identical(first$information, retry$information)

  rotated_release <- formal_glm_joint_test_release(
    fixture$point, fixture$compilation, fixture$keys, fixture$oracle,
    fixture$capabilities, label = "k3/binomial/unit/inference",
    epoch = "rotated-inference")
  rotated <- dsVertClient:::.dsvert_formal_glm_joint_inference_result(
    fixture$point, rotated_release, fixture$keys$pinset)
  expect_false(identical(first$public_release_sha256,
                         rotated$public_release_sha256))
  expect_false(identical(first$privacy_epoch_sha256,
                         rotated$privacy_epoch_sha256))
  expect_identical(rotated$root_rotation,
                   "new_visible_composed_release_never_blocks")
  expect_identical(rotated$privacy$history_can_deny_operation, FALSE)
})

test_that("rank, box and weight gates preserve the estimand without fallback", {
  fixture <- formal_glm_joint_test_fixture()
  rank_fixture <- formal_glm_joint_test_fixture(
    ridge = "0.0001", label = "rank-fixture")
  rank_oracle <- rank_fixture$oracle
  rank_oracle$information[,] <- 0
  rank_release <- formal_glm_joint_test_release(
    rank_fixture$point, rank_fixture$compilation, rank_fixture$keys,
    rank_oracle, rank_fixture$capabilities, label = "rank-deficient")
  rank_result <- dsVertClient:::.dsvert_formal_glm_joint_inference_result(
    rank_fixture$point, rank_release, rank_fixture$keys$pinset)
  expect_s3_class(rank_result$penalized_curvature,
                  "dsvert_formal_glm_inference_unavailable")
  expect_s3_class(rank_result$robust,
                  "dsvert_formal_glm_inference_unavailable")
  expect_match(rank_result$penalized_curvature$reason,
               "not stably positive definite")
  expect_null(rank_result$penalized_curvature$inverse)

  boundary_coefficients <- stats::setNames(
    c(0.01, rep(0, length(fixture$point$coefficients) - 1L)),
    names(fixture$point$coefficients))
  boundary_point_release <- formal_glm_joint_test_point_release(
    fixture$compilation, fixture$keys, boundary_coefficients,
    label = "boundary-point")
  boundary_point <- dsVertClient:::.dsvert_formal_glm_dp_result(
    fixture$compilation, boundary_point_release, fixture$keys$pinset)
  boundary_oracle <- dsVertClient:::.dsvert_formal_glm_joint_inference_oracle(
    boundary_point, fixture$compilation, fixture$data)
  boundary_release <- formal_glm_joint_test_release(
    boundary_point, fixture$compilation, fixture$keys, boundary_oracle,
    label = "boundary-inference")
  boundary <- dsVertClient:::.dsvert_formal_glm_joint_inference_result(
    boundary_point, boundary_release, fixture$keys$pinset)
  expect_s3_class(boundary$penalized_curvature,
                  "dsvert_formal_glm_inference_unavailable")
  expect_match(boundary$penalized_curvature$reason, "touches the signed")

  weighted <- formal_glm_joint_test_fixture(weights = "bounded")
  weighted_result <- dsVertClient:::.dsvert_formal_glm_joint_inference_result(
    weighted$point, weighted$release, weighted$keys$pinset)
  expect_s3_class(weighted_result$model_based,
                  "dsvert_formal_glm_inference_unavailable")
  expect_identical(weighted_result$penalized_curvature$status,
                   "available_penalized_curvature_inverse")
  expect_identical(weighted_result$robust$status,
                   "available_conditional_robust_hc0_covariance")
})

test_that("DP-noisy non-PSD meat is projected explicitly before sandwiching", {
  fixture <- formal_glm_joint_test_fixture()
  noisy <- fixture$oracle
  noisy$score_meat <- diag(c(-1, 1, 2))
  dimnames(noisy$score_meat) <- dimnames(fixture$oracle$score_meat)
  release <- formal_glm_joint_test_release(
    fixture$point, fixture$compilation, fixture$keys, noisy,
    fixture$capabilities, label = "non-psd-meat")
  result <- dsVertClient:::.dsvert_formal_glm_joint_inference_result(
    fixture$point, release, fixture$keys$pinset)

  expect_true(result$robust$score_meat_psd_projection$projection_applied)
  expect_true(all(result$robust$score_meat_psd_projection$
                    projected_eigenvalues >= 0))
  expect_true(min(eigen(result$score_meat, symmetric = TRUE,
                        only.values = TRUE)$values) >= -1e-12)
  expect_identical(result$robust$formula,
                   "H_inverse_J_psd_H_inverse")
})

test_that("PWL derivative ties and nested-model restrictions are explicit", {
  fixture <- formal_glm_joint_test_fixture()
  table <- fixture$compilation$unsigned_schema$link_surrogate
  internal_knot <- table$knots[[2L]]
  slope <- dsVertClient:::.dsvert_formal_glm_pwl_slope_at(
    internal_knot, table)
  expect_identical(dsVertClient:::.dsvert_glm_rat_cmp(
    slope, table$slopes[[1L]]), 0L)

  result <- dsVertClient:::.dsvert_formal_glm_joint_inference_result(
    fixture$point, fixture$release, fixture$keys$pinset)
  same_model <- dsVertClient:::.dsvert_formal_glm_joint_nested_compatibility(
    result, result)
  expect_identical(same_model$status,
                   "incompatible_for_likelihood_ratio")
  expect_true("not_strictly_nested_coefficient_space" %in%
                same_model$reasons)

  poisson <- formal_glm_joint_test_fixture(family = "poisson")
  poisson_result <- dsVertClient:::.dsvert_formal_glm_joint_inference_result(
    poisson$point, poisson$release, poisson$keys$pinset)
  incompatible <- dsVertClient:::.dsvert_formal_glm_joint_nested_compatibility(
    result, poisson_result)
  expect_identical(incompatible$status,
                   "incompatible_for_likelihood_ratio")
  expect_true("different_family" %in% incompatible$reasons)
  lr <- dsVertClient:::.dsvert_formal_glm_joint_inference_statistic(
    result, "likelihood_ratio")
  expect_s3_class(lr, "dsvert_formal_glm_inference_unavailable")
  expect_match(lr$reason, "same protected snapshot")
})

test_that("formal joint inference remains internal and has no client DSI path", {
  expect_false(".dsvert_formal_glm_joint_inference_result" %in%
                 getNamespaceExports("dsVertClient"))
  bodies <- paste(c(
    deparse(body(
      dsVertClient:::.dsvert_formal_glm_joint_inference_result)),
    deparse(body(
      dsVertClient:::.dsvert_formal_glm_joint_inference_statistic))),
    collapse = "\n")
  expect_false(grepl("datashield|aggregate|DSI", bodies,
                     ignore.case = TRUE))
  result <- formal_glm_joint_test_fixture()
  expect_false(result$point$production_ready)
})
