test_that("public LASSO post-processors reject unvalidated ds.glm inputs", {
  legacy <- structure(
    list(family = "gaussian", lambda = 0, covariance = diag(2),
         n_obs = 10, coefficients = c("(Intercept)" = 0, x = 0)),
    class = c("ds.glm", "list"))
  expect_error(ds.vertLASSOProximal(legacy, lambda = 0.05),
               "validated ds.vertDPGaussian")
  expect_error(ds.vertLASSOCV(legacy, lambda_grid = c(0.1, 0)),
               "validated ds.vertDPGaussian")
})

test_that("DP LASSO numerical core matches a normalized analytic oracle", {
  terms <- c("(Intercept)", "x")
  # x=(0,.5,1), y=.25+.5*x. The signed-bound normalized objective has
  # E[x_c^2]=1/6 and E[x_c*y_c]=1/12. At lambda=.05 the exact slope is .2
  # and the unpenalized intercept is .4.
  gram <- matrix(c(1, 0.5, 0.5, 5 / 12), 2L, 2L,
                 dimnames = list(terms, terms))
  cross <- c("(Intercept)" = 0.5, x = 1 / 3)
  result <- .dsvert_lasso_dp_solver(
    gram, cross, lambda = 0.05, max_iter = 2000L, tol = 1e-12,
    keep_intercept = TRUE)

  expect_equal(result$coefficients,
               c("(Intercept)" = 0.4, x = 0.2), tolerance = 1e-9)
  expect_true(result$converged)
  expect_true(result$kkt$satisfied)
  expect_lte(result$kkt$max_violation, result$kkt$tolerance)
  expect_true(result$curvature$strong_convexity)
  expect_true(result$curvature$solution_unique)
  expect_false(result$implicit_ridge)
  expect_false(result$client_psd_projection_applied)
})

.lasso_active_set_oracle <- function(gram, cross, lambda,
                                     keep_intercept) {
  terms <- rownames(gram)
  intercept <- match("(Intercept)", terms, nomatch = 0L)
  unpenalised <- if (keep_intercept && intercept > 0L) {
    intercept
  } else {
    integer()
  }
  penalised <- setdiff(seq_along(terms), unpenalised)
  states <- if (length(penalised)) {
    as.matrix(expand.grid(rep(list(c(-1, 0, 1)), length(penalised))))
  } else {
    matrix(numeric(), nrow = 1L, ncol = 0L)
  }
  candidates <- lapply(seq_len(nrow(states)), function(row) {
    signs <- numeric(length(terms))
    signs[penalised] <- states[row, ]
    active <- sort(c(unpenalised, penalised[signs[penalised] != 0]))
    beta <- numeric(length(terms))
    if (length(active)) {
      beta[active] <- solve(
        gram[active, active, drop = FALSE],
        cross[active] - lambda * signs[active])
    }
    gradient <- as.numeric(gram %*% beta - cross)
    active_penalised <- intersect(active, penalised)
    inactive_penalised <- setdiff(penalised, active_penalised)
    valid <- (!length(unpenalised) ||
                max(abs(gradient[unpenalised])) <= 1e-8) &&
      (!length(active_penalised) || all(
        sign(beta[active_penalised]) == signs[active_penalised])) &&
      (!length(inactive_penalised) || all(
        abs(gradient[inactive_penalised]) <= lambda + 1e-8))
    if (!valid) return(NULL)
    objective <- 0.5 * sum(beta * as.numeric(gram %*% beta)) -
      sum(cross * beta) + lambda * sum(abs(beta[penalised]))
    list(beta = stats::setNames(beta, terms), objective = objective)
  })
  candidates <- Filter(Negate(is.null), candidates)
  stopifnot(length(candidates) > 0L)
  candidates[[which.min(vapply(
    candidates, `[[`, numeric(1L), "objective"))]]
}

test_that("DP LASSO agrees with an exhaustive active-set oracle", {
  withr::local_seed(20260801)
  terms <- c("(Intercept)", "x1", "x2", "x3")
  for (keep_intercept in c(TRUE, FALSE)) {
    for (case in seq_len(30L)) {
      design <- matrix(stats::rnorm(24L), nrow = 6L, ncol = 4L)
      gram <- crossprod(design) / nrow(design) + diag(0.25, 4L)
      dimnames(gram) <- list(terms, terms)
      cross <- stats::setNames(stats::rnorm(4L), terms)
      lambda <- stats::runif(1L, min = 0.01, max = 0.8)
      oracle <- .lasso_active_set_oracle(
        gram, cross, lambda, keep_intercept)
      result <- .dsvert_lasso_dp_solver(
        gram, cross, lambda = lambda, max_iter = 10000L,
        tol = 1e-12, keep_intercept = keep_intercept)

      expect_equal(
        result$coefficients, oracle$beta, tolerance = 2e-8,
        info = paste("case", case, "keep_intercept", keep_intercept))
      expect_equal(
        result$objective_without_constant, oracle$objective,
        tolerance = 2e-8,
        info = paste("case", case, "keep_intercept", keep_intercept))
      expect_true(result$kkt$satisfied)
    }
  }
})

test_that("DP LASSO accepts certified PSD curvature without inventing ridge", {
  terms <- c("x1", "x2")
  gram <- matrix(1, 2L, 2L, dimnames = list(terms, terms))
  cross <- c(x1 = 1, x2 = 1)
  result <- .dsvert_lasso_dp_solver(
    gram, cross, lambda = 0.1, max_iter = 2000L, tol = 1e-12,
    keep_intercept = FALSE)

  expect_true(result$kkt$satisfied)
  expect_false(result$curvature$strong_convexity)
  expect_true(is.na(result$curvature$solution_unique))
  expect_match(result$curvature$warning, "neither uniqueness", fixed = TRUE)
  expect_false(result$implicit_ridge)
  expect_equal(sum(result$coefficients), 0.9, tolerance = 1e-9)

  non_psd <- gram
  non_psd[1, 2] <- non_psd[2, 1] <- 2
  condition <- tryCatch(
    .dsvert_lasso_dp_solver(
      non_psd, cross, lambda = 0.1, max_iter = 100L, tol = 1e-10,
      keep_intercept = FALSE),
    error = identity)
  expect_s3_class(condition, "non_identifiable")
  expect_identical(condition$reason, "non_psd_dp_lasso_information")
})

test_that("DP LASSO types an incompatible zero-curvature score as non-identifiable", {
  terms <- c("x_zero", "x_regular")
  gram <- diag(c(0, 1))
  dimnames(gram) <- list(terms, terms)
  cross <- c(x_zero = 0.2, x_regular = 0)

  condition <- tryCatch(
    .dsvert_lasso_dp_solver(
      gram, cross, lambda = 0.1, max_iter = 10L, tol = 1e-12,
      keep_intercept = FALSE),
    error = identity)

  expect_s3_class(condition, "non_identifiable")
  expect_identical(condition$reason, "zero_curvature_dp_lasso_score")
  expect_match(condition$message, "x_zero", fixed = TRUE)

  boundary <- .dsvert_lasso_dp_solver(
    gram, c(x_zero = 0.1, x_regular = 0), lambda = 0.1,
    max_iter = 10L, tol = 1e-12, keep_intercept = FALSE)
  expect_true(boundary$kkt$satisfied)
  expect_equal(boundary$coefficients[["x_zero"]], 0)
})

test_that("DP LASSO uses the separate real DP count without Gram averaging", {
  terms <- c("(Intercept)", "x")
  gram_sum <- matrix(c(3, 1.5, 1.5, 1.25), 2L, 2L,
                     dimnames = list(terms, terms))
  cross_sum <- c("(Intercept)" = 1.5, x = 1)
  moments <- .dsvert_lasso_dp_moments(
    gram_sum, cross_sum, outcome_square_projected = 0.875,
    n_obs = 3.25)

  expect_equal(moments$n_obs, 3.25)
  expect_equal(moments$gram, gram_sum / 3.25)
  expect_equal(moments$cross, cross_sum / 3.25)
  expect_equal(moments$outcome_square, 0.875 / 3.25)
  expect_false(isTRUE(all.equal(moments$n_obs, gram_sum[1, 1])))
  expect_match(moments$normalization, "without_count_Gram11_averaging",
               fixed = TRUE)
  expect_true(moments$augmented_certificate$positive_semidefinite)
  expect_false(moments$augmented_certificate$client_reprojection_applied)

  incompatible_cross <- cross_sum
  incompatible_cross[[2L]] <- 10
  condition <- tryCatch(
    .dsvert_lasso_dp_moments(
      gram_sum, incompatible_cross, outcome_square_projected = 0.875,
      n_obs = 3.25),
    error = identity)
  expect_s3_class(condition, "non_identifiable")
  expect_identical(
    condition$reason, "non_psd_dp_lasso_augmented_moments")
})

make_dp_gaussian_lasso_fit <- function(integrity = TRUE) {
  terms <- c("(Intercept)", "x")
  gram <- matrix(c(3, 1.5, 1.5, 1.25), 2L, 2L,
                 dimnames = list(terms, terms))
  fit <- list(
    status = "ok", family = "gaussian", n_obs = 3,
    ridge = 0,
    source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    signed_artifact = list(
      predictor_order = "x",
      predictors = list(x = list(lower = -2, upper = 2)),
      outcome = list(lower = 10, upper = 20),
      intercept = TRUE,
      implementation_state = "same_owner_materialized",
      cross_owner_state = "reserved_not_materialized"),
    sufficient_statistics_dp = list(
      gram_projected = gram,
      cross_projected = c("(Intercept)" = 1.5, x = 1),
      outcome_square_projected = 1),
    query_contract_sha256 = paste(rep("a", 64L), collapse = ""),
    release_contract_hash = paste(rep("b", 64L), collapse = ""),
    accuracy = list(confidence = 0.95),
    provenance_certificate = list(test_integrity = integrity))
  class(fit) <- c("ds.vertDPGaussian", "list")
  fit
}

test_that("historical iterative LASSO admits only the signed Gaussian path", {
  fit <- make_dp_gaussian_lasso_fit()
  glm_call <- NULL
  testthat::local_mocked_bindings(
    ds.validateDPGaussianCertificate = function(...) {
      list(
        integrity_valid = TRUE,
        authenticity = "session_transport_anchored",
        anchor_source = "online_session")
    },
    ds.vertGLM = function(formula, data = NULL, family = "gaussian",
                          dp_analysis_id = NULL, missing = NULL,
                          verbose = TRUE, datasources = NULL, ...) {
      glm_call <<- list(
        formula = formula, data = data, family = family,
        dp_analysis_id = dp_analysis_id, missing = missing,
        verbose = verbose, datasources = datasources)
      fit
    },
    .dsvert_quarantine_test_mode = function() FALSE,
    .package = "dsVertClient")

  expect_error(
    ds.vertLASSOIter(
      y ~ x, data = "cohort", family = "gaussian",
      datasources = list(site_a = structure(list(), class = "mock_connection"))),
    "dp_analysis_id is required")
  expect_null(glm_call)

  result <- ds.vertLASSOIter(
    y ~ x, data = "cohort", family = "gaussian",
    lambda = c(0.1, 0.05, 0), max_outer = 20L, tol = 1e-12,
    dp_analysis_id = "gaussian-primary", verbose = FALSE,
    datasources = list(site_a = structure(list(), class = "mock_connection")))

  expect_s3_class(result, "ds.vertLASSOIter")
  expect_s3_class(result, "ds.vertDPLASSOPath")
  expect_identical(result$family, "gaussian")
  expect_identical(result$input_provenance, "signed_dp_gaussian_capsule")
  expect_identical(result$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
  expect_true(all(vapply(result$path_certificates, function(value) {
    isTRUE(value$kkt$satisfied)
  }, logical(1L))))
  expect_identical(glm_call$data, "cohort")
  expect_identical(glm_call$family, "gaussian")
  expect_identical(glm_call$dp_analysis_id, "gaussian-primary")
  expect_identical(glm_call$missing, "complete_case_capsule")
  expect_false(glm_call$verbose)

  alias_result <- ds.vert.lasso_iter(
    y ~ x, data = "cohort", family = "gaussian", lambda = 0.1,
    dp_analysis_id = "gaussian-primary", verbose = FALSE,
    datasources = list(site_a = structure(list(), class = "mock_connection")))
  expect_s3_class(alias_result, "ds.vertLASSOIter")
  expect_identical(alias_result$frontdoor, "ds.vert.lasso_iter")
  expect_identical(alias_result$route,
                   "ds.vertLASSOIter(signed-gaussian-synopsis)")

  expect_error(
    ds.vertLASSOIter(
      y ~ x, data = "cohort", family = "binomial",
      dp_analysis_id = "gaussian-primary",
      datasources = list(site_a = structure(list(), class = "mock_connection"))),
    class = "dsvert_route_unavailable")
})

test_that("public proximal LASSO post-processes one certified DP release", {
  fit <- make_dp_gaussian_lasso_fit()
  testthat::local_mocked_bindings(
    ds.validateDPGaussianCertificate = function(value, ...) {
      list(
        integrity_valid = isTRUE(
          value$provenance_certificate$test_integrity),
        authenticity = "session_transport_anchored",
        anchor_source = "online_session")
    },
    .package = "dsVertClient")

  result <- ds.vertLASSOProximal(
    fit, lambda = 0.05, max_iter = 2000L, tol = 1e-12)

  expect_s3_class(result, "ds.vertDPLASSO")
  expect_equal(
    result$coefficients_normalized,
    c("(Intercept)" = 0.4, x = 0.2), tolerance = 1e-9)
  expect_equal(
    result$coefficients,
    c("(Intercept)" = 15, x = 0.5), tolerance = 1e-9)
  expect_true(result$kkt$satisfied)
  expect_identical(result$input_provenance,
                   "signed_dp_gaussian_capsule")
  expect_identical(result$additional_server_calls_after_capsule, 0L)
  expect_equal(result$additional_privacy_cost,
               c(epsilon = 0, delta = 0))
  expect_false(result$source_values_exposed)
  expect_false(result$intermediate_values_exposed)
  expect_false(result$coefficient_regions_available)
  printed <- capture.output(print(result))
  expect_true(any(grepl("extra DP cost = (0, 0)", printed, fixed = TRUE)))
})

test_that("historical LASSO frontdoor returns a certified DP path", {
  fit <- make_dp_gaussian_lasso_fit()
  testthat::local_mocked_bindings(
    ds.validateDPGaussianCertificate = function(value, ...) {
      list(
        integrity_valid = isTRUE(
          value$provenance_certificate$test_integrity),
        authenticity = "session_transport_anchored",
        anchor_source = "online_session")
    },
    .package = "dsVertClient")

  result <- ds.vertLASSO(
    fit, lambda_1 = 0.1, alpha_grid = c(1, 0.5, 0),
    max_iter = 2000L, tol = 1e-12)

  expect_s3_class(result, "ds.vertLASSO")
  expect_s3_class(result, "ds.vertDPLASSOPath")
  expect_equal(result$lambda_grid, c(0.1, 0.05, 0))
  expect_equal(result$original, result$paths[[3L]], tolerance = 1e-12)
  expect_true(all(vapply(result$path_certificates, function(value) {
    isTRUE(value$kkt$satisfied)
  }, logical(1L))))
  expect_identical(result$input_provenance, "signed_dp_gaussian_capsule")
  expect_identical(result$additional_server_calls_after_capsule, 0L)
  expect_equal(result$additional_privacy_cost,
               c(epsilon = 0, delta = 0))
  expect_false(result$source_values_exposed)
  expect_false(result$intermediate_values_exposed)
  expect_error(
    ds.vertLASSO(
      structure(list(coefficients = c(`(Intercept)` = 0, x = 0)),
                class = c("ds.glm", "list")),
      lambda_1 = 0.1),
    "validated ds.vertDPGaussian")
  fit$signed_artifact$cross_owner_state <- "exact_gc_to_joint_dp_vector_v1"
  expect_error(ds.vertLASSO(fit, lambda_1 = 0.1), "same-owner")
})

test_that("lowercase LASSO alias retains only the DP path", {
  fit <- make_dp_gaussian_lasso_fit()
  testthat::local_mocked_bindings(
    ds.validateDPGaussianCertificate = function(...) {
      list(integrity_valid = TRUE,
           authenticity = "session_transport_anchored")
    },
    .package = "dsVertClient")
  result <- ds.vert.lasso(fit, lambda_1 = 0.1, alpha_grid = 1)
  expect_s3_class(result, "ds.vertDPLASSOPath")
  expect_identical(result$frontdoor, "ds.vert.lasso")
  expect_identical(result$route, "ds.vertLASSO")
})

test_that("LASSO1Step is a certified DP path compatibility alias", {
  fit <- make_dp_gaussian_lasso_fit()
  testthat::local_mocked_bindings(
    ds.validateDPGaussianCertificate = function(...) {
      list(integrity_valid = TRUE,
           authenticity = "session_transport_anchored")
    },
    .package = "dsVertClient")
  result <- ds.vertLASSO1Step(
    fit, lambda = c(0.1, 0.05, 0), max_iter = 2000L, tol = 1e-12)
  expect_s3_class(result, "ds.vertLASSO1Step")
  expect_s3_class(result, "ds.vertDPLASSOPath")
  expect_equal(result$lambda, c(0.1, 0.05, 0))
  expect_equal(result$original, result$paths[[3L]], tolerance = 1e-12)
  expect_true(all(vapply(result$path_certificates, function(value) {
    isTRUE(value$kkt$satisfied)
  }, logical(1L))))
  expect_identical(result$additional_server_calls_after_capsule, 0L)

  alias_result <- ds.vert.lasso_1step(
    fit, lambda = c(0.1, 0.05, 0), max_iter = 2000L, tol = 1e-12)
  expect_s3_class(alias_result, "ds.vertLASSO1Step")
  expect_identical(alias_result$frontdoor, "ds.vert.lasso_1step")
  expect_identical(alias_result$route, "ds.vertLASSO1Step")
  expect_equal(alias_result$paths, result$paths, tolerance = 1e-12)
})

test_that("DP LASSO fails closed on certificate or design mismatch", {
  fit <- make_dp_gaussian_lasso_fit(integrity = FALSE)
  testthat::local_mocked_bindings(
    ds.validateDPGaussianCertificate = function(value, ...) {
      list(integrity_valid = isTRUE(
        value$provenance_certificate$test_integrity),
        authenticity = "session_transport_anchored")
    },
    .package = "dsVertClient")
  expect_error(
    ds.vertLASSOProximal(fit, lambda = 0.05),
    "failed integrity validation")

  fit$provenance_certificate$test_integrity <- TRUE
  colnames(fit$sufficient_statistics_dp$gram_projected) <- c("x", "z")
  expect_error(
    ds.vertLASSOProximal(fit, lambda = 0.05),
    "malformed|signed design terms")
})

test_that("DP LASSO warm starts are original-scale and never contact DSI", {
  fit <- make_dp_gaussian_lasso_fit()
  testthat::local_mocked_bindings(
    ds.validateDPGaussianCertificate = function(...) {
      list(integrity_valid = TRUE,
           authenticity = "session_transport_anchored")
    },
    .dsvert_aggregate_strict = function(...) {
      stop("unexpected DSI call", call. = FALSE)
    },
    .dsvert_fanout_by_site = function(...) {
      stop("unexpected DSI call", call. = FALSE)
    },
    .package = "dsVertClient")

  result <- ds.vertLASSOProximal(
    fit, lambda = 0.05,
    warm_start = c("(Intercept)" = 15, x = 0.5),
    max_iter = 2000L, tol = 1e-12)
  expect_true(result$kkt$satisfied)
  expect_equal(result$coefficients,
               c("(Intercept)" = 15, x = 0.5), tolerance = 1e-9)
  expect_error(
    ds.vertLASSOProximal(
      fit, lambda = 0.05,
      warm_start = c("(Intercept)" = 15, wrong = 0.5)),
    "must match the original-scale coefficient names")
})

test_that("public LASSOCV is honest deterministic DP pseudo-IC selection", {
  fit <- make_dp_gaussian_lasso_fit()
  testthat::local_mocked_bindings(
    ds.validateDPGaussianCertificate = function(...) {
      list(integrity_valid = TRUE,
           authenticity = "session_transport_anchored")
    },
    .dsvert_aggregate_strict = function(...) {
      stop("unexpected DSI call", call. = FALSE)
    },
    .dsvert_fanout_by_site = function(...) {
      stop("unexpected DSI call", call. = FALSE)
    },
    .package = "dsVertClient")

  result <- ds.vertLASSOCV(
    fit, lambda_grid = c(0.1, 0.05, 0), criterion = "BIC",
    se_threshold = 0.02)

  expect_s3_class(result, "ds.vertDPLASSOSelect")
  expect_true(result$selection_available)
  expect_identical(result$selection_method,
                   "DP_projected_pseudo_information_criterion")
  expect_false(result$classical_information_criterion)
  expect_false(result$cross_validation)
  expect_false(result$one_standard_error_rule)
  expect_length(result$paths, 3L)
  expect_true(all(vapply(
    result$path_certificates,
    function(value) isTRUE(value$kkt$satisfied), logical(1L))))
  expect_identical(result$additional_server_calls_after_capsule, 0L)
  expect_equal(result$additional_privacy_cost,
               c(epsilon = 0, delta = 0))
  expect_false(result$source_values_exposed)
  printed <- capture.output(print(result))
  expect_true(any(grepl("DP-projected pseudo-IC", printed, fixed = TRUE)))
  expect_true(any(grepl("extra DP cost = (0, 0)", printed, fixed = TRUE)))
})

test_that("DP pseudo-IC returns structured unavailability", {
  fit <- make_dp_gaussian_lasso_fit()
  fit$n_obs <- 1
  testthat::local_mocked_bindings(
    ds.validateDPGaussianCertificate = function(...) {
      list(integrity_valid = TRUE,
           authenticity = "session_transport_anchored")
    },
    .package = "dsVertClient")

  result <- ds.vertLASSOCV(
    fit, lambda_grid = c(0.1, 0), criterion = "AIC")
  expect_false(result$selection_available)
  expect_identical(result$selection_unavailable_reason,
                   "dp_effective_count_not_greater_than_one")
  expect_null(result$lambda.min)
  expect_null(result$beta.min)
  expect_true(any(grepl(
    "Selection unavailable", capture.output(print(result)), fixed = TRUE)))
})

test_that("DP LASSO requires transport- or caller-anchored provenance", {
  fit <- make_dp_gaussian_lasso_fit()
  trusted <- list(site_a = paste(rep("A", 43L), collapse = ""))
  testthat::local_mocked_bindings(
    ds.validateDPGaussianCertificate = function(value,
                                                  trusted_pinset = NULL) {
      list(
        integrity_valid = TRUE,
        authenticity = if (is.null(trusted_pinset)) {
          "unanchored"
        } else {
          "caller_anchored"
        })
    },
    .package = "dsVertClient")

  expect_error(
    ds.vertLASSOProximal(fit, lambda = 0.05),
    "not anchored")
  expect_error(
    ds.vertLASSOCV(fit, lambda_grid = c(0.1, 0)),
    "not anchored")
  expect_s3_class(
    ds.vertLASSOProximal(
      fit, lambda = 0.05, trusted_pinset = trusted),
    "ds.vertDPLASSO")
  expect_s3_class(
    ds.vertLASSOCV(
      fit, lambda_grid = c(0.1, 0), trusted_pinset = trusted),
    "ds.vertDPLASSOSelect")
})

test_that("DP LASSO coefficient scale transforms round-trip signed bounds", {
  artifact <- list(
    predictor_order = "x",
    predictors = list(x = list(lower = -2, upper = 2)),
    outcome = list(lower = 10, upper = 20),
    intercept = TRUE)
  theta <- c("(Intercept)" = 0.4, x = 0.2)
  original <- .dsvert_lasso_dp_original_scale(theta, artifact)

  expect_equal(original, c("(Intercept)" = 15, x = 0.5))
  expect_equal(.dsvert_lasso_dp_normalized_scale(original, artifact), theta)
  scales <- .dsvert_lasso_dp_artifact_scales(artifact)
  expect_equal(scales$original_penalty_weights, c(x = 0.4))

  artifact$intercept <- FALSE
  theta_no_intercept <- c(x = 0.2)
  original_no_intercept <- .dsvert_lasso_dp_original_scale(
    theta_no_intercept, artifact)
  expect_equal(original_no_intercept,
               c("(Intercept)" = 11, x = 0.5))
  expect_equal(
    .dsvert_lasso_dp_normalized_scale(original_no_intercept, artifact),
    theta_no_intercept)
  inconsistent <- original_no_intercept
  inconsistent[["(Intercept)"]] <- 0
  expect_error(
    .dsvert_lasso_dp_normalized_scale(inconsistent, artifact),
    "deterministic bound-transform offset", fixed = TRUE)
})

test_that("DP lambda_max conditions on the unpenalized intercept", {
  terms <- c("(Intercept)", "x")
  gram <- matrix(c(1, 0.5, 0.5, 5 / 12), 2L, 2L,
                 dimnames = list(terms, terms))
  cross <- c("(Intercept)" = 0.5, x = 1 / 3)
  maximum <- .dsvert_lasso_dp_lambda_max(
    gram, cross, keep_intercept = TRUE)

  expect_equal(maximum, 1 / 12, tolerance = 1e-14)
  expect_equal(
    .dsvert_lasso_dp_lambda_max(gram, cross, keep_intercept = FALSE),
    0.5)
  grid <- .dsvert_lasso_dp_default_lambda(
    gram, cross, keep_intercept = TRUE, length_out = 5L)
  expect_equal(grid[[1L]], maximum)
  expect_equal(tail(grid, 1L), maximum / 1000)
  expect_true(all(diff(grid) < 0))

  at_maximum <- .dsvert_lasso_dp_solver(
    gram, cross, lambda = maximum, max_iter = 2000L, tol = 1e-12,
    keep_intercept = TRUE)
  expect_equal(unname(at_maximum$coefficients[["x"]]), 0,
               tolerance = 1e-10)
  expect_true(at_maximum$kkt$satisfied)
})

test_that("DP pseudo-IC is deterministic heuristic selection, never CV", {
  terms <- c("(Intercept)", "x")
  gram_sum <- matrix(c(3, 1.5, 1.5, 1.25), 2L, 2L,
                     dimnames = list(terms, terms))
  cross_sum <- c("(Intercept)" = 1.5, x = 1)
  moments <- .dsvert_lasso_dp_moments(
    gram_sum, cross_sum, outcome_square_projected = 1,
    n_obs = 3)
  lambda <- c(0.1, 0.05, 0.01)
  solutions <- .dsvert_lasso_dp_lambda_path(
    moments$gram, moments$cross, lambda,
    max_iter = 5000L, tol = 1e-12, keep_intercept = TRUE)
  selected <- .dsvert_lasso_dp_pseudo_ic(
    moments, solutions, lambda, criterion = "BIC",
    parsimonious_delta = 0.25)

  expect_true(selected$available)
  expect_identical(selected$criterion, "DP_projected_pseudo_BIC")
  expect_false(selected$classical_information_criterion)
  expect_false(selected$cross_validation)
  expect_false(selected$one_standard_error_rule)
  expect_equal(selected$parsimonious_delta, 0.25)
  expect_equal(selected$lambda.parsimonious, selected$lambda.1se)
  expect_equal(
    selected$ic,
    3 * log(selected$rss_mean_dp_projected) + log(3) * selected$df)

  ebic <- .dsvert_lasso_dp_pseudo_ic(
    moments, solutions, lambda, criterion = "EBIC", ebic_gamma = 0.5)
  expected_extra <- vapply(ebic$active_penalised, function(active) {
    2 * 0.5 * lchoose(1, active)
  }, numeric(1L))
  expect_equal(
    ebic$ic,
    3 * log(ebic$rss_mean_dp_projected) + log(3) * ebic$df +
      expected_extra)
})

test_that("DP pseudo-IC reports unavailable counts and residual degeneracy", {
  terms <- c("(Intercept)", "x")
  gram_sum <- matrix(c(1, 0.5, 0.5, 0.5), 2L, 2L,
                     dimnames = list(terms, terms))
  cross_sum <- c("(Intercept)" = 0.5, x = 0.25)
  small_n <- .dsvert_lasso_dp_moments(
    gram_sum, cross_sum, outcome_square_projected = 0.5,
    n_obs = 1)
  solutions <- .dsvert_lasso_dp_lambda_path(
    small_n$gram, small_n$cross, lambda = 0.1,
    max_iter = 2000L, tol = 1e-10, keep_intercept = TRUE)
  unavailable <- .dsvert_lasso_dp_pseudo_ic(
    small_n, solutions, lambda = 0.1)
  expect_false(unavailable$available)
  expect_identical(
    unavailable$reason, "dp_effective_count_not_greater_than_one")
  expect_true(all(is.na(unavailable$ic)))

  perfect_terms <- c("(Intercept)", "x")
  perfect_gram <- matrix(c(3, 1.5, 1.5, 1.25), 2L, 2L,
                         dimnames = list(perfect_terms, perfect_terms))
  perfect_cross <- c("(Intercept)" = 1.5, x = 1)
  perfect <- .dsvert_lasso_dp_moments(
    perfect_gram, perfect_cross, outcome_square_projected = 0.875,
    n_obs = 3)
  perfect_solution <- .dsvert_lasso_dp_lambda_path(
    perfect$gram, perfect$cross, lambda = 0,
    max_iter = 5000L, tol = 1e-12, keep_intercept = TRUE)
  no_score <- .dsvert_lasso_dp_pseudo_ic(
    perfect, perfect_solution, lambda = 0)
  expect_false(no_score$available)
  expect_identical(
    no_score$reason, "non_positive_dp_projected_residual_mean_square")
  expect_null(no_score$lambda.min)
})
