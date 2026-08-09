.numeric_test_capability <- function(
    backend, available = TRUE, verified = TRUE, exact = TRUE,
    ring_bits = NA_integer_, frac_bits = NA_integer_,
    supported_ring_bits = integer(), max_ring_bits = ring_bits,
    max_frac_bits = frac_bits) {
  list(
    available = available,
    allowed = TRUE,
    e2e_verified = verified,
    canonical_encoding = verified,
    fail_closed_overflow = verified,
    runtime_bounds_enforced = verified,
    workload_adapter_e2e_verified = verified,
    public_scalar_mul_truncate_e2e_verified = verified,
    full_iteration_e2e_verified = verified,
    exact_truncation = exact,
    exact_comparison = exact,
    ring_bits = ring_bits,
    frac_bits = frac_bits,
    supported_ring_bits = supported_ring_bits,
    max_ring_bits = max_ring_bits,
    max_frac_bits = max_frac_bits,
    truncation_semantics = if (exact) "floor" else "local_probabilistic")
}

test_that("workload dimensions cannot overflow the intercept count", {
  expect_error(
    dsVertClient:::.dsvert_numeric_glm_workload(
      n_obs = 1L, n_predictors = .Machine$integer.max,
      family = "gaussian", max_iter = 1L),
    "too large to add the intercept")
})

test_that("large valid workload counters do not overflow R integers", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    n_obs = 1L, n_predictors = .Machine$integer.max - 1L,
    family = "gaussian", max_iter = .Machine$integer.max,
    compute_se = TRUE, compute_deviance = TRUE)

  expect_equal(workload$n_parameters, .Machine$integer.max)
  expect_true(is.finite(workload$operation_count))
  expect_true(is.finite(workload$truncation_count))
  expect_true(is.finite(workload$quantized_terms_per_path))
  expect_false(anyNA(unlist(workload[c(
    "n_parameters", "operation_count", "truncation_count",
    "quantized_terms_per_path")], use.names = FALSE)))
})

.numeric_test_policy <- function(id = "a", error_budget = 1e-3,
                                 eta = 2, predictor = 2,
                                 response = 2) {
  policy <- list(
    schema_version = 1L,
    policy_version = paste0("test-only-v1-", id),
    enabled = TRUE,
    workload = "glm",
    bounds = list(
      max_abs_predictor_input = max(1e3, predictor),
      max_abs_predictor = predictor,
      max_abs_response_input = list(
        gaussian = max(1e3, response), binomial = 1,
        poisson = max(1e3, response)),
      max_abs_response = list(
        gaussian = response, binomial = 1, poisson = response),
      max_abs_linear_predictor = list(
        gaussian = eta, binomial = eta, poisson = eta),
      max_abs_approximation_intermediate = list(
        gaussian = eta, binomial = 8, poisson = 1024),
      max_abs_offset = 1,
      max_abs_weight = 10,
      max_observations = 1e7,
      max_predictors = 10000L,
      max_iterations = 10000L,
      max_numeric_error = error_budget),
    approximation = list(
      gaussian = list(domain = NULL, max_abs_error = 0),
      binomial = list(domain = c(-8, 8), max_abs_error = 6.5e-6),
      poisson = list(domain = c(-5, 5), max_abs_error = 5e-10)),
    capabilities = list(
      ring63 = .numeric_test_capability(
        "ring63", ring_bits = 63L, frac_bits = 20L,
        supported_ring_bits = 63L, max_ring_bits = 63L,
        max_frac_bits = 20L),
      ring127 = .numeric_test_capability(
        "ring127", ring_bits = 127L, frac_bits = 50L,
        supported_ring_bits = 127L, max_ring_bits = 127L,
        max_frac_bits = 50L),
      exact_gc = .numeric_test_capability(
        "exact_gc", available = FALSE, verified = FALSE, exact = FALSE,
        supported_ring_bits = c(63L, 127L), max_ring_bits = 127L,
        max_frac_bits = 50L),
      multiprecision = .numeric_test_capability(
        "multiprecision", available = FALSE, verified = FALSE,
        exact = FALSE, max_ring_bits = 512L, max_frac_bits = 200L)))
  policy$policy_id <- dsVertClient:::.dsvert_numeric_policy_id(policy)
  policy
}

.numeric_test_rehash <- function(policy) {
  policy$policy_id <- dsVertClient:::.dsvert_numeric_policy_id(policy)
  policy
}

test_that("preflight selects the smallest jointly certified fast ring", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    n_obs = 100L, n_predictors = 2L, family = "gaussian", max_iter = 5L,
    compute_se = FALSE, compute_deviance = FALSE)
  policies <- list(site_a = .numeric_test_policy("a"),
                   site_b = .numeric_test_policy("b"))

  certificate <- dsVertClient:::.dsvert_numeric_preflight_from_policies(
    policies, workload, requested_backend = "auto", requested_ring = 63L)

  expect_s3_class(certificate, "dsvert_numeric_certificate")
  expect_identical(certificate$status, "preflight_eligible")
  expect_identical(certificate$effective_backend, "ring63")
  expect_identical(certificate$ring_bits, 63L)
  expect_identical(certificate$frac_bits, 20L)
  expect_lte(certificate$required_ring_bits, certificate$ring_bits)
  expect_lte(certificate$planned_total_numeric_error,
             certificate$numeric_error_budget)
  expect_true(certificate$bounds_custodian_owned)
  expect_true(certificate$capability_claims_validated)
  expect_true(certificate$workload_adapter_e2e_verified)
  expect_true(certificate$public_scalar_mul_truncate_e2e_verified)
  expect_true(certificate$full_iteration_e2e_verified)
  expect_false(certificate$capabilities_attested)
  expect_false(certificate$backend_e2e_verified)
  expect_false(certificate$numeric_error_bound_certified)
  expect_false(certificate$numerically_certified)
  expect_true(is.na(certificate$total_numeric_error_max))
  expect_true(certificate$estimator_unchanged)
  expect_identical(certificate$error_bound_scope,
                   "one_element_nonlinear_arithmetic_path")
  expect_true(is.na(certificate$aggregate_output_error_bound))
  expect_true(is.na(certificate$estimator_error_bound))
})

test_that("preflight promotes precision but never changes the estimator", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 2L, "gaussian", 5L, FALSE, FALSE)
  policies <- list(
    a = .numeric_test_policy("a", error_budget = 1e-12),
    b = .numeric_test_policy("b", error_budget = 1e-12))

  certificate <- dsVertClient:::.dsvert_numeric_preflight_from_policies(
    policies, workload, requested_backend = "auto", requested_ring = 63L)

  expect_identical(certificate$effective_backend, "ring127")
  expect_identical(certificate$requested_ring, 63L)
  expect_true(certificate$estimator_unchanged)
})

test_that("non-Gaussian GLM uses public domain bounds, not a phantom comparison", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 2L, "binomial", 5L, FALSE, FALSE)
  policies <- list(a = .numeric_test_policy("a", eta = 6),
                   b = .numeric_test_policy("b", eta = 6))

  certificate <- dsVertClient:::.dsvert_numeric_preflight_from_policies(
    policies, workload, requested_backend = "auto", requested_ring = 63L)

  expect_identical(certificate$comparison_count, 0)
  expect_identical(
    certificate$approximation_domain_guard,
    "custodian_input_bounds_plus_public_coefficient_bound")
  expect_identical(certificate$effective_backend, "ring127")

  policies$a$capabilities$ring127$exact_comparison <- FALSE
  policies$b$capabilities$ring127$exact_comparison <- FALSE
  policies <- lapply(policies, .numeric_test_rehash)
  certificate_without_comparison <-
    dsVertClient:::.dsvert_numeric_preflight_from_policies(
      policies, workload, requested_backend = "auto", requested_ring = 63L)
  expect_identical(certificate_without_comparison$effective_backend,
                   "ring127")
})

test_that("exact GC is selected only from simulated E2E-attested policies", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 2L, "gaussian", 5L, FALSE, FALSE)
  policies <- list(a = .numeric_test_policy("a"),
                   b = .numeric_test_policy("b"))
  for (site in names(policies)) {
    policies[[site]]$capabilities$ring63$exact_truncation <- FALSE
    policies[[site]]$capabilities$ring127$exact_truncation <- FALSE
    policies[[site]]$capabilities$exact_gc <- .numeric_test_capability(
      "exact_gc", available = TRUE, verified = TRUE, exact = TRUE,
      supported_ring_bits = c(63L, 127L), max_ring_bits = 127L,
      max_frac_bits = 50L)
    policies[[site]] <- .numeric_test_rehash(policies[[site]])
  }

  certificate <- dsVertClient:::.dsvert_numeric_preflight_from_policies(
    policies, workload, requested_backend = "auto", requested_ring = 63L)

  expect_identical(certificate$effective_backend, "exact_gc")
  expect_identical(certificate$truncation_semantics, "floor")
  expect_identical(certificate$ring_bits, 63L)
})

test_that("exact GC increases fractional precision before widening the ring", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 2L, "gaussian", 5L, FALSE, FALSE)
  policies <- list(
    a = .numeric_test_policy("a", error_budget = 1e-12),
    b = .numeric_test_policy("b", error_budget = 1e-12))
  for (site in names(policies)) {
    policies[[site]]$capabilities$ring63$exact_truncation <- FALSE
    policies[[site]]$capabilities$ring127$exact_truncation <- FALSE
    policies[[site]]$capabilities$ring127$allowed <- FALSE
    policies[[site]]$capabilities$exact_gc <- .numeric_test_capability(
      "exact_gc", available = TRUE, verified = TRUE, exact = TRUE,
      supported_ring_bits = c(63L, 127L), max_ring_bits = 127L,
      max_frac_bits = 50L)
    policies[[site]] <- .numeric_test_rehash(policies[[site]])
  }

  certificate <- dsVertClient:::.dsvert_numeric_preflight_from_policies(
    policies, workload, requested_backend = "auto", requested_ring = 63L)

  expect_identical(certificate$effective_backend, "exact_gc")
  expect_identical(certificate$ring_bits, 63L)
  expect_gt(certificate$frac_bits, 20L)
  expect_lte(certificate$planned_total_numeric_error,
             certificate$numeric_error_budget)
})

test_that("a separately attested multiprecision backend can widen above Ring127", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 2L, "gaussian", 5L, FALSE, FALSE)
  policies <- list(
    a = .numeric_test_policy(
      "a", predictor = 1e20, response = 1e20),
    b = .numeric_test_policy(
      "b", predictor = 1e20, response = 1e20))
  for (site in names(policies)) {
    policies[[site]]$capabilities$multiprecision <-
      .numeric_test_capability(
        "multiprecision", available = TRUE, verified = TRUE, exact = TRUE,
        max_ring_bits = 512L, max_frac_bits = 200L)
    policies[[site]] <- .numeric_test_rehash(policies[[site]])
  }

  certificate <- dsVertClient:::.dsvert_numeric_preflight_from_policies(
    policies, workload, requested_backend = "auto", requested_ring = 63L)

  expect_identical(certificate$effective_backend, "multiprecision")
  expect_gt(certificate$ring_bits, 127L)
  expect_lte(certificate$required_ring_bits, certificate$ring_bits)
})

test_that("dynamic exact GC itself plans signed rings through Ring4096", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 2L, "gaussian", 5L, FALSE, FALSE)
  policies <- list(
    a = .numeric_test_policy(
      "exact-wide-a", predictor = 1e20, response = 1e20),
    b = .numeric_test_policy(
      "exact-wide-b", predictor = 1e20, response = 1e20))
  for (site in names(policies)) {
    policies[[site]]$capabilities$exact_gc <- .numeric_test_capability(
      "exact_gc", available = TRUE, verified = TRUE, exact = TRUE,
      supported_ring_bits = 63L:4096L, max_ring_bits = 4096L,
      max_frac_bits = 4095L)
    policies[[site]] <- .numeric_test_rehash(policies[[site]])
  }

  certificate <- dsVertClient:::.dsvert_numeric_preflight_from_policies(
    policies, workload, requested_backend = "exact_gc",
    requested_ring = 127L)

  expect_identical(certificate$effective_backend, "exact_gc")
  expect_gt(certificate$ring_bits, 127L)
  expect_lte(certificate$ring_bits, 4096L)
  expect_lte(certificate$truncated_output_required_ring_bits,
             certificate$ring_bits)
  expect_true(certificate$product_headroom_verified)
})

test_that("exact GC precision is minimal after fast floors and monotone", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 2L, "gaussian", 5L, FALSE, FALSE)
  budgets <- c(1e-3, 1e-6, 1e-9, 1e-12, 1e-15)
  certificates <- lapply(seq_along(budgets), function(index) {
    policies <- list(
      a = .numeric_test_policy(
        paste0("monotone-a-", index), error_budget = budgets[[index]]),
      b = .numeric_test_policy(
        paste0("monotone-b-", index), error_budget = budgets[[index]]))
    for (site in names(policies)) {
      policies[[site]]$capabilities$exact_gc <- .numeric_test_capability(
        "exact_gc", available = TRUE, verified = TRUE, exact = TRUE,
        supported_ring_bits = 63L:4096L, max_ring_bits = 4096L,
        max_frac_bits = 4095L)
      policies[[site]] <- .numeric_test_rehash(policies[[site]])
    }
    dsVertClient:::.dsvert_numeric_preflight_from_policies(
      policies, workload, requested_backend = "exact_gc",
      requested_ring = 63L)
  })

  required <- vapply(certificates, `[[`, numeric(1L),
                     "required_frac_bits")
  planned <- vapply(certificates, `[[`, numeric(1L), "frac_bits")
  rings <- vapply(certificates, `[[`, numeric(1L), "ring_bits")
  expect_true(all(diff(required) >= 0))
  expect_true(all(diff(planned) >= 0))
  expect_true(all(diff(rings) >= 0))
  expect_true(all(planned >= required))
  expect_true(all(vapply(certificates, function(certificate) {
    certificate$planned_total_numeric_error <=
      certificate$numeric_error_budget
  }, logical(1L))))
  for (certificate in certificates) {
    fast_floor <- if (certificate$ring_bits == 63L) 20L else 50L
    expect_identical(
      certificate$frac_bits,
      as.integer(max(fast_floor, certificate$required_frac_bits)))
  }
})

test_that("the exact planner preserves the Ring127 f50 fast route", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 2L, "gaussian", 5L, FALSE, FALSE)
  policies <- list(a = .numeric_test_policy("f50-a"),
                   b = .numeric_test_policy("f50-b"))
  for (site in names(policies)) {
    policies[[site]]$capabilities$exact_gc <- .numeric_test_capability(
      "exact_gc", available = TRUE, verified = TRUE, exact = TRUE,
      supported_ring_bits = 127L:4096L, max_ring_bits = 4096L,
      max_frac_bits = 4095L)
    policies[[site]] <- .numeric_test_rehash(policies[[site]])
  }

  certificate <- dsVertClient:::.dsvert_numeric_preflight_from_policies(
    policies, workload, requested_backend = "exact_gc",
    requested_ring = 127L)

  expect_identical(certificate$ring_bits, 127L)
  expect_identical(certificate$frac_bits, 50L)
  expect_lt(certificate$required_frac_bits, certificate$frac_bits)
})

test_that("log-domain planning reaches the exact Ring4096 boundary", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 2L, "poisson", 5L, FALSE, FALSE)
  policies <- list(
    a = .numeric_test_policy(
      "ring4096-a", eta = 1399, predictor = 1e300, response = 1),
    b = .numeric_test_policy(
      "ring4096-b", eta = 1399, predictor = 1e300, response = 1))
  for (site in names(policies)) {
    policies[[site]]$approximation$poisson$domain <- c(-2000, 2000)
    policies[[site]]$capabilities$exact_gc <- .numeric_test_capability(
      "exact_gc", available = TRUE, verified = TRUE, exact = TRUE,
      supported_ring_bits = 127L:4096L, max_ring_bits = 4096L,
      max_frac_bits = 4095L)
    policies[[site]] <- .numeric_test_rehash(policies[[site]])
  }

  certificate <- dsVertClient:::.dsvert_numeric_preflight_from_policies(
    policies, workload, requested_backend = "exact_gc",
    requested_ring = 127L)

  expect_identical(certificate$ring_bits, 4096L)
  expect_identical(certificate$required_ring_bits, 4096L)
  expect_identical(certificate$frac_bits, 50L)
  expect_true(is.infinite(certificate$public_product_bound))
  expect_true(is.finite(certificate$public_product_log2_bound))
  expect_true(certificate$product_headroom_verified)
})

test_that("a public bound beyond Ring4096 fails with exact requirements", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 2L, "poisson", 5L, FALSE, FALSE)
  policies <- list(
    a = .numeric_test_policy(
      "over4096-a", eta = 1399.5, predictor = 1e300, response = 1),
    b = .numeric_test_policy(
      "over4096-b", eta = 1399.5, predictor = 1e300, response = 1))
  for (site in names(policies)) {
    policies[[site]]$approximation$poisson$domain <- c(-2000, 2000)
    policies[[site]]$capabilities$exact_gc <- .numeric_test_capability(
      "exact_gc", available = TRUE, verified = TRUE, exact = TRUE,
      supported_ring_bits = 127L:4096L, max_ring_bits = 4096L,
      max_frac_bits = 4095L)
    policies[[site]] <- .numeric_test_rehash(policies[[site]])
  }

  error <- tryCatch(
    dsVertClient:::.dsvert_numeric_preflight_from_policies(
      policies, workload, requested_backend = "exact_gc",
      requested_ring = 127L),
    error = identity)

  expect_s3_class(error, "dsvert_numeric_backend_unrepresentable")
  expect_s3_class(error, "numeric_backend_unrepresentable")
  expect_identical(error$certificate$status,
                   "numeric_backend_unrepresentable")
  expect_identical(error$certificate$failed_backend, "exact_gc")
  expect_identical(error$required_ring_bits, 4097L)
  expect_identical(error$certificate$failed_frac_bits, 50L)
  expect_true(error$required_frac_bits <= error$certificate$failed_frac_bits)
  expect_identical(error$numeric_error_budget,
                   error$certificate$numeric_error_budget)
})

test_that("unavailable multiprecision metadata is never selected or named", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 2L, "gaussian", 5L, FALSE, FALSE)
  policies <- list(a = .numeric_test_policy("no-multi-a"),
                   b = .numeric_test_policy("no-multi-b"))
  for (site in names(policies)) {
    policies[[site]]$capabilities$multiprecision$max_ring_bits <- 4096L
    policies[[site]]$capabilities$multiprecision$max_frac_bits <- 4095L
    policies[[site]] <- .numeric_test_rehash(policies[[site]])
  }

  error <- tryCatch(
    dsVertClient:::.dsvert_numeric_preflight_from_policies(
      policies, workload, requested_backend = "multiprecision",
      requested_ring = 127L),
    error = identity)

  expect_s3_class(error, "numeric_backend_unavailable")
  expect_false(inherits(error, "numeric_backend_unrepresentable"))
  expect_true(is.na(error$certificate$effective_backend))
  expect_true(is.na(error$certificate$failed_backend))
})

test_that("missing exact primitives produce a structured unavailable error", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 2L, "gaussian", 5L, FALSE, FALSE)
  policies <- list(a = .numeric_test_policy("a"),
                   b = .numeric_test_policy("b"))
  for (site in names(policies)) {
    policies[[site]]$capabilities$ring63$exact_truncation <- FALSE
    policies[[site]]$capabilities$ring127$exact_truncation <- FALSE
    policies[[site]] <- .numeric_test_rehash(policies[[site]])
  }

  error <- tryCatch(
    dsVertClient:::.dsvert_numeric_preflight_from_policies(
      policies, workload, requested_backend = "auto", requested_ring = 63L),
    error = identity)

  expect_s3_class(error, "numeric_backend_unavailable")
  expect_s3_class(error$certificate, "dsvert_numeric_certificate")
  expect_identical(error$certificate$status,
                   "numeric_backend_unavailable")
  expect_identical(error$reason, "no_certified_backend")
})

test_that("approximation-domain failure is distinct and fail-closed", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 2L, "binomial", 5L, FALSE, FALSE)
  policies <- list(a = .numeric_test_policy("a", eta = 9),
                   b = .numeric_test_policy("b", eta = 9))

  error <- tryCatch(
    dsVertClient:::.dsvert_numeric_preflight_from_policies(
      policies, workload, requested_backend = "auto", requested_ring = 63L),
    error = identity)

  expect_s3_class(error, "approximation_domain_failure")
  expect_identical(error$certificate$status,
                   "approximation_domain_failure")
  expect_equal(error$certificate$approximation_domain, c(-8, 8))
})

test_that("missing numeric policy always fails and obsolete selectors do nothing", {
  testthat::local_mocked_bindings(
    .dsvert_fetch_numeric_policies = function(...) {
      stop("endpoint missing", call. = FALSE)
    },
    .package = "dsVertClient")
  withr::local_options(list(dsvert.allow_unattested_numeric = TRUE))
  withr::local_envvar(c(DSVERT_ALLOW_UNATTESTED_NUMERIC = "true"))
  datasources <- list(site_a = list(), site_b = list())

  expect_error(
    ds.vertNumericPreflight(100L, 2L, datasources = datasources),
    class = "numeric_backend_unavailable")
  expect_error(
    ds.vertNumericPreflight(
      100L, 2L, datasources = datasources,
      allow_unattested_numeric = TRUE),
    "unused argument")
  expect_false("allow_unattested_numeric" %in%
                 names(formals(ds.vertNumericPreflight)))
})

test_that("GLM rejects a missing numeric policy before any data endpoint", {
  data_calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_fetch_numeric_policies = function(...) {
      stop("endpoint missing", call. = FALSE)
    },
    .dsvert_aggregate_strict = function(...) {
      data_calls <<- data_calls + 1L
      stop("protected endpoint must not run", call. = FALSE)
    },
    .package = "dsVertClient")

  expect_error(
    ds.vertGLM(y ~ x, data = "D", verbose = FALSE,
               datasources = list(site_a = list())),
    class = "numeric_backend_unavailable")
  expect_identical(data_calls, 0L)
})

test_that("malformed policy fails without a compatibility conversion", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 2L, "gaussian", 5L, FALSE, FALSE)
  malformed <- .numeric_test_policy("a")
  malformed$policy_id <- "not-a-hash"

  expect_error(
    dsVertClient:::.dsvert_numeric_preflight_from_policies(
      list(site = malformed), workload),
    class = "numeric_backend_unavailable")
})

test_that("a custodian-disabled numeric policy remains a structured refusal", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 2L, "gaussian", 5L, FALSE, FALSE)
  disabled <- .numeric_test_policy("a")
  disabled$enabled <- FALSE
  disabled <- .numeric_test_rehash(disabled)

  error <- tryCatch(
    dsVertClient:::.dsvert_numeric_preflight_from_policies(
      list(site = disabled), workload),
    error = identity)

  expect_s3_class(error, "numeric_backend_unavailable")
  expect_identical(error$reason, "custodian_numeric_policy_disabled")
})

test_that("non-identifiability has its own condition and no estimator fallback", {
  condition <- dsVertClient:::.dsvert_non_identifiable(
    "information matrix is singular")

  expect_s3_class(condition, "non_identifiable")
  expect_identical(condition$reason, "singular_information")
  expect_false(inherits(condition, "numeric_backend_unavailable"))

  glm_formals <- names(formals(ds.vertGLM))
  expect_true(all(c("numeric_backend", "dp_analysis_id") %in% glm_formals))
  expect_false("allow_unattested_numeric" %in% glm_formals)
  expect_identical(formals(ds.vertGLM)$dp_analysis_id, NULL)
  expect_identical(glm_formals[match("ring", glm_formals) + 1L],
                   "binomial_sigmoid_intervals")
  body_text <- paste(deparse(body(ds.vertGLM)), collapse = "\n")
  expect_false(grepl("MASS::ginv|qr.solve\\(|firth\\(",
                     body_text, ignore.case = TRUE))
  expect_match(body_text, ".dsvert_solve_identifiable", fixed = TRUE)
  expect_lt(regexpr(".dsvert_require_numeric_policies", body_text,
                    fixed = TRUE)[[1L]],
            regexpr("session_id <-", body_text, fixed = TRUE)[[1L]])
})

test_that("Ring127 planning accounts for the raw frac50 product", {
  physical_limit <- dsVertClient:::.dsvert_numeric_raw_product_limit(
    127L, 50L)
  expect_equal(physical_limit, 2^26)
  expect_gt(dsVertClient:::.dsvert_numeric_required_bits(
    public_bound = 1, product_bound = 2^26, frac_bits = 50L), 127L)

  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 1L, "gaussian", 1L, FALSE, FALSE)
  policies <- list(
    a = .numeric_test_policy("raw-a", eta = 2,
                             predictor = 2^13, response = 2^13),
    b = .numeric_test_policy("raw-b", eta = 2,
                             predictor = 2^13, response = 2^13))
  error <- tryCatch(
    dsVertClient:::.dsvert_numeric_preflight_from_policies(
      policies, workload, requested_backend = "ring127",
      requested_ring = 127L),
    error = identity)

  expect_s3_class(error, "numeric_backend_unavailable")
  expect_true("raw_product_capacity" %in%
                error$certificate$attempts$ring127$reasons)
  expect_false(error$certificate$attempts$ring127$
                 public_product_within_raw_ring)
})

test_that("exact GC accepts direct-wide products when the truncated output fits", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 1L, "gaussian", 1L, FALSE, FALSE)
  policies <- list(
    a = .numeric_test_policy("direct-wide-a", eta = 2,
                             predictor = 2^13, response = 2^13),
    b = .numeric_test_policy("direct-wide-b", eta = 2,
                             predictor = 2^13, response = 2^13))
  for (site in names(policies)) {
    policies[[site]]$capabilities$exact_gc <- .numeric_test_capability(
      "exact_gc", available = TRUE, verified = TRUE, exact = TRUE,
      supported_ring_bits = 127L, max_ring_bits = 127L,
      max_frac_bits = 50L)
    policies[[site]] <- .numeric_test_rehash(policies[[site]])
  }

  certificate <- dsVertClient:::.dsvert_numeric_preflight_from_policies(
    policies, workload, requested_backend = "exact_gc",
    requested_ring = 127L)

  expect_identical(certificate$effective_backend, "exact_gc")
  expect_identical(certificate$ring_bits, 127L)
  expect_false(certificate$public_product_within_raw_ring)
  expect_gt(certificate$raw_product_required_ring_bits,
            certificate$ring_bits)
  expect_lte(certificate$truncated_output_required_ring_bits,
             certificate$ring_bits)
  expect_identical(
    certificate$product_headroom_proof,
    "truncated_output_headroom_direct_wide")
  expect_true(certificate$product_headroom_verified)
  expect_false("raw_product_capacity" %in%
                 certificate$attempts$exact_gc$reasons)
})

test_that("ring planning includes raw product accumulators before truncation", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    n_obs = 3000000L, n_predictors = 1L, family = "gaussian",
    max_iter = 1L, compute_se = FALSE, compute_deviance = FALSE)
  policies <- list(a = .numeric_test_policy(
    "accum-a", eta = 1, predictor = 1, response = 1))

  certificate <- dsVertClient:::.dsvert_numeric_preflight_from_policies(
    policies, workload, requested_backend = "auto", requested_ring = 63L)

  # The individual product is tiny, but sum_i x_i * residual_i remains at
  # scale 2^(2f) until truncation and exceeds Ring63's 2^22 real-value limit.
  expect_gte(certificate$public_product_bound, 6000000)
  expect_true("raw_product_capacity" %in%
                certificate$attempts$ring63$reasons)
  expect_identical(certificate$effective_backend, "ring127")
})

test_that("policy hash authenticates the complete canonical policy body", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 1L, "gaussian", 1L, FALSE, FALSE)
  policy <- .numeric_test_policy("hash")
  policy$bounds$max_abs_predictor <- 3

  error <- tryCatch(
    dsVertClient:::.dsvert_numeric_preflight_from_policies(
      list(site = policy), workload),
    error = identity)
  expect_s3_class(error, "numeric_backend_unavailable")
  expect_identical(error$reason, "malformed_custodian_policy")
  expect_match(conditionMessage(error), "does not authenticate")
})

test_that("uncertified bounded backends have no compatibility fallback", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 2L, "gaussian", 2L, FALSE, FALSE)
  policies <- list(a = .numeric_test_policy("legacy-a"),
                   b = .numeric_test_policy("legacy-b"))
  for (site in names(policies)) {
    for (backend in c("ring63", "ring127")) {
      policies[[site]]$capabilities[[backend]]$e2e_verified <- FALSE
      policies[[site]]$capabilities[[backend]]$fail_closed_overflow <- FALSE
      policies[[site]]$capabilities[[backend]]$runtime_bounds_enforced <- FALSE
      policies[[site]]$capabilities[[backend]]$exact_truncation <- FALSE
      policies[[site]]$capabilities[[backend]]$exact_comparison <- FALSE
      policies[[site]]$capabilities[[backend]]$truncation_semantics <-
        "local_probabilistic"
    }
    policies[[site]] <- .numeric_test_rehash(policies[[site]])
  }

  expect_error(
    dsVertClient:::.dsvert_numeric_preflight_from_policies(
      policies, workload),
    class = "numeric_backend_unavailable")
  expect_error(
    dsVertClient:::.dsvert_numeric_preflight_from_policies(
      policies, workload, allow_legacy_unattested = TRUE),
    "unused argument")
})

test_that("execution attestations bind policy session dataset and variables", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 1L, "gaussian", 1L, FALSE, FALSE)
  policies <- list(site = .numeric_test_policy("binding"))
  certificate <- dsVertClient:::.dsvert_numeric_preflight_from_policies(
    policies, workload)
  sid <- "01234567-89ab-4def-8123-456789abcdef"
  binding <- dsVertClient:::.dsvert_numeric_attestation_binding(
    "glm_standardized_input", policies$site$policy_id, sid, "D",
    c("x", "y"), "gaussian", 63L, 100L)
  attestation <- list(
    schema_version = 1L,
    kind = "glm_standardized_input",
    policy_id = policies$site$policy_id,
    binding_id = binding,
    ring = 63L,
    n = 100L,
    checks = list(ieee_finite = TRUE, encoded_input_bounds = TRUE),
    runtime_input_bounds_enforced = TRUE,
    runtime_intermediate_bounds_enforced = FALSE,
    observed_extrema_released = FALSE,
    attestation_scope = "inputs only")

  validated <- dsVertClient:::.dsvert_numeric_validate_attestation(
    attestation, certificate, "site", "glm_standardized_input", sid,
    "D", c("x", "y"), "gaussian", 63L, 100L)
  certificate <- dsVertClient:::.dsvert_numeric_attach_attestations(
    certificate, validated, all_inputs = TRUE)
  certificate <- dsVertClient:::.dsvert_numeric_finalize_certificate(
    certificate, c("(Intercept)" = 0, x = 1), converged = TRUE)
  expect_true(certificate$runtime_input_bounds_attested)
  expect_true(certificate$result_postcondition_verified)
  expect_false(certificate$runtime_intermediate_bounds_attested)
  expect_false(certificate$numerically_certified)

  expect_error(
    dsVertClient:::.dsvert_numeric_validate_attestation(
      attestation, certificate, "site", "glm_standardized_input",
      sub("0", "1", sid), "D", c("x", "y"), "gaussian", 63L, 100L),
    class = "numeric_backend_unavailable")
})

test_that("final certification requires convergence and end-to-end error bounds", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 1L, "gaussian", 1L, FALSE, FALSE)
  certificate <- dsVertClient:::.dsvert_numeric_preflight_from_policies(
    list(site = .numeric_test_policy("final-error-bound")), workload)

  certificate$runtime_input_bounds_attested <- TRUE
  certificate$runtime_intermediate_bounds_attested <- TRUE
  expect_identical(certificate$estimability_status, "not_evaluated")
  without_bounds <- dsVertClient:::.dsvert_numeric_finalize_certificate(
    certificate, c("(Intercept)" = 0, x = 1), converged = TRUE)
  expect_false(without_bounds$numerically_certified)
  expect_false(without_bounds$numeric_error_bound_certified)

  certificate$aggregate_output_error_bound <- 1e-6
  certificate$estimator_error_bound <- 1e-5
  not_evaluated <- dsVertClient:::.dsvert_numeric_finalize_certificate(
    certificate, c("(Intercept)" = 0, x = 1), converged = TRUE)
  expect_false(not_evaluated$numerically_certified)
  expect_identical(not_evaluated$estimability_status, "not_evaluated")

  certificate <- dsVertClient:::.dsvert_numeric_mark_estimable(certificate)
  not_converged <- dsVertClient:::.dsvert_numeric_finalize_certificate(
    certificate, c("(Intercept)" = 0, x = 1), converged = FALSE)
  expect_false(not_converged$numerically_certified)

  over_budget <- certificate
  over_budget$estimator_error_bound <- over_budget$numeric_error_budget * 2
  over_budget <- dsVertClient:::.dsvert_numeric_finalize_certificate(
    over_budget, c("(Intercept)" = 0, x = 1), converged = TRUE)
  expect_false(over_budget$numerically_certified)

  certified <- dsVertClient:::.dsvert_numeric_finalize_certificate(
    certificate, c("(Intercept)" = 0, x = 1), converged = TRUE)
  expect_true(certified$numerically_certified)
  expect_true(certified$numeric_error_bound_certified)
  expect_true(certified$estimator_converged)
  expect_identical(certified$estimability_status, "estimable")
})

test_that("non-identifiability carries a failed numeric certificate", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 1L, "gaussian", 1L, FALSE, FALSE)
  certificate <- dsVertClient:::.dsvert_numeric_preflight_from_policies(
    list(site = .numeric_test_policy("non-identifiable-certificate")),
    workload)

  error <- tryCatch(
    dsVertClient:::.dsvert_solve_identifiable(
      matrix(c(1, 1, 1, 1), 2L), c(1, 1),
      context = "Protected model", reason = "rank_deficient_design",
      symmetric = TRUE, certificate = certificate),
    error = identity)

  expect_s3_class(error, "non_identifiable")
  expect_s3_class(error$certificate, "dsvert_numeric_certificate")
  expect_identical(error$certificate$status, "non_identifiable")
  expect_identical(error$certificate$estimability_status,
                   "non_identifiable")
  expect_identical(error$certificate$estimability_reason,
                   "rank_deficient_design")
  expect_true(error$certificate$estimability_evaluated)
  expect_false(error$certificate$numerically_certified)
})

test_that("public beta guard fails before eta leaves its declared domain", {
  workload <- dsVertClient:::.dsvert_numeric_glm_workload(
    100L, 1L, "binomial", 1L, FALSE, FALSE)
  policies <- list(site = .numeric_test_policy(
    "eta-guard", eta = 2, predictor = 2))
  certificate <- dsVertClient:::.dsvert_numeric_preflight_from_policies(
    policies, workload)

  expect_invisible(dsVertClient:::.dsvert_numeric_assert_eta_bound(
    beta = 0.5, intercept = 0.5, certificate = certificate))
  expect_error(
    dsVertClient:::.dsvert_numeric_assert_eta_bound(
      beta = 1, intercept = 0.5, certificate = certificate),
    class = "approximation_domain_failure")
})

test_that("every numeric certificate must carry finite eta bounds", {
  certificate <- list(
    status = "preflight_eligible",
    public_predictor_bound = NA_real_,
    public_linear_predictor_bound = NA_real_,
    public_offset_bound = NA_real_
  )
  expect_error(
    dsVertClient:::.dsvert_numeric_assert_eta_bound(
      beta = 1, intercept = 0, certificate = certificate),
    class = "numeric_backend_unavailable"
  )
})

test_that("integer workload values cannot overflow R integer serialization", {
  expect_error(
    dsVertClient:::.dsvert_numeric_glm_workload(
      3e9, 1L, "gaussian", 1L, FALSE, FALSE),
    "finite integer")
  expect_error(ds.vertNumericPreflight(
    100L, 1L, ring = 63.5,
    datasources = list(site = list())), "ring must be 63 or 127")
})
test_that("numeric-policy handshake uses the strict direct DSI dispatcher", {
  called <- FALSE
  policy <- list(schema_version = 1L)
  conns <- list(site = structure(
    list(sid = "numeric-policy-session"),
    class = c("DSLiteConnection", "list")))
  result <- testthat::with_mocked_bindings(
    .dsvert_fetch_numeric_policies(
      conns, .aggregate = DSI::datashield.aggregate),
    .dsvert_dsi_direct_aggregate = function(
        conns, expr, async, error, errors.print, ...) {
      called <<- TRUE
      list(site = policy)
    })

  expect_true(called)
  expect_identical(result, list(site = policy))
})
