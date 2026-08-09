test_that("effective GLM ring reflects the protocol actually executed", {
  effective_ring <- dsVertClient:::.dsvert_effective_glm_ring

  expect_identical(effective_ring("gaussian", 63L), 63L)
  expect_identical(effective_ring("gaussian", 127L), 127L)
  expect_identical(effective_ring("binomial", 63L), 127L)
  expect_identical(effective_ring("poisson", 63L), 127L)
})

test_that("Gaussian one-shot is used only for its exact contract", {
  eligible <- dsVertClient:::.k2_use_gaussian_oneshot
  base <- list(
    family = "gaussian", no_intercept = FALSE, weights_active = FALSE,
    gradient_only = FALSE, start = NULL, offset_active = FALSE,
    lambda = 0, ring = 63L, option_enabled = TRUE
  )

  expect_true(do.call(eligible, base))
  expect_false(do.call(eligible, utils::modifyList(base, list(lambda = 1e-8))))
  expect_false(do.call(eligible, utils::modifyList(base, list(offset_active = TRUE))))
  expect_false(do.call(eligible, utils::modifyList(base, list(ring = 127L))))
  expect_false(do.call(eligible, utils::modifyList(base, list(weights_active = TRUE))))
  expect_false(do.call(eligible, utils::modifyList(base, list(no_intercept = TRUE))))

  expect_error(
    ds.vertGLM(y ~ x, data = "D", family = "gaussian", offset = "baseline",
               datasources = list()),
    "Gaussian offsets are not implemented"
  )
})

test_that("final mean refresh cannot leave Gaussian or Poisson mu stale", {
  refresh <- dsVertClient:::.glm_refresh_final_mean
  datasources <- list(site_a = list(), site_b = list())
  aggregate_calls <- list()
  link_calls <- list()
  fake_aggregate <- function(conns, expr, ...) {
    aggregate_calls[[length(aggregate_calls) + 1L]] <<- expr
    invisible(list(ok = TRUE))
  }
  fake_link <- function(...) {
    link_calls[[length(link_calls) + 1L]] <<- list(...)
    invisible("secure_mu_share")
  }
  link_args <- list(
    n = 20L, datasources = datasources, dealer_ci = 2L,
    server_list = c("site_a", "site_b"),
    server_names = c("site_a", "site_b"), y_server = "site_a",
    nl = "site_b", transport_pks = list(), session_id = "session",
    .dsAgg = fake_aggregate, .sendBlob = function(...) NULL
  )

  expect_identical(
    refresh(
      family = "gaussian", datasources = datasources,
      server_list = c("site_a", "site_b"),
      server_names = c("site_a", "site_b"), session_id = "session",
      .dsAgg = fake_aggregate, link_args = link_args,
      .share_link = fake_link),
    "identity"
  )
  expect_length(aggregate_calls, 2L)
  expect_true(all(vapply(
    aggregate_calls,
    function(expr) identical(as.character(expr[[1L]]), "k2IdentityLinkDS"),
    logical(1L))))
  expect_length(link_calls, 0L)

  aggregate_calls <- list()
  expect_identical(
    refresh(
      family = "poisson", datasources = datasources,
      server_list = c("site_a", "site_b"),
      server_names = c("site_a", "site_b"), session_id = "session",
      .dsAgg = fake_aggregate, link_args = link_args,
      .share_link = fake_link),
    "exp"
  )
  expect_length(aggregate_calls, 0L)
  expect_length(link_calls, 1L)
  expect_identical(link_calls[[1L]]$family, "poisson")

  link_calls <- list()
  expect_identical(
    refresh(
      family = "binomial", datasources = datasources,
      server_list = c("site_a", "site_b"),
      server_names = c("site_a", "site_b"), session_id = "session",
      .dsAgg = fake_aggregate, link_args = link_args,
      .share_link = fake_link),
    "not_needed"
  )
  expect_length(link_calls, 0L)

  k2_body <- paste(deparse(body(dsVertClient:::.k2_strict_loop)),
                   collapse = "\n")
  k3_body <- paste(deparse(body(dsVertClient:::.k3_ring63_gradient_loop)),
                   collapse = "\n")
  expect_match(k2_body, ".glm_refresh_final_mean", fixed = TRUE)
  expect_match(k3_body, ".glm_refresh_final_mean", fixed = TRUE)
})

test_that("Gaussian fit statistics count dispersion and suppress invalid null models", {
  fit_statistics <- dsVertClient:::.dsvert_glm_fit_statistics
  ordinary <- fit_statistics(
    family = "gaussian", deviance = 120, n_obs = 100L,
    n_parameters = 3L, y_sd = 2, weights_active = FALSE,
    offset_active = FALSE, intercept_included = TRUE)

  expect_equal(ordinary$null_deviance, 99 * 4)
  expect_equal(ordinary$pseudo_r2, 1 - 120 / (99 * 4))
  expect_equal(
    ordinary$aic,
    100 * (log(2 * pi) + 1 + log(120 / 100)) + 2 * (3 + 1))
  expect_identical(ordinary$aic_type,
                   "gaussian_ml_including_dispersion")

  weighted <- fit_statistics(
    family = "gaussian", deviance = 120, n_obs = 100L,
    n_parameters = 3L, y_sd = 2, weights_active = TRUE,
    offset_active = FALSE, intercept_included = TRUE)
  expect_true(is.na(weighted$null_deviance))
  expect_true(is.na(weighted$pseudo_r2))
  expect_true(is.na(weighted$aic))

  offset <- fit_statistics(
    family = "gaussian", deviance = 120, n_obs = 100L,
    n_parameters = 3L, y_sd = 2, weights_active = FALSE,
    offset_active = TRUE, intercept_included = TRUE)
  expect_true(is.na(offset$null_deviance))
  expect_true(is.na(offset$pseudo_r2))

  no_intercept <- fit_statistics(
    family = "gaussian", deviance = 120, n_obs = 100L,
    n_parameters = 2L, y_sd = 2, weights_active = FALSE,
    offset_active = FALSE, intercept_included = FALSE)
  expect_true(is.na(no_intercept$null_deviance))
})

test_that("a real matching PSI manifest becomes cohort metadata", {
  manifest_hash <- paste(rep("a", 64L), collapse = "")
  status <- list(
    aligned = TRUE, n_common = 20L,
    manifests = list(
      site_a = list(version = 1L, hash = manifest_hash, n = 20L),
      site_b = list(version = 1L, hash = manifest_hash, n = 20L)))

  metadata <- dsVertClient:::.dsvert_glm_alignment_metadata(status, 20L)
  expect_true(metadata$alignment_attested)
  expect_identical(metadata$alignment_manifest_hash, manifest_hash)
  expect_identical(metadata$cohort_id, manifest_hash)

  mismatch <- dsVertClient:::.dsvert_glm_alignment_metadata(status, 19L)
  expect_false(mismatch$alignment_attested)
  expect_null(mismatch$alignment_manifest_hash)
  expect_null(mismatch$cohort_id)
})

test_that("Gaussian client-side inference uses t and F references", {
  coefficients <- c(`(Intercept)` = 2, x = 1)
  covariance <- diag(c(0.25, 0.0625))
  dimnames(covariance) <- list(names(coefficients), names(coefficients))
  fit <- list(
    coefficients = coefficients,
    std_errors = sqrt(diag(covariance)), covariance = covariance,
    family = "gaussian", n_obs = 12L, n_vars = 2L,
    converged = TRUE, lambda = 0)
  class(fit) <- c("ds.glm", "list")

  ci <- ds.vertConfint(fit)
  expect_equal(ci["x", "lower"], 1 - qt(0.975, df = 10) * 0.25)
  expect_identical(attr(ci, "distribution"), "t")
  expect_identical(attr(ci, "df"), 10L)

  wald <- ds.vertWald(fit, "x")
  expect_equal(wald$t, 4)
  expect_identical(wald$distribution, "t")
  expect_identical(wald$df, 10L)
  expect_equal(wald$p_value, 2 * pt(4, df = 10, lower.tail = FALSE))

  contrast <- ds.vertContrast(fit, K = matrix(c(0, 1), nrow = 1L))
  expect_equal(contrast$statistic, 16)
  expect_equal(contrast$wald_statistic, 16)
  expect_identical(contrast$distribution, "F")
  expect_identical(contrast$df, 1L)
  expect_identical(contrast$df_residual, 10L)
  expect_equal(contrast$p_value, pf(16, 1, 10, lower.tail = FALSE))
})

test_that("Poisson LR remains enabled only for refreshed canonical fits", {
  make_poisson_fit <- function(coefs, deviance) {
    covariance <- diag(rep(0.1, length(coefs)), nrow = length(coefs))
    dimnames(covariance) <- list(names(coefs), names(coefs))
    fit <- list(
      coefficients = coefs, std_errors = sqrt(diag(covariance)),
      covariance = covariance, family = "poisson", n_obs = 100L,
      n_vars = length(coefs), converged = TRUE, lambda = 0,
      deviance = deviance, deviance_type = "canonical_poisson",
      weights = NULL, offset = "log_exposure", cohort_id = "cohort-a")
    class(fit) <- c("ds.glm", "list")
    fit
  }

  reduced <- make_poisson_fit(c(`(Intercept)` = 0.1), 110)
  full <- make_poisson_fit(c(`(Intercept)` = 0.1, x = 0.2), 100)
  expect_equal(ds.vertLR(reduced, full)$statistic, 10)
})
