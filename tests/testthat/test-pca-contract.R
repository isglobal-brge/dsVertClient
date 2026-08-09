make_cor_result <- function(R, n = 100L) {
  p <- ncol(R)
  vars <- paste0("x", seq_len(p))
  dimnames(R) <- list(vars, vars)
  lower <- R - 0.01
  upper <- R + 0.01
  lower[lower < -1] <- -1
  upper[upper > 1] <- 1
  tolerance <- 256 * .Machine$double.eps * p
  out <- list(
    released = TRUE, source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE, legacy_exact_route_called = FALSE,
    history_gate = TRUE, request_limit = FALSE, operation_limit = TRUE,
    psd_projection_applied = TRUE, capsule_id = strrep("a", 64L),
    analysis_id = "cohort::site_a", correlation = R,
    source_artifact_family = "gaussian_models",
    estimand_missingness = "complete_case_joint", pca_eligible = TRUE,
    correlation_raw_complete_case = R,
    correlation_raw_pairwise = NULL, n_obs = n,
    complete_case_n = matrix(n, p, p, dimnames = dimnames(R)),
    pairwise_n = NULL, var_names = vars,
    epsilon = 1, delta = 2^-100, mechanism = "fixture",
    psd_projection = list(
      method = "fixture_psd_projection",
      numerical_tolerance = tolerance),
    correlation_95_enclosure_raw_estimand_around_projected_release = list(
      lower = lower, upper = upper),
    cross_owner_state = "reserved_not_materialized")
  class(out) <- c("ds.vertDPCor", "ds.cor", "list")
  out
}

test_that("PCA clamps only quantization-scale negative eigenvalues", {
  R <- matrix(c(1, 1 + 1e-14, 1 + 1e-14, 1), 2L)
  fit <- ds.vertPCA(cor_result = make_cor_result(R), verbose = FALSE)

  expect_equal(min(fit$eigenvalues), 0)
  expect_identical(fit$psd_diagnostic$tiny_negative_eigenvalues_clamped, 1L)
  expect_lt(abs(fit$psd_diagnostic$minimum_eigenvalue),
            fit$psd_diagnostic$tolerance)
})

test_that("PCA rejects materially indefinite or malformed correlations", {
  indefinite <- matrix(c(1, 1.2, 1.2, 1), 2L)
  expect_error(
    ds.vertPCA(cor_result = make_cor_result(indefinite), verbose = FALSE),
    "PSD certificate"
  )

  asymmetric <- matrix(c(1, 0.1, 0.2, 1), 2L)
  expect_error(
    ds.vertPCA(cor_result = make_cor_result(asymmetric), verbose = FALSE),
    "signed DP correlation matrix is invalid"
  )

  nonfinite <- diag(2)
  nonfinite[1, 2] <- nonfinite[2, 1] <- NA_real_
  expect_error(
    ds.vertPCA(cor_result = make_cor_result(nonfinite), verbose = FALSE),
    "signed DP correlation matrix is invalid"
  )
})

test_that("PCA component count is explicit and bounded", {
  cor <- make_cor_result(diag(3))
  expect_error(ds.vertPCA(cor_result = cor, n_components = 0,
                          verbose = FALSE), "n_components")
  expect_error(ds.vertPCA(cor_result = cor, n_components = 4,
                          verbose = FALSE), "n_components")
  expect_equal(ncol(ds.vertPCA(cor_result = cor, n_components = 2,
                               verbose = FALSE)$loadings), 2L)
  tied <- ds.vertPCA(cor_result = cor, n_components = 2, verbose = FALSE)
  expect_false(tied$loading_identifiability$individual_directions_identifiable)
  expect_true(
    tied$loading_identifiability$tied_or_unresolved_eigenspaces_identifiable)
})
