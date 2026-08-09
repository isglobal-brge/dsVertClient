test_that("identifiability solver fails closed with a stable condition", {
  err <- tryCatch(
    dsVertClient:::.dsvert_solve_identifiable(
      matrix(c(1, 1, 1, 1), 2L), c(1, 1),
      context = "Test model",
      reason = "rank_deficient_test_design",
      symmetric = TRUE),
    error = identity)

  expect_s3_class(err, "non_identifiable")
  expect_identical(err$reason, "rank_deficient_test_design")
  expect_match(conditionMessage(err), "No ridge, Firth, or pseudoinverse")
})

test_that("identifiability solver preserves a regular estimable system", {
  out <- dsVertClient:::.dsvert_solve_identifiable(
    diag(c(2, 4)), c(2, 8), context = "Test model", symmetric = TRUE)

  expect_equal(as.numeric(out), c(1, 2))
})

test_that("indefinite information is rejected even when algebraically invertible", {
  err <- tryCatch(
    dsVertClient:::.dsvert_solve_identifiable(
      diag(c(1, -1)), c(0, 0),
      context = "Indefinite model",
      reason = "indefinite_information",
      symmetric = TRUE),
    error = identity)

  expect_s3_class(err, "non_identifiable")
  expect_identical(err$reason, "indefinite_information")
  expect_match(conditionMessage(err), "not positive definite")
})

test_that("binomial separation and singular information have a typed cause", {
  err <- tryCatch(
    dsVertClient:::.dsvert_solve_identifiable(
      diag(c(1, 0)), c(0, 0),
      context = "Protected binomial model",
      reason = "separation_or_singular_information",
      symmetric = TRUE),
    error = identity)

  expect_s3_class(err, "non_identifiable")
  expect_identical(err$reason, "separation_or_singular_information")
})

test_that("empty contingency margins fail instead of using epsilon divisors", {
  err <- tryCatch(
    dsVertClient:::.dsvert_chisq_compute(
      matrix(c(4, 0, 6, 0), 2L),
      row_margins = c(10, 0), col_margins = c(4, 6), n = 10),
    error = identity)

  expect_s3_class(err, "non_identifiable")
  expect_identical(err$reason, "degenerate_contingency_margins")
})

test_that("GLMM and LASSO singular systems do not acquire a ridge fallback", {
  glmm_err <- tryCatch(
    dsVertClient:::.ds_glmm_safe_solve(
      matrix(c(1, 1, 1, 1), 2L), c(1, 1)),
    error = identity)
  expect_s3_class(glmm_err, "non_identifiable")
  expect_identical(glmm_err$reason,
                   "singular_glmm_fixed_effect_information")

  fit <- list(
    coefficients = c("(Intercept)" = 0, x = 0),
    covariance = diag(2),
    covariance_information = matrix(c(1, 1, 1, 1), 2L),
    n_obs = 20L,
    deviance = 10,
    family = "gaussian",
    lambda = 0)
  class(fit) <- c("ds.glm", "list")
  lasso_err <- tryCatch(
    ds.vertLASSOProximal(fit, lambda = 0.1),
    error = identity)
  expect_s3_class(lasso_err, "non_identifiable")
  expect_identical(lasso_err$reason, "rank_deficient_lasso_design")
})

test_that("regularization and imputation alternatives require explicit input", {
  expect_identical(formals(ds.vertGEE)$lambda, 0)
  expect_identical(formals(ds.vertMI)$lambda, 0)
  expect_identical(eval(formals(ds.vertMI)$intercept_only),
                   c("error", "aggregate"))
})

test_that("statistical routes contain no automatic estimator-changing solve", {
  src <- function(fun) paste(deparse(body(fun)), collapse = "\n")

  expect_match(src(ds.vertGLM), ".dsvert_solve_identifiable", fixed = TRUE)
  expect_false(grepl("qr.solve", src(ds.vertNBFullRegTheta), fixed = TRUE))
  expect_false(grepl("solve(fisher +", src(ds.vertNBFullRegTheta),
                     fixed = TRUE))

  lmm_src <- src(ds.vertLMM)
  expect_false(grepl("fit <- ols_ref", lmm_src, fixed = TRUE))
  expect_match(lmm_src, "No OLS fallback was applied", fixed = TRUE)

  cox_src <- src(ds.vertCoxDiscreteNonDisclosive)
  expect_false(grepl("ridge bumped", cox_src, fixed = TRUE))
  expect_match(cox_src, "information_step", fixed = TRUE)

  multinom_src <- src(ds.vertMultinomJointNewton)
  expect_false(grepl("0.1 * g_stacked", multinom_src, fixed = TRUE))
  expect_false(grepl("B_reg", multinom_src, fixed = TRUE))
  expect_match(multinom_src, "step_solver_history", fixed = TRUE)
})
