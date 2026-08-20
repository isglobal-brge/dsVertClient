# Tests for client-side wrapper functions (IPW shape, multinom shape and NB
# shape). Historical ds.glm LASSO sketches must stay unavailable.

library(testthat)

mock_fit <- function(coefs = c(`(Intercept)` = 1, x1 = 0.5, x2 = -0.3, x3 = 0.05)) {
  ses <- rep(0.1, length(coefs)); names(ses) <- names(coefs)
  fit <- list(
    coefficients = coefs,
    std_errors = ses,
    z_values = coefs / ses,
    p_values = 2 * pnorm(abs(coefs / ses), lower.tail = FALSE),
    family = "gaussian", n_obs = 100L, n_vars = length(coefs),
    lambda = 1e-4, deviance = 50, null_deviance = 75, pseudo_r2 = 0.33,
    aic = NA_real_, iterations = 5L, converged = TRUE)
  class(fit) <- c("ds.glm", "list")
  fit
}

# =============================================================================
# Legacy ds.vertLASSO boundary
# =============================================================================
test_that("LASSO rejects the retired ds.glm thresholding sketch", {
  fit <- mock_fit()
  expect_error(
    ds.vertLASSO(fit, lambda_1 = 0.3, alpha_grid = c(1, 0.5)),
    "validated ds.vertDPGaussian")
})

# =============================================================================
# ds.vertMultinom input validation
# =============================================================================
test_that("multinom requires >=2 non-reference classes", {
  expect_error(
    ds.vertMultinom(y ~ x1 + x2, data = "D", classes = c("A")),
    "Need at least 2")
})

# =============================================================================
# ds.vertNB class inheritance
# =============================================================================
test_that("ds.vertNB object inherits from ds.glm", {
  fit <- mock_fit()
  out <- c(unclass(fit), list(theta = 2.5, nb_correction = "placeholder"))
  class(out) <- c("ds.vertNB", "ds.glm", "list")
  expect_true(inherits(out, "ds.glm"))
  expect_true(inherits(out, "ds.vertNB"))
})

test_that("lowercase ds.vert aliases are exported", {
  aliases <- c(
    "ds.vert.align", "ds.vert.desc", "ds.vert.cor", "ds.vert.pca",
    "ds.vert.chisq", "ds.vert.fisher", "ds.vert.chisq_cross",
    "ds.vert.glm", "ds.vert.cox", "ds.vert.coxph", "ds.vert.nb",
    "ds.vert.multinom", "ds.vert.ordinal", "ds.vert.lmm",
    "ds.vert.gee", "ds.vert.glmm", "ds.vert.ipw", "ds.vert.mi",
    "ds.vert.lasso", "ds.vert.lasso_iter", "ds.vert.lasso_proximal",
    "ds.vert.lasso_1step", "ds.vert.lasso_cv", "ds.vert.lr",
    "ds.vert.confint", "ds.vert.wald", "ds.vert.contrast")
  expect_true(all(vapply(aliases, exists, logical(1),
                         envir = asNamespace("dsVertClient"),
                         inherits = FALSE)))
})

test_that("GLMM frontdoor exposes the validated PQL product method", {
  expect_identical(formals(ds.vert.glmm)$method, "pql")
  expect_false("ds.vert.glmer" %in% getNamespaceExports("dsVertClient"))
  expect_false("ds.vertGLMMLaplace" %in% getNamespaceExports("dsVertClient"))
})

test_that("frontdoors expose adaptive precision and method selectors", {
  expect_identical(formals(ds.vert.nb)$method, "accurate")
  expect_identical(formals(ds.vert.cox)$method,
                   quote(c("profile", "discrete")))
  expect_identical(formals(ds.vert.lasso_iter)$method,
                   quote(c("auto", "accurate", "fast")))
  expect_identical(formals(ds.vert.glm)$precision,
                   quote(c("auto", "high", "fast")))
  expect_identical(formals(ds.vert.gee)$precision,
                   quote(c("auto", "high", "fast")))
  expect_identical(formals(ds.vert.ipw)$precision,
                   quote(c("auto", "high", "fast")))
})

test_that("frontdoor precision helper only raises binomial precision safely", {
  ns <- asNamespace("dsVertClient")
  f <- get(".dsvert_apply_binomial_precision", envir = ns)
  expect_equal(f(list(family = "binomial"), precision = "auto")$
                 binomial_sigmoid_intervals, 150L)
  expect_null(f(list(family = "binomial"), precision = "fast")$
                binomial_sigmoid_intervals)
  expect_null(f(list(family = "gaussian"), precision = "auto")$
                binomial_sigmoid_intervals)
  expect_equal(f(list(), precision = "auto", force_binomial = TRUE)$
                 binomial_sigmoid_intervals, 150L)
})

test_that("frontdoor route metadata distinguishes every supported K mode", {
  x <- list(coefficients = c(a = 1))
  set_frontdoor <- dsVertClient:::.dsvert_set_frontdoor
  modes <- vapply(
    list(1L, 2L, 3L, 7L, NULL, NA_integer_, 0L),
    function(K) set_frontdoor(
      x, frontdoor = "ds.vert.demo", route = "backend", K = K)$k_mode,
    character(1L))

  expect_identical(modes,
                   c("K=1", "K=2", "K>=3", "K>=3",
                     "unknown", "unknown", "unknown"))
  out <- set_frontdoor(x, "ds.vert.demo", route = "backend", K = 3L)
  expect_equal(out$frontdoor, "ds.vert.demo")
  expect_equal(out$route, "backend")
})

test_that("ds.vert.coxph cannot dispatch to a non-Cox estimand", {
  calls <- list()
  local_mocked_bindings(
    ds.vert.cox = function(formula, data = NULL, method, ...) {
      calls[[length(calls) + 1L]] <<- list(method = method, data = data)
      list(coefficients = c(x = 0.25), frontdoor = "ds.vert.cox")
    },
    .package = "dsVertClient")

  out <- ds.vert.coxph(y ~ x, data = "D")
  expect_identical(out$frontdoor, "ds.vert.coxph")
  expect_identical(calls[[1L]]$method, "profile")
  expect_error(ds.vert.coxph(y ~ x, data = "D", method = "discrete"),
               "profile")
  expect_length(calls, 1L)
})

test_that("frontdoor aliases preserve backend outputs and route metadata", {
  fit <- mock_fit()
  expect_error(ds.vert.lasso(fit, lambda_1 = 0.3, alpha_grid = 1),
               "validated ds.vertDPGaussian")
})

# =============================================================================
# ds.vertLASSO1Step legacy boundary
# =============================================================================

test_that("LASSO1Step rejects the retired quadratic-surrogate ds.glm route", {
  expect_error(
    ds.vertLASSO1Step(mock_fit(), lambda = c(0.1, 0.05)),
    "validated ds.vertDPGaussian")
})

# =============================================================================
# ds.vertLASSOCV legacy boundary
# =============================================================================
test_that("LASSOCV rejects the retired ds.glm covariance route", {
  expect_error(ds.vertLASSOCV(mock_fit(), lambda_grid = c(0.1, 0)),
               "validated ds.vertDPGaussian")
})
