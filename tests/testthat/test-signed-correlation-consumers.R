.signed_consumer_conns <- function(k = 2L) {
  peers <- paste0("site_", letters[seq_len(k)])
  stats::setNames(lapply(peers, function(peer) {
    structure(list(peer = peer), class = "fake")
  }), peers)
}

.signed_consumer_lasso_fit <- function() {
  structure(list(
    coefficients = c("(Intercept)" = 0, x = 0.2, z = -0.1),
    x_means = c(x = 0, z = 0), x_sds = c(x = 1, z = 1),
    x_vars = list(site_a = "x", site_c = "z"),
    n_obs = 100, family = "binomial"),
    class = c("ds.glm", "list"))
}

test_that("binomial LASSO requires a signed correlation id before DSI", {
  conns <- .signed_consumer_conns(3L)
  glm_calls <- cor_calls <- 0L
  testthat::local_mocked_bindings(
    ds.vertGLM = function(...) {
      glm_calls <<- glm_calls + 1L
      stop("legacy GLM must not start", call. = FALSE)
    },
    ds.vertCor = function(...) {
      cor_calls <<- cor_calls + 1L
      stop("correlation must not start", call. = FALSE)
    },
    .package = "dsVertClient")

  condition <- tryCatch(
    ds.vertLASSOIter(
      y ~ x + z, data = "cohort", family = "binomial",
      lambda = 0.1, max_outer = 1L, inner_iter = 0L,
      verbose = FALSE, datasources = conns),
    dsvert_signed_analysis_required = function(condition) condition)
  expect_s3_class(condition, "dsvert_signed_analysis_required")
  expect_identical(condition$argument, "cor_analysis_id")
  expect_identical(condition$required_artifact_family, "gaussian_models")
  expect_match(condition$message, "cor_analysis_id")
  expect_identical(glm_calls, 0L)
  expect_identical(cor_calls, 0L)
})

test_that("binomial LASSO rejects clipped/raw design mismatch before DSI", {
  conns <- .signed_consumer_conns(3L)
  calls <- 0L
  testthat::local_mocked_bindings(
    ds.vertGLM = function(...) {
      calls <<- calls + 1L
      stop("raw score must not start", call. = FALSE)
    },
    ds.vertCor = function(...) {
      calls <<- calls + 1L
      stop("Gaussian capsule must not be misbound", call. = FALSE)
    },
    ds.vertDPCor = function(...) {
      calls <<- calls + 1L
      stop("pairwise fallback called", call. = FALSE)
    },
    .package = "dsVertClient")

  condition <- tryCatch(
    ds.vertLASSOIter(
      y ~ x + z, data = "cohort", family = "binomial",
      lambda = 0.1, max_outer = 1L, inner_iter = 0L,
      cor_analysis_id = "signed_lasso_gram", verbose = FALSE,
      datasources = conns),
    dsvert_signed_workload_unavailable = function(condition) condition)
  expect_s3_class(condition, "dsvert_signed_workload_unavailable")
  expect_identical(
    condition$required_artifact_family, "binomial_lasso_design_grams")
  expect_identical(
    condition$reason, "signed_clipped_design_not_bound_to_raw_score_mpc")
  expect_match(condition$message, "same clipping, normalization")
  expect_identical(calls, 0L)
})

test_that("non-correlation LASSO variants do not invent an artifact id", {
  conns <- .signed_consumer_conns(2L)
  fit <- .signed_consumer_lasso_fit()
  cor_calls <- 0L
  testthat::local_mocked_bindings(
    ds.vertGLM = function(...) fit,
    ds.vertLASSO1Step = function(...) list(
      paths = list(c("(Intercept)" = 0, x = 0, z = 0)),
      objective = 0),
    ds.vertCor = function(...) {
      cor_calls <<- cor_calls + 1L
      stop("unexpected correlation", call. = FALSE)
    },
    .package = "dsVertClient")
  result <- ds.vertLASSOIter(
    y ~ x + z, data = "cohort", family = "binomial",
    exact_non_gaussian = FALSE, lambda = 0.1, max_outer = 1L,
    inner_iter = 0L, verbose = FALSE, datasources = conns)
  expect_identical(cor_calls, 0L)
  expect_identical(result$cor_analysis_id, NULL)
  expect_identical(result$correlation_calls_per_fit, 0L)
})

test_that("multinomial signed-design gap fails before DSI for K2 and K3", {
  for (k in 2:3) {
    conns <- .signed_consumer_conns(k)
    calls <- 0L
    testthat::local_mocked_bindings(
      .ds_vertMultinomWarm = function(...) {
        calls <<- calls + 1L
        stop("warm start must not run", call. = FALSE)
      },
      .dsvert_aggregate_strict = function(...) {
        calls <<- calls + 1L
        stop("DSI must not run", call. = FALSE)
      },
      ds.vertCor = function(...) {
        calls <<- calls + 1L
        stop("Gaussian correlation is not a multinomial Gram contract",
             call. = FALSE)
      },
      .package = "dsVertClient")
    missing_condition <- tryCatch(
      ds.vertMultinomJointNewton(
        outcome ~ x + z, data = "cohort",
        levels = c("ref", "a", "b"),
        max_outer = 1L, verbose = FALSE, datasources = conns),
      dsvert_signed_analysis_required = function(condition) condition)
    expect_s3_class(missing_condition, "dsvert_signed_analysis_required")
    expect_identical(missing_condition$argument, "design_analysis_id")
    expect_identical(
      missing_condition$required_artifact_family,
      "multinomial_design_grams")
    expect_identical(calls, 0L)

    condition <- tryCatch(
      ds.vertMultinomJointNewton(
        outcome ~ x + z, data = "cohort",
        levels = c("ref", "a", "b"),
        design_analysis_id = "signed_multinomial_design",
        max_outer = 1L, verbose = FALSE, datasources = conns),
      dsvert_signed_workload_unavailable = function(condition) condition)
    expect_s3_class(condition, "dsvert_signed_workload_unavailable")
    expect_identical(condition$argument, "design_analysis_id")
    expect_identical(
      condition$required_artifact_family, "multinomial_design_grams")
    expect_match(condition$message, "raw design used by the score MPC")
    expect_identical(calls, 0L)
  }
})

test_that("signed design ids are explicit public formals", {
  expect_true("cor_analysis_id" %in% names(formals(ds.vertLASSOIter)))
  expect_gt(match("cor_analysis_id", names(formals(ds.vertLASSOIter))),
            match("datasources", names(formals(ds.vertLASSOIter))))
  for (fun in list(ds.vertMultinom, ds.vertMultinomJoint,
                   ds.vertMultinomJointNewton)) {
    expect_true("design_analysis_id" %in% names(formals(fun)))
    expect_gt(match("design_analysis_id", names(formals(fun))),
              match("datasources", names(formals(fun))))
  }
})

test_that("multinomial wrappers propagate design_analysis_id", {
  conns <- .signed_consumer_conns(2L)
  observed <- list()
  sentinel <- structure(list(ok = TRUE),
                        class = c("ds.vertMultinomJointNewton", "list"))
  testthat::local_mocked_bindings(
    ds.vertMultinomJointNewton = function(..., design_analysis_id) {
      observed[[length(observed) + 1L]] <<- design_analysis_id
      sentinel
    },
    .ds_gee_find_server_holding = function(...) "site_a",
    .package = "dsVertClient")
  first <- ds.vertMultinom(
    outcome ~ x, data = "cohort", classes = c("ref", "a", "b"),
    design_analysis_id = "signed_design_a", verbose = FALSE,
    datasources = conns)
  second <- ds.vertMultinomJoint(
    outcome ~ x, data = "cohort", levels = c("ref", "a", "b"),
    design_analysis_id = "signed_design_b", verbose = FALSE,
    datasources = conns)
  expect_identical(first, sentinel)
  expect_identical(second, sentinel)
  expect_identical(observed, list("signed_design_a", "signed_design_b"))
})

test_that("multinomial Gram helper has no legacy correlation or moment route", {
  globals <- codetools::findGlobals(
    .mnl_joint_xtx_over_n, merge = FALSE)$functions
  expect_length(intersect(globals, c(
    "ds.vertCor", "ds.vertDPCor", ".dsvert_aggregate_strict")), 0L)
  body_text <- paste(deparse(body(.mnl_joint_xtx_over_n)), collapse = "\n")
  for (forbidden in c(
      "dsvertLocalMomentsDS", "localCorDS", "numeric_pair_moments")) {
    expect_false(grepl(forbidden, body_text, fixed = TRUE), info = forbidden)
  }
})
