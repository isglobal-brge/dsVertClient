# These tests keep real signed-contract/schema validation. Only authenticated
# discovery and the already-tested protected release boundary are replaced.
.cross_frontdoor_trusted_fixture <- function(contract, policy, schema) {
  artifact <- .dsvert_dp_glm_grid_cross_workload_artifact(contract)
  list(manifest = list(workload = list(families = list(gaussian_models = list(
    artifacts = stats::setNames(list(artifact), contract$spec$analysis_id))))),
    policy = policy, schema = schema,
    datasources = list(peer = structure(list(), class = "mock_connection")))
}

test_that("all four promoted grid entries and the safe NB alias are exported", {
  for (name in c("dp_nb_grid", "dp_lasso_grid", "dp_multinomial_grid",
                 "dp_ordinal_grid", "ds.vertNB")) {
    expect_true(name %in% getNamespaceExports("dsVertClient"))
    expect_true(is.function(getExportedValue("dsVertClient", name)))
  }
  expect_identical(formals(ds.vertNB), formals(ds.vertNBFullRegTheta))
  expect_identical(body(ds.vertNB), body(ds.vertNBFullRegTheta))
})

test_that("K2 K3 K5 NB and categorical frontdoors use their signed cross release", {
  for (family in c("nb", "multinomial", "ordinal")) for (k in c(2L, 3L, 5L)) {
    f <- if (family == "nb") .nb_grid_cross_fixture("dsVertClient", peer_count = k) else
      .categorical_cross_client_fixture(family, peer_count = k)
    schema <- if (family == "nb") f$schema else f$schema_manifest
    trusted <- .cross_frontdoor_trusted_fixture(f$contract, f$policy, schema)
    spec <- f$contract$spec
    formula <- stats::as.formula(paste(spec$outcome$reference, "~",
      paste(unlist(spec$predictor_order), collapse = "+")))
    seen <- 0L
    run <- function(fun, ...) testthat::with_mocked_bindings(
      fun(formula, data = spec$dataset, analysis_id = spec$analysis_id, ...),
      .dsvert_dp_grid_cross_bootstrap = function(...) trusted,
      .dsvert_dp_likelihood_grid_cross_release = function(contract, policy, schema, ...) {
        seen <<- seen + 1L
        expect_true(.dsvert_dp_glm_grid_cross_equal(contract, f$contract))
        expect_identical(policy, f$policy)
        list(selected_candidate = 1L, certificate = "boundary-sentinel")
      }, .package = "dsVertClient")
    front <- switch(family, nb = ds.vertNBFullRegTheta,
      multinomial = ds.vertMultinom, ordinal = ds.vertOrdinal)
    expect_identical(run(front), list(selected_candidate = 1L, certificate = "boundary-sentinel"))
    if (family == "nb") {
      expect_identical(run(ds.vertNB), run(front))
      expect_identical(run(ds.vert.nb)$selected_candidate, run(front)$selected_candidate)
    }
    if (family == "multinomial") {
      expect_identical(run(front, classes = unlist(spec$class_order), reference = spec$reference), run(front))
    }
    if (family == "ordinal") expect_identical(run(front,
      levels_ordered = unlist(spec$class_order)), run(front))
    expect_gte(seen, 1L)
  }
})

test_that("cross frontdoor formula and semantic controls fail before release", {
  calls <- 0L
  f <- .categorical_cross_client_fixture("multinomial")
  trusted <- .cross_frontdoor_trusted_fixture(f$contract, f$policy, f$schema_manifest)
  testthat::local_mocked_bindings(
    .dsvert_dp_grid_cross_bootstrap = function(...) trusted,
    .dsvert_dp_likelihood_grid_cross_release = function(...) {
      calls <<- calls + 1L
      stop("private release must not occur")
    }, .package = "dsVertClient")
  for (formula in list(site_a$y ~ site_a$x + site_b$z + site_b$z,
                       site_a$y ~ x + site_b$z,
                       site_a$y ~ log(site_a$x) + site_b$z,
                       site_a$y ~ site_a$x * site_b$z,
                       site_a$y ~ site_a$x + site_b$z - 1,
                       site_a$y ~ site_a$x)) {
    expect_error(ds.vertMultinom(formula, "aligned", analysis_id = "cross_grid"))
  }
  formula <- site_a$y ~ site_a$x + site_b$z
  expect_error(ds.vertMultinom(formula, "wrong", analysis_id = "cross_grid"))
  expect_error(ds.vertMultinom(formula, "aligned", analysis_id = "wrong"))
  expect_error(ds.vertMultinom(formula, "aligned", analysis_id = "cross_grid", classes = c("B", "A", "C")), "exactly match")
  expect_error(ds.vertMultinom(formula, "aligned", analysis_id = "cross_grid", reference = "B"), "exactly match")
  expect_error(ds.vertMultinom(formula, "aligned", analysis_id = "cross_grid", max_iter = 2), "legacy controls")
  expect_error(ds.vertOrdinal(formula, "aligned", analysis_id = "cross_grid"))
  # A malformed embedded contract cannot be repaired by carrying a valid schema.
  broken <- f$contract
  broken$signatures[[1L]] <- NULL
  trusted$manifest$workload$families$gaussian_models$artifacts$cross_grid$signed_contract <-
    .dsvert_joint_dp_client_json(broken)
  expect_error(ds.vertMultinom(formula, "aligned", analysis_id = "cross_grid"),
    class = "dsvert_dp_public_failure")
  expect_identical(calls, 0L)
})

test_that("cross LASSO explicit extensions validate K2 K3 K5 and reuse one base", {
  for (family in c("binomial", "poisson")) for (k in c(2L, 3L, 5L)) {
    f <- .lasso_cross_fixture(family, peer_count = k, execution_profile = "piecewise_v2")
    trusted <- .cross_frontdoor_trusted_fixture(f$base_contract, f$policy, f$schema)
    formula <- stats::as.formula(paste(f$base_contract$spec$outcome$reference, "~",
      paste(unlist(f$base_contract$spec$predictor_order), collapse = "+")))
    releases <- 0L
    result <- testthat::with_mocked_bindings(ds.vertLASSOIter(formula, "cohort",
      family = family, analysis_id = f$spec$analysis_id, cross_contract = f$contract,
      base_contract = f$base_contract, lambda = 0),
      .dsvert_dp_grid_cross_bootstrap = function(...) trusted,
      dp_lasso_grid = function(covariates, signed_contract, base_contract, policy, schema_manifest, ...) {
        releases <<- releases + 1L
        expect_true(.dsvert_dp_glm_grid_cross_equal(signed_contract, f$contract))
        expect_true(.dsvert_dp_glm_grid_cross_equal(base_contract, f$base_contract))
        out <- .dsvert_dp_lasso_cross_postprocess(rep(100, f$base$coordinate_count), f$spec)
        out$paths <- out$normalized_paths
        out
      }, .package = "dsVertClient")
    expect_identical(releases, 1L)
    expect_equal(result$lambda, 0)
    expect_length(result$paths, 1L)
    expect_output(print(result), "signed finite candidates")
    expect_s3_class(result, "ds.vertDPLASSOPath")
    expect_equal(result$additional_privacy_cost, c(epsilon = 0, delta = 0))
    expect_identical(result$additional_server_calls_after_synopsis, 0L)
  }
})

test_that("LASSO contract, formula, family and penalty tampering never reaches release", {
  f <- .lasso_cross_fixture("poisson", execution_profile = "piecewise_v2")
  trusted <- .cross_frontdoor_trusted_fixture(f$base_contract, f$policy, f$schema)
  releases <- 0L
  testthat::local_mocked_bindings(.dsvert_dp_grid_cross_bootstrap = function(...) trusted,
    dp_lasso_grid = function(...) { releases <<- releases + 1L; stop("unexpected release") },
    .package = "dsVertClient")
  args <- list(formula = peer_a$y ~ peer_a$x + peer_b$z, data = "cohort", family = "poisson",
    analysis_id = f$spec$analysis_id, cross_contract = f$contract, base_contract = f$base_contract)
  for (change in list(list(cross_contract = NULL), list(base_contract = NULL),
      list(family = "gaussian"), list(family = "binomial"), list(data = "other"),
      list(analysis_id = "other"), list(lambda = 0.13), list(lambda = c(0, 0)),
      list(formula = peer_a$y ~ peer_a$x), list(max_outer = 3L),
      list(dp_analysis_id = "other"), list(cor_analysis_id = "other"))) {
    changed <- args
    changed[names(change)] <- change
    expect_error(do.call(ds.vertLASSOIter, changed))
  }
  for (key in c("cross_contract", "base_contract")) {
    changed <- args
    changed[[key]]$signatures[[1L]] <- NULL
    expect_error(do.call(ds.vertLASSOIter, changed), class = "dsvert_dp_public_failure")
  }
  changed <- args
  changed$cross_contract$spec$candidate_grid[[1L]]$lambda <- 100
  expect_error(do.call(ds.vertLASSOIter, changed), class = "dsvert_dp_public_failure")
  expect_identical(releases, 0L)
})

test_that("same-owner LASSO dispatch and supplied-argument NB semantics stay unchanged", {
  seen <- list()
  testthat::local_mocked_bindings(
    .dsvert_dp_lasso_grid_impl = function(formula, data_name, analysis_id,
        family, lambda, datasources, .aggregate) {
      list(paths = list(c(`(Intercept)` = 0, x = 1)), family = family,
           analysis_id = analysis_id, lambda = lambda)
    },
    .dsvert_dp_grid_cross_bootstrap = function(...) stop("unexpected cross discovery"),
    .dsvert_formal_nb2_frequency_adapter = function(explicit_arguments, ...) {
      seen[[length(seen) + 1L]] <<- explicit_arguments
      explicit_arguments
    }, .package = "dsVertClient")
  for (family in c("binomial", "poisson")) {
    direct <- ds.vertLASSOIter(y ~ x, "cohort", family = family,
      lambda = 0.1, analysis_id = "same-owner")
    alias <- ds.vert.lasso_iter(y ~ x, "cohort", family = family,
      lambda = 0.1, analysis_id = "same-owner")
    expect_identical(direct$paths, alias$paths)
    expect_identical(direct, ds.vertLASSOIter(y ~ x, "cohort", family = family,
      lambda = 0.1, analysis_id = "same-owner", cross_contract = NULL,
      base_contract = NULL))
    expect_identical(direct$called_via, "ds.vertLASSOIter_analysis_id")
  }
  expect_identical(ds.vertNB(y ~ 1, frequency = list()),
                   ds.vertNBFullRegTheta(y ~ 1, frequency = list()))
  expect_identical(ds.vertNB(y ~ 1, frequency = list(), theta = NULL),
                   ds.vertNBFullRegTheta(y ~ 1, frequency = list(), theta = NULL))
  expect_identical(seen[[1L]], c("formula", "frequency"))
  expect_identical(seen[[3L]], c("formula", "theta", "frequency"))
})

test_that("dotted LASSO forwards explicit contracts without injecting retired controls", {
  f <- .lasso_cross_fixture("binomial", execution_profile = "piecewise_v2")
  trusted <- .cross_frontdoor_trusted_fixture(f$base_contract, f$policy, f$schema)
  seen <- 0L
  testthat::local_mocked_bindings(
    .dsvert_dp_grid_cross_bootstrap = function(...) trusted,
    dp_lasso_grid = function(...) {
      seen <<- seen + 1L
      out <- .dsvert_dp_lasso_cross_postprocess(rep(100, f$base$coordinate_count), f$spec)
      out$paths <- out$normalized_paths
      out
    }, .package = "dsVertClient")
  result <- ds.vert.lasso_iter(peer_a$y ~ peer_a$x + peer_b$z,
    data = "cohort", family = "binomial", analysis_id = f$spec$analysis_id,
    cross_contract = f$contract, base_contract = f$base_contract)
  expect_identical(seen, 1L)
  expect_identical(result$called_via, "ds.vertLASSOIter_cross_contract")
  expect_length(result$paths, 2L)
})
