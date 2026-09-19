test_that("both categorical cross-owner validators require both authentic signatures", {
  for (family in c("multinomial", "ordinal")) {
    fixture <- .categorical_cross_client_fixture(family)
    registration <- fixture$registration
    validate <- function(contract = fixture$contract,
                         schema = fixture$schema_manifest) {
      registration$validate(contract, fixture$policy, schema)
    }
    contract <- validate()
    expect_identical(contract, .dsvert_joint_dp_client_canonical(fixture$contract))
    expect_identical(contract$artifact$implementation_state,
                     "cross_owner_exact_gc_materialized")
    expect_identical(contract$artifact$cross_owner_state,
                     "exact_gc_to_joint_dp_vector_v1")
    expect_false(registration$release_enabled)
    wire <- jsonlite::fromJSON(.dsvert_joint_dp_client_json(contract),
                              simplifyVector = FALSE)
    expect_identical(validate(wire), contract)
    schema_wire <- jsonlite::fromJSON(
      .dsvert_joint_dp_client_json(fixture$schema_manifest), simplifyVector = FALSE)
    expect_identical(validate(wire, schema_wire), contract)
    for (peer in names(fixture$policy$peer_pinset)) {
      bad <- fixture$contract
      bad$signatures[[peer]] <- NULL
      .categorical_cross_client_reject(function() validate(bad))
      bad <- fixture$contract
      bad$signatures[[peer]] <- strrep("A", 86L)
      .categorical_cross_client_reject(function() validate(bad))
      schema <- fixture$schema_manifest
      schema$signatures[[peer]] <- NULL
      .categorical_cross_client_reject(function() validate(schema = schema))
    }
  }
})

test_that("categorical contracts reject signed tampering in every binding", {
  mutations <- list(
    function(x) { x$spec$predictor_order <- rev(x$spec$predictor_order); x },
    function(x) { x$spec$predictors[["site_b$z"]]$upper <- 100; x },
    function(x) { x$spec$beta_encoded[[1L]][[1L]] <- "1"; x },
    function(x) { x$spec$computation_peers <- list("site_a", "site_a"); x },
    function(x) { x$spec$numeric_grid_bits <- 7; x },
    function(x) { x$spec$sensitivity$raw_l1_sensitivity <- 1; x },
    function(x) { x$spec$sensitivity$raw_l2_sensitivity <- 1; x },
    function(x) { x$spec$sensitivity$maximum_coordinates[[1L]] <- 1; x },
    function(x) { x$spec$alignment$public_patient_dependent_hash <- TRUE; x },
    function(x) { x$spec$numeric_contract$arithmetic_width_bits <- 128; x },
    function(x) { x$spec$numeric_contract$profile_sha256 <- strrep("0", 64); x },
    function(x) { x$spec$numeric_contract$certified_uniform_error <- "0"; x },
    function(x) { x$spec$numeric_contract$nonlinear_fraction_bits <- 64; x },
    function(x) { x$spec$numeric_contract$profile_degree <- 32; x },
    function(x) { x$artifact$implementation_state <- "same_owner_materialized"; x },
    function(x) { x$artifact$cross_owner_state <- "reserved_not_materialized"; x },
    function(x) { x$artifact$result_evidence_required <- FALSE; x },
    function(x) { x$artifact$transcript$row_batch_size <- 100; x },
    function(x) { x$artifact$unexpected <- "patient-loss"; x },
    function(x) { x$source_contract$recipients <- list("site_a"); x },
    function(x) { x$source_contract$private_layout$blocks[[1L]]$value_fraction_bits <- 0; x },
    function(x) { x$source_contract$private_layout$padding_validity <- 1; x })
  for (family in c("multinomial", "ordinal")) {
    fixture <- .categorical_cross_client_fixture(family)
    for (mutate in mutations) {
      bad <- fixture$sign(mutate(fixture$contract))
      .categorical_cross_client_reject(function() {
        fixture$registration$validate(
          bad, fixture$policy, fixture$schema_manifest)
      })
    }
  }
})

test_that("multinomial DP postprocessing converts units and breaks ties canonically", {
  fixture <- .categorical_cross_client_fixture()
  fit <- .dsvert_dp_multinomial_grid_cross_postprocess(
    fixture$contract, c(100, 90))
  expect_identical(fit$selected_candidate, 2L)
  expect_equal(fit$coefficients, structure(c(0.5, 0.5, 0.1, -0.5, -0.5, -0.1),
    dim = c(3L, 2L), dimnames = list(c("(Intercept)", "site_a$x", "site_b$z"),
                                    c("B", "C"))))
  expect_equal(fit$selected_dp_negative_log_likelihood, 90 / 256)
  expect_identical(fit$coordinate_scale, 256)
  expect_null(fit$standard_errors)
  expect_false(fit$source_values_exposed)
  expect_identical(.dsvert_dp_multinomial_grid_cross_postprocess(
    fixture$contract, c(90, 90))$selected_candidate, 1L)
  for (losses in list(c(1), c(-1, 0), c(0, NA), c(0, Inf), c(0, 0.5),
                     c(0, 2^53), list(1, 2))) {
    .categorical_cross_client_reject(function() {
      .dsvert_dp_multinomial_grid_cross_postprocess(fixture$contract, losses)
    })
  }
})

test_that("new family frontdoors validate then fail closed without a producer", {
  for (family in c("multinomial", "ordinal")) {
    fixture <- .categorical_cross_client_fixture(family)
    frontdoor <- if (family == "multinomial") dp_multinomial_grid else dp_ordinal_grid
    formulas <- list(site_a$y ~ site_a$x + site_b$z,
      site_a$y ~ site_b$z + site_a$x, y ~ x + z,
      site_a$y ~ log(site_a$x) + site_b$z,
      site_a$y ~ site_a$x * site_b$z,
      site_a$y ~ site_a$x + site_a$x + site_b$z,
      site_a$y ~ site_a$x)
    for (formula in formulas) {
      .categorical_cross_client_reject(function() frontdoor(
        formula, "aligned", "cross_grid", fixture$contract,
        fixture$policy, fixture$schema_manifest,
        datasources = list(not_a_connection = "never_resolved")))
    }
    datasource_evaluated <- FALSE
    .categorical_cross_client_reject(function() frontdoor(
      formulas[[1L]], "aligned", "cross_grid", fixture$contract,
      fixture$policy, fixture$schema_manifest, datasources = {
        datasource_evaluated <<- TRUE
        stop("must not contact a datasource")
      }))
    expect_false(datasource_evaluated)
    expect_silent(.dsvert_dp_categorical_grid_cross_formula(
      formulas[[1L]], fixture$contract$spec))
    expect_silent(.dsvert_dp_categorical_grid_cross_formula(
      formulas[[2L]], fixture$contract$spec))
    for (formula in formulas[-c(1L, 2L)]) {
      .categorical_cross_client_reject(function() {
        .dsvert_dp_categorical_grid_cross_formula(formula, fixture$contract$spec)
      })
    }
  }
})

test_that("existing same-owner family frontdoors retain their grid dispatch", {
  seen <- list()
  capture <- function(family) function(formula, data_name, analysis_id,
      server = NULL, datasources = NULL, .aggregate) {
    seen[[family]] <<- list(formula = formula, data = data_name,
      analysis_id = analysis_id, datasources = datasources)
    family
  }
  testthat::local_mocked_bindings(
    .dsvert_dp_multinom_grid_impl = capture("multinomial"),
    .dsvert_dp_ordinal_grid_impl = capture("ordinal"),
    .dsvert_dp_multinomial_grid_cross_release = function(...) {
      stop("cross-owner release must not be dispatched")
    },
    .dsvert_dp_ordinal_grid_cross_release = function(...) {
      stop("cross-owner release must not be dispatched")
    },
    .package = "dsVertClient")
  formula <- y ~ x + z
  connections <- list(site_a = structure(list(), class = "mock_connection"))
  expect_identical(ds.vertMultinom(formula, data = "cohort",
    analysis_id = "grid", datasources = connections), "multinomial")
  expect_identical(ds.vertOrdinal(formula, data = "cohort",
    analysis_id = "grid", datasources = connections), "ordinal")
  expected <- list(formula = formula, data = "cohort", analysis_id = "grid",
                   datasources = connections)
  expect_identical(seen$multinomial, expected)
  expect_identical(seen$ordinal, expected)
})

test_that("a nonfirst multinomial reference retains signed class associations", {
  fixture <- .categorical_cross_client_fixture()
  raw <- fixture$raw
  raw$reference <- "B"
  spec <- fixture$registration$spec(raw, fixture$policy, fixture$schema)
  contract <- fixture$contract
  contract$spec <- spec
  contract$artifact <- fixture$registration$artifact(spec)
  contract$source_contract <- fixture$registration$source_contract(
    spec, contract$artifact)
  contract <- fixture$sign(contract)
  validated <- fixture$registration$validate(
    contract, fixture$policy, fixture$schema_manifest)
  fit <- fixture$registration$postprocess(validated, c(100, 90))
  expect_identical(fit$reference, "B")
  expect_identical(colnames(fit$coefficients), c("A", "C"))
  expect_equal(fit$coefficients["site_a$x", ], c(A = 0.5, C = -0.5))
})

test_that("categorical sensitivity caps retain both adjacency multipliers", {
  for (family in c("multinomial", "ordinal")) {
    for (bits in c(8L, 16L, 18L)) {
      fixture <- .categorical_cross_client_fixture(family, bits)
      spec <- fixture$contract$spec
      candidates <- if (family == "multinomial") spec$beta_grid else spec$candidate_grid
      addition <- .dsvert_dp_categorical_grid_cross_sensitivity(
        candidates, family, 3, 2, bits, 20, "add_remove_patient")
      replacement <- .dsvert_dp_categorical_grid_cross_sensitivity(
        candidates, family, 3, 2, bits, 20, "replace_one_fixed_cohort")
      caps <- vapply(addition$candidate_bounds, `[[`, numeric(1L), "per_patient_cap")
      numeric_contract <- .dsvert_dp_categorical_grid_cross_numeric(family)
      expect_identical(numeric_contract$nonlinear_fraction_bits, 16)
      expect_identical(numeric_contract$arithmetic_width_bits, 32)
      expect_identical(numeric_contract$profile_pieces, 64)
      expect_identical(numeric_contract$profile_degree, 2)
      for (bound in addition$candidate_bounds) {
        expect_equal(bound$profile_error_bound,
                     as.numeric(numeric_contract$certified_uniform_error))
        expect_equal(bound$output_rounding_error_bound,
                     if (bits < 16) 2^(-bits - 1) else 0)
        expect_gte(bound$per_patient_cap / 2^bits,
          bound$exact_loss_bound + 2 * bound$profile_error_bound +
            bound$output_rounding_error_bound)
      }
      expect_equal(unlist(addition$maximum_coordinates), 20 * caps)
      expect_equal(addition$raw_l1_sensitivity, sum(caps))
      expect_gte(addition$raw_l2_sensitivity, sqrt(sum(caps^2)))
      expect_equal(replacement$raw_l1_sensitivity, 2 * addition$raw_l1_sensitivity)
      expect_equal(replacement$raw_l2_sensitivity, 2 * addition$raw_l2_sensitivity)
      expect_equal(addition$natural_l1_sensitivity,
                   addition$raw_l1_sensitivity / 2^bits)
      extreme <- if (family == "multinomial") list(c(8, 8, 0, 8, 8, 0)) else
        list(list(beta = c(0, 8, 0), thresholds = c(-8, 8)))
      .categorical_cross_client_reject(function() {
        .dsvert_dp_categorical_grid_cross_sensitivity(
          extreme, family, 3, 2, 18, 2^31 - 1, "add_remove_patient")
      })
    }
  }
})
