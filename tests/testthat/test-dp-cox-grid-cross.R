
test_that("Cox client accepts only a literal additive Surv formula", {
  value <- .dsvert_dp_cox_grid_cross_formula(
    Surv(site_a$time, site_a$event) ~ site_a$x + site_b$z)
  expect_identical(value$predictors, c("site_a$x", "site_b$z"))
  for (formula in list(y ~ x, Surv(t, e) ~ log(x), Surv(t, e) ~ x * z,
      Surv(t, e) ~ x + strata(z), Surv(t, e) ~ x - z,
      Surv(start, t, e) ~ x, Surv(t, e) ~ x + x)) {
    expect_error(.dsvert_dp_cox_grid_cross_formula(formula),
                 class = "dsvert_dp_public_failure")
  }
})

test_that("Cox cross-owner contracts require both real custodian signatures", {
  fixture <- .cox_cross_client_fixture()
  validate <- function(value) .dsvert_dp_cox_grid_cross_contract_validate(
    value, fixture$policy, fixture$schema_manifest)
  expect_identical(validate(fixture$contract),
                   .dsvert_joint_dp_client_canonical(fixture$contract))
  for (peer in names(fixture$contract$signatures)) {
    changed <- fixture$contract
    changed$signatures[[peer]] <- NULL
    expect_error(validate(changed), class = "dsvert_dp_public_failure")
    changed <- fixture$contract
    changed$signatures[[peer]] <- strrep("A", 86)
    expect_error(validate(changed), class = "dsvert_dp_public_failure")
  }
  for (change in list(
      function(x) { x$spec$ties <- "efron"; x },
      function(x) { x$spec$intercept <- TRUE; x },
      function(x) { x$spec$time_semantics <- "binned"; x },
      function(x) { x$spec$numeric_contract <- list(); x },
      function(x) { x$spec$sensitivity$raw_l1_sensitivity <- 0; x },
      function(x) { x$artifact$implementation_state <- "same_owner_materialized"; x },
      function(x) { x$source_contract$purpose <- "generic_nonlinear_rpc"; x })) {
    expect_error(validate(fixture$sign(change(fixture$contract))),
                 class = "dsvert_dp_public_failure")
  }
})

test_that("Cox client canonical argmin rescales slopes and omits inference", {
  fixture <- .cox_cross_client_fixture()
  spec <- fixture$contract$spec
  coordinates <- rep(0, length(spec$beta_grid))
  moment <- .dsvert_dp_cox_grid_cross_moment(coordinates, spec)
  expect_identical(moment$selected_candidate, 1L)
  expect_equal(unname(moment$coefficients), unlist(spec$beta_grid[[1L]]))
  for (bad in list(-coordinates - 1, coordinates + 0.5,
                   rep(Inf, length(coordinates)), coordinates[-1L],
                   unlist(spec$sensitivity$maximum_coordinates) + 1)) {
    expect_error(.dsvert_dp_cox_grid_cross_moment(bad, spec),
                 class = "dsvert_dp_public_failure")
  }
  release <- function(...) c(fixture[c("contract", "policy", "schema_manifest")],
                             list(coordinates = coordinates))
  result <- .dsvert_dp_cox_grid_cross_impl(
    Surv(site_a$time, site_a$event) ~ site_a$x + site_b$z,
    "aligned", "cox_grid", list(site_a = NULL, site_b = NULL), release)
  expect_identical(result$selected_candidate, 1L)
  expect_null(result$std_errors)
  expect_null(result$baseline_hazard)
  expect_false(result$production_ready)
  expect_false(result$protected_optimizer_called)
  expect_error(.dsvert_dp_cox_grid_cross_impl(
    Surv(site_a$time, site_a$event) ~ site_a$x + site_b$wrong,
    "aligned", "cox_grid", list(site_a = NULL, site_b = NULL), release),
    class = "dsvert_dp_public_failure")
})

test_that("production Cox client is closed and exposes no test evaluator", {
  expect_identical(names(formals(dp_cox_grid)),
    c("formula", "data", "analysis_id", "datasources"))
  expect_error(dp_cox_grid(Surv(site_a$time, site_a$event) ~ site_a$x + site_b$z,
    "aligned", "cox_grid", list(site_a = NULL, site_b = NULL)),
    class = "dsvert_dp_public_failure")
  expect_false(.dsvert_dp_cox_grid_cross_client_register()$runtime_enabled)
  expect_identical(.dsvert_dp_cox_grid_cross_client_register()$entry, dp_cox_grid)
  expect_false("dp_cox_grid" %in% getNamespaceExports("dsVertClient"))
})

test_that("Cox server and client canonical contracts agree", {
  skip_if_not_installed("dsVert")
  server <- asNamespace("dsVert")
  if (!exists(".dsvert_dp_cox_grid_cross_register", server, inherits=FALSE)) {
    skip("companion server lacks additive Cox registration")
  }
  f <- .cox_cross_client_fixture()
  validated <- get(".dsvert_dp_cox_grid_cross_contract_validate", server)(
    f$contract, f$policy, f$schema_manifest)
  expect_identical(.dsvert_joint_dp_client_json(validated),
    .dsvert_joint_dp_client_json(.dsvert_dp_cox_grid_cross_contract_validate(
      f$contract,f$policy,f$schema_manifest)))
})


test_that("Cox measured admission is enforced even with both valid signatures", {
  admission <- .dsvert_dp_cox_grid_cross_admission(2, 1)
  n <- admission$maximum_rows; j <- admission$maximum_candidates
  grid <- lapply(seq_len(j), function(i) c(i / 100, 0))
  f <- .cox_cross_client_fixture(capacity = n, beta_grid = grid)
  validate <- function(value) .dsvert_dp_cox_grid_cross_contract_validate(
    value, f$policy, f$schema_manifest)
  expect_silent(validate(f$contract))
  expect_identical(f$contract$spec$resource_admission, admission)
  expect_error(.cox_cross_client_fixture(capacity = n + 1), class = "dsvert_dp_public_failure")
  expect_error(.cox_cross_client_fixture(beta_grid = c(grid, list(c(1, 0)))),
    class = "dsvert_dp_public_failure")
  for (field in names(admission)) {
    changed <- f$contract
    changed$spec$resource_admission[[field]] <- NULL
    expect_error(validate(f$sign(changed)), class = "dsvert_dp_public_failure")
  }
  changed <- f$contract
  changed$spec$resource_admission$maximum_rows <- n + 1
  expect_error(validate(f$sign(changed)), class = "dsvert_dp_public_failure")
})


test_that("Cox signed layout pins the terminal share conversion for fusion", {
  f <- .cox_cross_client_fixture()
  layout <- f$contract$source_contract$private_layout
  expect_equal(layout$kernel_output_ring_bits, 64)
  expect_equal(layout$joint_dp_source_ring_bits, 128)
  expect_identical(layout$output_conversion,
    "mod64_reconstruct_cast_remask_in_authenticated_fusion_v1")
  changed <- f$contract
  changed$source_contract$private_layout$kernel_output_ring_bits <- 128
  expect_error(.dsvert_dp_cox_grid_cross_contract_validate(
    f$sign(changed), f$policy, f$schema_manifest), class = "dsvert_dp_public_failure")
})

test_that("Cox K-owner source contracts retain exactly two computation authorities", {
  for (owners in c(2L, 3L, 5L)) {
    f <- .cox_cross_client_fixture(owners = owners)
    validate <- function(value) .dsvert_dp_cox_grid_cross_contract_validate(
      value, f$policy, f$schema_manifest)
    spec <- validate(f$contract)$spec
    expect_length(spec$participating_peers, owners)
    expect_identical(unlist(spec$computation_peers, use.names = FALSE), c("site_a", "site_b"))
    source <- f$contract$source_contract
    expect_identical(source$source_peers, spec$participating_peers)
    expect_identical(source$recipients, spec$computation_peers)
    layout <- .dsvert_dp_cox_grid_cross_layout(spec)
    expect_identical(layout$version, if (owners == 2L)
      "dsvert-cox-cross-private-source-layout-v1" else "dsvert-cox-cross-private-source-layout-v2")
    expect_identical(layout$partial_predictors, if (owners == 2L)
      "two_owners_exact_f100_additive_ring128_v1" else "all_source_owners_exact_f100_additive_ring128_v1")
    expect_identical(layout$release_prefix_source_rule,
      "all_zero_until_authenticated_result_injection_v1")
    expect_identical(spec$numeric_contract$eta_rounding, "one_rne_after_owner_sum_v1")
    for (peer in names(f$keys)) {
      changed <- f$contract; changed$signatures[[peer]] <- NULL
      expect_error(validate(changed), class = "dsvert_dp_public_failure")
    }
    for (mutate in list(
      function(x) { x$spec$participating_peers <- x$spec$participating_peers[-1]; x },
      function(x) { x$spec$computation_peers <- rev(x$spec$computation_peers); x },
      function(x) { x$source_contract$source_peers <- x$source_contract$source_peers[-1]; x },
      function(x) { x$source_contract$recipients <- x$source_contract$recipients[1]; x })) {
      expect_error(validate(f$sign(mutate(f$contract))), class = "dsvert_dp_public_failure")
    }
    for (compute in list("site_a", c("site_a", "site_a"), c("site_a", "absent"),
                         c("site_a", "site_b", "site_c"))) {
      policy <- f$policy; policy$designated_noise_peers <- compute
      expect_error(.dsvert_dp_cox_grid_cross_spec(f$raw, policy, f$schema),
        class = "dsvert_dp_public_failure")
    }
    if (owners > 2L) {
      changed <- f$contract
      changed$source_contract$private_layout$partial_predictors <-
        "two_owners_exact_f100_additive_ring128_v1"
      expect_error(validate(f$sign(changed)), class = "dsvert_dp_public_failure")
      raw <- f$raw; raw$predictor_order <- raw$predictor_order[-length(raw$predictor_order)]
      raw$beta_grid <- lapply(raw$beta_grid, function(beta) beta[-length(beta)])
      expect_error(.dsvert_dp_cox_grid_cross_spec(raw, f$policy, f$schema),
        class = "dsvert_dp_public_failure")
    }
  }
})
