.dsvert_retired_frontdoor_cases <- list(
  legacy_glm = c("ds.vertGLM", "ds.vert.glm"),
  cox = c(
    "ds.vertCox", "ds.vertCoxDiscreteNonDisclosive",
    "ds.vertCoxProfileNonDisclosive", "ds.vert.cox", "ds.vert.coxph"),
  negative_binomial = c("ds.vertNBFullRegTheta", "ds.vert.nb"),
  multinomial = c(
    "ds.vertMultinom", "ds.vertMultinomJoint",
    "ds.vertMultinomJointNewton", "ds.vert.multinom"),
  ordinal = c(
    "ds.vertOrdinal", "ds.vertOrdinalJointNewton", "ds.vert.ordinal"),
  lmm = c("ds.vertLMM", "ds.vertLMM.k3", "ds.vert.lmm"),
  gee = c("ds.vertGEE", "ds.vert.gee"),
  glmm = c("ds.vertGLMM", "ds.vert.glmm"),
  ipw = c("ds.vertIPW", "ds.vert.ipw"),
  mi = c("ds.vertMI", "ds.vert.mi")
)

test_that("every retired analytic frontdoor fails locally before DSI", {
  remote_calls <- new.env(parent = emptyenv())
  remote_calls$n <- 0L
  reached_remote <- function(...) {
    remote_calls$n <- remote_calls$n + 1L
    stop("test observed a remote operation", call. = FALSE)
  }
  testthat::local_mocked_bindings(
    .dsvert_quarantine_test_mode = function() FALSE,
    .dsvert_aggregate_strict = reached_remote,
    .dsvert_fanout_by_site = reached_remote,
    .dsvert_assign_strict = reached_remote,
    .dsvert_cleanup_best_effort = reached_remote,
    .dsvert_datasources = reached_remote,
    .package = "dsVertClient")

  contracts <- get(
    ".DSVERT_RETIRED_REMOTE_ROUTES", envir = asNamespace("dsVertClient"),
    inherits = FALSE)
  for (route in names(.dsvert_retired_frontdoor_cases)) {
    contract <- contracts[[route]]
    for (frontdoor in .dsvert_retired_frontdoor_cases[[route]]) {
      condition <- tryCatch(
        do.call(getExportedValue("dsVertClient", frontdoor), list()),
        dsvert_route_unavailable = identity)
      expect_s3_class(condition, "dsvert_route_unavailable")
      expect_identical(condition$code, "dsvert_route_unavailable")
      expect_identical(condition$method, contract$method)
      expect_identical(condition$state, contract$state)
      expect_identical(condition$replacement, contract$replacement)
      expect_match(
        condition$message, "^[[]dsvert_route_unavailable:v1[]] ")
      expect_match(condition$message, "unavailable before DSI", fixed = TRUE)
      expect_identical(remote_calls$n, 0L)
    }
  }
})

test_that("the quarantined Cox frontdoor has no unreachable exact-profile fallback", {
  body_text <- paste(deparse(body(ds.vertCox)), collapse = "\n")
  expect_match(
    body_text,
    '.dsvert_block_retired_remote_route("cox")',
    fixed = TRUE)
  expect_false(grepl("ds.vertCoxProfileNonDisclosive", body_text, fixed = TRUE))
})

test_that("unregistered internal routes are blocked before DSI", {
  testthat::local_mocked_bindings(
    .dsvert_quarantine_test_mode = function() FALSE,
    .dsvert_aggregate_strict = function(...) {
      stop("test observed a remote operation", call. = FALSE)
    },
    .package = "dsVertClient")
  routes <- list(
    legacy_joint_dp_capsule = ".dsvert_joint_dp_capsule_status_impl",
    formal_finalizer_handoff = ".dsvert_relay_formal_finalizer_handoff_v1",
    formal_glm_control = ".dsvert_relay_formal_glm_control_v1")
  for (route in names(routes)) {
    impl <- get(routes[[route]], envir = asNamespace("dsVertClient"),
                inherits = FALSE)
    condition <- tryCatch(do.call(impl, list()),
                          dsvert_route_unavailable = identity)
    expect_s3_class(condition, "dsvert_route_unavailable")
    expected_state <- if (identical(route, "legacy_joint_dp_capsule")) {
      "lifetime_admission_route_removed"
    } else {
      "unregistered_source_route_removed"
    }
    expect_identical(condition$state, expected_state, info = route)
  }
})

test_that("the explicit Gaussian DP adapter remains reachable", {
  observed <- NULL
  testthat::local_mocked_bindings(
    .dsvert_quarantine_test_mode = function() FALSE,
    .dsvert_dp_gaussian_glm_adapter = function(...) {
      observed <<- list(...)
      structure(list(ok = TRUE), class = "dp_adapter_sentinel")
    },
    .package = "dsVertClient")

  value <- ds.vertGLM(y ~ x, dp_analysis_id = "analysis-1")
  expect_s3_class(value, "dp_adapter_sentinel")
  expect_identical(observed$analysis_id, "analysis-1")
  expect_identical(observed$family, "gaussian")
})

test_that("generic exact-GLM helpers are production-blocked before DSI", {
  testthat::local_mocked_bindings(
    .dsvert_quarantine_test_mode = function() FALSE,
    .package = "dsVertClient")

  softplus <- get(
    ".ring127_glm_softplus_manifest_mul",
    envir = asNamespace("dsVertClient"), inherits = FALSE)
  expect_error(
    do.call(softplus, list()), class = "dsvert_route_unavailable")

  vecmul <- get(
    ".dsvert_exact_gc_vecmul_run",
    envir = asNamespace("dsVertClient"), inherits = FALSE)
  text <- paste(deparse(body(vecmul)), collapse = "\n")
  branch <- regexpr("if (is.null(input_manifests))", text, fixed = TRUE)[[1L]]
  guard <- regexpr(
    '.dsvert_block_retired_remote_route("legacy_glm")', text,
    fixed = TRUE)[[1L]]
  remote <- regexpr('name = "exactGCVecmulBindInputsDS"', text,
                    fixed = TRUE)[[1L]]
  expect_gt(branch, 0L)
  expect_gt(guard, branch)
  expect_gt(remote, guard)
})
