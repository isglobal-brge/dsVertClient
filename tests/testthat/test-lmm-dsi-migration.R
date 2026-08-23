.lmm_client_source <- function(name) {
  readLines(.dsvert_client_source_file(name), warn = FALSE)
}

test_that("LMM orchestrators have no direct DSI aggregate surface", {
  files <- c("ds.vertLMM.R", "ds.vertLMM.closed_form.R", "ds.vertLMM.k3.R")
  source <- paste(unlist(lapply(files, .lmm_client_source), use.names = FALSE),
                  collapse = "\n")

  expect_false(grepl("DSI::datashield.aggregate", source, fixed = TRUE))
  expect_false(grepl("DSI::datashield.assign", source, fixed = TRUE))
  expect_false(grepl("DSI::datashield.errors", source, fixed = TRUE))
  expect_true(grepl(".dsvert_lmm_aggregate_strict", source, fixed = TRUE))
  expect_true(grepl(".dsvert_fanout_by_site", source, fixed = TRUE))
  expect_true(grepl(".dsvert_cleanup_best_effort", source, fixed = TRUE))
})

test_that("LMM help distinguishes the signed moment route from quarantined K>=3", {
  main <- paste(.lmm_client_source("ds.vertLMM.R"), collapse = "\n")
  k3 <- paste(.lmm_client_source("ds.vertLMM.k3.R"), collapse = "\n")

  expect_match(main, "signed random-intercept", fixed = TRUE)
  expect_match(main, "method-of-moments", fixed = TRUE)
  expect_match(main, "slopes, ML/REML", fixed = TRUE)
  expect_match(main, "@return A \\code{ds.vertLMM}", fixed = TRUE)
  expect_match(k3, "raises a typed \\code{dsvert_route_unavailable}",
               fixed = TRUE)
  expect_match(k3, "before any DSI call", fixed = TRUE)
  expect_match(k3, "public frontdoor fails locally", fixed = TRUE)
  expect_match(k3, "@return No fitted object", fixed = TRUE)
  expect_false(grepl("individual observations are not revealed to the client",
                     main, fixed = TRUE))
  expect_false(grepl("analyst client never", k3, fixed = TRUE))
})

test_that("the K>=3 implementation accepts K=4 before data discovery", {
  discovered_K <- NULL
  local_mocked_bindings(
    .dsvert_fanout_by_site = function(datasources, ...) {
      discovered_K <<- length(datasources)
      stop("K4_REACHED_COLUMN_DISCOVERY", call. = FALSE)
    },
    .package = "dsVertClient")
  conns <- stats::setNames(rep(list(structure(
    list(), class = "mock_connection")), 4L), paste0("site", 1:4))

  expect_error(
    .ds_vertLMM_k3_impl(
      y ~ x, data = "D", cluster_col = "cluster",
      datasources = conns, verbose = FALSE),
    "K4_REACHED_COLUMN_DISCOVERY", fixed = TRUE)
  expect_identical(discovered_K, 4L)
})

test_that("K>=3 LMM diagnostics never label the route as exactly K=3", {
  fit <- list(
    n_clusters = 2L, cluster_sizes = c(4L, 5L), rho_hat = 0.2,
    sigma_b2 = 0.5, sigma2 = 1.5,
    coefficients = c(`(Intercept)` = 1, x = 0.25))
  class(fit) <- c("ds.vertLMM.k3", "list")

  printed <- capture.output(print(fit))
  expect_match(printed[[1L]], "K>=3", fixed = TRUE)
  expect_false(any(grepl("K=3", printed, fixed = TRUE)))
  source <- paste(.lmm_client_source("ds.vertLMM.k3.R"), collapse = "\n")
  expect_match(source, "requires K>=3 connections", fixed = TRUE)
  expect_false(grepl("DataSHIELD K=3 connections", source, fixed = TRUE))
})

test_that("LMM aggregate helper delegates once to the strict contract", {
  calls <- list()
  local_mocked_bindings(
    .dsvert_aggregate_strict = function(conns, expr, operation, ...) {
      calls[[length(calls) + 1L]] <<- list(
        conns = conns, expr = expr, operation = operation)
      stats::setNames(list(list(ok = TRUE)), names(conns))
    },
    .package = "dsVertClient")

  conns <- list(site_a = structure(list(), class = "mock_connection"))
  result <- .dsvert_lmm_aggregate_strict(
    conns, call(name = "dsvertClusterSizesDS", data_name = "D",
                cluster_col = "id"))

  expect_length(calls, 1L)
  expect_identical(names(calls[[1L]]$conns), "site_a")
  expect_identical(calls[[1L]]$operation, "LMM dsvertClusterSizesDS")
  expect_identical(result$site_a$ok, TRUE)
})

test_that("LMM strict path maps a real DSLite multi-site phase", {
  skip_if_not_installed("DSLite")
  method <- paste0("dsvertLMMTransportProbe", Sys.getpid())
  server <- DSLite::newDSLiteServer(tables = list(
    site_a = data.frame(x = 1), site_b = data.frame(x = 2)))
  server$aggregateMethod(method, function(ordinal) {
    list(ordinal = as.integer(ordinal), accepted = TRUE)
  })
  object <- paste0("dsvert_lmm_transport_server_", Sys.getpid())
  assign(object, server, envir = .GlobalEnv)
  on.exit(rm(list = object, envir = .GlobalEnv), add = TRUE)

  builder <- DSI::newDSLoginBuilder()
  builder$append(
    server = "site_a", url = object, table = "site_a",
    driver = "DSLiteDriver")
  builder$append(
    server = "site_b", url = object, table = "site_b",
    driver = "DSLiteDriver")
  conns <- DSI::datashield.login(builder$build(), assign = FALSE)
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  expressions <- list(
    site_a = call(name = method, ordinal = 1L),
    site_b = call(name = method, ordinal = 2L))
  result <- .dsvert_lmm_aggregate_strict(conns, expressions)

  expect_named(result, c("site_a", "site_b"))
  expect_identical(
    unname(vapply(result, `[[`, integer(1L), "ordinal")), 1:2)
  expect_true(all(vapply(result, `[[`, logical(1L), "accepted")))
})

test_that("LMM optional fallbacks never hide identity or session failures", {
  peer <- .dsvert_client_peer_not_recognized_condition(
    "site_a", strrep("a", 64L), strrep("b", 64L))
  poison <- .dsvert_dsi_poisoned_session_condition("site_a")

  peer_seen <- tryCatch(
    .dsvert_lmm_transport_fallback(peer, "fallback"),
    dsvert_peer_not_recognized = identity)
  poison_seen <- tryCatch(
    .dsvert_lmm_transport_fallback(poison, "fallback"),
    dsvert_dsi_poisoned_session = identity)

  expect_s3_class(peer_seen, "dsvert_peer_not_recognized")
  expect_s3_class(poison_seen, "dsvert_dsi_poisoned_session")
  expect_identical(
    .dsvert_lmm_transport_fallback(simpleError("optional"), "fallback"),
    "fallback")
})

test_that("Beaver policy discovery preserves typed transport conditions", {
  conns <- list(site_a = structure(list(), class = "mock_connection"))
  peer <- .dsvert_client_peer_not_recognized_condition(
    "site_a", strrep("c", 64L), strrep("d", 64L))
  poison <- .dsvert_dsi_poisoned_session_condition("site_a")

  peer_seen <- tryCatch(
    .beaver_fetch_policies(
      conns, 1L, function(...) stop(peer),
      session_id = paste0("lmm-peer-", .dsvert_uuid4())),
    dsvert_peer_not_recognized = identity)
  poison_seen <- tryCatch(
    .beaver_fetch_policies(
      conns, 1L, function(...) stop(poison),
      session_id = paste0("lmm-poison-", .dsvert_uuid4())),
    dsvert_dsi_poisoned_session = identity)

  expect_s3_class(peer_seen, "dsvert_peer_not_recognized")
  expect_s3_class(poison_seen, "dsvert_dsi_poisoned_session")

  ordinary <- tryCatch(
    .beaver_fetch_policies(
      conns, 1L, function(...) stop("ordinary"),
      session_id = paste0("lmm-ordinary-", .dsvert_uuid4())),
    error = identity)
  expect_s3_class(ordinary, "error")
  expect_false(inherits(ordinary, "dsvert_peer_not_recognized"))
  expect_false(inherits(ordinary, "dsvert_dsi_poisoned_session"))
  expect_match(conditionMessage(ordinary),
               "Could not verify the server IKNP policy", fixed = TRUE)
})

test_that("K=2 closed-form fit returns coherent inference and prints", {
  coefficient_names <- c("(Intercept)", "x")
  information <- matrix(c(4, 1, 1, 3), 2L, 2L,
                        dimnames = list(coefficient_names,
                                        coefficient_names))
  coefficients <- c("(Intercept)" = 1.25, x = -0.5)

  fit <- .dsvert_lmm_closed_form_fit(
    list(XtX = information), coefficients, sigma2 = 2)
  expected_covariance <- 2 * solve(information)
  expect_identical(names(fit), c("coefficients", "covariance", "std_errors"))
  expect_equal(fit$coefficients, coefficients)
  expect_equal(fit$covariance, expected_covariance)
  expect_equal(fit$std_errors, sqrt(diag(expected_covariance)))
  expect_identical(names(fit$std_errors), coefficient_names)

  result <- c(fit, list(
    sigma2 = 2, sigma_b2 = 0.5, icc = 0.2,
    n_clusters = 2L, cluster_sizes = c(5L, 5L),
    converged = TRUE, iterations = 1L))
  class(result) <- c("ds.vertLMM", "list")
  printed <- capture.output(returned <- print(result))
  expect_identical(returned, result)
  expect_match(paste(printed, collapse = "\n"),
               "Fixed effects:", fixed = TRUE)
  expect_false(any(grepl("NA", printed, fixed = TRUE)))
})
