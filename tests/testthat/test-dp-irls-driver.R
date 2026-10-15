.irls_driver_fixture <- function() {
  f <- jsonlite::fromJSON(paste(readLines(system.file("dp-irls-fixtures",
    "synthetic-chains.json", package = "dsVertClient"), warn = FALSE), collapse = ""),
    simplifyVector = FALSE)[[1L]]
  list(f = f, conns = stats::setNames(rep(list(structure(list(), class = "mock")),
    length(f$pins)), names(f$pins)), metadata = list(version = "dsvert-dp-irls-metadata-v1",
      analysis_id = "smoke", root = f$stages[[1L]]$body$root,
      signed_contract = f$signed_contract, stages = list()))
}

test_that("IRLS public driver rejects analyst controls before DSI", {
  calls <- 0L
  testthat::local_mocked_bindings(.dsvert_dp_irls_client_call = function(...) {
    calls <<- calls + 1L; stop("must not reach DSI")
  }, .package = "dsVertClient")
  f <- .irls_driver_fixture()
  for (name in c("beta", "predecessor", "epsilon", "T", "family", "seed", "bounds")) {
    args <- c(list(analysis_id = "smoke", datasources = f$conns), setNames(list("x"), name))
    expect_error(do.call(ds.vertDPIRLS, args), "unused argument")
  }
  expect_error(ds.vertDPIRLS("bad selector", f$conns), "identifier")
  withr::local_options(dsvert.dp_irls.trusted_pins = NULL)
  expect_error(ds.vertDPIRLS("smoke", f$conns), "trusted_pins")
  expect_identical(calls, 0L)
})

test_that("IRLS metadata authenticates pins and rejects unrequested openings", {
  f <- .irls_driver_fixture(); value <- f$metadata
  testthat::local_mocked_bindings(.dsvert_dp_irls_client_call = function(...) {
    stats::setNames(rep(list(value), 2L), names(f$conns))
  }, .package = "dsVertClient")
  expect_identical(.dsvert_dp_irls_client_metadata("smoke", f$conns,
    f$f$pins)$metadata, value)
  value$exact_score <- list("0", "0")
  expect_error(.dsvert_dp_irls_client_metadata("smoke", f$conns, f$f$pins), "metadata")
  value <- f$metadata; value$signed_contract$signatures[[1L]] <- strrep("A", 86L)
  expect_error(.dsvert_dp_irls_client_metadata("smoke", f$conns, f$f$pins), "signature")
})

test_that("IRLS signed translation preserves negatives and flags finalizer endpoints", {
  f <- .irls_driver_fixture(); contract <- f$f$signed_contract$contract
  stage <- f$f$stages[[1L]]
  stage$body$payload$g <- list("0", "15")
  stage$body$payload$H <- list("32", "16", "16")
  result <- .dsvert_dp_irls_client_decode(stage, contract)
  expect_true(result$endpoint_hit)
  expect_identical(result$released$g, list(list(n = "-16", d = "1"), list(n = "-1", d = "1")))
  expect_identical(result$released$H[[1L]], list(n = "16", d = "1"))
  stage$body$payload <- list(g = list("15", "17"), H = list("17", "15", "17"))
  expect_false(.dsvert_dp_irls_client_decode(stage, contract)$endpoint_hit)
})

test_that("IRLS publication requires both matching signed chains", {
  f <- .irls_driver_fixture()
  value <- list(version = "dsvert-dp-irls-publication-v1", analysis_id = "smoke",
    root = f$metadata$root, stages = f$f$stages[1L])
  right <- value
  testthat::local_mocked_bindings(.dsvert_dp_irls_client_call = function(...) {
    stats::setNames(list(value, right), names(f$conns))
  }, .package = "dsVertClient")
  expect_identical(.dsvert_dp_irls_client_publication("smoke", f$conns,
    f$metadata, f$f$pins), value)
  right$stages <- list()
  expect_error(.dsvert_dp_irls_client_publication("smoke", f$conns,
    f$metadata, f$f$pins), "both authorities")
  right <- value; right$stages[[1L]]$body$beta[[1L]]$n <- "1"
  expect_error(.dsvert_dp_irls_client_publication("smoke", f$conns,
    f$metadata, f$f$pins), "signature")
})
