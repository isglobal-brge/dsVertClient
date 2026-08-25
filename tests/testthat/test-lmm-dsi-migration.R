.lmm_client_source <- function(name) {
  readLines(.dsvert_client_source_file(name), warn = FALSE)
}

test_that("LMM compatibility frontdoors retain no legacy remote implementation", {
  files <- c("ds.vertLMM.R", "ds.vertLMM.k3.R")
  source <- paste(unlist(lapply(files, .lmm_client_source), use.names = FALSE),
                  collapse = "\n")

  expect_false(grepl("DSI::datashield.aggregate", source, fixed = TRUE))
  expect_false(grepl("DSI::datashield.assign", source, fixed = TRUE))
  expect_false(grepl("dsvertClusterSizesDS", source, fixed = TRUE))
  expect_false(grepl(".ds_vertLMM_k3_impl", source, fixed = TRUE))
})

test_that("LMM help keeps both historical names on the signed finite-grid route", {
  main <- paste(.lmm_client_source("ds.vertLMM.R"), collapse = "\n")
  k3 <- paste(.lmm_client_source("ds.vertLMM.k3.R"), collapse = "\n")

  expect_match(main, "signed random-intercept", fixed = TRUE)
  expect_match(main, "method-of-moments", fixed = TRUE)
  expect_match(main, "unrestricted profile optimisation", fixed = TRUE)
  expect_match(main, "signed ML or REML", fixed = TRUE)
  expect_match(main, "@return A \\code{ds.vertLMM}", fixed = TRUE)
  expect_match(k3, "signed random-intercept", fixed = TRUE)
  expect_match(k3, "ML or REML profile", fixed = TRUE)
  expect_match(k3, "additive fixed-effect formula", fixed = TRUE)
  expect_match(k3, "at least three DataSHIELD connections", fixed = TRUE)
})

test_that("the K>=3 compatibility printer describes only its signed finite grid", {
  fit <- list(
    status = "ok", cluster_count = 2L, n_obs = 9L,
    sigma_b2 = 0.5, sigma2 = 1.5, icc = 0.2,
    coefficients = c(`(Intercept)` = 1, x = 0.25))
  class(fit) <- c("ds.vertLMM.k3", "list")

  printed <- capture.output(returned <- print(fit))
  expect_identical(returned, fit)
  expect_match(printed[[1L]], "K>=3", fixed = TRUE)
  expect_match(paste(printed, collapse = "\n"), "fixed signed grid",
               fixed = TRUE)
  expect_false(any(grepl("unrestricted", printed, fixed = TRUE)))
})
