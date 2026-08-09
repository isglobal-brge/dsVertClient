test_that("epi 2x2 measures match standard closed-form estimators", {
  tab <- matrix(
    c(90, 60, 10, 40), nrow = 2L,
    dimnames = list(exposure = c("unexposed", "exposed"),
                    outcome = c("no", "event")))

  fit <- ds.vertEpi2x2(tab, exposed = "exposed", event = "event")

  expect_equal(fit$risks[["exposed"]], 0.4)
  expect_equal(fit$risks[["unexposed"]], 0.1)
  expect_equal(fit$risks[["population"]], 0.25)
  expect_equal(fit$risk_difference$estimate, 0.3)
  expect_equal(fit$risk_ratio$estimate, 4)
  expect_equal(fit$odds_ratio$estimate, 6)
  expect_equal(fit$attributable_fraction_exposed$estimate, 0.75)
  expect_equal(fit$population_attributable_fraction$estimate, 0.6)
  expect_equal(fit$nnt$estimate, 1 / 0.3)
  expect_equal(fit$number_needed$estimate, 1 / 0.3)
  expect_identical(fit$number_needed$direction, "harm")
  expect_true(all(is.finite(unlist(fit$risk_ratio[c("lower", "upper")]))))
  expect_identical(fit$correction, 0)
  expect_identical(fit$method_status, "promoted")
})

test_that("epi 2x2 accepts an existing disclosure-checked table object", {
  x <- list(
    observed = matrix(c(70, 50, 30, 50), 2L,
                      dimnames = list(c("no", "yes"), c("control", "case"))),
    disclosure_guard = list(passed = TRUE, threshold = 5L))
  class(x) <- c("ds.vertChisq", "list")

  fit <- ds.vertEpi2x2(x, exposed = "yes", event = "case")

  expect_identical(fit$disclosure_guard, x$disclosure_guard)
  expect_identical(fit$input_provenance, "disclosure_checked_table")
  expect_equal(fit$n, 200)
})

test_that("epi 2x2 handles zero cells only with an explicit reported correction", {
  tab <- matrix(c(20, 10, 0, 10), 2L,
                dimnames = list(c("u", "e"), c("no", "event")))

  corrected <- ds.vertEpi2x2(
    tab, exposed = "e", event = "event", zero_correction = "if_zero")
  expect_equal(corrected$correction, 0.5)
  expect_match(corrected$notes, "Haldane-Anscombe")
  expect_true(is.finite(corrected$odds_ratio$estimate))

  uncorrected <- ds.vertEpi2x2(tab, exposed = "e", event = "event")
  expect_equal(uncorrected$correction, 0)
  expect_true(is.infinite(uncorrected$odds_ratio$estimate))
  expect_identical(uncorrected$quality$status, "degraded")
  expect_true(uncorrected$quality$metrics$zero_cell_boundary)
})

test_that("epi 2x2 validates orientation and aggregate counts", {
  expect_error(ds.vertEpi2x2(matrix(1:6, 2L)), "exactly 2 x 2")
  expect_error(ds.vertEpi2x2(matrix(c(1, 2, -1, 3), 2L)),
               "non-negative")
  expect_error(ds.vertEpi2x2(matrix(c(1, 2, 3, 4.1), 2L)),
               "whole-number")

  tab <- matrix(1:4, 2L,
                dimnames = list(c("u", "e"), c("no", "event")))
  expect_error(ds.vertEpi2x2(tab, exposed = "missing"), "Unknown exposed")
  expect_error(ds.vertEpi2x2(tab, event = "missing"), "Unknown event")
})
