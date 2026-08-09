test_that("plain additive GLM formulas have an explicit design contract", {
  parsed <- dsVertClient:::.dsvert_plain_formula(y ~ x + z)
  expect_identical(parsed$response, "y")
  expect_identical(parsed$predictors, c("x", "z"))
  expect_true(parsed$intercept)

  no_intercept <- dsVertClient:::.dsvert_plain_formula("y ~ x + z - 1")
  expect_false(no_intercept$intercept)
  expect_identical(no_intercept$predictors, c("x", "z"))
})

test_that("unsupported formula designs fail before any server call", {
  unsupported <- list(
    y ~ log(x),
    y ~ I(x^2),
    y ~ factor(group),
    y ~ x * z,
    y ~ x:z,
    y ~ .,
    cbind(y1, y2) ~ x
  )
  for (formula in unsupported) {
    expect_error(
      dsVertClient:::.dsvert_plain_formula(formula),
      "pre-existing numeric column|Only additive"
    )
  }
  expect_error(dsVertClient:::.dsvert_plain_formula("y ~~ x"),
               "Invalid formula")
  expect_error(dsVertClient:::.dsvert_plain_formula(y ~ x | z),
               "Only additive")
})

test_that("ds.vertGLM rejects a misleading transformed formula eagerly", {
  expect_error(
    ds.vertGLM(y ~ log(x), data = "D", datasources = list()),
    "Only additive pre-existing numeric columns"
  )
})
