test_that("ds.vertOrdinalJointNewton strict path is the only route", {
  expect_error(
    ds.vertOrdinalJointNewton(
      y ~ x,
      levels_ordered = c("low", "mid", "high"),
      datasources = list()),
    "requires at least two servers")
})
