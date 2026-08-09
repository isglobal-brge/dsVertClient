test_that("multinomial joint output cannot expose warm-start inference as final", {
  warm <- list(
    fits = list(class_a = list(coefficients = c(x = 1))),
    coefficients = matrix(1, 1, 1),
    coefficients_ovr = matrix(2, 1, 1),
    std_errors = matrix(3, 1, 1),
    family = "warm",
    n_obs = 12L
  )

  out <- dsVertClient:::.dsvert_label_joint_warm_start(warm, "multinomial")

  expect_null(out$fits)
  expect_null(out$coefficients)
  expect_null(out$std_errors)
  expect_equal(out$warm_start$std_errors, warm$std_errors)
  expect_equal(out$n_obs, 12L)
  expect_identical(out$inference_status,
                   "unavailable_for_joint_estimator")
})

test_that("ordinal joint output segregates every warm-start inference field", {
  warm <- list(
    fits = list(), thresholds = 1:2, thresholds_ovr = 2:3,
    beta = matrix(1), beta_po = c(x = 1), covariance_po = matrix(1),
    po_test = list(p_value = 0.1), std_errors_po = c(x = 1),
    joint_mle = list(covariance = matrix(1)), family = "warm",
    levels = letters[1:3], n_obs = 20L
  )

  out <- dsVertClient:::.dsvert_label_joint_warm_start(warm, "ordinal")

  expect_null(out$beta_po)
  expect_null(out$covariance_po)
  expect_null(out$po_test)
  expect_null(out$joint_mle)
  expect_equal(out$warm_start$po_test, warm$po_test)
  expect_equal(out$levels, warm$levels)
  expect_match(out$inference_reason, "not valid for the final joint")
})

test_that("joint output contract rejects unknown estimators", {
  expect_error(
    dsVertClient:::.dsvert_label_joint_warm_start(list(), "other"),
    "unknown joint estimator"
  )
})
