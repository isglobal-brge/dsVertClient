test_that("NB full_reg legacy is redirected unless explicitly allowed", {
  fn <- get("ds.vertNBFullRegTheta", envir = asNamespace("dsVertClient"))
  expect_true("allow_disclosive_legacy" %in% names(formals(fn)))

  src <- paste(deparse(body(fn)), collapse = "\n")
  expect_true(any(grepl("variant <- \"full_reg_nd\"", src, fixed = TRUE)))
  expect_true(any(grepl("per-patient non-label eta", src, fixed = TRUE)))
  expect_true(any(grepl("allow_disclosive_legacy = TRUE", src, fixed = TRUE)))
})

test_that("multinomial joint legacy wrapper dispatches to JointNewton", {
  fn <- get("ds.vertMultinomJoint", envir = asNamespace("dsVertClient"))
  expect_true("allow_legacy_ovr" %in% names(formals(fn)))

  src <- paste(deparse(body(fn)), collapse = "\n")
  expect_true(any(grepl("ds.vertMultinomJointNewton", src, fixed = TRUE)))
  expect_true(any(grepl("allow_legacy_ovr", src, fixed = TRUE)))
  expect_true(any(grepl("superseded", src, fixed = TRUE)))
})
