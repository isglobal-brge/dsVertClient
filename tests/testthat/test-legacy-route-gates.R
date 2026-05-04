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

test_that("user-facing Cox and categorical wrappers default to paper-safe routes", {
  cox_fn <- get("ds.vertCox", envir = asNamespace("dsVertClient"))
  expect_true("method" %in% names(formals(cox_fn)))
  cox_src <- paste(deparse(body(cox_fn)), collapse = "\n")
  expect_true(any(grepl("ds.vertCoxProfileNonDisclosive", cox_src,
                        fixed = TRUE)))
  expect_true(any(grepl(".ds.vertCoxLegacyRank", cox_src, fixed = TRUE)))

  mn_fn <- get("ds.vertMultinom", envir = asNamespace("dsVertClient"))
  expect_equal(formals(mn_fn)$method,
               as.call(list(as.name("c"), "joint", "warm")))
  mn_src <- paste(deparse(body(mn_fn)), collapse = "\n")
  expect_true(any(grepl("ds.vertMultinomJointNewton", mn_src, fixed = TRUE)))
  expect_true(any(grepl("paper-safe joint softmax", mn_src, fixed = TRUE)))

  ord_fn <- get("ds.vertOrdinal", envir = asNamespace("dsVertClient"))
  expect_equal(formals(ord_fn)$method,
               as.call(list(as.name("c"), "joint", "warm")))
  ord_src <- paste(deparse(body(ord_fn)), collapse = "\n")
  expect_true(any(grepl("ds.vertOrdinalJointNewton", ord_src, fixed = TRUE)))
  expect_true(any(grepl("paper-safe proportional odds", ord_src,
                        fixed = TRUE)))
})
