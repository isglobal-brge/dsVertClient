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

test_that("multinomial legacy OVR and GLMM legacy EM fail closed", {
  old_mn <- getOption("dsvert.allow_multinom_legacy_ovr", NULL)
  old_glmm <- getOption("dsvert.allow_glmm_legacy_em", NULL)
  on.exit({
    if (is.null(old_mn)) {
      options(dsvert.allow_multinom_legacy_ovr = NULL)
    } else {
      options(dsvert.allow_multinom_legacy_ovr = old_mn)
    }
    if (is.null(old_glmm)) {
      options(dsvert.allow_glmm_legacy_em = NULL)
    } else {
      options(dsvert.allow_glmm_legacy_em = old_glmm)
    }
  }, add = TRUE)
  options(dsvert.allow_multinom_legacy_ovr = FALSE,
          dsvert.allow_glmm_legacy_em = FALSE)

  expect_error(
    ds.vertMultinomJoint(y ~ x, data = "D", levels = c("A", "B", "C"),
                         allow_legacy_ovr = TRUE, datasources = list()),
    "disabled by default")

  expect_error(
    ds.vertGLMM(y ~ x, data = "D", cluster_col = "id", method = "em",
                datasources = list()),
    "disabled by default")
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

test_that("GEE AR1 fails closed instead of falling back to independence", {
  expect_error(
    ds.vertGEE(y ~ x, data = "D", family = "gaussian", corstr = "ar1",
               datasources = list()),
    "not implemented as a true working-correlation estimator")
})
