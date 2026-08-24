test_that("discarded routes are absent from the exported API", {
  ns <- asNamespace("dsVertClient")

  expect_false(exists("ds.vertCox.k3", envir = ns, inherits = FALSE))
  expect_false("method" %in% names(formals(get("ds.vertCox", ns))))
  expect_false("allow_disclosive_legacy" %in%
                 names(formals(get("ds.vertNBFullRegTheta", ns))))
  expect_false("allow_legacy_ovr" %in%
                 names(formals(get("ds.vertMultinomJoint", ns))))
  expect_false("method" %in% names(formals(get("ds.vertGLMM", ns))))
  expect_false("ds.vert.glmer" %in% getNamespaceExports(ns))
  expect_false("ds.vertGLMMLaplace" %in% getNamespaceExports(ns))
  expect_false("ds.vertNB" %in% getNamespaceExports(ns))
  expect_false("ds.vertNBMoMTheta" %in% getNamespaceExports(ns))
  expect_false("method" %in% names(formals(get("ds.vertMultinom", ns))))
  expect_false("method" %in% names(formals(get("ds.vertOrdinal", ns))))
  expect_false("allow_warm_diagnostic" %in%
                 names(formals(get("ds.vertMultinom", ns))))
  expect_false("allow_warm_diagnostic" %in%
                 names(formals(get("ds.vertOrdinal", ns))))
})

test_that("discarded route arguments fail as removed, not gated", {
  expect_error(
    ds.vertCox(as.formula("Surv(time, event) ~ x1"), data = "D",
               method = "legacy_rank", datasources = list()),
    "unused argument")

  expect_error(
    ds.vertNBFullRegTheta(y ~ x, data = "D", variant = "full_reg",
                          datasources = list()),
    "variant must")

  expect_error(
    ds.vertMultinomJoint(y ~ x, data = "D", levels = c("A", "B", "C"),
                         allow_legacy_ovr = TRUE, datasources = list()),
    "unused argument")

  expect_error(
    ds.vertGLMM(y ~ x, data = "D", cluster_col = "id", method = "em",
                datasources = list()),
    "unused argument")

  expect_error(
    ds.vertMultinom(y ~ x, data = "D", classes = c("A", "B", "C"),
                    method = "warm", datasources = list()),
    "unused argument")

  expect_error(
    ds.vertOrdinal(y ~ x, data = "D", levels_ordered = c("low", "mid", "hi"),
                   method = "warm", datasources = list()),
    "unused argument")
})

test_that("user-facing wrappers dispatch only to product routes", {
  cox_src <- paste(deparse(body(ds.vertCox)), collapse = "\n")
  expect_true(any(grepl(".dsvert_block_retired_remote_route", cox_src,
                        fixed = TRUE)))
  expect_false(any(grepl(".ds.vertCoxLegacyRank", cox_src, fixed = TRUE)))

  mn_src <- paste(deparse(body(ds.vertMultinom)), collapse = "\n")
  expect_true(any(grepl("ds.vertMultinomJointNewton", mn_src, fixed = TRUE)))
  expect_false(any(grepl(".ds_vertMultinomWarm", mn_src, fixed = TRUE)))

  mnj_src <- paste(deparse(body(ds.vertMultinomJoint)), collapse = "\n")
  expect_true(any(grepl("ds.vertMultinomJointNewton", mnj_src,
                        fixed = TRUE)))
  expect_false(any(grepl("allow_legacy_ovr", mnj_src, fixed = TRUE)))

  ord_src <- paste(deparse(body(ds.vertOrdinal)), collapse = "\n")
  expect_true(any(grepl("ds.vertOrdinalJointNewton", ord_src, fixed = TRUE)))
  expect_false(any(grepl(".ds_vertOrdinalWarm", ord_src, fixed = TRUE)))

  glmm_src <- paste(deparse(body(ds.vertGLMM)), collapse = "\n")
  expect_true(any(grepl("ds.vertDPLMM", glmm_src, fixed = TRUE)))
  expect_false(any(grepl(".ds_glmm_pql_aggregate_loop", glmm_src,
                         fixed = TRUE)))

  lmm_src <- paste(deparse(body(ds.vertLMM)), collapse = "\n")
  expect_true(any(grepl("ds.vertDPLMM", lmm_src, fixed = TRUE)))
  expect_false(any(grepl(".dsvert_lmm_aggregate_strict", lmm_src,
                         fixed = TRUE)))
  expect_false(any(grepl("dsvertClusterSizesDS", lmm_src, fixed = TRUE)))
  expect_false(any(grepl(".ds_vertLMM_k3_impl", lmm_src, fixed = TRUE)))
  print_lmm_src <- paste(deparse(getS3method("print", "ds.vertLMM")),
                         collapse = "\n")
  expect_false(any(grepl("std_errors", print_lmm_src, fixed = TRUE)))
})

test_that("retired GLMM PQL implementation is not retained", {
  glmm_source <- paste(
    readLines(.dsvert_client_source_file("ds.vertGLMM.R"), warn = FALSE),
    collapse = "\n")

  expect_false(grepl(".ds_glmm_pql_aggregate_loop", glmm_source,
                     fixed = TRUE))
  expect_false(grepl(".ds_glmm_pql_solve_components", glmm_source,
                     fixed = TRUE))
  expect_false(grepl("k2SetOffsetDS", glmm_source, fixed = TRUE))
})

test_that("internal warm starts remain internal helpers for joint methods", {
  ns <- asNamespace("dsVertClient")
  expect_false(".ds_vertMultinomWarm" %in% getNamespaceExports(ns))
  expect_false(".ds_vertOrdinalWarm" %in% getNamespaceExports(ns))

  mnj_src <- paste(deparse(body(ds.vertMultinomJointNewton)), collapse = "\n")
  expect_true(any(grepl(".ds_vertMultinomWarm", mnj_src, fixed = TRUE)))

  ordj_src <- paste(deparse(body(ds.vertOrdinalJointNewton)), collapse = "\n")
  expect_true(any(grepl(".ord_joint_secure_fit", ordj_src, fixed = TRUE)))
  secure_src <- paste(deparse(get(".ord_joint_secure_fit", ns)), collapse = "\n")
  expect_true(any(grepl(".ds_vertOrdinalWarm", secure_src, fixed = TRUE)))
})

test_that("clustered GEE remains unavailable before DSI", {
  expect_error(
    ds.vertGEE(y ~ x, data = "D", family = "gaussian", corstr = "ar1",
               datasources = list()),
    class = "dsvert_route_unavailable")
  src <- paste(deparse(body(ds.vertGEE)), collapse = "\n")
  expect_true(any(grepl(".dsvert_block_retired_remote_route", src,
                        fixed = TRUE)))
  expect_false(any(grepl(".ds_gee_secure_", src, fixed = TRUE)))
})
