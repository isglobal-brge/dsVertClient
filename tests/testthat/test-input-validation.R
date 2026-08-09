# Tests for input validation in client functions (no server needed)

# =============================================================================
# ds.vertCor input validation
# =============================================================================

test_that("ds.vertCor rejects non-character data_name", {
  conns <- list(a = structure(list(), class = "fake"))
  expect_error(ds.vertCor(
    123, list(a = c("x", "y")), analysis_id = "D::a",
    datasources = conns), "data_name must be one non-empty capsule identifier")
  expect_error(ds.vertCor(
    c("a", "b"), list(a = c("x", "y")), analysis_id = "D::a",
    datasources = conns), "data_name must be one non-empty capsule identifier")
})

test_that("ds.vertCor rejects non-list variables", {
  conns <- list(a = structure(list(), class = "fake"))
  expect_error(ds.vertCor(
    "D", 123, analysis_id = "D::a", datasources = conns),
    "variables must be NULL")
})

test_that("ds.vertCor rejects unnamed list variables", {
  conns <- list(a = structure(list(), class = "fake"))
  expect_error(ds.vertCor(
    "D", list(c("x"), c("y")), analysis_id = "D::a",
    datasources = conns), "variables must be a character vector or a named")
})

test_that("ds.vertCor accepts explicit cross-owner request syntax", {
  conns <- list(
    a = structure(list(), class = "fake"),
    b = structure(list(), class = "fake"))
  request <- .dsvert_dp_cor_request(
    list(a = "x", b = "y"), server = NULL, datasources = conns)
  expect_identical(request$variables, c("x", "y"))
  expect_identical(request$owner_map, c(x = "a", y = "b"))
  expect_identical(request$owner_count, 2L)
})

# =============================================================================
# ds.vertPCA input validation
# =============================================================================

test_that("ds.vertPCA rejects non-character data_name", {
  expect_error(ds.vertPCA(
    123, list(a = c("x", "y")), analysis_id = "D::a"),
    "data_name must be one non-empty capsule identifier")
})

test_that("ds.vertPCA rejects non-list variables", {
  conns <- list(a = structure(list(), class = "fake"))
  expect_error(ds.vertPCA(
    "D", 123, analysis_id = "D::a", datasources = conns),
    "variables must be NULL")
})

# =============================================================================
# ds.vertGLM input validation
# =============================================================================

test_that("ds.vertGLM rejects non-character data_name", {
  expect_error(ds.vertGLM(y ~ x, data = 123, datasources = list()),
               "data must be a single character string")
})

test_that("ds.vertGLM rejects non-list x_vars", {
  expect_error(ds.vertGLM(y ~ x, data = "D", x_vars = 123,
                          datasources = list()),
               "x_vars must be NULL")
})

test_that("ds.vertGLM requires an explicit supported missing-data policy", {
  expect_identical(formals(ds.vertGLM)$missing, "fail")
  expect_identical(formals(ds.vertGLM)$lambda, 0)
  expect_error(
    ds.vertGLM("D", "y", list(s1 = "x", s2 = "z"),
               missing = "silent_mean_impute", datasources = list()),
    "missing must be one of 'fail' or 'mean_impute'",
    fixed = TRUE)
})

test_that("ds.vertGLM rejects invalid numerical and protocol controls eagerly", {
  expect_error(ds.vertGLM(y ~ x, data = "D", max_iter = -1,
                          datasources = list()), "max_iter")
  expect_error(ds.vertGLM(y ~ x, data = "D", tol = 0,
                          datasources = list()), "tol")
  expect_error(ds.vertGLM(y ~ x, data = "D", lambda = -1,
                          datasources = list()), "lambda")
  expect_error(ds.vertGLM(y ~ x, data = "D", ring = 64,
                          datasources = list()), "ring")
  expect_error(ds.vertGLM(y ~ x, data = "D", eta_privacy = "plaintext",
                          datasources = list()), "eta_privacy")
  expect_error(ds.vertGLM(y ~ x, data = "D", compute_se = NA,
                          datasources = list()), "compute_se")
  expect_error(ds.vertGLM(y ~ x, data = "D", offset = c("a", "b"),
                          datasources = list()), "offset")
})

# =============================================================================
# ds.psiAlign input validation
# =============================================================================

test_that("ds.psiAlign rejects non-character data_name", {
  expect_error(ds.psiAlign(123, "id"),
               "data_name must be a single non-empty character string")
})

test_that("ds.psiAlign rejects non-character id_col", {
  expect_error(ds.psiAlign("D", 123),
               "id_col must be a single non-empty character string")
})

test_that("ds.psiAlign rejects non-character newobj", {
  expect_error(ds.psiAlign("D", "id", newobj = 123),
               "newobj must be a single non-empty character string")
})

test_that("ds.vertCox legacy rank route is removed", {
  expect_false("method" %in% names(formals(ds.vertCox)))
  expect_error(
    ds.vertCox(as.formula("Surv(time, event) ~ x1"), data = "D",
               method = "legacy_rank", datasources = list()),
    "unused argument")
})

test_that("Beaver preprocessing selection negotiates server policy", {
  mode <- getFromNamespace(".beaver_preprocessing_mode", "dsVertClient")
  old <- getOption("dsvert.beaver_preprocessing")
  on.exit(options(dsvert.beaver_preprocessing = old), add = TRUE)

  options(dsvert.beaver_preprocessing = "auto")
  expect_identical(mode("vecmul", 1L, 1L, 63L), "iknp")

  options(dsvert.beaver_preprocessing = "dealer")
  expect_error(mode("vecmul", 1L, 1L, 63L), "Invalid")

  for (value in c("ot", "iknp")) {
    options(dsvert.beaver_preprocessing = value)
    expect_identical(mode("vecmul", 1L, 1L, 63L), "iknp")
  }

  dsAgg_requires_iknp <- function(conns, expr) {
    list(list(supported = "iknp",
              allowed = "iknp",
              preferred = "iknp",
              minimum = "iknp",
              requires_iknp = TRUE))
  }
  options(dsvert.beaver_preprocessing = "auto")
  expect_identical(
    mode("vecmul", 1L, 1L, 63L,
         datasources = list(s1 = NULL, s2 = NULL), party_conns = 1:2,
         .dsAgg = dsAgg_requires_iknp, session_id = "test"),
    "iknp")

  options(dsvert.beaver_preprocessing = "direct_ot")
  expect_error(mode("vecmul", 1L, 1L, 63L), "Invalid")
})

test_that("client IKNP orchestration has no dealer or KOS downgrade path", {
  expect_false(exists(".dealer_prepare_vecmul",
                      envir = asNamespace("dsVertClient"), inherits = FALSE))
  expect_false(exists(".dealer_prepare_grad",
                      envir = asNamespace("dsVertClient"), inherits = FALSE))
  expect_false(exists(".dealer_prepare_spline",
                      envir = asNamespace("dsVertClient"), inherits = FALSE))
  implementation <- paste(deparse(body(.iknp_beaver_direction)),
                          collapse = "\n")
  expect_false(grepl("iknp_kos_check", implementation, fixed = TRUE))
  expect_match(implementation, "KOS consistency-check opener missing")
  expect_match(implementation, "kos_check = ext\\$kos_check")
})
