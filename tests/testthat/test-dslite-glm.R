# DSLite Integration Tests: ds.vertGLM with Different Server Counts
#
# These tests run the FULL distributed GLM pipeline (transport keys, PSI
# alignment, L-BFGS, deviance) using DSLite and compare against local lm()/glm().
# Requires: DSLite, dsVert (with dsvert-mpc binary), dsVertClient.

skip_on_cran()
skip_if_not_installed("DSLite")
skip_if_not_installed("dsVert")
skip_if_not(dsVert::mpcAvailable(), "dsvert-mpc binary not available")

# =============================================================================
# Helper: setup and run a full DSLite GLM scenario
# =============================================================================

setup_dslite_glm <- function(K, family = "gaussian", n = 200, seed = 42,
                              eta_privacy = "auto", ...) {
  set.seed(seed)

  # Generate data: 2 predictors on label server + 1 per non-label server
  n_predictors <- K + 1
  X <- matrix(rnorm(n * n_predictors), n, n_predictors)
  colnames(X) <- paste0("x", seq_len(n_predictors))
  true_beta <- seq(0.5, by = 0.3, length.out = n_predictors)

  eta <- as.vector(X %*% true_beta)
  if (family == "gaussian") {
    y <- eta + rnorm(n, 0, 0.5)
  } else if (family == "binomial") {
    prob <- 1 / (1 + exp(-eta))
    y <- rbinom(n, 1, prob)
  } else if (family == "poisson") {
    y <- rpois(n, exp(pmin(eta, 5)))
  }

  patient_ids <- paste0("P", sprintf("%04d", seq_len(n)))

  # Partition variables across K servers
  server_names <- paste0("s", LETTERS[seq_len(K)])  # sA, sB, sC, sD
  x_vars <- list()
  x_vars[[server_names[1]]] <- c("x1", "x2")  # label server gets 2 predictors
  for (k in 2:K) {
    x_vars[[server_names[k]]] <- paste0("x", k + 1)  # 1 predictor each
  }

  # Create data frames per server (shuffled independently for PSI)
  tables <- list()
  for (k in seq_len(K)) {
    order_k <- sample(n)
    vars_k <- x_vars[[server_names[k]]]
    df_k <- data.frame(patient_id = patient_ids[order_k],
                       stringsAsFactors = FALSE)
    for (v in vars_k) {
      col_idx <- as.integer(sub("x", "", v))
      df_k[[v]] <- X[order_k, col_idx]
    }
    if (k == 1) df_k[["outcome"]] <- y[order_k]
    tables[[server_names[k]]] <- df_k
  }

  # Local ground truth
  pooled <- as.data.frame(X)
  pooled$outcome <- as.numeric(y)
  if (family == "gaussian") {
    local_fit <- lm(outcome ~ . - 1, data = pooled)
  } else {
    local_fit <- glm(outcome ~ . - 1, data = pooled, family = family)
  }

  # Setup DSLite server
  dslite_server <- DSLite::newDSLiteServer(tables = tables)
  dslite_server$config(DSLite::defaultDSConfiguration(include = c("dsVert")))
  assign("dslite_server", dslite_server, envir = globalenv())

  builder <- DSI::newDSLoginBuilder()
  for (sname in server_names) {
    builder$append(server = sname, url = "dslite_server",
                   table = sname, driver = "DSLiteDriver")
  }
  conns <- DSI::datashield.login(builder$build(), assign = TRUE, symbol = "D")

  # PSI alignment
  dsVertClient::ds.psiAlign("D", "patient_id", "D_aligned", datasources = conns)

  # Run distributed GLM
  distributed_fit <- dsVertClient::ds.vertGLM(
    "D_aligned", "outcome", x_vars,
    y_server = server_names[1],
    family = family,
    lambda = 1e-6,
    tol = 1e-6,
    max_iter = 200,
    verbose = FALSE,
    datasources = conns,
    eta_privacy = eta_privacy,
    ...
  )

  list(
    local = local_fit,
    distributed = distributed_fit,
    conns = conns,
    K = K,
    family = family,
    server_names = server_names
  )
}

# Helper to clean up after each test
cleanup <- function(result) {
  tryCatch(DSI::datashield.logout(result$conns), error = function(e) NULL)
  if (exists("dslite_server", envir = globalenv()))
    rm("dslite_server", envir = globalenv())
}

# Helper to compare coefficients (distributed includes intercept, local doesn't)
compare_coefs <- function(result, tolerance) {
  local_coefs <- coef(result$local)
  dist_coefs <- coef(result$distributed)
  dist_coefs_no_int <- dist_coefs[names(dist_coefs) != "(Intercept)"]

  for (v in names(local_coefs)) {
    expect_equal(
      unname(dist_coefs_no_int[[v]]),
      unname(local_coefs[[v]]),
      tolerance = tolerance,
      label = paste("coefficient", v)
    )
  }
}

# =============================================================================
# K=2, gaussian → transport mode
# =============================================================================

test_that("K=2 gaussian: auto selects transport, matches lm()", {
  result <- setup_dslite_glm(K = 2, family = "gaussian")
  on.exit(cleanup(result))

  expect_equal(result$distributed$eta_privacy, "transport")
  expect_true(result$distributed$converged)
  compare_coefs(result, tolerance = 0.05)
})

# =============================================================================
# K=2, binomial → he_link mode
# =============================================================================

test_that("K=2 binomial: auto selects he_link, matches glm()", {
  result <- setup_dslite_glm(K = 2, family = "binomial")
  on.exit(cleanup(result))

  expect_equal(result$distributed$eta_privacy, "he_link")
  expect_true(result$distributed$converged)
  compare_coefs(result, tolerance = 0.15)
})

# =============================================================================
# K=2, binomial with forced transport (user opt-in)
# =============================================================================

test_that("K=2 binomial: forced transport works", {
  result <- setup_dslite_glm(K = 2, family = "binomial", eta_privacy = "transport")
  on.exit(cleanup(result))

  expect_equal(result$distributed$eta_privacy, "transport")
  expect_true(result$distributed$converged)
  compare_coefs(result, tolerance = 0.15)
})

# =============================================================================
# K=3, gaussian → secure_agg mode
# =============================================================================

test_that("K=3 gaussian: auto selects secure_agg, matches lm()", {
  result <- setup_dslite_glm(K = 3, family = "gaussian")
  on.exit(cleanup(result))

  expect_equal(result$distributed$eta_privacy, "secure_agg")
  expect_true(result$distributed$converged)
  compare_coefs(result, tolerance = 0.05)
})

# =============================================================================
# K=3, binomial → secure_agg mode
# =============================================================================

test_that("K=3 binomial: auto selects secure_agg, matches glm()", {
  result <- setup_dslite_glm(K = 3, family = "binomial")
  on.exit(cleanup(result))

  expect_equal(result$distributed$eta_privacy, "secure_agg")
  expect_true(result$distributed$converged)
  compare_coefs(result, tolerance = 0.15)
})

# =============================================================================
# K=4, gaussian → secure_agg mode
# =============================================================================

test_that("K=4 gaussian: auto selects secure_agg, matches lm()", {
  result <- setup_dslite_glm(K = 4, family = "gaussian")
  on.exit(cleanup(result))

  expect_equal(result$distributed$eta_privacy, "secure_agg")
  expect_true(result$distributed$converged)
  compare_coefs(result, tolerance = 0.05)
})

# =============================================================================
# Model diagnostics (deviance, pseudo_r2, AIC)
# =============================================================================

test_that("Distributed GLM returns valid deviance and diagnostics", {
  result <- setup_dslite_glm(K = 3, family = "gaussian")
  on.exit(cleanup(result))

  expect_true(result$distributed$deviance > 0)
  expect_true(result$distributed$null_deviance > result$distributed$deviance)
  expect_true(result$distributed$pseudo_r2 > 0)
  expect_true(result$distributed$pseudo_r2 < 1)
  expect_true(is.finite(result$distributed$aic))
  expect_equal(result$distributed$n_obs, 200)
})
