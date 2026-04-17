#' @title Generalized Linear Model for Vertically Partitioned Data
#' @description Fits a GLM across vertically partitioned data using Ring63
#'   Beaver MPC with DCF wide spline for the link function. The system
#'   auto-detects which server holds each variable. Only p-dimensional
#'   aggregate gradients are revealed to the client per iteration.
#'   No observation-level data is ever disclosed.
#'
#' @param formula A formula (e.g. \code{npreg ~ age + bmi + glu}) or
#'   character string (\code{"npreg ~ age + bmi + glu"}). Can also be
#'   a data_name string for backward compatibility.
#' @param data Character string. Name of the (aligned) data frame on
#'   each server.
#' @param x_vars Optional. Character vector of predictor names, or a
#'   named list mapping server names to variable vectors. If NULL and
#'   formula is used, extracted from the formula. If NULL and no formula,
#'   all available columns (minus y and IDs) are used.
#' @param y_server Character string. Name of the server holding the response.
#' @param family Character string. GLM family: "gaussian", "binomial",
#'   or "poisson". Default is "gaussian".
#' @param max_iter Integer. Maximum L-BFGS iterations. Default is 100.
#' @param tol Numeric. Convergence tolerance on coefficient change.
#'   Default is 1e-4.
#' @param lambda Numeric. L2 regularization parameter. Default is 1e-4.
#' @param log_n Integer. Legacy parameter (ignored). Kept for backward
#'   compatibility.
#' @param verbose Logical. Print progress messages. Default is TRUE.
#' @param datasources DataSHIELD connection object or list of connections.
#' @param eta_privacy Character. \code{"auto"} (default) selects
#'   \code{"k2_beaver"} for K=2 or \code{"secure_agg"} for K>=3.
#'
#' @return A list with class "ds.glm" containing:
#'   \itemize{
#'     \item \code{coefficients}: Named coefficient vector (original scale)
#'     \item \code{std_errors}: Standard errors (finite-difference Hessian)
#'     \item \code{z_values}: z-statistics (coef / SE)
#'     \item \code{p_values}: Two-sided p-values
#'     \item \code{iterations}: Number of iterations
#'     \item \code{converged}: Logical
#'     \item \code{family}: Family used
#'     \item \code{n_obs}: Number of observations
#'     \item \code{deviance}: Residual sum of squares
#'     \item \code{pseudo_r2}: 1 - deviance/null_deviance
#'   }
#'
#' @details
#' \subsection{Protocol}{
#' All computation uses Ring63 fixed-point arithmetic with Beaver MPC:
#' \enumerate{
#'   \item \strong{Transport keys}: X25519 keypairs on all servers (~0.5s)
#'   \item \strong{Standardize}: Each server standardizes its features
#'   \item \strong{Input sharing}: Features split into additive Ring63 shares
#'     between 2 DCF parties. Non-DCF servers contribute shares.
#'   \item \strong{L-BFGS loop}: Per iteration:
#'     \itemize{
#'       \item Compute eta shares (Ring63 matrix-vector)
#'       \item DCF wide spline for sigmoid/exp (binomial/Poisson) or
#'         identity link (Gaussian)
#'       \item Beaver matvec for gradient (server-generated triples)
#'       \item Client aggregates Ring63 shares -> p gradient scalars
#'       \item L-BFGS quasi-Newton update
#'     }
#'   \item \strong{SE}: p+1 gradient evaluations (finite-difference Hessian)
#'   \item \strong{Deviance}: Beaver dot-product for residual sum of squares
#'   \item \strong{Unstandardize}: Coefficients + SE via Jacobian transform
#' }
#' }
#'
#' \subsection{Security}{
#' No observation-level data is disclosed. The client sees only p-dimensional
#' aggregate gradients per iteration. Beaver triples are generated server-side
#' (never seen by client). Dealer rotation for K>=4 ensures the analyst must
#' compromise (K-1)/K servers to extract data.
#' }
#'
#' @references
#' van Kesteren, E.J. et al. (2019). Privacy-preserving generalized linear
#' models using distributed block coordinate descent. arXiv:1911.03183.
#'
#' Mouchet, C. et al. (2021). "Multiparty Homomorphic Encryption from
#' Ring-Learning-With-Errors". \emph{Proceedings on Privacy Enhancing Technologies}.
#'
#' @seealso \code{\link{ds.vertCor}} for correlation analysis,
#'   \code{\link{ds.vertPCA}} for PCA analysis
#'
#' @examples
#' \dontrun{
#' # Simplest: formula interface (auto-detects everything)
#' model <- ds.vertGLM(npreg ~ age + bmi + glu + bp + skin,
#'                      data = "DA", family = "gaussian")
#'
#' # String formula
#' model <- ds.vertGLM("diabetes ~ age + bmi + glu",
#'                      data = "DA", family = "binomial")
#'
#' # Auto-detect all features
#' model <- ds.vertGLM("DA", "npreg", family = "poisson")
#'
#' # Manual server mapping (legacy)
#' model <- ds.vertGLM("DA", "npreg",
#'   list(s1 = c("age", "bmi"), s2 = c("glu", "bp")),
#'   y_server = "s2", family = "gaussian")
#' }
#'
#' @importFrom DSI datashield.aggregate datashield.connections_find
#' @export
ds.vertGLM <- function(formula, data = NULL, x_vars = NULL, y_server = NULL,
                       family = "gaussian", max_iter = 100, tol = 1e-4,
                       lambda = 1e-4, log_n = 12,
                       verbose = TRUE, datasources = NULL,
                       eta_privacy = "auto",
                       # Legacy positional args for backward compatibility
                       data_name = NULL, y_var = NULL) {
  call_matched <- match.call()

  # ===========================================================================
  # Parse formula or legacy arguments
  # ===========================================================================
  if (!missing(formula)) {
    if (inherits(formula, "formula")) {
      # R formula object: npreg ~ age + bmi + ped
      f_terms <- terms(formula)
      y_var <- as.character(attr(f_terms, "variables")[[2]])
      x_from_formula <- attr(f_terms, "term.labels")
      if (is.null(x_vars)) x_vars <- x_from_formula
      data_name <- data
    } else if (is.character(formula) && grepl("~", formula)) {
      # String formula: "npreg ~ age + bmi + ped"
      f <- as.formula(formula)
      f_terms <- terms(f)
      y_var <- as.character(attr(f_terms, "variables")[[2]])
      x_from_formula <- attr(f_terms, "term.labels")
      if (is.null(x_vars)) x_vars <- x_from_formula
      data_name <- data
    } else if (is.character(formula) && !grepl("~", formula)) {
      # Legacy: first arg is data_name (backward compat)
      data_name <- formula
      y_var <- data
      # x_vars stays as passed
    }
  }
  if (is.null(data_name) && !is.null(data)) data_name <- data

  # ===========================================================================
  # Input Validation + Smart Auto-Detection
  # ===========================================================================

  # data can be: "DA" (same name on all servers) or list(s1="tableA", s2="tableB")
  # For now we require all servers use the same data frame name (standard in DataSHIELD).
  if (is.list(data_name) && !is.null(names(data_name))) {
    # Named list: validate and extract. For future use.
    # Currently DataSHIELD requires the same assign name on all servers.
    stop("Named list for data is not yet supported. Use ds.psiAlign() to align data first, ",
         "then pass the aligned name (e.g. 'DA').", call. = FALSE)
  }
  if (!is.character(data_name) || length(data_name) != 1)
    stop("data must be a single character string (the name of the aligned data frame on all servers).",
         call. = FALSE)
  if (!is.character(y_var) || length(y_var) != 1)
    stop("y_var must be a single character string", call. = FALSE)
  if (!family %in% c("gaussian", "binomial", "poisson"))
    stop("family must be 'gaussian', 'binomial', or 'poisson'",
         call. = FALSE)

  if (is.null(datasources))
    datasources <- DSI::datashield.connections_find()
  if (length(datasources) == 0)
    stop("No DataSHIELD connections found", call. = FALSE)

  # Auto-detect: query servers for their columns and map variables automatically
  if (is.null(x_vars) || is.character(x_vars)) {
    user_x_vars <- x_vars  # NULL = use all available, character = specific vars
    if (verbose) message("[Auto-detect] Querying server columns...")
    col_results <- DSI::datashield.aggregate(datasources,
      call("dsvertColNamesDS", data_name = data_name))
    server_names <- names(datasources)

    # Build column map: which server has which variable (exclude IDs)
    col_map <- list()
    all_available <- character(0)
    for (srv in server_names) {
      cols <- setdiff(col_results[[srv]]$columns, c("id", "patient_id"))
      col_map[[srv]] <- cols
      all_available <- c(all_available, cols)
    }
    all_available <- unique(all_available)

    # Validate: all requested variables must exist somewhere
    if (!is.null(user_x_vars)) {
      missing <- setdiff(user_x_vars, all_available)
      if (length(missing) > 0)
        stop("Variables not found on any server: ", paste(missing, collapse = ", "),
             "\n  Available: ", paste(all_available, collapse = ", "), call. = FALSE)
    }
    if (!y_var %in% all_available)
      stop("Response variable '", y_var, "' not found on any server.\n",
           "  Available columns per server:\n",
           paste(sapply(server_names, function(s)
             paste0("    ", s, ": ", paste(col_map[[s]], collapse = ", "))),
             collapse = "\n"), call. = FALSE)

    # Find y_server automatically
    if (is.null(y_server)) {
      y_servers <- server_names[sapply(server_names, function(s) y_var %in% col_map[[s]])]
      y_server <- y_servers[1]
      if (length(y_servers) > 1 && verbose)
        message("  Note: '", y_var, "' found on multiple servers (",
                paste(y_servers, collapse = ", "), "). Using: ", y_server)
      if (verbose) message("  y_var '", y_var, "' found on: ", y_server)
    }

    # Build x_vars automatically
    x_vars <- list()
    for (srv in server_names) {
      feats <- col_map[[srv]]
      # Always exclude y_var from features (even on non-label servers)
      feats <- setdiff(feats, y_var)
      # If user specified specific vars, filter to those
      if (!is.null(user_x_vars))
        feats <- intersect(feats, user_x_vars)
      x_vars[[srv]] <- feats
    }

    # Remove response-only vars from feature lists
    # Keep server in x_vars even if it has 0 features (label server)
    if (verbose) {
      for (srv in server_names) {
        role <- if (srv == y_server) " (label)" else ""
        if (length(x_vars[[srv]]) > 0)
          message("  ", srv, role, ": ", paste(x_vars[[srv]], collapse = ", "))
        else
          message("  ", srv, role, ": (response only)")
      }
    }
  }

  if (!is.list(x_vars) || is.null(names(x_vars)))
    stop("x_vars must be a named list mapping server names to variable vectors",
         call. = FALSE)
  if (is.null(y_server))
    stop("y_server must be specified: the server holding '", y_var, "'",
         call. = FALSE)
  if (!y_server %in% names(x_vars))
    stop("y_server '", y_server, "' must be in x_vars", call. = FALSE)

  n_partitions_check <- length(x_vars)
  non_label_count <- n_partitions_check - 1

  # Route to the appropriate protocol:
  #   K=2: Beaver MPC with wide spline DCF (non-disclosive, ~1e-3 precision)
  #   K>=3: Secure aggregation with Enc(r) + L-BFGS (~8e-3 precision)
  if (eta_privacy == "auto") {
    if (non_label_count >= 2) {
      eta_privacy <- "secure_agg"
    } else if (non_label_count >= 1) {
      eta_privacy <- "k2_beaver"
    } else {
      stop("Need at least 2 servers (1 label + 1 non-label)", call. = FALSE)
    }
  }

  use_secure_agg <- (eta_privacy == "secure_agg")
  use_k2_beaver <- (eta_privacy == "k2_beaver")

  if (use_secure_agg && non_label_count < 2)
    stop("secure_agg requires >= 3 servers", call. = FALSE)
  if (use_k2_beaver && non_label_count != 1)
    stop("K=2 mode requires exactly 2 servers", call. = FALSE)

  # Adaptive log_n for large n: max_slots = 2^(log_n-1)
  # n ≤ 4096 → log_n=13, n ≤ 8192 → log_n=14, n ≤ 16384 → log_n=15
  # Note: this is set here before n_obs is known; adjusted later after getObsCountDS

  # ===========================================================================
  # Setup (datasources already resolved above for auto-detect)
  # ===========================================================================

  server_names <- names(datasources)
  if (!all(names(x_vars) %in% server_names)) {
    missing <- setdiff(names(x_vars), server_names)
    stop("Unknown server(s) in x_vars: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  server_list <- names(x_vars)
  non_label_servers <- setdiff(server_list, y_server)
  n_partitions <- length(x_vars)

  # Get observation count (lightweight sync call before helpers are available)
  first_conn <- which(server_names == server_list[1])
  count_result <- DSI::datashield.aggregate(
    conns = datasources[first_conn],
    expr = call("getObsCountDS", data_name)
  )
  if (is.list(count_result) && length(count_result) == 1)
    count_result <- count_result[[1]]
  n_obs <- count_result$n_obs

  # Adaptive log_n: ensure max_slots >= n_obs
  max_slots <- 2^(log_n - 1)
  if (n_obs > max_slots) {
    new_log_n <- ceiling(log2(n_obs)) + 1
    if (new_log_n > 15) {
      warning(sprintf("n_obs=%d exceeds max for log_n=15 (%d slots). Using log_n=15 with observation batching.",
                       n_obs, 2^14))
      new_log_n <- 15L
    }
    if (verbose) message(sprintf("  [Adaptive] log_n bumped %d→%d for n=%d observations", log_n, new_log_n, n_obs))
    log_n <- as.integer(new_log_n)
  }

  n_vars_total <- sum(sapply(x_vars, length))

  if (verbose) {
    message(sprintf("=== Encrypted-Label BCD-IRLS for %s GLM ===", family))
    message(sprintf("Observations: %d, Variables: %d, Partitions: %d",
                    n_obs, n_vars_total, n_partitions))
    message(sprintf("Label server: %s (holds '%s')", y_server, y_var))
    if (length(non_label_servers) > 0)
      message(sprintf("Non-label servers: %s",
                      paste(non_label_servers, collapse = ", ")))
  }

  # Generate session_id for all protocol phases (crypto state isolation)
  session_id <- local({
    hex <- sample(c(0:9, letters[1:6]), 32, replace = TRUE)
    hex[13] <- "4"  # UUID v4
    hex[17] <- sample(c("8","9","a","b"), 1)  # variant 1
    paste0(
      paste(hex[1:8], collapse = ""), "-",
      paste(hex[9:12], collapse = ""), "-",
      paste(hex[13:16], collapse = ""), "-",
      paste(hex[17:20], collapse = ""), "-",
      paste(hex[21:32], collapse = "")
    )
  })

  # Guaranteed cleanup on exit (even if error occurs mid-protocol).
  # Uses DSI::datashield.aggregate directly (not .dsAgg) so cleanup works
  # even if .glm_mpc_setup() fails before returning the closure.
  on.exit({
    for (.srv in server_list) {
      .ci <- which(server_names == .srv)
      tryCatch(
        DSI::datashield.aggregate(conns = datasources[.ci],
          expr = call("mpcCleanupDS", session_id = session_id)),
        error = function(e) NULL)
      tryCatch(
        DSI::datashield.aggregate(conns = datasources[.ci],
          expr = call("mpcGcDS")),
        error = function(e) NULL)
    }
    .dsvert_reset_chunk_size()
  }, add = TRUE)

  # ===========================================================================
  # Phase 0-1: Transport key setup + standardize features
  #   (delegated to .glm_mpc_setup in ds.vertGLM.setup.R)
  # ===========================================================================
  setup <- .glm_mpc_setup(
    datasources      = datasources,
    server_names     = server_names,
    server_list      = server_list,
    non_label_servers = non_label_servers,
    y_server         = y_server,
    y_var            = y_var,
    x_vars           = x_vars,
    data_name        = data_name,
    family           = family,
    session_id       = session_id,
    verbose          = verbose
  )

  # Unpack setup results
  transport_pks <- setup$transport_pks
  x_means       <- setup$x_means
  x_sds         <- setup$x_sds
  y_mean        <- setup$y_mean
  y_sd          <- setup$y_sd
  std_data      <- setup$std_data
  standardize_y <- setup$standardize_y
  .dsAgg        <- setup$.dsAgg
  .sendBlob     <- setup$.sendBlob

  # ===========================================================================
  # Phase 3: Iterative Ring63 Beaver (on standardized scale)
  # ===========================================================================
  label_intercept <- !standardize_y

  coordinator <- y_server
  coordinator_conn <- which(server_names == coordinator)

  betas <- list()
  for (server in server_list) {
    p <- length(x_vars[[server]])
    if (server == coordinator && label_intercept && !use_k2_beaver) p <- p + 1
    betas[[server]] <- rep(0, p)
  }

  converged <- FALSE
  final_iter <- 0

  # Encrypted eta blobs (opaque to client): server_name -> base64url string
  encrypted_etas <- list()

  if (use_k2_beaver) {
    # =========================================================================
    # Phase 3 (K=2 Strict v2): Chebyshev Beaver MPC
    # =========================================================================
    # Secure polynomial evaluation of sigmoid/exp on secret shares.
    # Neither party sees eta, mu, residuals, or weights in plaintext.
    # Only p_k gradient scalars + 2 intercept scalars revealed per iteration.

    if (verbose) message("\n[Phase 3] BCD iterations (K=2 strict Chebyshev Beaver)...")

    nl <- non_label_servers[1]
    nl_conn <- which(server_names == nl)

    loop_result <- .k2_strict_loop(
      datasources = datasources,
      server_names = server_names,
      server_list = server_list,
      coordinator = coordinator,
      coordinator_conn = coordinator_conn,
      non_label_servers = non_label_servers,
      nl = nl,
      nl_conn = nl_conn,
      x_vars = x_vars,
      y_var = y_var,
      std_data = std_data,
      transport_pks = transport_pks,
      session_id = session_id,
      family = family,
      lambda = lambda,
      max_iter = max_iter,
      tol = tol,
      n_obs = n_obs,
      verbose = verbose,
      .dsAgg = .dsAgg,
      .sendBlob = .sendBlob
    )

    betas <- loop_result$betas
    converged <- loop_result$converged
    final_iter <- loop_result$iterations
    # Store the loop's intercept for destandardization
    k2_loop_intercept <- loop_result$intercept

  } else if (use_secure_agg) {
    # ALL families: Ring63 Beaver gradient (Gaussian=identity link, others=DCF wide spline)
    if (verbose) message("\n[Phase 3] Ring63 Beaver Gradient (K=",
                         length(server_list), " servers, family=", family, ")...")
    k3_result <- .k3_ring63_gradient_loop(
      datasources = datasources, server_list = server_list,
      server_names = server_names, x_vars = x_vars,
      coordinator = coordinator, coordinator_conn = coordinator_conn,
      non_label_servers = non_label_servers, transport_pks = transport_pks,
      std_data = std_data, y_var = y_var, family = family,
      betas = betas, n_obs = n_obs, lambda = lambda,
      session_id = session_id,
      max_iter = max_iter, tol = tol, verbose = verbose,
      label_intercept = label_intercept,
      .dsAgg = .dsAgg, .sendBlob = .sendBlob)
    betas <- k3_result$betas
    converged <- k3_result$converged
    final_iter <- k3_result$final_iter

  } else {
    stop("Unsupported eta_privacy mode: ", eta_privacy, call. = FALSE)
  }
  # ===========================================================================
  # Phase 4: Unstandardize coefficients
  # ===========================================================================
  all_coefs_std <- numeric()
  all_x_means <- numeric()
  all_x_sds <- numeric()
  all_names <- character()
  beta_0_from_label <- 0
  # For K=2 beaver (wide spline), use the intercept from the loop
  if (use_k2_beaver && exists("k2_loop_intercept")) {
    beta_0_from_label <- k2_loop_intercept
  }

  for (server in server_list) {
    server_beta <- betas[[server]]
    if (server == coordinator && label_intercept && !use_k2_beaver) {
      beta_0_from_label <- server_beta[1]
      server_beta <- server_beta[-1]
    }
    all_coefs_std <- c(all_coefs_std, server_beta)
    all_x_means <- c(all_x_means, x_means[[server]])
    all_x_sds <- c(all_x_sds, x_sds[[server]])
    all_names <- c(all_names, x_vars[[server]])
  }

  if (standardize_y && !is.null(y_sd)) {
    all_coefs_orig <- all_coefs_std * y_sd / all_x_sds
    intercept <- y_mean - sum(all_coefs_orig * all_x_means)
  } else {
    all_coefs_orig <- all_coefs_std / all_x_sds
    intercept <- beta_0_from_label - sum(all_coefs_orig * all_x_means)
  }

  all_coefs <- c(intercept, all_coefs_orig)
  names(all_coefs) <- c("(Intercept)", all_names)
  n_vars_total <- length(all_coefs)

  # ===========================================================================
  # Phase 5: Deviance (server-side, on ORIGINAL-scale data)
  # ===========================================================================
  # Deviance is computed on the coordinator (label server) using stored
  # eta values. The client relays the final encrypted etas from non-label
  # servers so the coordinator can reconstruct eta_total server-side.
  # The client NEVER sees the n-length eta_total vector.

  if (use_secure_agg) {
    # Secure deviance already computed in the Ring63 loop via Beaver dot-product.
    # No individual η values are revealed. Only the scalar Σ(mu-y)² is returned.
    if (!is.null(k3_result$deviance)) {
      if (verbose) message(sprintf("\n[Phase 5] Secure deviance (Beaver): %.4f", k3_result$deviance))
    }
  } else if (use_k2_beaver) {
    # K=2: secure deviance already computed in loop
    if (!is.null(loop_result$deviance)) {
      if (verbose) message(sprintf("\n[Phase 5] Secure deviance (K=2 Beaver): %.4f", loop_result$deviance))
    }
  }

  # Deviance: from secure Beaver Σr² (available in both K=2 and K≥3)
  deviance <- NA; null_deviance <- NA
  if (use_secure_agg && exists("k3_result") && !is.null(k3_result$deviance)) {
    deviance <- k3_result$deviance
  } else if (use_k2_beaver && exists("loop_result") && !is.null(loop_result$deviance)) {
    deviance <- loop_result$deviance
  }
  # Gaussian deviance is computed in standardized space — destandardize
  if (family == "gaussian" && !is.null(y_sd) && !is.na(deviance)) {
    deviance <- deviance * y_sd^2
  }
  # Null deviance: for Gaussian = Σ(y-ȳ)² = (n-1)*var(y). In std space, var(y_std)=1.
  null_deviance <- if (family == "gaussian") (n_obs - 1) * (y_sd %||% 1)^2 else NA
  pseudo_r2 <- if (!is.na(null_deviance) && null_deviance > 0) 1 - (deviance / null_deviance) else NA
  aic <- if (!is.na(deviance)) deviance + 2 * n_vars_total else NA

  if (verbose)
    message(sprintf("\nDeviance: %.4f, Null deviance: %.4f, Pseudo R2: %.4f",
                    deviance, null_deviance, pseudo_r2))

  # ===========================================================================
  # Standard Errors + P-values (exact, via finite-difference Hessian)
  # Computed in STANDARDIZED space, then destandardized via Jacobian.
  # ===========================================================================
  inv_H <- NULL
  if (use_secure_agg && exists("k3_result")) inv_H <- k3_result$inv_hessian
  if (use_k2_beaver && exists("loop_result")) inv_H <- loop_result$inv_hessian

  std_errors <- rep(NA, n_vars_total)
  z_values <- rep(NA, n_vars_total)
  p_values <- rep(NA, n_vars_total)
  covariance <- NULL

  if (!is.null(inv_H) && !is.null(attr(inv_H, "raw_hessian"))) {
    H_raw <- attr(inv_H, "raw_hessian")
    # Fisher = n × (Hessian - λI) where Hessian = X_std^T W X_std / n + λI
    H_adj <- H_raw - lambda * diag(nrow(H_raw))
    fisher_std <- n_obs * H_adj
    cov_std <- tryCatch(solve(fisher_std), error = function(e) NULL)

    if (!is.null(cov_std)) {
      # Destandardize: construct Jacobian J where θ_orig = J × θ_std + const
      p_feat <- length(all_x_sds)
      J <- diag(p_feat + 1)  # (intercept + features)

      if (standardize_y && !is.null(y_sd)) {
        # Gaussian: β_orig_j = β_std_j × y_sd / x_sd_j
        for (jj in seq_len(p_feat)) {
          J[jj + 1, jj + 1] <- y_sd / all_x_sds[jj]
          J[1, jj + 1] <- -y_sd * all_x_means[jj] / all_x_sds[jj]
        }
        J[1, 1] <- y_sd
      } else {
        # Binomial/Poisson: β_orig_j = β_std_j / x_sd_j
        for (jj in seq_len(p_feat)) {
          J[jj + 1, jj + 1] <- 1.0 / all_x_sds[jj]
          J[1, jj + 1] <- -all_x_means[jj] / all_x_sds[jj]
        }
      }

      # Transform covariance to original space
      cov_orig <- J %*% cov_std %*% t(J)
      dimnames(cov_orig) <- list(names(all_coefs), names(all_coefs))
      se_orig <- sqrt(pmax(diag(cov_orig), 0))
      std_errors <- se_orig
      names(std_errors) <- names(all_coefs)

      z_values <- all_coefs / std_errors
      z_values[!is.finite(z_values)] <- NA
      p_values <- 2 * stats::pnorm(-abs(z_values))
      names(z_values) <- names(p_values) <- names(all_coefs)
    }
    # Expose the full covariance matrix (original scale) so downstream
    # client-side inference (multi-coef Wald, GEE sandwich, CI on linear
    # contrasts) can reuse it without another MPC round.
    if (exists("cov_orig", inherits = FALSE)) covariance <- cov_orig

    if (verbose && any(!is.na(std_errors))) {
      message("\nCoefficients:")
      message(sprintf("  %-15s %10s %10s %10s %10s", "", "Estimate", "Std.Error", "z value", "Pr(>|z|)"))
      for (nm in names(all_coefs)) {
        sig <- if (!is.na(p_values[nm]) && p_values[nm] < 0.001) "***"
               else if (!is.na(p_values[nm]) && p_values[nm] < 0.01) "**"
               else if (!is.na(p_values[nm]) && p_values[nm] < 0.05) "*"
               else ""
        message(sprintf("  %-15s %10.4f %10.4f %10.3f %10.6f %s",
          nm, all_coefs[nm], std_errors[nm], z_values[nm], p_values[nm], sig))
      }
    }
  }

  # ===========================================================================
  # Assemble Result
  # ===========================================================================
  result <- list(
    coefficients = all_coefs,
    std_errors = std_errors,
    z_values = z_values,
    p_values = p_values,
    covariance = covariance,
    iterations = final_iter,
    converged = converged,
    family = family,
    n_obs = n_obs,
    n_vars = n_vars_total,
    lambda = lambda,
    deviance = deviance,
    null_deviance = null_deviance,
    pseudo_r2 = pseudo_r2,
    aic = aic,
    y_server = y_server,
    eta_privacy = eta_privacy,
    call = call_matched
  )

  # Cleanup handled by on.exit()

  class(result) <- c("ds.glm", "list")
  return(result)
}

#' @title Print Method for ds.glm Objects
#' @description Prints a summary of GLM results.
#' @param x A ds.glm object
#' @param ... Additional arguments (ignored)
#' @export
print.ds.glm <- function(x, ...) {
  cat("\nVertically Partitioned GLM (Block Coordinate Descent)\n")
  cat("=======================================================\n\n")

  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat("Family:", x$family, "\n")
  cat("Observations:", x$n_obs, "\n")
  cat("Predictors:", x$n_vars, "\n")
  cat("Regularization (lambda):", x$lambda, "\n")
  if (!is.null(x$y_server))
    cat("Label server:", x$y_server, "\n")
  cat("Iterations:", x$iterations, "\n")
  cat("Converged:", x$converged, "\n\n")

  cat("Coefficients:\n")
  print(round(x$coefficients, 6))

  invisible(x)
}

#' @title Summary Method for ds.glm Objects
#' @description Prints detailed summary including deviance and fit statistics.
#' @param object A ds.glm object
#' @param ... Additional arguments (ignored)
#' @export
summary.ds.glm <- function(object, ...) {
  cat("\nVertically Partitioned GLM - Summary\n")
  cat("====================================\n\n")

  cat("Call:\n")
  print(object$call)
  cat("\n")

  cat("Family:", object$family, "\n")
  cat("Observations:", object$n_obs, "\n")
  cat("Predictors:", object$n_vars, "\n")
  cat("Regularization (lambda):", object$lambda, "\n")
  if (!is.null(object$y_server))
    cat("Label server:", object$y_server, "\n")
  cat("\n")

  cat("Convergence:\n")
  cat("  Iterations:", object$iterations, "\n")
  cat("  Converged:", object$converged, "\n\n")

  cat("Deviance:\n")
  cat("  Null deviance:    ", sprintf("%.4f", object$null_deviance),
      " on", object$n_obs - 1, "degrees of freedom\n")
  cat("  Residual deviance:", sprintf("%.4f", object$deviance),
      " on", object$n_obs - object$n_vars, "degrees of freedom\n\n")

  cat("Model Fit:\n")
  cat("  Pseudo R-squared (McFadden):", sprintf("%.4f", object$pseudo_r2), "\n")
  cat("  AIC:", sprintf("%.4f", object$aic), "\n\n")

  cat("Coefficients:\n")
  coef_df <- data.frame(Estimate = object$coefficients)
  print(round(coef_df, 6))

  invisible(object)
}

#' @title Coefficients Method for ds.glm Objects
#' @description Extract coefficients from a ds.glm object.
#' @param object A ds.glm object
#' @param ... Additional arguments (ignored)
#' @return Named numeric vector of coefficients
#' @export
coef.ds.glm <- function(object, ...) {
  object$coefficients
}
