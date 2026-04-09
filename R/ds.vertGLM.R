#' @title Generalized Linear Model for Vertically Partitioned Data
#' @description Client-side function that fits a Generalized Linear Model
#'   across vertically partitioned data using Block Coordinate Descent with
#'   encrypted labels. The response variable y only needs to exist on ONE
#'   server (the "label server"). Non-label servers compute gradient updates
#'   using y encrypted under the MHE collective public key, and only the
#'   aggregated p_k-length gradient is revealed via threshold decryption.
#'
#' @param data_name Character string. Name of the (aligned) data frame on
#'   each server.
#' @param y_var Character string. Name of the response variable (must exist
#'   on the label server specified by \code{y_server}).
#' @param x_vars A named list where each name corresponds to a server name
#'   and each element is a character vector of predictor variable names
#'   from that server.
#' @param y_server Character string. Name of the server holding the response
#'   variable. This server uses plaintext IRLS; all other servers use the
#'   encrypted gradient protocol.
#' @param family Character string. GLM family: "gaussian", "binomial",
#'   or "poisson". Default is "gaussian".
#' @param max_iter Integer. Maximum number of BCD iterations. Default is 100.
#' @param tol Numeric. Convergence tolerance on coefficient change.
#'   Default is 1e-4 (accounts for CKKS approximation noise).
#' @param lambda Numeric. L2 regularization parameter. Default is 1e-4.
#' @param log_n Integer. CKKS ring dimension parameter (12, 13, or 14).
#'   Default is 12 (2048 slots, supports up to 2048 observations).
#' @param log_scale Integer. CKKS scale parameter. Default is 40.
#' @param verbose Logical. Print progress messages. Default is TRUE.
#' @param datasources DataSHIELD connection object or list of connections.
#'   If NULL, uses all available connections.
#' @param eta_privacy Character string. Privacy mode for linear predictor
#'   aggregation: \code{"auto"} (default) selects the appropriate mode based
#'   on K and family — secure aggregation for K>=3, strict HE-Link for K=2
#'   Gaussian, and policy-dependent for K=2 binomial/Poisson (strict by
#'   default, pragmatic if server admin enables it). Other values:
#'   \code{"secure_agg"} (pairwise PRG masks, K>=3 only),
#'   \code{"he_link"} (encrypted link function, K=2, requires log_n>=14 for
#'   non-Gaussian), \code{"k2_mpc"} (Gauss-Seidel BCD-IRLS, K=2 pragmatic),
#'   \code{"transport"} (standard secure routing).
#' @param topology Character string. Seed derivation topology for secure
#'   aggregation: \code{"pairwise"} (default, O(K-1) seeds per server) or
#'   \code{"ring"} (O(2) seeds per server for K>=4). For K=3, ring and
#'   pairwise are identical. Only used when \code{eta_privacy = "secure_agg"}.
#' @param reuse_mhe Logical. If TRUE, reuse cached MHE context (keys, parameters)
#'   from a previous analysis sharing the same CKKS parameters. Saves key
#'   generation time but uses fresh transport keys for forward secrecy.
#'   Default is FALSE.
#'
#' @return A list with class "ds.glm" containing:
#'   \itemize{
#'     \item \code{coefficients}: Named vector of coefficient estimates
#'       (on original scale, including intercept)
#'     \item \code{iterations}: Number of iterations until convergence
#'     \item \code{converged}: Logical indicating convergence
#'     \item \code{family}: Family used
#'     \item \code{n_obs}: Number of observations
#'     \item \code{n_vars}: Number of predictor variables (including intercept)
#'     \item \code{lambda}: Regularization parameter used
#'     \item \code{deviance}: Residual deviance of the fitted model
#'     \item \code{null_deviance}: Null deviance (intercept-only model)
#'     \item \code{pseudo_r2}: McFadden's pseudo R-squared
#'     \item \code{aic}: Akaike Information Criterion
#'     \item \code{y_server}: Name of the label server
#'     \item \code{call}: The matched call
#'   }
#'
#' @details
#' \subsection{Feature Standardization}{
#' Features are automatically standardized (centered and scaled) on each
#' server before BCD to ensure fast convergence. For Gaussian family, the
#' response is also standardized. Coefficients are transformed back to the
#' original scale after convergence, and an intercept is computed.
#' }
#'
#' \subsection{Encrypted-Label BCD-IRLS Protocol}{
#' The response variable y resides on a single "label server". Non-label
#' servers never see y in plaintext. The protocol proceeds as:
#'
#' \enumerate{
#'   \item \strong{MHE Key Setup}: All servers generate key shares and
#'     combine them into a Collective Public Key (CPK) with Galois keys.
#'   \item \strong{Standardize}: Each server standardizes its features.
#'   \item \strong{Encrypt y}: The label server encrypts (standardized) y
#'     under the CPK and distributes the ciphertext to non-label servers.
#'   \item \strong{BCD Loop}: For each iteration, each server updates
#'     its block of coefficients on the standardized scale.
#'   \item \strong{Unstandardize}: Coefficients are transformed back to the
#'     original scale and an intercept is computed.
#'   \item \strong{Deviance}: Computed on the label server using
#'     plaintext y and the final linear predictor (original scale).
#' }
#' }
#'
#' \subsection{K=2 Encrypted Link-Function Mode (HE-Link)}{
#' With exactly two servers, secure aggregation cannot hide the non-label
#' server's eta contribution. dsVert automatically switches to encrypted
#' link-function mode, which evaluates the inverse link under CKKS:
#' \itemize{
#'   \item \strong{Gaussian}: Identity link (mu = eta); no polynomial
#'     evaluation needed, uses log_n = 12 and exact Newton block solve.
#'   \item \strong{Binomial}: Degree-7 sigmoid polynomial on [-8, 8];
#'     requires log_n >= 14 for multiplicative depth.
#'   \item \strong{Poisson}: Degree-7 Chebyshev exp polynomial on [-3, 3]
#'     with per-server eta clipping; requires log_n >= 14.
#' }
#' Binomial and Poisson use gradient descent with a warm-start step-size
#' schedule. The converged coefficient error is dominated by the polynomial
#' approximation gap rather than optimizer choice. Heavy CKKS operations
#' at log_n = 14 are dispatched asynchronously with short-poll result
#' retrieval to prevent reverse-proxy timeouts.
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
#' x_vars <- list(
#'   server1 = c("age", "bmi"),
#'   server2 = c("glucose"),
#'   server3 = c("cholesterol", "heart_rate")
#' )
#'
#' # Gaussian GLM (bp on server2)
#' model <- ds.vertGLM("D_aligned", "bp", x_vars,
#'                      y_server = "server2", family = "gaussian")
#' print(model)
#' }
#'
#' @importFrom DSI datashield.aggregate datashield.connections_find
#' @export
ds.vertGLM <- function(data_name, y_var, x_vars, y_server = NULL,
                       family = "gaussian", max_iter = 100, tol = 1e-4,
                       lambda = 1e-4, log_n = 12, log_scale = 40,
                       verbose = TRUE, datasources = NULL,
                       eta_privacy = "auto",
                       topology = "pairwise",
                       reuse_mhe = FALSE) {
  call_matched <- match.call()

  # ===========================================================================
  # Input Validation
  # ===========================================================================
  if (!is.character(data_name) || length(data_name) != 1)
    stop("data_name must be a single character string", call. = FALSE)
  if (!is.character(y_var) || length(y_var) != 1)
    stop("y_var must be a single character string", call. = FALSE)
  if (!is.list(x_vars) || is.null(names(x_vars)))
    stop("x_vars must be a named list mapping server names to variable vectors",
         call. = FALSE)
  if (!family %in% c("gaussian", "binomial", "poisson"))
    stop("family must be 'gaussian', 'binomial', or 'poisson'",
         call. = FALSE)
  if (is.null(y_server))
    stop("y_server must be specified: the server holding '", y_var, "'",
         call. = FALSE)
  if (!y_server %in% names(x_vars))
    stop("y_server '", y_server, "' must be in x_vars", call. = FALSE)

  # Validate topology parameter
  if (!topology %in% c("pairwise", "ring"))
    stop("topology must be 'pairwise' or 'ring'", call. = FALSE)

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

  # Ring63 DCF: NO CKKS infrastructure needed (no CPK, Galois, RLK, encrypted y)
  generate_rlk <- FALSE
  # Ring63 for ALL families (Gaussian=identity link, others=DCF spline)
  skip_ckks <- use_secure_agg

  # Adaptive log_n for large n: max_slots = 2^(log_n-1)
  # n ≤ 4096 → log_n=13, n ≤ 8192 → log_n=14, n ≤ 16384 → log_n=15
  # Note: this is set here before n_obs is known; adjusted later after getObsCountDS

  # ===========================================================================
  # Setup
  # ===========================================================================
  if (is.null(datasources))
    datasources <- DSI::datashield.connections_find()
  if (length(datasources) == 0)
    stop("No DataSHIELD connections found", call. = FALSE)

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
  # even if .glm_mhe_setup() fails before returning the closure.
  on.exit({
    for (.srv in server_list) {
      .ci <- which(server_names == .srv)
      tryCatch(
        DSI::datashield.aggregate(conns = datasources[.ci],
          expr = call("mheCleanupDS", session_id = session_id)),
        error = function(e) NULL)
      tryCatch(
        DSI::datashield.aggregate(conns = datasources[.ci],
          expr = call("mheGcDS")),
        error = function(e) NULL)
    }
    .dsvert_reset_chunk_size()
  }, add = TRUE)

  # ===========================================================================
  # Phase 0-2: MHE key setup, standardize, encrypt y
  #   (delegated to .glm_mhe_setup in ds.vertGLM.setup.R)
  # ===========================================================================
  setup <- .glm_mhe_setup(
    datasources      = datasources,
    server_names     = server_names,
    server_list      = server_list,
    non_label_servers = non_label_servers,
    y_server         = y_server,
    y_var            = y_var,
    x_vars           = x_vars,
    data_name        = data_name,
    family           = family,
    n_obs            = n_obs,
    log_n            = log_n,
    log_scale        = log_scale,
    generate_rlk     = generate_rlk,
    use_secure_agg   = use_secure_agg,
    use_k2_beaver    = use_k2_beaver,
    reuse_mhe        = reuse_mhe,
    session_id       = session_id,
    verbose          = verbose,
    skip_ckks        = skip_ckks
  )

  # Unpack setup results
  transport_pks <- setup$transport_pks
  cpk           <- setup$cpk
  x_means       <- setup$x_means
  x_sds         <- setup$x_sds
  y_mean        <- setup$y_mean
  y_sd          <- setup$y_sd
  std_data      <- setup$std_data
  standardize_y <- setup$standardize_y
  .dsAgg        <- setup$.dsAgg
  .sendBlob     <- setup$.sendBlob

  # ===========================================================================
  # Phase 3: BCD Iterations (on standardized scale)
  # ===========================================================================
  label_intercept <- !standardize_y

  # Coordinator = label server (has y, runs IRLS)
  coordinator <- y_server
  coordinator_conn <- which(server_names == coordinator)
  coordinator_pk <- if (length(transport_pks) > 0) transport_pks[[coordinator]] else NULL

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
      log_n = log_n, log_scale = log_scale, session_id = session_id,
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
    null_deviance <- n_obs
    if (!is.null(y_sd)) null_deviance <- n_obs * y_sd^2
  } else if (use_k2_beaver && exists("loop_result") && !is.null(loop_result$deviance)) {
    deviance <- loop_result$deviance
    null_deviance <- n_obs
    if (!is.null(y_sd)) null_deviance <- n_obs * y_sd^2
  }
  pseudo_r2 <- 1 - (deviance / null_deviance)
  aic <- deviance + 2 * n_vars_total

  if (verbose)
    message(sprintf("\nDeviance: %.4f, Null deviance: %.4f, Pseudo R2: %.4f",
                    deviance, null_deviance, pseudo_r2))

  # ===========================================================================
  # Assemble Result
  # ===========================================================================
  result <- list(
    coefficients = all_coefs,
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
    topology = topology,
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
