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
#'   "poisson", "Gamma", or "inverse.gaussian". Default is "gaussian".
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
                       verbose = TRUE, datasources = NULL) {
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
  if (!family %in% c("gaussian", "binomial", "poisson", "Gamma", "inverse.gaussian"))
    stop("family must be 'gaussian', 'binomial', 'poisson', 'Gamma', or 'inverse.gaussian'",
         call. = FALSE)
  if (is.null(y_server))
    stop("y_server must be specified: the server holding '", y_var, "'",
         call. = FALSE)
  if (!y_server %in% names(x_vars))
    stop("y_server '", y_server, "' must be in x_vars", call. = FALSE)

  # ===========================================================================
  # Helpers
  # ===========================================================================

  .b64url_to_b64 <- function(x) {
    x <- gsub("-", "+", x, fixed = TRUE)
    x <- gsub("_", "/", x, fixed = TRUE)
    pad <- (4 - nchar(x) %% 4) %% 4
    if (pad > 0) x <- paste0(x, paste(rep("=", pad), collapse = ""))
    x
  }

  .findMheTool <- function() {
    os <- .Platform$OS.type
    binary_name <- if (os == "windows") "mhe-tool.exe" else "mhe-tool"
    arch <- Sys.info()["machine"]
    if (Sys.info()["sysname"] == "Darwin") {
      subdir <- if (arch == "arm64") "darwin-arm64" else "darwin-amd64"
    } else if (os == "windows") {
      subdir <- "windows-amd64"
    } else {
      subdir <- "linux-amd64"
    }
    bin_path <- system.file("bin", subdir, binary_name, package = "dsVertClient")
    if (bin_path != "" && file.exists(bin_path)) return(bin_path)
    bin_path <- system.file("bin", subdir, binary_name, package = "dsVert")
    if (bin_path != "" && file.exists(bin_path)) return(bin_path)
    env_path <- Sys.getenv("DSVERT_MHE_TOOL")
    if (env_path != "" && file.exists(env_path)) return(env_path)
    stop("mhe-tool binary not found for platform ", subdir, call. = FALSE)
  }

  .mheFuseLocal <- function(ciphertext, decryption_shares, log_n, log_scale,
                             num_slots = 0) {
    bin_path <- .findMheTool()
    ct_std <- .b64url_to_b64(ciphertext)
    shares_std <- sapply(decryption_shares, .b64url_to_b64, USE.NAMES = FALSE)
    input <- list(
      ciphertext = ct_std,
      decryption_shares = as.list(shares_std),
      num_slots = as.integer(num_slots),
      log_n = as.integer(log_n),
      log_scale = as.integer(log_scale)
    )
    input_file <- tempfile(fileext = ".json")
    output_file <- tempfile(fileext = ".json")
    on.exit({
      if (file.exists(input_file)) unlink(input_file)
      if (file.exists(output_file)) unlink(output_file)
    })
    jsonlite::write_json(input, input_file, auto_unbox = TRUE)
    status <- system2(bin_path, "mhe-fuse",
                      stdin = input_file, stdout = output_file,
                      stderr = output_file)
    output <- jsonlite::read_json(output_file, simplifyVector = TRUE)
    if (!is.null(output$error) && nzchar(output$error))
      stop("mhe-fuse error: ", output$error, call. = FALSE)
    output$value
  }

  # --- GLM link/variance helpers ---
  # These compute the mean, IRLS weights, and variance correction factor
  # on the client side. The client broadcasts mu and w to all servers at
  # each BCD iteration so they can compute their local contributions.

  .compute_mu <- function(eta, family) {
    # Inverse link: converts linear predictor eta to mean mu.
    # Clipping prevents numerical overflow in exp() and division by 0.
    # - eta in [-20, 20] keeps exp() in [~2e-9, ~5e8]
    # - mu clamped to [1e-10, 1-1e-10] for binomial avoids log(0)
    if (family == "gaussian") {
      eta
    } else if (family == "binomial") {
      eta <- pmax(pmin(eta, 20), -20)
      pmax(pmin(1 / (1 + exp(-eta)), 1 - 1e-10), 1e-10)
    } else {
      pmax(exp(pmin(eta, 20)), 1e-10)
    }
  }

  .compute_w <- function(mu, family) {
    # IRLS weights: w = 1 / (V(mu) * (d_eta/d_mu)^2)
    # For canonical links, simplifies to:
    #   gaussian: w = 1, binomial: w = mu(1-mu), poisson: w = mu
    # For non-canonical links with log link:
    #   Gamma (V=mu^2): w = 1/(mu^2 * 1/mu^2) = 1
    #   inverse.gaussian (V=mu^3): w = 1/(mu^3 * 1/mu^2) = 1/mu
    if (family == "gaussian") {
      rep(1, length(mu))
    } else if (family == "binomial") {
      mu * (1 - mu)
    } else if (family == "poisson") {
      mu
    } else if (family == "Gamma") {
      rep(1, length(mu))
    } else if (family == "inverse.gaussian") {
      1 / mu
    }
  }

  .compute_v <- function(mu, family) {
    # Variance correction factor for the encrypted gradient. For canonical
    # links, the gradient is X^T(y - mu) and v = NULL (no correction).
    # For non-canonical links, the gradient needs: X^T * v * (y - mu)
    # where v = (d_mu/d_eta) / V(mu).
    if (family %in% c("gaussian", "binomial", "poisson")) {
      NULL
    } else if (family == "Gamma") {
      1 / mu
    } else if (family == "inverse.gaussian") {
      1 / (mu^2)
    }
  }

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
  chunk_size <- 10000

  # Get observation count
  first_conn <- which(server_names == server_list[1])
  count_result <- DSI::datashield.aggregate(
    conns = datasources[first_conn],
    expr = call("getObsCountDS", data_name)
  )
  if (is.list(count_result) && length(count_result) == 1)
    count_result <- count_result[[1]]
  n_obs <- count_result$n_obs

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

  # ===========================================================================
  # Phase 0: MHE Key Setup (only needed if non-label servers exist)
  # ===========================================================================
  if (length(non_label_servers) > 0) {
    if (verbose) message("\n[Phase 0] MHE key setup...")

    conn_idx <- which(server_names == server_list[1])
    result0 <- DSI::datashield.aggregate(
      conns = datasources[conn_idx],
      expr = call("mheInitDS",
                  party_id = 0L, crp = NULL, gkg_seed = NULL,
                  num_obs = as.integer(n_obs),
                  log_n = as.integer(log_n),
                  log_scale = as.integer(log_scale))
    )
    if (is.list(result0)) result0 <- result0[[1]]
    crp <- result0$crp
    gkg_seed <- result0$gkg_seed

    pk_shares <- list()
    gkg_shares <- list()
    pk_shares[[server_list[1]]] <- result0$public_key_share
    gkg_shares[[server_list[1]]] <- result0$galois_key_shares

    for (i in 2:length(server_list)) {
      server <- server_list[i]
      conn_idx <- which(server_names == server)
      result <- DSI::datashield.aggregate(
        conns = datasources[conn_idx],
        expr = call("mheInitDS",
                    party_id = as.integer(i - 1), crp = crp,
                    gkg_seed = gkg_seed,
                    num_obs = as.integer(n_obs),
                    log_n = as.integer(log_n),
                    log_scale = as.integer(log_scale))
      )
      if (is.list(result)) result <- result[[1]]
      pk_shares[[server]] <- result$public_key_share
      gkg_shares[[server]] <- result$galois_key_shares
    }

    conn_idx <- which(server_names == server_list[1])
    gkg_shares_ordered <- lapply(server_list, function(s) gkg_shares[[s]])
    combined <- DSI::datashield.aggregate(
      conns = datasources[conn_idx],
      expr = call("mheCombineDS",
                  public_key_shares = unlist(pk_shares, use.names = FALSE),
                  crp = crp,
                  galois_key_shares = gkg_shares_ordered,
                  gkg_seed = gkg_seed,
                  num_obs = as.integer(n_obs),
                  log_n = as.integer(log_n),
                  log_scale = as.integer(log_scale))
    )
    if (is.list(combined)) combined <- combined[[1]]
    cpk <- combined$collective_public_key
    galois_keys <- combined$galois_keys

    for (server in server_list[-1]) {
      conn_idx <- which(server_names == server)
      DSI::datashield.aggregate(
        conns = datasources[conn_idx],
        expr = call("mheStoreCPKDS", cpk = cpk, galois_keys = galois_keys)
      )
    }
    if (verbose) message("  Key setup complete")
  }

  # ===========================================================================
  # Phase 1: Standardize features (and y for Gaussian) for fast BCD
  # ===========================================================================
  if (verbose) message("\n[Phase 1] Standardizing features...")

  std_data <- paste0(data_name, "_std")
  standardize_y <- (family == "gaussian")

  x_means <- list()
  x_sds <- list()
  y_mean <- NULL
  y_sd <- NULL

  for (server in server_list) {
    conn_idx <- which(server_names == server)
    y_arg <- if (server == y_server && standardize_y) y_var else NULL

    std_result <- DSI::datashield.aggregate(
      conns = datasources[conn_idx],
      expr = call("glmStandardizeDS",
                  data_name = data_name,
                  output_name = std_data,
                  x_vars = x_vars[[server]],
                  y_var = y_arg)
    )
    if (is.list(std_result) && length(std_result) == 1)
      std_result <- std_result[[1]]

    x_means[[server]] <- std_result$x_means
    x_sds[[server]] <- std_result$x_sds

    if (!is.null(std_result$y_mean)) {
      y_mean <- std_result$y_mean
      y_sd <- std_result$y_sd
    }
  }
  if (verbose) message("  Features standardized on all servers")

  # ===========================================================================
  # Phase 2: Encrypt y and distribute (only if non-label servers exist)
  # ===========================================================================
  if (length(non_label_servers) > 0) {
    if (verbose) message("\n[Phase 2] Encrypting response variable...")

    # Encrypt STANDARDIZED y (Gaussian) or raw y (non-Gaussian)
    enc_data <- if (standardize_y) std_data else data_name

    conn_idx <- which(server_names == y_server)
    enc_result <- DSI::datashield.aggregate(
      conns = datasources[conn_idx],
      expr = call("mheEncryptRawDS",
                  data_name = enc_data, y_var = y_var)
    )
    if (is.list(enc_result)) enc_result <- enc_result[[1]]
    ct_y <- enc_result$encrypted_y
    if (verbose) message("  Encrypted y on ", y_server)

    for (server in non_label_servers) {
      conn_idx <- which(server_names == server)
      n_chars <- nchar(ct_y)
      n_chunks <- ceiling(n_chars / chunk_size)

      for (ch in seq_len(n_chunks)) {
        start <- (ch - 1) * chunk_size + 1
        end <- min(ch * chunk_size, n_chars)
        DSI::datashield.aggregate(
          conns = datasources[conn_idx],
          expr = call("mheStoreEncChunkDS",
                      col_index = 1L,
                      chunk_index = as.integer(ch),
                      chunk = substr(ct_y, start, end))
        )
      }
      DSI::datashield.aggregate(
        conns = datasources[conn_idx],
        expr = call("mheAssembleEncColumnDS",
                    col_index = 1L, n_chunks = as.integer(n_chunks))
      )
      if (verbose)
        message("  ct_y transferred to ", server, " (", n_chunks, " chunks)")
    }
  }

  # ===========================================================================
  # Phase 3: BCD Iterations (on standardized scale)
  # ===========================================================================
  # Standardization ensures all features have comparable scale, which is
  # critical for BCD convergence: without it, blocks with larger-magnitude
  # features would dominate the linear predictor and convergence would stall.
  #
  # For Gaussian family, both X and y are standardized, so the intercept
  # is implicitly zero on the standardized scale (recovered in Phase 4).
  # For non-Gaussian families, y is NOT standardized (the link function is
  # nonlinear), so the label server needs an explicit intercept column.
  label_intercept <- !standardize_y

  betas <- list()
  etas <- list()
  for (server in server_list) {
    p <- length(x_vars[[server]])
    if (server == y_server && label_intercept) p <- p + 1
    betas[[server]] <- rep(0, p)
    etas[[server]] <- rep(0, n_obs)
  }

  converged <- FALSE
  final_iter <- 0

  if (verbose) message("\n[Phase 3] BCD iterations...")

  for (iter in seq_len(max_iter)) {
    betas_old <- betas
    max_diff <- 0

    for (server in server_list) {
      conn_idx <- which(server_names == server)
      vars <- x_vars[[server]]

      if (server == y_server) {
        # LABEL SERVER: has y in plaintext, so standard IRLS is sufficient.
        # No encryption needed -- data stays on this server.
        # eta_other = sum of linear predictor contributions from other servers.
        other_servers <- setdiff(server_list, server)
        if (length(other_servers) > 0) {
          eta_other <- Reduce(`+`, etas[other_servers])
        } else {
          eta_other <- rep(0, n_obs)
        }

        result <- DSI::datashield.aggregate(
          conns = datasources[conn_idx],
          expr = call("glmPartialFitDS",
                      std_data, y_var, vars,
                      eta_other, betas[[server]],
                      family, lambda, label_intercept)
        )
        if (is.list(result) && length(result) == 1) result <- result[[1]]

        betas[[server]] <- result$beta
        etas[[server]] <- result$eta

      } else {
        # NON-LABEL SERVER: does NOT have y. Must use encrypted gradient
        # protocol. Steps:
        #   1. Client computes mu, w from current eta (broadcast to server)
        #   2. Server computes g_k = X_k^T(ct_y - mu) homomorphically
        #   3. Threshold decryption reveals the p_k-length gradient
        #   4. Server solves BCD block update using plaintext gradient

        # Compute IRLS quantities from the current aggregate linear predictor.
        # These are broadcast to the server (they contain no individual y values).
        eta_total <- Reduce(`+`, etas)
        mu <- .compute_mu(eta_total, family)
        w <- .compute_w(mu, family)
        v <- .compute_v(mu, family)

        # Step 2: Server computes encrypted gradient using stored ct_y and local X_k.
        # The result is p_k encrypted scalars (one per feature on this server).
        grad_result <- DSI::datashield.aggregate(
          conns = datasources[conn_idx],
          expr = call("mheGLMGradientDS",
                      data_name = std_data,
                      x_vars = vars,
                      mu = mu,
                      v = v,
                      num_obs = as.integer(n_obs))
        )
        if (is.list(grad_result) && length(grad_result) == 1)
          grad_result <- grad_result[[1]]

        enc_gradients <- grad_result$encrypted_gradients
        ct_hashes <- grad_result$ct_hashes
        p_k <- length(vars)

        # Protocol Firewall: authorize gradient ciphertexts on non-producing servers.
        # The producing server (current non-label server) already has them registered.
        if (!is.null(ct_hashes) && length(ct_hashes) > 0) {
          for (auth_server in server_list) {
            if (auth_server != server) {
              auth_conn <- which(server_names == auth_server)
              DSI::datashield.aggregate(
                conns = datasources[auth_conn],
                expr = call("mheAuthorizeCTDS",
                            ct_hashes = ct_hashes,
                            op_type = "glm-gradient")
              )
            }
          }
        }

        # Step 3: Threshold decrypt each gradient component.
        # Each encrypted gradient[j] is a single CKKS ciphertext. We collect
        # partial decryption shares from ALL K servers, then fuse them on the
        # client to recover the plaintext gradient scalar. This reveals only
        # the aggregate gradient (safe), not individual observations of y.
        gradient <- numeric(p_k)
        for (j in seq_len(p_k)) {
          ct_grad <- enc_gradients[j]

          partial_shares <- list()
          for (dec_server in server_list) {
            dec_conn <- which(server_names == dec_server)

            ct_str <- ct_grad
            n_ct_chars <- nchar(ct_str)
            n_ct_chunks <- ceiling(n_ct_chars / chunk_size)

            for (ch in seq_len(n_ct_chunks)) {
              start <- (ch - 1) * chunk_size + 1
              end <- min(ch * chunk_size, n_ct_chars)
              DSI::datashield.aggregate(
                conns = datasources[dec_conn],
                expr = call("mheStoreCTChunkDS",
                            chunk_index = as.integer(ch),
                            chunk = substr(ct_str, start, end))
              )
            }

            pd <- DSI::datashield.aggregate(
              conns = datasources[dec_conn],
              expr = call("mhePartialDecryptDS",
                          n_chunks = as.integer(n_ct_chunks))
            )
            if (is.list(pd)) pd <- pd[[1]]
            partial_shares[[dec_server]] <- pd$decryption_share
          }

          gradient[j] <- .mheFuseLocal(
            ct_grad,
            unlist(partial_shares, use.names = FALSE),
            log_n = log_n, log_scale = log_scale,
            num_slots = 0
          )
        }

        # Step 4: BCD block update using the decrypted gradient.
        # The server solves: beta_new = (X^TWX + λI)^{-1}(X^TWX*beta + g_k)
        # This happens on the server using plaintext X_k and the gradient.
        solve_result <- DSI::datashield.aggregate(
          conns = datasources[conn_idx],
          expr = call("glmBlockSolveDS",
                      data_name = std_data,
                      x_vars = vars,
                      w = w,
                      beta_current = betas[[server]],
                      gradient = gradient,
                      lambda = lambda)
        )
        if (is.list(solve_result) && length(solve_result) == 1)
          solve_result <- solve_result[[1]]

        betas[[server]] <- solve_result$beta
        etas[[server]] <- solve_result$eta
      }

      diff <- sum(abs(betas[[server]] - betas_old[[server]]))
      max_diff <- max(max_diff, diff)
    }

    final_iter <- iter

    if (max_diff < tol) {
      converged <- TRUE
      if (verbose)
        message(sprintf("  Converged after %d iterations (diff = %.2e)",
                        iter, max_diff))
      break
    }

    if (verbose && iter %% 5 == 0)
      message(sprintf("  Iteration %d: max diff = %.2e", iter, max_diff))
  }

  if (!converged && verbose)
    warning(sprintf("Did not converge after %d iterations (diff = %.2e)",
                    max_iter, max_diff))

  # ===========================================================================
  # Phase 4: Unstandardize coefficients
  # ===========================================================================
  # BCD produced coefficients on the standardized scale: y_std = X_std * beta_std.
  # To recover original-scale coefficients: y = X * beta_orig + beta_0.
  #
  # Since x_std = (x - x_mean) / x_sd and y_std = (y - y_mean) / y_sd:
  #   Gaussian:     beta_orig[j] = beta_std[j] * (y_sd / x_sd[j])
  #                 beta_0 = y_mean - sum(beta_orig * x_mean)
  #   Non-Gaussian: beta_orig[j] = beta_std[j] / x_sd[j]  (y not standardized)
  #                 beta_0 = beta_0_std - sum(beta_orig * x_mean)
  #                 where beta_0_std comes from the label server's intercept column

  all_coefs_std <- numeric()
  all_x_means <- numeric()
  all_x_sds <- numeric()
  all_names <- character()
  beta_0_from_label <- 0

  for (server in server_list) {
    server_beta <- betas[[server]]
    if (server == y_server && label_intercept) {
      # First element is the intercept from the label server's IRLS
      beta_0_from_label <- server_beta[1]
      server_beta <- server_beta[-1]  # remaining are feature coefficients
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
  # Phase 5: Deviance (label server only, on ORIGINAL-scale data)
  # ===========================================================================
  # Compute eta_total on original scale
  eta_total_orig <- Reduce(`+`, etas)
  if (standardize_y && !is.null(y_sd)) {
    # Transform standardized eta back to original scale
    eta_total_orig <- eta_total_orig * y_sd + y_mean
  }

  conn_idx <- which(server_names == y_server)
  deviance_result <- DSI::datashield.aggregate(
    conns = datasources[conn_idx],
    expr = call("glmDevianceDS", data_name, y_var, eta_total_orig, family)
  )
  if (is.list(deviance_result) && length(deviance_result) == 1)
    deviance_result <- deviance_result[[1]]

  deviance <- deviance_result$deviance
  null_deviance <- deviance_result$null_deviance
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
    call = call_matched
  )

  # Clean up cryptographic state on all servers
  for (server in server_list) {
    conn_idx <- which(server_names == server)
    tryCatch(
      DSI::datashield.aggregate(conns = datasources[conn_idx],
                                expr = call("mheCleanupDS")),
      error = function(e) NULL
    )
  }

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
