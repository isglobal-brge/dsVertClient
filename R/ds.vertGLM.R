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
  if (!family %in% c("gaussian", "binomial", "poisson"))
    stop("family must be 'gaussian', 'binomial', or 'poisson'",
         call. = FALSE)
  if (is.null(y_server))
    stop("y_server must be specified: the server holding '", y_var, "'",
         call. = FALSE)
  if (!y_server %in% names(x_vars))
    stop("y_server '", y_server, "' must be in x_vars", call. = FALSE)

  # ===========================================================================
  # Helpers
  # ===========================================================================

  chunk_size <- 10000  # DataSHIELD parser-safe chunk size

  # Send CT chunks to a server (reusable helper)
  .sendCTChunks <- function(ct_str, conn_idx) {
    n_chars_ct <- nchar(ct_str)
    n_ct_chunks <- ceiling(n_chars_ct / chunk_size)
    for (ch in seq_len(n_ct_chunks)) {
      start <- (ch - 1) * chunk_size + 1
      end <- min(ch * chunk_size, n_chars_ct)
      DSI::datashield.aggregate(
        conns = datasources[conn_idx],
        expr = call("mheStoreCTChunkDS",
                    chunk_index = as.integer(ch),
                    chunk = substr(ct_str, start, end))
      )
    }
    n_ct_chunks
  }

  # Relay a wrapped share to the fusion server in chunks
  .sendWrappedShare <- function(wrapped_share, party_id, fusion_conn_idx) {
    n_chars <- nchar(wrapped_share)
    n_chunks <- ceiling(n_chars / chunk_size)
    for (ch in seq_len(n_chunks)) {
      start <- (ch - 1) * chunk_size + 1
      end <- min(ch * chunk_size, n_chars)
      DSI::datashield.aggregate(
        conns = datasources[fusion_conn_idx],
        expr = call("mheStoreWrappedShareDS",
                    party_id = as.integer(party_id),
                    share_data = substr(wrapped_share, start, end))
      )
    }
  }

  # Send a blob to a server, chunking if necessary
  .sendBlob <- function(blob, key, conn_idx) {
    n_chars <- nchar(blob)
    n_chunks <- ceiling(n_chars / chunk_size)
    if (n_chunks == 1) {
      DSI::datashield.aggregate(
        conns = datasources[conn_idx],
        expr = call("mheStoreBlobDS", key = key, chunk = blob)
      )
    } else {
      for (ch in seq_len(n_chunks)) {
        start <- (ch - 1) * chunk_size + 1
        end <- min(ch * chunk_size, n_chars)
        DSI::datashield.aggregate(
          conns = datasources[conn_idx],
          expr = call("mheStoreBlobDS",
                      key = key,
                      chunk = substr(blob, start, end),
                      chunk_index = as.integer(ch),
                      n_chunks = as.integer(n_chunks))
        )
      }
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
  # Phase 0: MHE Key Setup + Transport Keys (only needed if non-label servers exist)
  # ===========================================================================
  transport_pks <- list()  # Collect transport PKs for secure routing + share-wrapping

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
    transport_pks[[server_list[1]]] <- result0$transport_pk

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
      transport_pks[[server]] <- result$transport_pk
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

    # Distribute transport keys for share-wrapping and secure routing
    fusion_server <- server_list[1]  # Party 0
    fusion_conn_idx <- which(server_names == fusion_server)
    non_fusion_servers <- setdiff(server_list, fusion_server)

    for (server in server_list) {
      conn_idx <- which(server_names == server)
      tk_map <- list(fusion = transport_pks[[fusion_server]])
      for (s in server_list) {
        tk_map[[s]] <- transport_pks[[s]]
      }
      DSI::datashield.aggregate(
        conns = datasources[conn_idx],
        expr = call("mheStoreTransportKeysDS", transport_keys = tk_map)
      )
    }

    if (verbose) message("  Key setup + transport keys complete")
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
  # Phase 3: BCD Iterations with Secure Routing (on standardized scale)
  # ===========================================================================
  # Secure Routing Protocol: individual-level vectors (eta, mu, w, v) are
  # transport-encrypted end-to-end between the coordinator (label server)
  # and non-label servers. The client only handles:
  #   - Betas: p_k-length coefficient vectors (safe aggregate statistics)
  #   - Opaque encrypted blobs: transport-encrypted vectors the client cannot read
  #
  # BCD variant: Jacobi-style. All non-label servers receive the same mu/w/v
  # snapshot per iteration (computed by coordinator at iteration start).
  label_intercept <- !standardize_y

  # Coordinator = label server (has y, runs IRLS)
  coordinator <- y_server
  coordinator_conn <- which(server_names == coordinator)
  coordinator_pk <- if (length(transport_pks) > 0) transport_pks[[coordinator]] else NULL

  # Non-label server PKs for the coordinator to encrypt (mu, w, v) under
  non_label_pk_map <- list()
  for (s in non_label_servers) {
    non_label_pk_map[[s]] <- transport_pks[[s]]
  }

  betas <- list()
  for (server in server_list) {
    p <- length(x_vars[[server]])
    if (server == coordinator && label_intercept) p <- p + 1
    betas[[server]] <- rep(0, p)
  }

  # Encrypted eta blobs (opaque to client): server_name -> base64url string
  encrypted_etas <- list()

  converged <- FALSE
  final_iter <- 0

  if (verbose) message("\n[Phase 3] BCD iterations (secure routing)...")

  for (iter in seq_len(max_iter)) {
    betas_old <- betas
    max_diff <- 0

    # --- Phase 3a: Coordinator step (label server) ---
    # The coordinator receives encrypted etas from non-label servers,
    # decrypts them, runs IRLS, computes mu/w/v, and encrypts (mu, w, v)
    # under each non-label server's transport PK.

    # Send encrypted eta blobs to coordinator (if any, from previous iteration)
    eta_blob_keys <- character(0)
    if (length(encrypted_etas) > 0) {
      for (s in names(encrypted_etas)) {
        blob <- encrypted_etas[[s]]
        if (!is.null(blob) && nzchar(blob)) {
          key <- paste0("eta_", s)
          .sendBlob(blob, key, coordinator_conn)
          eta_blob_keys <- c(eta_blob_keys, key)
        }
      }
    }

    coord_result <- DSI::datashield.aggregate(
      conns = datasources[coordinator_conn],
      expr = call("glmCoordinatorStepDS",
                  data_name = std_data,
                  y_var = y_var,
                  x_vars = x_vars[[coordinator]],
                  encrypted_eta_blobs = NULL,
                  eta_blob_keys = if (length(eta_blob_keys) > 0) eta_blob_keys else NULL,
                  non_label_pks = non_label_pk_map,
                  family = family,
                  beta_current = betas[[coordinator]],
                  lambda = lambda,
                  intercept = label_intercept,
                  n_obs = as.integer(n_obs))
    )
    if (is.list(coord_result) && length(coord_result) == 1)
      coord_result <- coord_result[[1]]

    betas[[coordinator]] <- coord_result$beta
    mwv_blobs <- coord_result$encrypted_blobs  # Named list: server -> blob

    diff_coord <- sum(abs(betas[[coordinator]] - betas_old[[coordinator]]))
    max_diff <- max(max_diff, diff_coord)

    # --- Phase 3b: Non-label servers (secure gradient + share-wrapped decrypt + block solve) ---
    encrypted_etas <- list()  # Reset for this iteration

    for (server in non_label_servers) {
      conn_idx <- which(server_names == server)
      vars <- x_vars[[server]]
      p_k <- length(vars)

      # Step 1: Send encrypted (mu, w, v) blob to non-label server
      mwv_blob <- mwv_blobs[[server]]
      .sendBlob(mwv_blob, "mwv", conn_idx)

      # Step 2: Server computes encrypted gradient using decrypted mu/v + ct_y
      grad_result <- DSI::datashield.aggregate(
        conns = datasources[conn_idx],
        expr = call("glmSecureGradientDS",
                    data_name = std_data,
                    x_vars = vars,
                    encrypted_mwv = NULL,
                    num_obs = as.integer(n_obs))
      )
      if (is.list(grad_result) && length(grad_result) == 1)
        grad_result <- grad_result[[1]]

      enc_gradients <- grad_result$encrypted_gradients
      ct_hashes <- grad_result$ct_hashes

      # Protocol Firewall: authorize gradient ciphertexts on non-producing servers
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

      # Step 3: Share-wrapped threshold decryption of each gradient component
      gradient <- numeric(p_k)
      for (j in seq_len(p_k)) {
        ct_grad <- enc_gradients[j]

        # Collect wrapped shares from non-fusion servers
        for (nf_server in non_fusion_servers) {
          nf_conn <- which(server_names == nf_server)
          nf_party_id <- which(server_list == nf_server) - 1

          n_ct_chunks <- .sendCTChunks(ct_grad, nf_conn)
          pd <- DSI::datashield.aggregate(
            conns = datasources[nf_conn],
            expr = call("mhePartialDecryptWrappedDS",
                        n_chunks = as.integer(n_ct_chunks))
          )
          if (is.list(pd)) pd <- pd[[1]]
          .sendWrappedShare(pd$wrapped_share, nf_party_id, fusion_conn_idx)
        }

        # Send CT to fusion server and fuse
        n_ct_chunks <- .sendCTChunks(ct_grad, fusion_conn_idx)
        fuse_result <- DSI::datashield.aggregate(
          conns = datasources[fusion_conn_idx],
          expr = call("mheFuseServerDS",
                      n_parties = as.integer(length(server_list)),
                      n_ct_chunks = as.integer(n_ct_chunks),
                      num_slots = 0L)
        )
        if (is.list(fuse_result)) fuse_result <- fuse_result[[1]]
        gradient[j] <- fuse_result$value
      }

      # Step 4: BCD block solve with encrypted eta output
      # Send mwv blob again for the block solve (w needed for IRLS weights)
      .sendBlob(mwv_blob, "mwv", conn_idx)

      solve_result <- DSI::datashield.aggregate(
        conns = datasources[conn_idx],
        expr = call("glmSecureBlockSolveDS",
                    data_name = std_data,
                    x_vars = vars,
                    encrypted_mwv = NULL,
                    beta_current = betas[[server]],
                    gradient = gradient,
                    lambda = lambda,
                    coordinator_pk = coordinator_pk)
      )
      if (is.list(solve_result) && length(solve_result) == 1)
        solve_result <- solve_result[[1]]

      betas[[server]] <- solve_result$beta
      encrypted_etas[[server]] <- solve_result$encrypted_eta

      diff_server <- sum(abs(betas[[server]] - betas_old[[server]]))
      max_diff <- max(max_diff, diff_server)
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
  all_coefs_std <- numeric()
  all_x_means <- numeric()
  all_x_sds <- numeric()
  all_names <- character()
  beta_0_from_label <- 0

  for (server in server_list) {
    server_beta <- betas[[server]]
    if (server == coordinator && label_intercept) {
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

  # Send final encrypted etas to coordinator
  eta_blob_keys_final <- character(0)
  if (length(encrypted_etas) > 0) {
    for (s in names(encrypted_etas)) {
      blob <- encrypted_etas[[s]]
      if (!is.null(blob) && nzchar(blob)) {
        key <- paste0("eta_", s)
        .sendBlob(blob, key, coordinator_conn)
        eta_blob_keys_final <- c(eta_blob_keys_final, key)
      }
    }
  }

  deviance_result <- DSI::datashield.aggregate(
    conns = datasources[coordinator_conn],
    expr = call("glmSecureDevianceDS",
                data_name = data_name,
                y_var = y_var,
                encrypted_eta_blobs = NULL,
                eta_blob_keys = if (length(eta_blob_keys_final) > 0) eta_blob_keys_final else NULL,
                family = family,
                y_sd = if (!is.null(y_sd)) y_sd else NULL,
                y_mean = if (!is.null(y_mean)) y_mean else NULL)
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
