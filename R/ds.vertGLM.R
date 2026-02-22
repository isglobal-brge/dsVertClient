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
#'   aggregation: \code{"auto"} (default) selects HE-Link for K=2 binomial,
#'   \code{"transport"} uses standard secure routing, \code{"he_link"} forces
#'   homomorphic sigmoid evaluation (requires K=2, binomial, log_n >= 14).
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
                       verbose = TRUE, datasources = NULL,
                       eta_privacy = "auto") {
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

  # Resolve eta_privacy mode
  if (!eta_privacy %in% c("auto", "transport", "he_link", "secure_agg"))
    stop("eta_privacy must be 'auto', 'transport', 'he_link', or 'secure_agg'",
         call. = FALSE)

  n_partitions_check <- length(x_vars)
  non_label_count <- n_partitions_check - 1

  if (eta_privacy == "auto") {
    if (non_label_count >= 2) {
      eta_privacy <- "secure_agg"    # K>=3: pairwise PRG masks
    } else if (non_label_count == 1 && family %in% c("binomial", "poisson")) {
      eta_privacy <- "he_link"       # K=2 nonlinear: homomorphic sigmoid
      if (log_n < 14) {
        log_n <- 14L
      }
    } else {
      eta_privacy <- "transport"     # K=2 gaussian or K=1
    }
  }

  use_he_link <- (eta_privacy == "he_link")
  use_secure_agg <- (eta_privacy == "secure_agg")

  # Policy enforcement
  if (use_secure_agg && non_label_count < 2)
    stop("secure_agg requires >= 3 servers (>= 2 non-label). ",
         "For K=2, use 'he_link' (binomial) or 'transport'.", call. = FALSE)
  if (use_he_link) {
    if (non_label_count != 1)
      stop("HE-Link mode requires exactly 2 servers (1 label + 1 non-label). ",
           "Got ", n_partitions_check, " partitions.", call. = FALSE)
    if (!family %in% c("binomial"))
      stop("HE-Link mode currently supports binomial family only. ",
           "Got '", family, "'.", call. = FALSE)
    if (log_n < 14)
      stop("HE-Link mode requires log_n >= 14 (need 4 multiplicative levels). ",
           "Got log_n = ", log_n, ".", call. = FALSE)
  }

  # ===========================================================================
  # Helpers
  # ===========================================================================

  chunk_size <- 100000  # DataSHIELD parser-safe chunk size (100KB)

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
                    chunk = substr(ct_str, start, end),
                    session_id = session_id)
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
                    share_data = substr(wrapped_share, start, end),
                    session_id = session_id)
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
        expr = call("mheStoreBlobDS", key = key, chunk = blob,
                    session_id = session_id)
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
                      n_chunks = as.integer(n_chunks),
                      session_id = session_id)
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
  chunk_size <- 100000

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
                  log_scale = as.integer(log_scale),
                  generate_rlk = use_he_link,
                  session_id = session_id)
    )
    if (is.list(result0)) result0 <- result0[[1]]
    crp <- result0$crp
    gkg_seed <- result0$gkg_seed

    pk_shares <- list()
    gkg_shares <- list()
    rlk_r1_shares <- list()
    pk_shares[[server_list[1]]] <- result0$public_key_share
    gkg_shares[[server_list[1]]] <- result0$galois_key_shares
    transport_pks[[server_list[1]]] <- result0$transport_pk
    if (use_he_link && !is.null(result0$rlk_round1_share)) {
      rlk_r1_shares[[server_list[1]]] <- result0$rlk_round1_share
    }

    for (i in 2:length(server_list)) {
      server <- server_list[i]
      conn_idx <- which(server_names == server)

      # CRP can be several MB at log_n=14; send via blob storage
      .sendBlob(crp, "crp", conn_idx)
      .sendBlob(gkg_seed, "gkg_seed", conn_idx)

      result <- DSI::datashield.aggregate(
        conns = datasources[conn_idx],
        expr = call("mheInitDS",
                    party_id = as.integer(i - 1),
                    from_storage = TRUE,
                    num_obs = as.integer(n_obs),
                    log_n = as.integer(log_n),
                    log_scale = as.integer(log_scale),
                    generate_rlk = use_he_link,
                    session_id = session_id)
      )
      if (is.list(result)) result <- result[[1]]
      pk_shares[[server]] <- result$public_key_share
      gkg_shares[[server]] <- result$galois_key_shares
      transport_pks[[server]] <- result$transport_pk
      if (use_he_link && !is.null(result$rlk_round1_share)) {
        rlk_r1_shares[[server]] <- result$rlk_round1_share
      }
    }

    # Two-round RLK generation protocol (HE-Link mode only)
    if (use_he_link && length(rlk_r1_shares) > 0) {
      if (verbose) message("  Generating collective relinearization key (2-round protocol)...")

      # Round 1 aggregation: send all R1 shares to coordinator (party 0)
      combine_conn <- which(server_names == server_list[1])
      for (i in seq_along(server_list)) {
        .sendBlob(rlk_r1_shares[[server_list[i]]],
                  paste0("rlk_r1_", i - 1), combine_conn)
      }

      agg_r1_result <- DSI::datashield.aggregate(
        conns = datasources[combine_conn],
        expr = call("mheRLKAggregateR1DS",
                    from_storage = TRUE,
                    n_parties = as.integer(length(server_list)),
                    session_id = session_id)
      )
      if (is.list(agg_r1_result)) agg_r1_result <- agg_r1_result[[1]]
      agg_r1 <- agg_r1_result$aggregated_round1

      # Distribute aggregated R1 to non-coordinator servers for round 2
      for (i in 2:length(server_list)) {
        srv_conn <- which(server_names == server_list[i])
        .sendBlob(agg_r1, "rlk_agg_r1", srv_conn)
      }

      # Round 2: each server generates its R2 share
      rlk_r2_shares <- list()
      for (i in seq_along(server_list)) {
        srv_conn <- which(server_names == server_list[i])
        r2_result <- DSI::datashield.aggregate(
          conns = datasources[srv_conn],
          expr = call("mheRLKRound2DS",
                      from_storage = (i > 1),
                      session_id = session_id)
        )
        if (is.list(r2_result)) r2_result <- r2_result[[1]]
        rlk_r2_shares[[server_list[i]]] <- r2_result$rlk_round2_share
      }

      # Store aggregated R1 and R2 shares on coordinator for mheCombineDS
      .sendBlob(agg_r1, "rlk_agg_r1", combine_conn)
      for (i in seq_along(server_list)) {
        .sendBlob(rlk_r2_shares[[server_list[i]]],
                  paste0("rlk_r2_", i - 1), combine_conn)
      }
    }

    # Store PK shares, CRP, GKG seed, and GKG shares on the combining server
    # via chunked blob storage. Cryptographic objects can be several MB each,
    # exceeding R's expression parser stack limit if passed as call arguments.
    conn_idx <- which(server_names == server_list[1])

    for (i in seq_along(server_list)) {
      .sendBlob(pk_shares[[server_list[i]]], paste0("pk_", i - 1), conn_idx)
    }
    .sendBlob(crp, "crp", conn_idx)
    .sendBlob(gkg_seed, "gkg_seed", conn_idx)

    n_gkg_shares <- length(gkg_shares[[server_list[1]]])
    for (i in seq_along(server_list)) {
      shares <- gkg_shares[[server_list[i]]]
      for (j in seq_along(shares)) {
        .sendBlob(shares[j], paste0("gkg_", i - 1, "_", j - 1), conn_idx)
      }
    }

    combined <- DSI::datashield.aggregate(
      conns = datasources[conn_idx],
      expr = call("mheCombineDS",
                  from_storage = TRUE,
                  n_parties = as.integer(length(server_list)),
                  n_gkg_shares = as.integer(n_gkg_shares),
                  num_obs = as.integer(n_obs),
                  log_n = as.integer(log_n),
                  log_scale = as.integer(log_scale),
                  session_id = session_id)
    )
    if (is.list(combined)) combined <- combined[[1]]
    cpk <- combined$collective_public_key
    galois_keys <- combined$galois_keys
    relin_key <- combined$relin_key  # Non-NULL when RLK was generated

    # Distribute CPK, Galois keys, and RLK to other servers via chunked blob storage
    for (server in server_list[-1]) {
      srv_conn <- which(server_names == server)
      .sendBlob(cpk, "cpk", srv_conn)
      if (!is.null(galois_keys)) {
        for (gk_i in seq_along(galois_keys)) {
          .sendBlob(galois_keys[gk_i], paste0("gk_", gk_i - 1), srv_conn)
        }
      }
      if (!is.null(relin_key) && nzchar(relin_key)) {
        .sendBlob(relin_key, "rk", srv_conn)
      }
      DSI::datashield.aggregate(
        conns = datasources[srv_conn],
        expr = call("mheStoreCPKDS", from_storage = TRUE,
                    session_id = session_id)
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
        expr = call("mheStoreTransportKeysDS", transport_keys = tk_map,
                    session_id = session_id)
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
                  y_var = y_arg,
                  session_id = session_id)
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
                  data_name = enc_data, y_var = y_var,
                  session_id = session_id)
    )
    if (is.list(enc_result)) enc_result <- enc_result[[1]]
    ct_y <- enc_result$encrypted_y
    if (verbose) message("  Encrypted y on ", y_server)

    # In HE-Link mode, the label server also needs ct_y stored locally
    # (for computing gradient with encrypted mu: ct_y - ct_mu)
    if (use_he_link) {
      label_conn <- which(server_names == y_server)
      DSI::datashield.aggregate(
        conns = datasources[label_conn],
        expr = call("mheStoreEncYDS", enc_y = ct_y,
                    session_id = session_id)
      )
      if (verbose) message("  ct_y stored locally on label server (HE-Link mode)")
    }

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
                      chunk = substr(ct_y, start, end),
                      session_id = session_id)
        )
      }
      DSI::datashield.aggregate(
        conns = datasources[conn_idx],
        expr = call("mheAssembleEncColumnDS",
                    col_index = 1L, n_chunks = as.integer(n_chunks),
                    session_id = session_id)
      )
      if (verbose)
        message("  ct_y transferred to ", server, " (", n_chunks, " chunks)")
    }
  }

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
    if (server == coordinator && label_intercept && !use_he_link) p <- p + 1
    betas[[server]] <- rep(0, p)
  }

  converged <- FALSE
  final_iter <- 0

  # Encrypted eta blobs (opaque to client): server_name -> base64url string
  encrypted_etas <- list()

  if (use_secure_agg) {
    # =========================================================================
    # Phase 3 (Secure Aggregation): BCD with pairwise PRG masks (K>=3)
    # =========================================================================
    # Secure Aggregation Protocol: non-label servers add correlated masks to
    # their eta vectors before sending to coordinator. Masks cancel when summed,
    # so coordinator only sees aggregate eta. Individual per-server eta values
    # are never visible.
    #
    # Uses Jacobi (synchronous) BCD: all non-label servers receive the same
    # mu/w snapshot per iteration (computed by coordinator at iteration start).

    non_label_pk_map <- list()
    for (s in non_label_servers) {
      non_label_pk_map[[s]] <- transport_pks[[s]]
    }

    # Sort non-label server names for canonical ordering in seed derivation
    nonlabel_sorted <- sort(non_label_servers)

    # Fusion server and non-fusion servers (for threshold decryption)
    fusion_server <- server_list[1]  # Party 0
    fusion_conn_idx <- which(server_names == fusion_server)
    non_fusion_servers <- setdiff(server_list, fusion_server)

    if (verbose) message("\n[Phase 3] BCD iterations (secure aggregation, K=",
                         n_partitions, " servers)...")

    # Initialize FSM on coordinator
    DSI::datashield.aggregate(
      conns = datasources[coordinator_conn],
      expr = call("glmFSMInitDS",
                  session_id = session_id,
                  n_nonlabel = as.integer(non_label_count),
                  mode = "secure_agg")
    )

    # Initialize secure aggregation on each non-label server
    for (server in non_label_servers) {
      conn_idx <- which(server_names == server)
      DSI::datashield.aggregate(
        conns = datasources[conn_idx],
        expr = call("glmSecureAggInitDS",
                    self_name = server,
                    session_id = session_id,
                    nonlabel_names = nonlabel_sorted,
                    scale_bits = 20L)
      )
    }

    if (verbose) message("  Secure aggregation initialized on all non-label servers")

    for (iter in seq_len(max_iter)) {
      betas_old <- betas
      max_diff <- 0

      # --- Phase 3a: FSM check + Coordinator step (label server) ---
      # Send masked eta blobs to coordinator (from previous iteration)
      eta_blob_keys <- character(0)
      if (length(encrypted_etas) > 0) {
        for (s in names(encrypted_etas)) {
          blob <- encrypted_etas[[s]]
          if (!is.null(blob) && nzchar(blob)) {
            key <- paste0("eta_", s)
            .sendBlob(blob, key, coordinator_conn)
            eta_blob_keys <- c(eta_blob_keys, key)

            # FSM: register eta receipt
            DSI::datashield.aggregate(
              conns = datasources[coordinator_conn],
              expr = call("glmFSMCheckDS",
                          session_id = session_id,
                          action = "receive_eta",
                          server_name = s)
            )
          }
        }
      }

      # FSM: authorize coordinator step
      DSI::datashield.aggregate(
        conns = datasources[coordinator_conn],
        expr = call("glmFSMCheckDS",
                    session_id = session_id,
                    action = "coordinator_step",
                    iteration = as.integer(iter))
      )

      coord_result <- DSI::datashield.aggregate(
        conns = datasources[coordinator_conn],
        expr = call("glmSecureAggCoordinatorStepDS",
                    data_name = std_data,
                    y_var = y_var,
                    x_vars = x_vars[[coordinator]],
                    eta_blob_keys = if (length(eta_blob_keys) > 0) eta_blob_keys else NULL,
                    non_label_pks = non_label_pk_map,
                    family = family,
                    beta_current = betas[[coordinator]],
                    lambda = lambda,
                    intercept = label_intercept,
                    n_obs = as.integer(n_obs),
                    scale_bits = 20L,
                    session_id = session_id)
      )
      if (is.list(coord_result) && length(coord_result) == 1)
        coord_result <- coord_result[[1]]

      betas[[coordinator]] <- coord_result$beta
      mwv_blobs <- coord_result$encrypted_blobs

      diff_coord <- sum(abs(betas[[coordinator]] - betas_old[[coordinator]]))
      max_diff <- max(max_diff, diff_coord)

      # FSM: distribute_mwv
      DSI::datashield.aggregate(
        conns = datasources[coordinator_conn],
        expr = call("glmFSMCheckDS",
                    session_id = session_id,
                    action = "distribute_mwv")
      )

      # --- Phase 3b: Non-label servers (gradient + decrypt + masked block solve) ---
      encrypted_etas <- list()  # Reset for this iteration

      for (server in non_label_servers) {
        conn_idx <- which(server_names == server)
        vars <- x_vars[[server]]
        p_k <- length(vars)

        # Step 1: Send encrypted (mu, w) blob
        mwv_blob <- mwv_blobs[[server]]
        .sendBlob(mwv_blob, "mwv", conn_idx)

        # Step 2: Compute encrypted gradient (REUSE existing function)
        grad_result <- DSI::datashield.aggregate(
          conns = datasources[conn_idx],
          expr = call("glmSecureGradientDS",
                      data_name = std_data,
                      x_vars = vars,
                      encrypted_mwv = NULL,
                      num_obs = as.integer(n_obs),
                      session_id = session_id)
        )
        if (is.list(grad_result) && length(grad_result) == 1)
          grad_result <- grad_result[[1]]

        enc_gradients <- grad_result$encrypted_gradients
        ct_hashes <- grad_result$ct_hashes

        # Protocol Firewall: authorize gradient CTs (REUSE existing pattern)
        if (!is.null(ct_hashes) && length(ct_hashes) > 0) {
          ct_hashes_blob <- paste(ct_hashes, collapse = ",")
          for (auth_server in server_list) {
            if (auth_server != server) {
              auth_conn <- which(server_names == auth_server)
              .sendBlob(ct_hashes_blob, "ct_hashes", auth_conn)
              DSI::datashield.aggregate(
                conns = datasources[auth_conn],
                expr = call("mheAuthorizeCTDS",
                            op_type = "glm-gradient",
                            from_storage = TRUE,
                            session_id = session_id)
              )
            }
          }
        }

        # Step 3: Threshold decrypt each gradient component (REUSE)
        gradient <- numeric(p_k)
        for (j in seq_len(p_k)) {
          ct_grad <- enc_gradients[j]

          for (nf_server in non_fusion_servers) {
            nf_conn <- which(server_names == nf_server)
            nf_party_id <- which(server_list == nf_server) - 1

            n_ct_chunks <- .sendCTChunks(ct_grad, nf_conn)
            pd <- DSI::datashield.aggregate(
              conns = datasources[nf_conn],
              expr = call("mhePartialDecryptWrappedDS",
                          n_chunks = as.integer(n_ct_chunks),
                          session_id = session_id)
            )
            if (is.list(pd)) pd <- pd[[1]]
            .sendWrappedShare(pd$wrapped_share, nf_party_id, fusion_conn_idx)
          }

          n_ct_chunks <- .sendCTChunks(ct_grad, fusion_conn_idx)
          fuse_result <- DSI::datashield.aggregate(
            conns = datasources[fusion_conn_idx],
            expr = call("mheFuseServerDS",
                        n_parties = as.integer(length(server_list)),
                        n_ct_chunks = as.integer(n_ct_chunks),
                        num_slots = 0L,
                        session_id = session_id)
          )
          if (is.list(fuse_result)) fuse_result <- fuse_result[[1]]
          gradient[j] <- fuse_result$value
        }

        # Step 4: Send (mu, w) blob again for block solve (w needed for IRLS)
        .sendBlob(mwv_blob, "mwv", conn_idx)

        # Step 5: Masked block solve (secure aggregation version)
        solve_result <- DSI::datashield.aggregate(
          conns = datasources[conn_idx],
          expr = call("glmSecureAggBlockSolveDS",
                      data_name = std_data,
                      x_vars = vars,
                      beta_current = betas[[server]],
                      gradient = gradient,
                      lambda = lambda,
                      coordinator_pk = coordinator_pk,
                      iteration = as.integer(iter),
                      session_id = session_id)
        )
        if (is.list(solve_result) && length(solve_result) == 1)
          solve_result <- solve_result[[1]]

        betas[[server]] <- solve_result$beta
        encrypted_etas[[server]] <- solve_result$encrypted_masked_eta

        diff_server <- sum(abs(betas[[server]] - betas_old[[server]]))
        max_diff <- max(max_diff, diff_server)

        # FSM: block complete
        DSI::datashield.aggregate(
          conns = datasources[coordinator_conn],
          expr = call("glmFSMCheckDS",
                      session_id = session_id,
                      action = "block_complete")
        )
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

  } else if (use_he_link) {
    # =========================================================================
    # Phase 3 (HE-Link): BCD with homomorphic link function
    # =========================================================================
    # Instead of transporting plaintext (mu, w, v) between servers, we compute
    # mu = sigmoid(eta_total) entirely in the encrypted domain using polynomial
    # approximation. This prevents the K=2 eta privacy leak.
    #
    # Each iteration:
    #   1. Each server encrypts eta_k = X_k * beta_k under CPK
    #   2. Coordinator: ct_eta_total = ct_add(ct_eta_label, ct_eta_nonlabel)
    #   3. Coordinator: ct_mu = eval_poly(sigmoid_coeffs, ct_eta_total, rlk)
    #   4. Distribute ct_mu to non-label server
    #   5. Both servers: encrypted_gradient = X_k^T (ct_y - ct_mu)
    #   6. Threshold decrypt gradient scalars
    #   7. Both servers: GD block update with gradient

    if (verbose) message("\n[Phase 3] BCD iterations (HE-Link, homomorphic sigmoid)...")

    # Non-label server PKs for coordinator
    non_label_pk_map <- list()
    for (s in non_label_servers) {
      non_label_pk_map[[s]] <- transport_pks[[s]]
    }

    # Fusion server and non-fusion servers (for threshold decryption)
    fusion_server <- server_list[1]  # Party 0
    fusion_conn_idx <- which(server_names == fusion_server)
    non_fusion_servers <- setdiff(server_list, fusion_server)

    for (iter in seq_len(max_iter)) {
      betas_old <- betas
      max_diff <- 0

      # --- Step 1: Each server encrypts eta_k = X_k * beta_k ---
      for (server in server_list) {
        conn_idx <- which(server_names == server)
        party_idx <- which(server_list == server) - 1

        enc_eta_result <- DSI::datashield.aggregate(
          conns = datasources[conn_idx],
          expr = call("glmHEEncryptEtaDS",
                      data_name = std_data,
                      x_vars = x_vars[[server]],
                      beta = betas[[server]],
                      session_id = session_id)
        )
        if (is.list(enc_eta_result)) enc_eta_result <- enc_eta_result[[1]]

        # Send encrypted eta to coordinator via blob storage
        .sendBlob(enc_eta_result$encrypted_eta,
                  paste0("ct_eta_", party_idx), coordinator_conn)
      }

      # --- Step 2-3: Coordinator computes ct_mu = sigmoid(ct_eta_total) ---
      link_result <- DSI::datashield.aggregate(
        conns = datasources[coordinator_conn],
        expr = call("glmHELinkStepDS",
                    from_storage = TRUE,
                    n_parties = as.integer(length(server_list)),
                    session_id = session_id)
      )
      if (is.list(link_result)) link_result <- link_result[[1]]
      ct_mu <- link_result$ct_mu

      # --- Step 4: Distribute ct_mu to non-label servers ---
      for (server in non_label_servers) {
        conn_idx <- which(server_names == server)
        .sendBlob(ct_mu, "ct_mu", conn_idx)
      }

      # --- Step 5-7: Each server computes gradient and updates beta ---
      for (server in server_list) {
        conn_idx <- which(server_names == server)
        vars <- x_vars[[server]]
        p_k <- length(vars)

        # Compute encrypted gradient: X_k^T (ct_y - ct_mu)
        grad_result <- DSI::datashield.aggregate(
          conns = datasources[conn_idx],
          expr = call("glmHEGradientEncDS",
                      data_name = std_data,
                      x_vars = vars,
                      num_obs = as.integer(n_obs),
                      from_storage = (server != coordinator),
                      session_id = session_id)
        )
        if (is.list(grad_result)) grad_result <- grad_result[[1]]

        enc_gradients <- grad_result$encrypted_gradients
        ct_hashes <- grad_result$ct_hashes

        # Protocol Firewall: authorize gradient CTs on all other servers
        if (!is.null(ct_hashes) && length(ct_hashes) > 0) {
          ct_hashes_blob <- paste(ct_hashes, collapse = ",")
          for (auth_server in server_list) {
            if (auth_server != server) {
              auth_conn <- which(server_names == auth_server)
              .sendBlob(ct_hashes_blob, "ct_hashes", auth_conn)
              DSI::datashield.aggregate(
                conns = datasources[auth_conn],
                expr = call("mheAuthorizeCTDS",
                            op_type = "he-link-gradient",
                            from_storage = TRUE,
                            session_id = session_id)
              )
            }
          }
        }

        # Threshold decrypt each gradient component (share-wrapped)
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
                          n_chunks = as.integer(n_ct_chunks),
                          session_id = session_id)
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
                        num_slots = 0L,
                        session_id = session_id)
          )
          if (is.list(fuse_result)) fuse_result <- fuse_result[[1]]
          gradient[j] <- fuse_result$value
        }

        # GD block update (gradient descent with adaptive step size)
        # Use alpha = 1/(1 + iter/10) schedule: starts aggressive, decays
        he_alpha <- 1.0 / (1.0 + iter / 10.0)
        update_result <- DSI::datashield.aggregate(
          conns = datasources[conn_idx],
          expr = call("glmHEBlockUpdateDS",
                      beta_current = betas[[server]],
                      gradient = gradient,
                      alpha = he_alpha,
                      lambda = lambda,
                      n_obs = as.integer(n_obs),
                      session_id = session_id)
        )
        if (is.list(update_result)) update_result <- update_result[[1]]

        betas[[server]] <- update_result$beta

        diff_server <- max(abs(betas[[server]] - betas_old[[server]]))
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

      if (verbose && iter %% 10 == 0)
        message(sprintf("  Iteration %d: max diff = %.2e", iter, max_diff))
    }

    if (!converged && verbose)
      warning(sprintf("Did not converge after %d iterations (diff = %.2e)",
                      max_iter, max_diff))
  } else {
    # =========================================================================
    # Phase 3 (Secure Routing): standard transport-encrypted BCD-IRLS
    # =========================================================================
    # Secure Routing Protocol: individual-level vectors (eta, mu, w, v) are
    # transport-encrypted end-to-end between the coordinator (label server)
    # and non-label servers. The client only handles:
    #   - Betas: p_k-length coefficient vectors (safe aggregate statistics)
    #   - Opaque encrypted blobs: transport-encrypted vectors the client cannot read
    #
    # BCD variant: Jacobi-style. All non-label servers receive the same mu/w/v
    # snapshot per iteration (computed by coordinator at iteration start).

    # Non-label server PKs for the coordinator to encrypt (mu, w, v) under
    non_label_pk_map <- list()
    for (s in non_label_servers) {
      non_label_pk_map[[s]] <- transport_pks[[s]]
    }

    # Fusion server and non-fusion servers (for threshold decryption)
    fusion_server <- server_list[1]  # Party 0
    fusion_conn_idx <- which(server_names == fusion_server)
    non_fusion_servers <- setdiff(server_list, fusion_server)

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
                    n_obs = as.integer(n_obs),
                    session_id = session_id)
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
                      num_obs = as.integer(n_obs),
                      session_id = session_id)
        )
        if (is.list(grad_result) && length(grad_result) == 1)
          grad_result <- grad_result[[1]]

        enc_gradients <- grad_result$encrypted_gradients
        ct_hashes <- grad_result$ct_hashes

        # Protocol Firewall: authorize gradient ciphertexts on non-producing servers.
        # ct_hashes sent via blob storage for robustness with many variables.
        if (!is.null(ct_hashes) && length(ct_hashes) > 0) {
          ct_hashes_blob <- paste(ct_hashes, collapse = ",")
          for (auth_server in server_list) {
            if (auth_server != server) {
              auth_conn <- which(server_names == auth_server)
              .sendBlob(ct_hashes_blob, "ct_hashes", auth_conn)
              DSI::datashield.aggregate(
                conns = datasources[auth_conn],
                expr = call("mheAuthorizeCTDS",
                            op_type = "glm-gradient",
                            from_storage = TRUE,
                            session_id = session_id)
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
                          n_chunks = as.integer(n_ct_chunks),
                          session_id = session_id)
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
                        num_slots = 0L,
                        session_id = session_id)
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
                      coordinator_pk = coordinator_pk,
                      session_id = session_id)
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
  }

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
    if (server == coordinator && label_intercept && !use_he_link) {
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
    # Secure Aggregation deviance: one-time unmasked eta transmission.
    # Same pattern as HE-Link deviance: coordinator stores eta_label locally,
    # non-label servers send UNMASKED eta (one-time, transport-encrypted).
    if (verbose) message("\n[Phase 5] Computing deviance (one-time secure routing)...")

    # FSM: transition to deviance
    DSI::datashield.aggregate(
      conns = datasources[coordinator_conn],
      expr = call("glmFSMCheckDS",
                  session_id = session_id,
                  action = "deviance")
    )

    # Coordinator: store eta_label locally
    DSI::datashield.aggregate(
      conns = datasources[coordinator_conn],
      expr = call("glmSecureAggPrepDevianceDS",
                  data_name = std_data,
                  x_vars = x_vars[[coordinator]],
                  beta = betas[[coordinator]],
                  coordinator_pk = NULL,
                  session_id = session_id)
    )

    # Non-label servers: compute and transport-encrypt eta (UNMASKED)
    eta_blob_keys_final <- character(0)
    for (server in non_label_servers) {
      conn_idx <- which(server_names == server)
      prep_result <- DSI::datashield.aggregate(
        conns = datasources[conn_idx],
        expr = call("glmSecureAggPrepDevianceDS",
                    data_name = std_data,
                    x_vars = x_vars[[server]],
                    beta = betas[[server]],
                    coordinator_pk = coordinator_pk,
                    session_id = session_id)
      )
      if (is.list(prep_result)) prep_result <- prep_result[[1]]

      if (!is.null(prep_result$encrypted_eta)) {
        key <- paste0("eta_", server)
        .sendBlob(prep_result$encrypted_eta, key, coordinator_conn)
        eta_blob_keys_final <- c(eta_blob_keys_final, key)
      }
    }
  } else if (use_he_link) {
    # HE-Link deviance: one-time secure-routing deviance computation.
    # Each server computes its final eta = X*beta on the standardized scale.
    # The coordinator stores its own eta locally; non-label servers
    # transport-encrypt their eta to the coordinator. This briefly reveals
    # eta_nonlabel to the coordinator for the FINAL iteration only.
    if (verbose) message("\n[Phase 5] Computing deviance (one-time secure routing)...")

    # Coordinator: store eta_label locally
    DSI::datashield.aggregate(
      conns = datasources[coordinator_conn],
      expr = call("glmHEPrepDevianceDS",
                  data_name = std_data,
                  x_vars = x_vars[[coordinator]],
                  beta = betas[[coordinator]],
                  coordinator_pk = NULL,
                  session_id = session_id)
    )

    # Non-label servers: compute and transport-encrypt eta
    eta_blob_keys_final <- character(0)
    for (server in non_label_servers) {
      conn_idx <- which(server_names == server)
      prep_result <- DSI::datashield.aggregate(
        conns = datasources[conn_idx],
        expr = call("glmHEPrepDevianceDS",
                    data_name = std_data,
                    x_vars = x_vars[[server]],
                    beta = betas[[server]],
                    coordinator_pk = coordinator_pk,
                    session_id = session_id)
      )
      if (is.list(prep_result)) prep_result <- prep_result[[1]]

      if (!is.null(prep_result$encrypted_eta)) {
        key <- paste0("eta_", server)
        .sendBlob(prep_result$encrypted_eta, key, coordinator_conn)
        eta_blob_keys_final <- c(eta_blob_keys_final, key)
      }
    }
  } else {
    # Secure Routing deviance: use stored encrypted etas from last iteration
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
                y_mean = if (!is.null(y_mean)) y_mean else NULL,
                session_id = session_id)
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
    eta_privacy = eta_privacy,
    call = call_matched
  )

  # Clean up cryptographic state on all servers
  for (server in server_list) {
    conn_idx <- which(server_names == server)
    tryCatch(
      DSI::datashield.aggregate(conns = datasources[conn_idx],
                                expr = call("mheCleanupDS",
                                            session_id = session_id)),
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
