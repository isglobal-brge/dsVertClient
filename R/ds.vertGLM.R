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
#'   aggregation: \code{"auto"} (default) selects HE-Link for K=2 (all
#'   families), \code{"transport"} uses standard secure routing,
#'   \code{"he_link"} forces encrypted link-function mode (requires K=2,
#'   log_n >= 14).
#' @param topology Character string. Seed derivation topology for secure
#'   aggregation: \code{"pairwise"} (default, O(K-1) seeds per server) or
#'   \code{"ring"} (O(2) seeds per server for K>=4). For K=3, ring and
#'   pairwise are identical. Only used when \code{eta_privacy = "secure_agg"}.
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

  # Resolve eta_privacy mode
  if (!eta_privacy %in% c("auto", "transport", "he_link", "secure_agg", "k2_mpc"))
    stop("eta_privacy must be 'auto', 'transport', 'he_link', 'secure_agg', or 'k2_mpc'",
         call. = FALSE)

  # Validate topology parameter
  if (!topology %in% c("pairwise", "ring"))
    stop("topology must be 'pairwise' or 'ring'", call. = FALSE)

  n_partitions_check <- length(x_vars)
  non_label_count <- n_partitions_check - 1

  if (eta_privacy == "auto") {
    if (non_label_count >= 2) {
      eta_privacy <- "secure_agg"    # K>=3: pairwise PRG masks
    } else if (non_label_count == 1) {
      if (family == "gaussian") {
        eta_privacy <- "he_link"     # K=2 Gaussian: identity link, exact via HE
      } else {
        eta_privacy <- "k2_mpc"      # K=2 Binomial/Poisson: MPC piecewise exact
      }
    } else {
      eta_privacy <- "transport"     # K=1 (single server, no privacy concern)
    }
  }

  use_he_link <- (eta_privacy == "he_link")
  use_secure_agg <- (eta_privacy == "secure_agg")
  use_k2_mpc <- (eta_privacy == "k2_mpc")

  # Policy enforcement
  if (use_secure_agg && non_label_count < 2)
    stop("secure_agg requires >= 3 servers (>= 2 non-label). ",
         "For K=2, use 'he_link' (binomial) or 'transport'.", call. = FALSE)
  if (use_he_link) {
    if (non_label_count != 1)
      stop("HE-Link mode requires exactly 2 servers (1 label + 1 non-label). ",
           "Got ", n_partitions_check, " partitions.", call. = FALSE)
    if (family != "gaussian" && log_n < 14)
      stop("HE-Link mode with ", family, " requires log_n >= 14 ",
           "(need 4+ multiplicative levels for polynomial + gradient). ",
           "Got log_n = ", log_n, ".", call. = FALSE)
  }

  # Gaussian identity link doesn't need polynomial evaluation → no RLK
  generate_rlk <- use_he_link && family != "gaussian"

  # Binomial/Poisson HE-Link uses gradient descent → needs more iterations
  # Only auto-increase from the default (100); respect explicit user values
  if (use_he_link && family %in% c("binomial", "poisson") && max_iter == 100) {
    max_iter <- 300L
    if (verbose) message("  Auto-increasing max_iter to 300 for HE-Link GD mode")
  }

  # ===========================================================================
  # Helpers
  # ===========================================================================

  # Async-safe datashield.aggregate: submits the call asynchronously and polls
  # for the result, preventing reverse-proxy (nginx) read-timeout on heavy
  # computations (e.g., CKKS operations at log_n >= 14). Each poll is a short
  # HTTP GET request (< 5 s), so no connection is held open for minutes.
  # Falls back to standard synchronous call for non-Opal backends or when
  # operations are expected to be fast (log_n <= 13).
  # Async polling wrapper: prevents reverse-proxy (nginx) read-timeout on
  # heavy CKKS operations by submitting the command asynchronously and polling
  # for completion with short HTTP GETs. Only activates for known heavy
  # operations when log_n >= 14 (where single CKKS calls can exceed typical
  # 60s proxy timeouts). Chunk transfers and lightweight calls always use
  # standard synchronous dispatch.
  .need_async <- (use_he_link && log_n >= 14)

  .heavy_fns <- c("mheInitDS", "mheCombineDS", "mheEncryptRawDS",
    "glmHEEncryptEtaDS", "glmHELinkStepDS", "glmHEGradientEncDS",
    "mhePartialDecryptBatchWrappedDS", "mheFuseBatchDS",
    "mheRLKAggregateR1DS", "mheRLKRound2DS",
    "mhePartialDecryptWrappedDS", "mheFuseServerDS")

  .dsAgg <- function(conns, expr, ...) {
    if (.need_async && length(conns) == 1 && is.call(expr) &&
        as.character(expr[[1]]) %in% .heavy_fns) {
      conn <- conns[[1]]
      if (inherits(conn, "OpalConnection")) {
        opal <- conn@opal
        script <- asNamespace("DSOpal")$.deparse(expr)
        cmd_id <- opalr::opal.post(opal, "datashield", "session",
          opal$rid, "aggregate", query = list(async = "true"),
          body = script, contentType = "application/x-rscript",
          acceptType = "application/octet-stream")
        repeat {
          Sys.sleep(2)
          cmd <- opalr::opal.command(opal, cmd_id, wait = FALSE)
          if (cmd$status %in% c("COMPLETED", "FAILED", "CANCELED"))
            break
        }
        if (cmd$status != "COMPLETED") {
          msg <- cmd$error %||% cmd$status
          opalr::opal.command_rm(opal, cmd_id)
          stop("Server-side command failed on '", names(conns)[1],
               "': ", msg, call. = FALSE)
        }
        result <- opalr::opal.command_result(opal, cmd_id)
        opalr::opal.command_rm(opal, cmd_id)
        return(setNames(list(result), names(conns)[1]))
      }
    }
    # Retry once on transient 500 errors (Opal OrientDB/session bugs)
    tryCatch(
      DSI::datashield.aggregate(conns = conns, expr = expr, ...),
      error = function(e) {
        msg <- conditionMessage(e)
        if (grepl("500|NullPointer|Internal Server Error", msg)) {
          Sys.sleep(2)
          DSI::datashield.aggregate(conns = conns, expr = expr, ...)
        } else {
          stop(e)
        }
      }
    )
  }

  # Send CT chunks to a server (reusable helper, adaptive chunk size)
  .sendCTChunks <- function(ct_str, conn_idx) {
    .dsvert_adaptive_send(ct_str, function(chunk_str, chunk_idx, n_chunks) {
      .dsAgg(
        conns = datasources[conn_idx],
        expr = call("mheStoreCTChunkDS",
                    chunk_index = chunk_idx,
                    chunk = chunk_str,
                    session_id = session_id)
      )
    })
  }

  # Relay a wrapped share to the fusion server in chunks (adaptive chunk size)
  .sendWrappedShare <- function(wrapped_share, party_id, fusion_conn_idx) {
    chunk_size <- .dsvert_get_chunk_size()
    n_chars <- nchar(wrapped_share)
    n_chunks <- ceiling(n_chars / chunk_size)
    for (ch in seq_len(n_chunks)) {
      start <- (ch - 1) * chunk_size + 1
      end <- min(ch * chunk_size, n_chars)
      .dsAgg(
        conns = datasources[fusion_conn_idx],
        expr = call("mheStoreWrappedShareDS",
                    party_id = as.integer(party_id),
                    share_data = substr(wrapped_share, start, end),
                    session_id = session_id)
      )
    }
  }

  # Send a blob to a server with adaptive chunking and fallback
  .sendBlob <- function(blob, key, conn_idx) {
    .dsvert_adaptive_send(blob, function(chunk_str, chunk_idx, n_chunks) {
      if (n_chunks == 1L) {
        .dsAgg(
          conns = datasources[conn_idx],
          expr = call("mheStoreBlobDS", key = key, chunk = chunk_str,
                      session_id = session_id)
        )
      } else {
        .dsAgg(
          conns = datasources[conn_idx],
          expr = call("mheStoreBlobDS",
                      key = key,
                      chunk = chunk_str,
                      chunk_index = chunk_idx,
                      n_chunks = n_chunks,
                      session_id = session_id)
        )
      }
    })
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

  # Get observation count
  first_conn <- which(server_names == server_list[1])
  count_result <- .dsAgg(
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

  # Guaranteed cleanup on exit (even if error occurs mid-protocol)
  on.exit({
    for (.srv in server_list) {
      .ci <- which(server_names == .srv)
      tryCatch(
        .dsAgg(conns = datasources[.ci],
          expr = call("mheCleanupDS", session_id = session_id)),
        error = function(e) NULL)
    }
  }, add = TRUE)

  # ===========================================================================
  # Phase 0: MHE Key Setup + Transport Keys (only needed if non-label servers exist)
  # ===========================================================================
  transport_pks <- list()  # Collect transport PKs for secure routing + share-wrapping

  if (length(non_label_servers) > 0) {
    if (verbose) message("\n[Phase 0] MHE key setup...")

    mhe_reused <- FALSE

    # --- MHE Context Reuse Check ---
    # NOTE: Use '_' as separator — Opal's R expression parser rejects many
    # special characters ('|', '/', etc.) inside string arguments passed via
    # datashield.aggregate(). Only alphanumeric, '_', ',', ':', '-' are safe.
    context_id <- paste0(
      "peers_", paste(sort(server_list), collapse = ","),
      "_logn_", log_n,
      "_logscale_", log_scale,
      "_numobs_", n_obs,
      "_rlk_", tolower(as.character(generate_rlk))
    )

    if (isTRUE(reuse_mhe)) {
      reuse_results <- list()
      all_reusable <- TRUE
      for (server in server_list) {
        conn_idx <- which(server_names == server)
        res <- .dsAgg(
          conns = datasources[conn_idx],
          expr = call("mheReuseContextDS",
                      context_id = context_id,
                      session_id = session_id)
        )
        if (is.list(res) && !is.null(res[[1]])) res <- res[[1]]
        reuse_results[[server]] <- res
        if (!isTRUE(res$reusable)) {
          all_reusable <- FALSE
          break
        }
      }

      if (all_reusable) {
        mhe_reused <- TRUE
        for (server in server_list) {
          transport_pks[[server]] <- reuse_results[[server]]$transport_pk
        }
        if (verbose) message("  MHE keys reused from cached context (skipped)")
      }
    }

    if (!mhe_reused) {
      conn_idx <- which(server_names == server_list[1])
      result0 <- .dsAgg(
        conns = datasources[conn_idx],
        expr = call("mheInitDS",
                    party_id = 0L, crp = NULL, gkg_seed = NULL,
                    num_obs = as.integer(n_obs),
                    log_n = as.integer(log_n),
                    log_scale = as.integer(log_scale),
                    generate_rlk = generate_rlk,
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
      if (generate_rlk && !is.null(result0$rlk_round1_share)) {
        rlk_r1_shares[[server_list[1]]] <- result0$rlk_round1_share
      }

      # Send CRP, GKG seed, and party_id to all non-party0 servers (sequential
      # blob transfer, then parallel mheInitDS calls).
      other_servers <- server_list[-1]
      other_conn_idxs <- integer(length(other_servers))
      for (i in seq_along(other_servers)) {
        server <- other_servers[i]
        ci <- which(server_names == server)
        other_conn_idxs[i] <- ci
        .sendBlob(crp, "crp", ci)
        .sendBlob(gkg_seed, "gkg_seed", ci)
        .sendBlob(as.character(i), "party_id", ci)
      }

      # Launch mheInitDS on all non-party0 servers in parallel
      init_expr <- call("mheInitDS",
                        party_id = 0L,
                        from_storage = TRUE,
                        num_obs = as.integer(n_obs),
                        log_n = as.integer(log_n),
                        log_scale = as.integer(log_scale),
                        generate_rlk = generate_rlk,
                        session_id = session_id)

      init_results <- .dsAgg(
        conns = datasources[other_servers],
        expr = init_expr
      )

      for (i in seq_along(other_servers)) {
        server <- other_servers[i]
        result <- init_results[[server]]
        pk_shares[[server]] <- result$public_key_share
        gkg_shares[[server]] <- result$galois_key_shares
        transport_pks[[server]] <- result$transport_pk
        if (generate_rlk && !is.null(result$rlk_round1_share)) {
          rlk_r1_shares[[server]] <- result$rlk_round1_share
        }
      }

      # Two-round RLK generation protocol (binomial/Poisson HE-Link only)
      if (generate_rlk && length(rlk_r1_shares) > 0) {
        if (verbose) message("  Generating collective relinearization key (2-round protocol)...")

        # Round 1 aggregation: send all R1 shares to coordinator (party 0)
        combine_conn <- which(server_names == server_list[1])
        for (i in seq_along(server_list)) {
          .sendBlob(rlk_r1_shares[[server_list[i]]],
                    paste0("rlk_r1_", i - 1), combine_conn)
        }

        agg_r1_result <- .dsAgg(
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
          r2_result <- .dsAgg(
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

      combined <- .dsAgg(
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

      # Distribute CPK, Galois keys, and RLK to other servers (blob send then parallel store)
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
      }
      # Store CPK on all non-party0 servers in parallel
      .dsAgg(
        conns = datasources[server_list[-1]],
        expr = call("mheStoreCPKDS", from_storage = TRUE,
                    session_id = session_id)
      )

      # Register context for future reuse (parallel across all servers)
      .dsAgg(
        conns = datasources[server_list],
        expr = call("mheRegisterContextDS",
                    context_id = context_id,
                    session_id = session_id)
      )
    }  # end if (!mhe_reused)

    # Distribute transport keys for share-wrapping and secure routing
    # (always needed - fresh transport keys even when MHE keys are reused)
    fusion_server <- server_list[1]  # Party 0
    fusion_conn_idx <- which(server_names == fusion_server)
    non_fusion_servers <- setdiff(server_list, fusion_server)

    tk_map <- list(fusion = transport_pks[[fusion_server]])
    for (s in server_list) tk_map[[s]] <- transport_pks[[s]]
    # Distribute transport keys to all servers in parallel (same key map)
    .dsAgg(
      conns = datasources[server_list],
      expr = call("mheStoreTransportKeysDS", transport_keys = tk_map,
                  session_id = session_id)
    )

    if (verbose) message("  Key setup + transport keys complete")

    # =========================================================================
    # Phase 0.5: Peer Manifest Consensus (when use_secure_agg)
    # =========================================================================
    # Cross-validate that all servers see the same peer set and configuration.
    # Prevents phantom peer injection by a malicious client.
    if (use_secure_agg) {
      if (verbose) message("\n[Phase 0.5] Peer manifest consensus...")

      # Build canonical manifest (client-side helper)
      .buildManifest <- function() {
        peers <- sort(server_list)
        pk_sorted <- transport_pks[sort(names(transport_pks))]
        manifest <- list(
          job_id = session_id,
          method = family,
          non_label_servers = sort(non_label_servers),
          peers = peers,
          topology = topology,
          transport_pks = pk_sorted,
          y_server = y_server
        )
        jsonlite::toJSON(manifest, auto_unbox = TRUE)
      }
      manifest_json <- .buildManifest()

      # B64url-encode the manifest JSON to avoid Opal parser issues with
      # special characters ('{', '}', '+', '/', '=') in string arguments.
      manifest_b64 <- jsonlite::base64_enc(charToRaw(as.character(manifest_json)))
      manifest_b64 <- gsub("[\r\n]", "", manifest_b64)  # strip MIME line breaks
      manifest_b64 <- gsub("+", "-", manifest_b64, fixed = TRUE)
      manifest_b64 <- gsub("/", "_", manifest_b64, fixed = TRUE)
      manifest_b64 <- gsub("=", "", manifest_b64, fixed = TRUE)

      # Each server: store manifest + get encrypted hashes
      server_enc_hashes <- list()
      for (server in server_list) {
        conn_idx <- which(server_names == server)
        store_result <- .dsAgg(
          conns = datasources[conn_idx],
          expr = call("peerManifestStoreDS",
                      manifest_json = manifest_b64,
                      session_id = session_id)
        )
        if (is.list(store_result) && length(store_result) == 1)
          store_result <- store_result[[1]]
        server_enc_hashes[[server]] <- store_result
      }

      # Relay each server's encrypted hash to each peer via blob storage
      for (sender in server_list) {
        hashes <- server_enc_hashes[[sender]]
        for (recipient in server_list) {
          if (recipient == sender) next
          hash_blob <- hashes[[recipient]]
          if (is.null(hash_blob)) next
          recipient_conn <- which(server_names == recipient)
          blob_key <- paste0("manifest_hash_", sender)
          .sendBlob(hash_blob, blob_key, recipient_conn)
        }
      }

      # Each server: validate every peer's hash
      for (server in server_list) {
        conn_idx <- which(server_names == server)
        for (peer in server_list) {
          if (peer == server) next
          .dsAgg(
            conns = datasources[conn_idx],
            expr = call("peerManifestValidateDS",
                        peer_name = peer,
                        session_id = session_id)
          )
        }
      }

      if (verbose) message("  Manifest consensus validated across all servers")
    }
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

    std_result <- .dsAgg(
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
  if (length(non_label_servers) > 0 && !use_k2_mpc) {
    if (verbose) message("\n[Phase 2] Encrypting response variable...")

    # Encrypt STANDARDIZED y (Gaussian) or raw y (non-Gaussian)
    enc_data <- if (standardize_y) std_data else data_name

    conn_idx <- which(server_names == y_server)
    enc_result <- .dsAgg(
      conns = datasources[conn_idx],
      expr = call("mheEncryptRawDS",
                  data_name = enc_data, y_var = y_var,
                  store_local = use_he_link,
                  session_id = session_id)
    )
    if (is.list(enc_result)) enc_result <- enc_result[[1]]
    ct_y <- enc_result$encrypted_y
    if (verbose) message("  Encrypted y on ", y_server)
    if (use_he_link && verbose)
      message("  ct_y stored locally on label server (HE-Link mode)")

    for (server in non_label_servers) {
      conn_idx <- which(server_names == server)
      n_chunks_sent <- .dsvert_adaptive_send(ct_y, function(chunk_str, chunk_idx, n_chunks) {
        .dsAgg(
          conns = datasources[conn_idx],
          expr = call("mheStoreEncChunkDS",
                      col_index = 1L,
                      chunk_index = chunk_idx,
                      chunk = chunk_str,
                      session_id = session_id)
        )
      })
      .dsAgg(
        conns = datasources[conn_idx],
        expr = call("mheAssembleEncColumnDS",
                    col_index = 1L, n_chunks = as.integer(n_chunks_sent),
                    session_id = session_id)
      )
      if (verbose)
        message("  ct_y transferred to ", server, " (", n_chunks_sent, " chunks)")
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
    if (server == coordinator && label_intercept && !use_he_link && !use_k2_mpc) p <- p + 1
    betas[[server]] <- rep(0, p)
  }

  converged <- FALSE
  final_iter <- 0

  # Encrypted eta blobs (opaque to client): server_name -> base64url string
  encrypted_etas <- list()

  if (use_k2_mpc) {
    # =========================================================================
    # Phase 3 (K=2 MPC): Secure Beaver-Triple BCD with Gradient Descent
    # =========================================================================
    # FULLY SECURE: neither party sees eta_total, mu, residual, or peer data.
    # Only p_k gradient scalars revealed per server per iteration (same as HE-Link).
    #
    # Protocol per iteration:
    #   A. Split eta → additive shares of eta_total
    #   B. Beaver polynomial eval → shares of mu (degree-21 for binomial)
    #   C. Residual shares + local gradient (plaintext × own share)
    #   D. Cross-gradient via Beaver (X_k × peer_residual on shares)
    #   E. Reconstruct p_k gradient scalars → GD update

    if (verbose) message("\n[Phase 3] BCD iterations (K=2 secure Beaver-triple MPC)...")

    nl <- non_label_servers[1]
    nl_conn <- which(server_names == nl)
    frac_bits <- 20L

    # Polynomial setup
    poly_degree <- if (family == "poisson") 15L else 21L
    poly_info <- dsVert:::.callMheTool("mpc-get-poly-coeffs", list(
      family = family, degree = poly_degree))
    poly_coeffs <- poly_info$coefficients
    if (verbose) message(sprintf("  Degree-%d polynomial, max approx error: %.2e",
                                  poly_degree, poly_info$max_error))

    # p_k for each server (for cross-gradient triple allocation)
    p_coord <- length(x_vars[[coordinator]])
    p_nl <- length(x_vars[[nl]])
    gd_alpha <- 0.5  # aggressive start, decay after 20 iters

    for (iter in seq_len(max_iter)) {
      betas_old <- betas

      # === A: Split eta → shares ===
      for (server in server_list) {
        ci <- which(server_names == server)
        peer <- setdiff(server_list, server)
        r <- .dsAgg(datasources[ci], call("k2MpcSplitEtaDS",
          data_name=std_data, x_vars=x_vars[[server]], beta=betas[[server]],
          peer_pk=transport_pks[[peer]], frac_bits=frac_bits,
          session_id=session_id))
        if (is.list(r)) r <- r[[1]]
        # Relay peer's share
        peer_ci <- which(server_names == peer)
        .sendBlob(r$peer_share_enc, "mpc_peer_eta_share", peer_ci)
      }
      # Combine shares on both servers
      for (server in server_list) {
        ci <- which(server_names == server)
        .dsAgg(datasources[ci], call("k2MpcSecureStepDS", step="combine_eta",
               session_id=session_id))
      }

      # === B: Generate + distribute Beaver triples ===
      n_poly_triples <- poly_degree * n_obs
      n_cross_coord <- n_obs * p_coord   # for coordinator's cross-gradient
      n_cross_nl <- n_obs * p_nl         # for non-label's cross-gradient
      n_total <- n_poly_triples + n_cross_coord + n_cross_nl

      # Client generates triples (trusted dealer)
      tu <- runif(n_total, -1, 1)
      tv <- runif(n_total, -1, 1)
      tw <- tu * tv
      tu0 <- runif(n_total, -5, 5); tu1 <- tu - tu0
      tv0 <- runif(n_total, -5, 5); tv1 <- tv - tv0
      tw0 <- runif(n_total, -5, 5); tw1 <- tw - tw0

      # Send to both servers
      .b64url_to_b64 <- function(x) {
        x <- gsub("-", "+", gsub("_", "/", x, fixed=TRUE), fixed=TRUE)
        pad <- nchar(x) %% 4
        if (pad == 2) x <- paste0(x, "==")
        if (pad == 3) x <- paste0(x, "=")
        x
      }
      .b64_to_b64url <- function(x) {
        gsub("+", "-", gsub("/", "_", gsub("=+$", "", x, perl=TRUE), fixed=TRUE), fixed=TRUE)
      }

      for (server in server_list) {
        ci <- which(server_names == server)
        is_coord <- (server == coordinator)
        us <- if (is_coord) tu0 else tu1
        vs <- if (is_coord) tv0 else tv1
        ws <- if (is_coord) tw0 else tw1
        pk_b64 <- .b64url_to_b64(transport_pks[[server]])
        sealed <- dsVert:::.callMheTool("transport-encrypt-vectors", list(
          vectors = list(u=us, v=vs, w=ws), recipient_pk=pk_b64))
        .sendBlob(.b64_to_b64url(sealed$sealed), "beaver_triples", ci)
        .dsAgg(datasources[ci], call("k2MpcSecureStepDS", step="store_triples",
               session_id=session_id))
      }

      # === C: Beaver polynomial evaluation ===
      tri_offset <- 0
      for (k in seq_len(poly_degree)) {
        a_key <- if (k == 1) "secure_eta_share" else paste0("secure_pow", k)
        b_key <- "secure_eta_share"
        result_key <- paste0("secure_pow", k + 1)
        idx <- (tri_offset + 1):(tri_offset + n_obs)

        # Open on both servers (pass triple offset + count, not values)
        open_results <- list()
        for (server in server_list) {
          ci <- which(server_names == server)
          peer <- setdiff(server_list, server)
          r <- .dsAgg(datasources[ci], call("k2MpcSecureStepDS", step="beaver_open",
            a_share_key=a_key, b_share_key=b_key,
            u_vals=tri_offset, v_vals=n_obs,
            peer_pk=transport_pks[[peer]], frac_bits=frac_bits,
            session_id=session_id))
          if (is.list(r)) r <- r[[1]]
          open_results[[server]] <- r
        }
        .sendBlob(open_results[[coordinator]]$peer_de_enc, "beaver_peer_de", nl_conn)
        .sendBlob(open_results[[nl]]$peer_de_enc, "beaver_peer_de", coordinator_conn)

        for (server in server_list) {
          ci <- which(server_names == server)
          is_coord <- (server == coordinator)
          .dsAgg(datasources[ci], call("k2MpcSecureStepDS", step="beaver_close",
            result_key=result_key,
            u_vals=tri_offset, v_vals=n_obs,
            party_id=if(is_coord) 0L else 1L,
            frac_bits=frac_bits, session_id=session_id))
        }
        tri_offset <- tri_offset + n_obs
      }

      # Polynomial assembly
      power_keys <- c("secure_eta_share", paste0("secure_pow", 2:(poly_degree+1)))
      for (server in server_list) {
        ci <- which(server_names == server)
        is_coord <- (server == coordinator)
        .dsAgg(datasources[ci], call("k2MpcSecureStepDS", step="poly_eval",
          power_keys=power_keys, coefficients=poly_coeffs,
          party_id=if(is_coord) 0L else 1L,
          frac_bits=frac_bits, session_id=session_id))
      }

      # === D: Secure gradient computation ===
      # D1: Prepare gradient (residual share + local grad + share private inputs)
      prep_results <- list()
      for (server in server_list) {
        ci <- which(server_names == server)
        peer <- setdiff(server_list, server)
        r <- .dsAgg(datasources[ci], call("k2MpcSecureStepDS", step="prepare_gradient",
          data_name=std_data, x_vars=x_vars[[server]],
          y_var=if(server==coordinator) y_var else NULL,
          role=if(server==coordinator) "label" else "nonlabel",
          peer_pk=transport_pks[[peer]],
          frac_bits=frac_bits, session_id=session_id))
        if (is.list(r)) r <- r[[1]]
        prep_results[[server]] <- r
      }

      # D2: Relay shared inputs (both FixedPoint, via transport-encrypt from server)
      .sendBlob(prep_results[[coordinator]]$x_share_for_peer_enc, "peer_x_share", nl_conn)
      .sendBlob(prep_results[[coordinator]]$r_share_for_peer_enc, "peer_r_share", nl_conn)
      .sendBlob(prep_results[[nl]]$x_share_for_peer_enc, "peer_x_share", coordinator_conn)
      .sendBlob(prep_results[[nl]]$r_share_for_peer_enc, "peer_r_share", coordinator_conn)

      for (server in server_list) {
        ci <- which(server_names == server)
        .dsAgg(datasources[ci], call("k2MpcSecureStepDS", step="receive_peer_shares",
               session_id=session_id))
      }

      # D3-D5: Sequential cross-gradient rounds (one per target server)
      # Each round computes one server's cross-gradient term.
      # Round for target T: both servers do Beaver(X_T_shares, resid_peer_T)
      # where peer_T contributed its residual share and T shared its X columns.
      # The Beaver size = n_obs * p_T (T's feature count).
      cross_grad_shares <- list()  # target → cross-gradient share from peer

      for (target in server_list) {
        peer_of_target <- setdiff(server_list, target)
        is_target_coord <- (target == coordinator)
        # p_target includes intercept for coordinator
        p_target <- if (is_target_coord) p_coord else p_nl
        cross_off <- n_poly_triples + (if(is_target_coord) 0 else n_cross_coord)
        cross_cnt <- n_obs * p_target

        # Both servers do Beaver open (same size: n_obs * p_target)
        open_r <- list()
        for (server in server_list) {
          ci <- which(server_names == server)
          server_role <- if (server == target) "target" else "peer"
          r <- .dsAgg(datasources[ci], call("k2MpcSecureStepDS",
            step="cross_gradient_open",
            u_vals=cross_off, v_vals=cross_cnt,
            role=server_role,
            peer_pk=transport_pks[[setdiff(server_list, server)]],
            frac_bits=frac_bits, session_id=session_id))
          if (is.list(r)) r <- r[[1]]
          open_r[[server]] <- r
        }

        # Relay DEs
        .sendBlob(open_r[[coordinator]]$peer_de_enc, "cross_beaver_peer_de", nl_conn)
        .sendBlob(open_r[[nl]]$peer_de_enc, "cross_beaver_peer_de", coordinator_conn)

        # Both do Beaver close (with role: target stores own share, peer encrypts)
        close_r <- list()
        for (server in server_list) {
          ci <- which(server_names == server)
          is_coord <- (server == coordinator)
          server_role <- if (server == target) "target" else "peer"
          r <- .dsAgg(datasources[ci], call("k2MpcSecureStepDS",
            step="cross_gradient_close",
            u_vals=cross_off, v_vals=cross_cnt,
            party_id=if(is_coord) 0L else 1L,
            role=server_role,
            peer_pk=transport_pks[[setdiff(server_list, server)]],
            frac_bits=frac_bits, session_id=session_id))
          if (is.list(r)) r <- r[[1]]
          close_r[[server]] <- r
        }

        # peer's close result has the encrypted share for the target
        cross_grad_shares[[target]] <- close_r[[peer_of_target]]
      }

      # Relay cross-gradient shares to the target servers
      for (target in server_list) {
        target_ci <- which(server_names == target)
        .sendBlob(cross_grad_shares[[target]]$cross_for_peer_enc,
                  "cross_gradient_from_peer", target_ci)
      }

      # D6: Combine all gradient parts + GD update
      for (server in server_list) {
        ci <- which(server_names == server)
        r <- .dsAgg(datasources[ci], call("k2MpcSecureStepDS",
          step="combine_gradient",
          frac_bits=frac_bits, session_id=session_id))
        if (is.list(r)) r <- r[[1]]

        gradient <- r$gradient
        # Learning rate with decay (aggressive start for fast convergence)
        lr <- if (iter <= 20) gd_alpha else gd_alpha / (1 + (iter - 20) / 15)
        grad_update <- gradient / n_obs - lambda * betas[[server]]
        betas[[server]] <- betas[[server]] + lr * grad_update
      }

      # Check convergence
      max_diff <- 0
      for (server in server_list) {
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

  } else if (use_secure_agg) {
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
    .dsAgg(
      conns = datasources[coordinator_conn],
      expr = call("glmFSMInitDS",
                  session_id = session_id,
                  n_nonlabel = as.integer(non_label_count),
                  mode = "secure_agg")
    )

    # Initialize secure aggregation on each non-label server
    for (server in non_label_servers) {
      conn_idx <- which(server_names == server)
      .dsAgg(
        conns = datasources[conn_idx],
        expr = call("glmSecureAggInitDS",
                    self_name = server,
                    session_id = session_id,
                    nonlabel_names = nonlabel_sorted,
                    scale_bits = 20L,
                    topology = topology)
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
            .dsAgg(
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
      .dsAgg(
        conns = datasources[coordinator_conn],
        expr = call("glmFSMCheckDS",
                    session_id = session_id,
                    action = "coordinator_step",
                    iteration = as.integer(iter))
      )

      coord_result <- .dsAgg(
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
      .dsAgg(
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
        grad_result <- .dsAgg(
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
              .dsAgg(
                conns = datasources[auth_conn],
                expr = call("mheAuthorizeCTDS",
                            op_type = "glm-gradient",
                            from_storage = TRUE,
                            session_id = session_id)
              )
            }
          }
        }

        # Step 3: Batched threshold decryption of gradient components
        # Send all p_k gradient CTs as batch blobs, then one round-trip per
        # server instead of O(p_k * K) individual round-trips.
        gradient <- numeric(p_k)

        # Send all gradient CTs to each non-fusion server
        for (nf_server in non_fusion_servers) {
          nf_conn <- which(server_names == nf_server)
          nf_party_id <- which(server_list == nf_server) - 1
          for (j in seq_len(p_k)) {
            .sendBlob(enc_gradients[j], paste0("ct_batch_", j), nf_conn)
          }
          pd <- .dsAgg(
            conns = datasources[nf_conn],
            expr = call("mhePartialDecryptBatchWrappedDS",
                        n_cts = as.integer(p_k),
                        session_id = session_id)
          )
          if (is.list(pd)) pd <- pd[[1]]
          for (j in seq_len(p_k)) {
            share_key <- paste0("wrapped_share_", nf_party_id, "_ct_", j)
            .sendBlob(pd$wrapped_shares[j], share_key, fusion_conn_idx)
          }
        }

        # Send all gradient CTs to fusion server and batch fuse
        for (j in seq_len(p_k)) {
          .sendBlob(enc_gradients[j], paste0("ct_batch_", j), fusion_conn_idx)
        }
        fuse_result <- .dsAgg(
          conns = datasources[fusion_conn_idx],
          expr = call("mheFuseBatchDS",
                      n_cts = as.integer(p_k),
                      n_parties = as.integer(length(server_list)),
                      num_slots = 0L,
                      session_id = session_id)
        )
        if (is.list(fuse_result)) fuse_result <- fuse_result[[1]]
        gradient <- fuse_result$values

        # Step 4: Send (mu, w) blob again for block solve (w needed for IRLS)
        .sendBlob(mwv_blob, "mwv", conn_idx)

        # Step 5: Masked block solve (secure aggregation version)
        solve_result <- .dsAgg(
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
        .dsAgg(
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
    # Phase 3 (HE-Link): BCD with encrypted link function (K=2)
    # =========================================================================
    # Instead of transporting plaintext (mu, w, v) between servers, we compute
    # mu = g^{-1}(eta_total) entirely in the encrypted domain. This prevents
    # the K=2 eta privacy leak.
    #
    # Family-specific pathways:
    #   - Gaussian: identity link, ct_mu = ct_eta_total (no polynomial, no RLK)
    #              → exact block solve with unit weights
    #   - Binomial: sigmoid polynomial under CKKS → GD block update
    #   - Poisson:  exp polynomial under CKKS with eta clipping → GD block update
    #
    # Each iteration:
    #   1. Each server encrypts eta_k = X_k * beta_k under CPK (Poisson: clip)
    #   2. Coordinator: ct_eta_total = ct_add(ct_eta_label, ct_eta_nonlabel)
    #   3. Coordinator: ct_mu = poly(ct_eta_total) or ct_eta_total (Gaussian)
    #   4. Distribute ct_mu to non-label server
    #   5. Both servers: encrypted_gradient = X_k^T (ct_y - ct_mu)
    #   6. Threshold decrypt gradient scalars
    #   7. Block update: exact solve (Gaussian) or GD (binomial/Poisson)

    if (verbose) {
      link_desc <- switch(family,
        gaussian = "encrypted residual, identity link",
        binomial = "homomorphic sigmoid",
        poisson  = "homomorphic exp"
      )
      message(sprintf("\n[Phase 3] BCD iterations (HE-Link, %s)...", link_desc))
    }

    # Poisson: degree-7 Chebyshev approximation to exp(x) on [-3, 3]
    he_link_poly <- NULL
    he_link_clip_radius <- NULL
    if (family == "poisson") {
      he_link_poly <- c(9.989069e-01, 9.994338e-01, 5.043165e-01, 1.675787e-01,
                        3.908437e-02, 7.945614e-03, 1.865407e-03, 2.576334e-04)
      he_link_clip_radius <- 1.5  # Per-server clip → total in [-3, 3]
    }

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

      # --- Step 1: Each server encrypts eta_k = X_k * beta_k (parallel) ---
      enc_eta_results <- list()
      for (server in server_list) {
        conn_idx <- which(server_names == server)
        enc_eta_results[[server]] <- .dsAgg(
          conns = datasources[conn_idx],
          expr = call("glmHEEncryptEtaDS",
                      data_name = std_data,
                      x_vars = x_vars[[server]],
                      beta = betas[[server]],
                      clip_radius = he_link_clip_radius,
                      session_id = session_id)
        )
      }
      # Send encrypted etas to coordinator
      for (server in server_list) {
        party_idx <- which(server_list == server) - 1
        enc_eta_result <- enc_eta_results[[server]]
        if (is.list(enc_eta_result)) enc_eta_result <- enc_eta_result[[1]]
        .sendBlob(enc_eta_result$encrypted_eta,
                  paste0("ct_eta_", party_idx), coordinator_conn)
      }

      # --- Step 2-3: Coordinator computes ct_mu = g^{-1}(ct_eta_total) ---
      if (family == "gaussian") {
        # Gaussian: identity link → ct_mu = ct_eta_total (no polynomial)
        link_result <- .dsAgg(
          conns = datasources[coordinator_conn],
          expr = call("glmHELinkStepDS",
                      from_storage = TRUE,
                      n_parties = as.integer(length(server_list)),
                      skip_poly = TRUE,
                      session_id = session_id)
        )
      } else if (family == "poisson") {
        # Poisson: exp polynomial on ct_eta_total
        link_result <- .dsAgg(
          conns = datasources[coordinator_conn],
          expr = call("glmHELinkStepDS",
                      from_storage = TRUE,
                      n_parties = as.integer(length(server_list)),
                      poly_coefficients = he_link_poly,
                      session_id = session_id)
        )
      } else {
        # Binomial: default sigmoid polynomial
        link_result <- .dsAgg(
          conns = datasources[coordinator_conn],
          expr = call("glmHELinkStepDS",
                      from_storage = TRUE,
                      n_parties = as.integer(length(server_list)),
                      session_id = session_id)
        )
      }
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
        grad_result <- .dsAgg(
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
              .dsAgg(
                conns = datasources[auth_conn],
                expr = call("mheAuthorizeCTDS",
                            op_type = "he-link-gradient",
                            from_storage = TRUE,
                            session_id = session_id)
              )
            }
          }
        }

        # Batched threshold decryption of all p_k gradient components
        # Send all gradient CTs to each non-fusion server in one batch
        gradient <- numeric(p_k)
        for (nf_server in non_fusion_servers) {
          nf_conn <- which(server_names == nf_server)
          nf_party_id <- which(server_list == nf_server) - 1
          for (j in seq_len(p_k)) {
            .sendBlob(enc_gradients[j], paste0("ct_batch_", j), nf_conn)
          }
          pd <- .dsAgg(
            conns = datasources[nf_conn],
            expr = call("mhePartialDecryptBatchWrappedDS",
                        n_cts = as.integer(p_k),
                        session_id = session_id)
          )
          if (is.list(pd)) pd <- pd[[1]]
          for (j in seq_len(p_k)) {
            share_key <- paste0("wrapped_share_", nf_party_id, "_ct_", j)
            .sendBlob(pd$wrapped_shares[j], share_key, fusion_conn_idx)
          }
        }

        # Send all gradient CTs to fusion server and batch fuse
        for (j in seq_len(p_k)) {
          .sendBlob(enc_gradients[j], paste0("ct_batch_", j), fusion_conn_idx)
        }
        fuse_result <- .dsAgg(
          conns = datasources[fusion_conn_idx],
          expr = call("mheFuseBatchDS",
                      n_cts = as.integer(p_k),
                      n_parties = as.integer(length(server_list)),
                      num_slots = 0L,
                      session_id = session_id)
        )
        if (is.list(fuse_result)) fuse_result <- fuse_result[[1]]
        gradient <- fuse_result$values

        # Block update: family-specific strategy
        if (family == "gaussian") {
          # Exact Newton step with unit weights (identity link, w=1)
          update_result <- .dsAgg(
            conns = datasources[conn_idx],
            expr = call("glmBlockSolveDS",
                        data_name = std_data,
                        x_vars = vars,
                        w = NULL,
                        beta_current = betas[[server]],
                        gradient = gradient,
                        lambda = lambda,
                        session_id = session_id)
          )
          if (is.list(update_result)) update_result <- update_result[[1]]
          betas[[server]] <- update_result$beta
        } else {
          # Plain GD for binomial/Poisson (client-side, no server call)
          # Nesterov momentum was benchmarked but does not improve the
          # converged error, which is dominated by the polynomial sigmoid
          # approximation gap (see k2_optimizer_benchmarks.md).
          if (family == "binomial") {
            he_alpha <- if (iter <= 50) 0.3 else 0.3 / (1.0 + (iter - 50) / 20.0)
          } else {
            he_alpha <- if (iter <= 50) 0.2 else 0.2 / (1.0 + (iter - 50) / 20.0)
          }
          grad_update <- gradient / n_obs - lambda * betas[[server]]
          betas[[server]] <- betas[[server]] + he_alpha * grad_update
        }

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

      # Periodic server-side GC to prevent memory buildup on long runs
      if (iter %% 20 == 0) {
        for (server in server_list) {
          tryCatch(.dsAgg(conns = datasources[which(server_names == server)],
            expr = call("mheGcDS")), error = function(e) NULL)
        }
      }
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

      coord_result <- .dsAgg(
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
        grad_result <- .dsAgg(
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
              .dsAgg(
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
            pd <- .dsAgg(
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
          fuse_result <- .dsAgg(
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

        solve_result <- .dsAgg(
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
    if (server == coordinator && label_intercept && !use_he_link && !use_k2_mpc) {
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
    .dsAgg(
      conns = datasources[coordinator_conn],
      expr = call("glmFSMCheckDS",
                  session_id = session_id,
                  action = "deviance")
    )

    # Coordinator: store eta_label locally
    .dsAgg(
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
      prep_result <- .dsAgg(
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
  } else if (use_he_link || use_k2_mpc) {
    # HE-Link / K2-MPC deviance: one-time secure-routing deviance computation.
    # Each server computes its final eta = X*beta on the standardized scale.
    # The coordinator stores its own eta locally; non-label servers
    # transport-encrypt their eta to the coordinator. This briefly reveals
    # eta_nonlabel to the coordinator for the FINAL iteration only.
    if (verbose) message("\n[Phase 5] Computing deviance (one-time secure routing)...")

    # Coordinator: store eta_label locally
    .dsAgg(
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
      prep_result <- .dsAgg(
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

  deviance_result <- .dsAgg(
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
