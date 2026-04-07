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

  # RLK not needed (no CKKS polynomial evaluation in secure_agg or K=2 Beaver)
  generate_rlk <- FALSE

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
  .need_async <- FALSE  # sync dispatch for all current paths

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
      tryCatch(
        .dsAgg(conns = datasources[.ci],
          expr = call("mheGcDS")),
        error = function(e) NULL)
    }
    .dsvert_reset_chunk_size()
  }, add = TRUE)

  # ===========================================================================
  # Phase 0: MHE Key Setup + Transport Keys (only needed if non-label servers exist)
  # ===========================================================================
  transport_pks <- list()  # Collect transport PKs for secure routing + share-wrapping

  if (use_k2_beaver && length(non_label_servers) > 0) {
    # K=2 wide spline: generate transport keys on SERVERS via lightweight mheInitDS.
    # This creates the X25519 keypair server-side (needed for k2ReceiveShareDS decrypt).
    if (verbose) message("\n[Phase 0] Transport key setup (K=2 wide spline)...")
    crp_k2 <- NULL
    gkg_k2 <- NULL
    for (server in server_list) {
      conn_idx <- which(server_names == server)
      party_id <- as.integer(which(server_list == server) - 1L)
      tk_result <- .dsAgg(
        conns = datasources[conn_idx],
        expr = call("mheInitDS",
                    party_id = party_id,
                    crp = crp_k2, gkg_seed = gkg_k2,
                    num_obs = as.integer(n_obs),
                    log_n = 14L, log_scale = 40L,
                    generate_rlk = FALSE,
                    session_id = session_id)
      )
      if (is.list(tk_result)) tk_result <- tk_result[[1]]
      transport_pks[[server]] <- tk_result$transport_pk
      if (is.null(crp_k2)) {
        crp_k2 <- tk_result$crp
        gkg_k2 <- tk_result$gkg_seed
      }
    }
    # Store transport keys on each server (for peer PK lookup in k2ShareInputDS)
    for (server in server_list) {
      conn_idx <- which(server_names == server)
      pk_sorted <- transport_pks[sort(names(transport_pks))]
      .dsAgg(
        conns = datasources[conn_idx],
        expr = call("mheStoreTransportKeysDS",
                    transport_keys = pk_sorted,
                    session_id = session_id)
      )
    }
    if (verbose) message("  Transport keys exchanged for ", length(server_list), " servers")
  }

  if (length(non_label_servers) > 0 && !use_k2_beaver) {
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
      # ---- Coin-Tossing CRP (distributed randomness) ----
      if (verbose) message("  Coin-tossing CRP (distributed randomness)...")

      commit_results <- DSI::datashield.aggregate(
        conns = datasources[server_list],
        expr = call("mheCoinTossCommitDS", session_id = session_id)
      )
      commitments <- character(length(server_list))
      for (i in seq_along(server_list)) {
        r <- commit_results[[server_list[i]]]
        if (is.list(r) && length(r) == 1 && is.list(r[[1]])) r <- r[[1]]
        commitments[i] <- r$commitment
      }

      reveal_results <- DSI::datashield.aggregate(
        conns = datasources[server_list],
        expr = call("mheCoinTossRevealDS", session_id = session_id)
      )
      contributions <- character(length(server_list))
      for (i in seq_along(server_list)) {
        r <- reveal_results[[server_list[i]]]
        if (is.list(r) && length(r) == 1 && is.list(r[[1]])) r <- r[[1]]
        contributions[i] <- r$contribution
      }

      for (i in seq_along(server_list)) {
        ci <- which(server_names == server_list[i])
        DSI::datashield.aggregate(
          conns = datasources[ci],
          expr = call("mheCoinTossDeriveCRPDS",
                      contributions = unname(contributions),
                      commitments = unname(commitments),
                      log_n = as.integer(log_n),
                      log_scale = as.integer(log_scale),
                      party_id = as.integer(i - 1L),
                      session_id = session_id)
        )
      }

      if (verbose) message("  CRP derived from ", length(server_list), "-party coin-toss")

      # ---- Key Setup (all from coin-toss CRP) ----
      conn_idx <- which(server_names == server_list[1])
      result0 <- .dsAgg(
        conns = datasources[conn_idx],
        expr = call("mheInitDS",
                    party_id = 0L, from_storage = TRUE,
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
  if (length(non_label_servers) > 0 && !use_k2_beaver) {
    if (verbose) message("\n[Phase 2] Encrypting response variable...")

    # Encrypt STANDARDIZED y (Gaussian) or raw y (non-Gaussian)
    enc_data <- if (standardize_y) std_data else data_name

    conn_idx <- which(server_names == y_server)
    enc_result <- .dsAgg(
      conns = datasources[conn_idx],
      expr = call("mheEncryptRawDS",
                  data_name = enc_data, y_var = y_var,
                  store_local = FALSE,
                  session_id = session_id)
    )
    if (is.list(enc_result)) enc_result <- enc_result[[1]]
    ct_y <- enc_result$encrypted_y
    if (verbose) message("  Encrypted y on ", y_server)
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

  # K=2 MPC: uses Beaver-triple protocol with L-BFGS optimizer (below in Phase 3)

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
    k3_result <- .k3_secure_agg_loop(
      datasources = datasources, server_list = server_list,
      server_names = server_names, x_vars = x_vars,
      coordinator = coordinator, coordinator_conn = coordinator_conn,
      non_label_servers = non_label_servers, transport_pks = transport_pks,
      std_data = std_data, y_var = y_var, family = family,
      betas = betas, n_obs = n_obs, lambda = lambda,
      log_n = log_n, log_scale = log_scale, session_id = session_id,
      max_iter = max_iter, tol = tol, verbose = verbose,
      topology = topology, label_intercept = label_intercept,
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
  } else if (use_k2_beaver) {
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
