#' @title Privacy-Preserving Correlation for Vertically Partitioned Data
#' @description Computes the Pearson correlation matrix for vertically partitioned
#'   data using Multiparty Homomorphic Encryption (MHE) with threshold decryption.
#'
#' @param data_name Character string. Name of the (aligned) data frame on
#'   each server.
#' @param variables A named list where each name corresponds to a server name
#'   and each element is a character vector of variable names from that server.
#' @param log_n Integer. CKKS ring dimension parameter (12, 13, or 14).
#'   Controls the number of slots: N/2 = 2^(logN-1).
#'   Default is 12 (2048 slots, fast). Use 13 or 14 for larger datasets.
#' @param log_scale Integer. CKKS scale parameter controlling precision.
#'   Default is 40 (approximately 12 decimal digits of precision).
#' @param datasources DataSHIELD connection object or list of connections.
#'   If NULL, uses all available connections.
#'
#' @return A list with class \code{"ds.cor"} containing:
#'   \itemize{
#'     \item \code{correlation}: The full correlation matrix (p x p)
#'     \item \code{var_names}: Variable names in order
#'     \item \code{n_obs}: Number of observations
#'     \item \code{method}: \code{"MHE-CKKS-Threshold"} indicating the method used
#'     \item \code{servers}: Names of servers involved
#'     \item \code{local_correlations}: List of within-server correlation matrices
#'     \item \code{cross_correlations}: List of cross-server correlation matrices
#'     \item \code{mhe_params}: List with \code{log_n} and \code{log_scale} used
#'   }
#'
#' @details
#' \subsection{Threshold MHE Protocol (6 phases)}{
#' This function implements a full multiparty homomorphic encryption protocol
#' with threshold decryption:
#'
#' \enumerate{
#'   \item \strong{Key Generation}: Each server generates its own secret key share
#'     and a public key share. Party 0 also generates the Common Reference
#'     Polynomial (CRP) shared by all parties.
#'   \item \strong{Key Combination}: Public key shares are combined into a
#'     Collective Public Key (CPK). Data encrypted under the CPK can only be
#'     decrypted with ALL servers cooperating.
#'   \item \strong{Encryption}: Each server standardizes its columns (Z-scores)
#'     and encrypts them column-by-column under the CPK.
#'   \item \strong{Local Correlation}: Within-server correlations are computed in
#'     plaintext (no encryption needed for data that stays on-server).
#'   \item \strong{Cross-Server Correlation}: For each pair of servers (A, B),
#'     server A receives Enc(Z_B) and computes the element-wise product
#'     Z_A * Enc(Z_B) homomorphically. The result is still encrypted and
#'     requires threshold decryption. Each server produces a partial decryption
#'     share, and the client fuses all shares to recover the inner product.
#'   \item \strong{Assembly}: The client assembles the full p x p correlation
#'     matrix from local correlations (diagonal blocks) and cross-server
#'     correlations (off-diagonal blocks).
#' }
#' }
#'
#' \subsection{Security Guarantees}{
#' \describe{
#'   \item{Client privacy}{The client (researcher) CANNOT decrypt any individual
#'     data. It only receives partial decryption shares that are useless alone.
#'     Only the final aggregate statistics (correlation coefficients) are revealed
#'     after fusing ALL shares.}
#'   \item{Server privacy}{Each server's raw data never leaves the server. Other
#'     servers only see encrypted columns (opaque ciphertexts).}
#'   \item{Collusion resistance}{Even K-1 colluding servers cannot decrypt
#'     without the K-th server's key share. Full decryption requires cooperation
#'     from ALL K servers.}
#' }
#' }
#'
#' @section Performance Notes:
#' \itemize{
#'   \item \code{log_n = 12}: Up to 2048 observations, fastest
#'   \item \code{log_n = 13}: Up to 4096 observations
#'   \item \code{log_n = 14}: Up to 8192 observations, slowest
#'   \item Precision: approximately 10^-3 to 10^-4 error due to CKKS approximation
#'   \item The dominant cost is the threshold decryption loop (Phase 5), which
#'     requires one round-trip per server per cross-correlation element.
#' }
#'
#' @seealso \code{\link{ds.vertPCA}} for PCA analysis built on this function
#'
#' @examples
#' \dontrun{
#' # Connect to Opal/DataSHIELD servers
#' connections <- DSI::datashield.login(builder$build())
#'
#' # Align records across servers using ECDH-PSI
#' ds.psiAlign("D", "patient_id", "D_aligned", datasources = connections)
#'
#' # Define which variables are on which server
#' vars <- list(
#'   hospital_A = c("age", "bmi"),
#'   hospital_B = c("glucose", "systolic_bp")
#' )
#'
#' # Compute privacy-preserving correlation
#' result <- ds.vertCor("D_aligned", vars, datasources = connections)
#' print(result)
#' }
#'
#' @references
#' Mouchet, C. et al. (2021). "Multiparty Homomorphic Encryption from
#' Ring-Learning-With-Errors". \emph{Proceedings on Privacy Enhancing Technologies (PETS)}.
#'
#' Cheon, J.H. et al. (2017). "Homomorphic Encryption for Arithmetic of
#' Approximate Numbers". \emph{ASIACRYPT 2017}.
#'
#' @importFrom DSI datashield.aggregate datashield.connections_find
#' @export
ds.vertCor <- function(data_name, variables, log_n = 12, log_scale = 40,
                       datasources = NULL, reuse_mhe = FALSE) {

  # ===========================================================================
  # Input Validation
  # ===========================================================================
  if (!is.character(data_name) || length(data_name) != 1) {
    stop("data_name must be a single character string", call. = FALSE)
  }
  if (!is.list(variables) || is.null(names(variables))) {
    stop("variables must be a named list mapping server names to variable vectors",
         call. = FALSE)
  }
  if (length(variables) < 2) {
    stop("At least 2 servers required for cross-server correlation. ",
         "For single-server correlation, use standard cor() via datashield.aggregate.",
         call. = FALSE)
  }

  # Generate session_id for MHE protocol state isolation
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

  # Guaranteed cleanup: if the function errors at any point, clean up server-side
  # cryptographic state. This prevents session leaks that accumulate memory.
  on.exit({
    for (.srv in server_list) {
      .ci <- which(server_names == .srv)
      tryCatch(
        DSI::datashield.aggregate(conns = datasources[.ci],
          expr = call("mheCleanupDS", session_id = session_id)),
        error = function(e) NULL)
    }
    .dsvert_reset_chunk_size()
  }, add = TRUE)

  # Helpers -------------------------------------------------------------------

  # Send CT chunks to a server (adaptive chunk size)
  .sendCTChunks <- function(ct_str, conn_idx) {
    .dsvert_adaptive_send(ct_str, function(chunk_str, chunk_idx, n_chunks) {
      DSI::datashield.aggregate(
        conns = datasources[conn_idx],
        expr = call("mheStoreCTChunkDS",
                    chunk_index = chunk_idx,
                    chunk = chunk_str,
                    session_id = session_id)
      )
    })
  }

  # Store a large base64url blob on a server with adaptive chunking and fallback
  .storeLargeBlob <- function(key, data, conn_idx) {
    .dsvert_adaptive_send(data, function(chunk_str, chunk_idx, n_chunks) {
      if (n_chunks == 1L) {
        DSI::datashield.aggregate(
          conns = datasources[conn_idx],
          expr = call("mheStoreBlobDS", key = key, chunk = chunk_str,
                      session_id = session_id)
        )
      } else {
        DSI::datashield.aggregate(
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

  # Relay a wrapped share to the fusion server in chunks (adaptive chunk size)
  .sendWrappedShare <- function(wrapped_share, party_id, fusion_conn_idx) {
    chunk_size <- .dsvert_get_chunk_size()
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

  # ===========================================================================
  # Setup
  # ===========================================================================
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  if (length(datasources) == 0) {
    stop("No DataSHIELD connections found", call. = FALSE)
  }

  server_names <- names(datasources)
  server_list <- names(variables)

  missing_servers <- setdiff(server_list, server_names)
  if (length(missing_servers) > 0) {
    stop("Servers not found in connections: ",
         paste(missing_servers, collapse = ", "), call. = FALSE)
  }

  message("=== MHE Correlation with Threshold Decryption ===")
  message("Servers: ", paste(server_list, collapse = ", "))

  # ===========================================================================
  # Phase 1: Key Generation
  # ===========================================================================
  # Each server generates an RLWE secret key share sk_k and a public key share.
  # Party 0 (first server) additionally generates:
  #   - CRP (Common Reference Polynomial): shared randomness needed by all
  #     parties to construct compatible public key shares.
  #   - GKG seed: shared seed for Galois Key Generation, enabling all parties
  #     to independently generate matching Galois rotation keys (needed for
  #     ciphertext slot rotations in inner-product computation).
  # Get observation count (always needed, regardless of key reuse)
  n_obs <- NULL
  for (server in server_list) {
    conn_idx <- which(server_names == server)
    call_expr <- call("mheGetObsDS",
                      data_name = data_name,
                      variables = variables[[server]],
                      session_id = session_id)
    count <- DSI::datashield.aggregate(conns = datasources[conn_idx], expr = call_expr)
    if (is.list(count)) count <- count[[1]]

    if (is.null(n_obs)) {
      n_obs <- count
    } else if (count != n_obs) {
      stop("Observation count mismatch: ", server, " has ", count,
           ", expected ", n_obs, ". Ensure data is aligned.", call. = FALSE)
    }
  }
  message("  Observations: ", n_obs)

  first_server <- server_list[1]
  transport_pks <- list()
  mhe_reused <- FALSE

  # --- MHE Context Reuse Check ---
  # Compute a canonical context_id from peer set + params. If all servers
  # have a cached session with matching context_id, skip key generation
  # and combination entirely (only fresh transport keys are generated).
  # NOTE: Use '_' as separator — Opal's R expression parser rejects many
  # special characters ('|', '/', etc.) inside string arguments passed via
  # datashield.aggregate(). Only alphanumeric, '_', ',', ':', '-' are safe.
  context_id <- paste0(
    "peers_", paste(sort(server_list), collapse = ","),
    "_logn_", log_n,
    "_logscale_", log_scale,
    "_numobs_", n_obs,
    "_rlk_false"
  )

  if (isTRUE(reuse_mhe)) {
    reuse_results <- list()
    all_reusable <- TRUE
    for (server in server_list) {
      conn_idx <- which(server_names == server)
      res <- DSI::datashield.aggregate(
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
      message("\n[Phase 1-2] MHE keys reused from cached context (skipped)")
    }
  }

  if (!mhe_reused) {
    # ===========================================================================
    # Phase 1: Key Generation
    # ===========================================================================
    message("\n[Phase 1] Key generation on each server...")

    # Party 0 generates CRP and GKG seed
    conn_idx <- which(server_names == first_server)

    call_expr <- call("mheInitDS",
                      party_id = 0L,
                      crp = NULL,
                      gkg_seed = NULL,
                      num_obs = as.integer(n_obs),
                      log_n = as.integer(log_n),
                      log_scale = as.integer(log_scale),
                      session_id = session_id)

    result0 <- DSI::datashield.aggregate(conns = datasources[conn_idx], expr = call_expr)
    if (is.list(result0)) result0 <- result0[[1]]
    crp <- result0$crp
    gkg_seed <- result0$gkg_seed
    message("  ", first_server, ": Party 0 initialized (CRP generated)")

    # Other parties initialize with CRP and shared GKG seed
    pk_shares <- list()
    gkg_shares <- list()
    pk_shares[[first_server]] <- result0$public_key_share
    gkg_shares[[first_server]] <- result0$galois_key_shares
    transport_pks[[first_server]] <- result0$transport_pk

    # Send CRP, GKG seed, and party_id to all non-party0 servers (sequential
    # blob transfer, then parallel mheInitDS calls).
    other_servers <- server_list[-1]
    other_conn_idxs <- integer(length(other_servers))
    for (i in seq_along(other_servers)) {
      server <- other_servers[i]
      ci <- which(server_names == server)
      other_conn_idxs[i] <- ci
      .storeLargeBlob("crp", crp, ci)
      .storeLargeBlob("gkg_seed", gkg_seed, ci)
      .storeLargeBlob("party_id", as.character(i), ci)
    }

    # Launch mheInitDS on all non-party0 servers in parallel (same expression,
    # each server reads its own party_id from blob storage).
    init_expr <- call("mheInitDS",
                      party_id = 0L,
                      from_storage = TRUE,
                      num_obs = as.integer(n_obs),
                      log_n = as.integer(log_n),
                      log_scale = as.integer(log_scale),
                      session_id = session_id)

    init_results <- DSI::datashield.aggregate(
      conns = datasources[other_servers],
      expr = init_expr
    )

    for (i in seq_along(other_servers)) {
      server <- other_servers[i]
      result <- init_results[[server]]
      pk_shares[[server]] <- result$public_key_share
      gkg_shares[[server]] <- result$galois_key_shares
      transport_pks[[server]] <- result$transport_pk
      message("  ", server, ": Party ", i, " initialized")
    }

    # ===========================================================================
    # Phase 2: Combine public key shares into CPK
    # ===========================================================================
    # Public key shares are added together to form the Collective Public Key (CPK).
    # The CPK has a critical property: data encrypted under it can ONLY be
    # decrypted when ALL K servers cooperate (K-of-K threshold). No subset of
    # K-1 servers (or the client alone) can decrypt. The combining is done on
    # one server; the resulting CPK + Galois keys are distributed to all others.
    message("\n[Phase 2] Combining public key shares...")

    conn_idx <- which(server_names == first_server)

    # Store PK shares, CRP, GKG seed, and GKG shares on the combining server
    # via chunked blob storage. Cryptographic objects can be several MB each,
    # exceeding R's expression parser stack limit if passed as call arguments.
    for (i in seq_along(server_list)) {
      .storeLargeBlob(paste0("pk_", i - 1), pk_shares[[server_list[i]]], conn_idx)
    }
    .storeLargeBlob("crp", crp, conn_idx)
    .storeLargeBlob("gkg_seed", gkg_seed, conn_idx)

    n_gkg_shares <- length(gkg_shares[[server_list[1]]])
    for (i in seq_along(server_list)) {
      shares <- gkg_shares[[server_list[i]]]
      for (j in seq_along(shares)) {
        .storeLargeBlob(paste0("gkg_", i - 1, "_", j - 1), shares[j], conn_idx)
      }
    }

    call_expr <- call("mheCombineDS",
                      from_storage = TRUE,
                      n_parties = as.integer(length(server_list)),
                      n_gkg_shares = as.integer(n_gkg_shares),
                      num_obs = as.integer(n_obs),
                      log_n = as.integer(log_n),
                      log_scale = as.integer(log_scale),
                      session_id = session_id)

    combined <- DSI::datashield.aggregate(conns = datasources[conn_idx], expr = call_expr)
    if (is.list(combined)) combined <- combined[[1]]
    cpk <- combined$collective_public_key
    galois_keys <- combined$galois_keys
    message("  CPK created on ", first_server)

    # Distribute CPK and Galois keys to other servers via chunked blob storage
    other_srv_conns <- integer(length(server_list) - 1)
    for (si in seq_along(server_list[-1])) {
      server <- server_list[-1][si]
      srv_conn <- which(server_names == server)
      other_srv_conns[si] <- srv_conn
      .storeLargeBlob("cpk", cpk, srv_conn)
      if (!is.null(galois_keys)) {
        for (gk_i in seq_along(galois_keys)) {
          .storeLargeBlob(paste0("gk_", gk_i - 1), galois_keys[gk_i], srv_conn)
        }
      }
    }
    # Store CPK on all non-party0 servers in parallel
    DSI::datashield.aggregate(
      conns = datasources[server_list[-1]],
      expr = call("mheStoreCPKDS", from_storage = TRUE,
                  session_id = session_id)
    )
    for (server in server_list[-1]) message("  CPK stored on ", server)

    # Register context for future reuse (parallel across all servers)
    DSI::datashield.aggregate(
      conns = datasources[server_list],
      expr = call("mheRegisterContextDS",
                  context_id = context_id,
                  session_id = session_id)
    )
  }  # end if (!mhe_reused)

  # Distribute transport keys for share-wrapping protocol (always needed -
  # fresh transport keys are generated even when MHE keys are reused).
  # The "fusion" key identifies the fusion server (Party 0) whose transport PK
  # is used by non-fusion servers to wrap their partial decryption shares.
  fusion_server <- first_server
  fusion_conn_idx <- which(server_names == fusion_server)

  for (server in server_list) {
    conn_idx <- which(server_names == server)
    tk_map <- list(fusion = transport_pks[[fusion_server]])
    # Also store all server PKs (needed for GLM secure routing)
    for (s in server_list) {
      tk_map[[s]] <- transport_pks[[s]]
    }
    DSI::datashield.aggregate(
      conns = datasources[conn_idx],
      expr = call("mheStoreTransportKeysDS", transport_keys = tk_map,
                  session_id = session_id)
    )
  }
  message("  Transport keys distributed for share-wrapping")

  # ===========================================================================
  # Phase 3: Encrypt data on each server
  # ===========================================================================
  message("\n[Phase 3] Encrypting data on each server...")

  server_encrypted <- list()
  for (server in server_list) {
    conn_idx <- which(server_names == server)

    call_expr <- call("mheEncryptLocalDS",
                      data_name = data_name,
                      variables = variables[[server]],
                      session_id = session_id)

    result <- DSI::datashield.aggregate(conns = datasources[conn_idx], expr = call_expr)
    if (is.list(result)) result <- result[[1]]
    server_encrypted[[server]] <- result
    message("  ", server, ": ", length(result$encrypted_columns), " variables encrypted")
  }

  # ===========================================================================
  # Phase 4: Compute local correlations (no MHE needed)
  # ===========================================================================
  message("\n[Phase 4] Computing local correlations...")

  local_cors <- list()
  for (server in server_list) {
    conn_idx <- which(server_names == server)
    call_expr <- call("localCorDS",
                      data_name = data_name,
                      variables = variables[[server]],
                      session_id = session_id)

    result <- DSI::datashield.aggregate(conns = datasources[conn_idx], expr = call_expr)
    if (is.list(result)) result <- result[[1]]
    local_cors[[server]] <- result$correlation
    message("  ", server, ": ", nrow(result$correlation), "x", ncol(result$correlation))
  }

  # ===========================================================================
  # Phase 5: Cross-server correlations with SHARE-WRAPPED THRESHOLD DECRYPTION
  # ===========================================================================
  # This is the core privacy-preserving computation. For each server pair
  # (A, B), server A computes z_A * Enc(z_B) homomorphically -- the result
  # is still encrypted. To recover the inner product, we use share-wrapped
  # threshold decryption: each non-fusion server computes a partial decryption
  # share, wraps it under the fusion server's transport PK, and the client
  # relays the opaque wrapped shares to the fusion server. The fusion server
  # unwraps all shares, computes its own share, and fuses them to recover
  # the plaintext. The client NEVER sees raw shares or unsanitized plaintext.
  message("\n[Phase 5] Computing cross-server correlations (share-wrapped threshold decryption)...")

  cross_cors <- list()
  non_fusion_servers <- setdiff(server_list, fusion_server)

  for (i in 1:(length(server_list) - 1)) {
    for (j in (i + 1):length(server_list)) {
      server_A <- server_list[i]
      server_B <- server_list[j]

      message("  ", server_A, " x ", server_B, ":")

      # Transfer encrypted columns from server_B to server_A using chunking
      conn_idx_A <- which(server_names == server_A)
      enc_cols_B <- server_encrypted[[server_B]]$encrypted_columns

      for (col_k in seq_along(enc_cols_B)) {
        col_str <- enc_cols_B[[col_k]]
        n_chunks_sent <- .dsvert_adaptive_send(col_str, function(chunk_str, chunk_idx, n_chunks) {
          DSI::datashield.aggregate(
            conns = datasources[conn_idx_A],
            expr = call("mheStoreEncChunkDS",
                        col_index = as.integer(col_k),
                        chunk_index = chunk_idx,
                        chunk = chunk_str,
                        session_id = session_id)
          )
        })

        # Assemble the column from chunks
        assemble_expr <- call("mheAssembleEncColumnDS",
                              col_index = as.integer(col_k),
                              n_chunks = as.integer(n_chunks_sent),
                              session_id = session_id)
        DSI::datashield.aggregate(conns = datasources[conn_idx_A], expr = assemble_expr)
      }

      # Compute encrypted cross-product on server_A
      call_expr <- call("mheCrossProductEncDS",
                        data_name = data_name,
                        variables = variables[[server_A]],
                        n_enc_cols = as.integer(length(enc_cols_B)),
                        n_obs = as.integer(n_obs),
                        session_id = session_id)

      enc_result <- DSI::datashield.aggregate(conns = datasources[conn_idx_A], expr = call_expr)
      if (is.list(enc_result)) enc_result <- enc_result[[1]]

      n_A <- length(variables[[server_A]])
      n_B <- length(variables[[server_B]])
      message("    Encrypted cross-product: ", n_A, "x", n_B)

      # Protocol Firewall: authorize ciphertexts on non-producing servers.
      # ct_hashes scales with p_A * p_B; sent via blob storage to handle
      # cases where many variables produce a large hash vector.
      ct_hashes <- enc_result$ct_hashes
      if (!is.null(ct_hashes) && length(ct_hashes) > 0) {
        ct_hashes_blob <- paste(ct_hashes, collapse = ",")
        for (server in server_list) {
          if (server != server_A) {
            conn_idx <- which(server_names == server)
            .storeLargeBlob("ct_hashes", ct_hashes_blob, conn_idx)
            DSI::datashield.aggregate(
              conns = datasources[conn_idx],
              expr = call("mheAuthorizeCTDS",
                          op_type = "cross-product",
                          from_storage = TRUE,
                          session_id = session_id))
          }
        }
      }

      # Batched share-wrapped threshold decryption.
      # Instead of one round-trip per CT per server, we send ALL CTs for this
      # server pair as blobs, then call batch decrypt once per server.
      er <- enc_result$encrypted_results
      n_cts <- n_A * n_B
      all_cts <- character(n_cts)

      for (row_i in 1:n_A) {
        for (col_j in 1:n_B) {
          idx <- (row_i - 1) * n_B + col_j
          if (is.matrix(er)) {
            all_cts[idx] <- er[row_i, col_j]
          } else if (is.list(er) && length(er) >= row_i) {
            all_cts[idx] <- er[[row_i]][col_j]
          } else {
            all_cts[idx] <- er[[idx]]
          }
        }
      }

      # Step 1: Send all CTs to each non-fusion server as batch blobs,
      # then call mhePartialDecryptBatchWrappedDS once per server.
      for (nf_server in non_fusion_servers) {
        nf_conn <- which(server_names == nf_server)
        nf_party_id <- which(server_list == nf_server) - 1

        for (ct_i in seq_len(n_cts)) {
          .storeLargeBlob(paste0("ct_batch_", ct_i), all_cts[ct_i], nf_conn)
        }

        pd <- DSI::datashield.aggregate(
          conns = datasources[nf_conn],
          expr = call("mhePartialDecryptBatchWrappedDS",
                      n_cts = as.integer(n_cts),
                      session_id = session_id)
        )
        if (is.list(pd)) pd <- pd[[1]]

        # Relay each wrapped share to fusion server
        for (ct_i in seq_len(n_cts)) {
          share_key <- paste0("wrapped_share_", nf_party_id, "_ct_", ct_i)
          .storeLargeBlob(share_key, pd$wrapped_shares[ct_i], fusion_conn_idx)
        }
      }

      # Step 2: Send all CTs to fusion server, then fuse all at once.
      for (ct_i in seq_len(n_cts)) {
        .storeLargeBlob(paste0("ct_batch_", ct_i), all_cts[ct_i], fusion_conn_idx)
      }

      fuse_result <- DSI::datashield.aggregate(
        conns = datasources[fusion_conn_idx],
        expr = call("mheFuseBatchDS",
                    n_cts = as.integer(n_cts),
                    n_parties = as.integer(length(server_list)),
                    num_slots = as.integer(n_obs),
                    session_id = session_id)
      )
      if (is.list(fuse_result)) fuse_result <- fuse_result[[1]]

      result_matrix <- matrix(fuse_result$values, nrow = n_A, ncol = n_B, byrow = TRUE)

      # The homomorphic computation yields sum(z_A * z_B) -- the raw inner
      # product of standardized columns. Dividing by (n-1) converts to
      # Pearson correlation (same as cor() on standardized data).
      cross_cor <- result_matrix / (n_obs - 1)
      rownames(cross_cor) <- variables[[server_A]]
      colnames(cross_cor) <- variables[[server_B]]
      cross_cors[[paste(server_A, server_B, sep = "_")]] <- cross_cor

      message("    Threshold decryption complete (server-side fusion)")
    }
  }

  # ===========================================================================
  # Phase 6: Assemble correlation matrix
  # ===========================================================================
  message("\n[Phase 6] Assembling correlation matrix...")

  all_vars <- unlist(lapply(server_list, function(s) variables[[s]]))
  p <- length(all_vars)
  full_cor <- matrix(NA, nrow = p, ncol = p)
  rownames(full_cor) <- all_vars
  colnames(full_cor) <- all_vars

  # Fill diagonal blocks (local correlations)
  start_idx <- 1
  for (server in server_list) {
    server_vars <- variables[[server]]
    n_vars <- length(server_vars)
    end_idx <- start_idx + n_vars - 1

    full_cor[start_idx:end_idx, start_idx:end_idx] <- local_cors[[server]]
    start_idx <- end_idx + 1
  }

  # Fill off-diagonal blocks (cross-server correlations)
  for (i in 1:(length(server_list) - 1)) {
    for (j in (i + 1):length(server_list)) {
      key <- paste(server_list[i], server_list[j], sep = "_")
      if (!is.null(cross_cors[[key]])) {
        vars_i <- variables[[server_list[i]]]
        vars_j <- variables[[server_list[j]]]
        full_cor[vars_i, vars_j] <- cross_cors[[key]]
        full_cor[vars_j, vars_i] <- t(cross_cors[[key]])
      }
    }
  }

  message("Correlation matrix: ", p, "x", p)

  result <- list(
    correlation = full_cor,
    var_names = all_vars,
    n_obs = n_obs,
    method = "MHE-CKKS-Threshold",
    servers = server_list,
    local_correlations = local_cors,
    cross_correlations = cross_cors,
    mhe_params = list(log_n = log_n, log_scale = log_scale)
  )

  # Cleanup is handled by on.exit() — runs whether we succeed or fail

  class(result) <- c("ds.cor", "list")
  return(result)
}

#' @title Print Method for ds.cor Objects
#' @description Prints a summary of correlation results.
#' @param x A ds.cor object
#' @param digits Number of digits to display
#' @param ... Additional arguments (ignored)
#' @export
print.ds.cor <- function(x, digits = 4, ...) {
  cat("Correlation Matrix (", x$method, ")\n", sep = "")
  cat("=====================================\n\n")
  cat("Observations:", x$n_obs, "\n")
  cat("Variables:", length(x$var_names), "\n")
  if (!is.null(x$servers)) {
    cat("Servers:", paste(x$servers, collapse = ", "), "\n")
  }
  cat("\n")
  print(round(x$correlation, digits))
  invisible(x)
}
