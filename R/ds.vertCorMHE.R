#' @title Vertical Correlation with Full MHE (Threshold Decryption)
#' @description Computes correlation matrix across vertically partitioned data
#'   using full Multiparty Homomorphic Encryption. The client CANNOT decrypt
#'   individual data - only the aggregated result after ALL servers cooperate.
#'
#' @param data_name Character. Name of the data frame on each server.
#' @param variables Named list. Variable names for each server.
#' @param log_n Integer. CKKS parameter (default 12 for speed)
#' @param log_scale Integer. CKKS scale parameter (default 40)
#' @param datasources DSConnection objects
#'
#' @return List with correlation matrix
#' @export
ds.vertCorMHE <- function(data_name, variables, log_n = 12, log_scale = 40,
                          datasources = NULL) {

  # Helper: Convert base64url to standard base64 (with padding)
  .b64url_to_b64 <- function(x) {
    x <- gsub("-", "+", x, fixed = TRUE)
    x <- gsub("_", "/", x, fixed = TRUE)
    pad <- (4 - nchar(x) %% 4) %% 4
    if (pad > 0) x <- paste0(x, paste(rep("=", pad), collapse = ""))
    x
  }

  # Helper: Fuse partial decryptions locally (on client)
  .mheFuseLocal <- function(ciphertext, decryption_shares, log_n, log_scale, num_slots = 0) {
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

    status <- system2(bin_path, "mhe-fuse", stdin = input_file,
                      stdout = output_file, stderr = output_file)

    output <- jsonlite::read_json(output_file, simplifyVector = TRUE)
    if (!is.null(output$error) && nzchar(output$error)) {
      stop("mhe-fuse error: ", output$error)
    }

    output$value
  }

  # Find mhe-tool binary for client-side operations
  .findMheTool <- function() {
    # Check package installation
    bin_path <- system.file("mhe-tool", "mhe-tool", package = "dsVertClient")
    if (bin_path != "" && file.exists(bin_path)) return(bin_path)

    # Check dsVert package
    bin_path <- system.file("mhe-tool", "mhe-tool", package = "dsVert")
    if (bin_path != "" && file.exists(bin_path)) return(bin_path)

    # Check environment variable
    env_path <- Sys.getenv("DSVERT_MHE_TOOL")
    if (env_path != "" && file.exists(env_path)) return(env_path)

    # Check common locations
    common_paths <- c(
      file.path(getwd(), "mhe-tool"),
      "/usr/local/bin/mhe-tool",
      "~/.local/bin/mhe-tool"
    )
    for (p in common_paths) {
      if (file.exists(p)) return(p)
    }

    stop("mhe-tool binary not found. Install dsVert or set DSVERT_MHE_TOOL env var.")
  }

  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }

  server_names <- names(datasources)
  server_list <- names(variables)

  if (!all(server_list %in% server_names)) {
    stop("Variable list contains servers not in datasources")
  }

  message("=== MHE Correlation with Threshold Decryption ===")
  message("Servers: ", paste(server_list, collapse = ", "))

  # =========================================================================
  # Phase 1: Each server initializes and generates key shares
  # =========================================================================
  message("\n[Phase 1] Key generation on each server...")

  # Get observation count first
  n_obs <- NULL
  for (server in server_list) {
    conn_idx <- which(server_names == server)
    call_expr <- call("mheGetObsDS",
                      data_name = data_name,
                      variables = variables[[server]])
    count <- DSI::datashield.aggregate(conns = datasources[conn_idx], expr = call_expr)
    if (is.list(count)) count <- count[[1]]

    if (is.null(n_obs)) {
      n_obs <- count
    } else if (count != n_obs) {
      stop("Observation count mismatch: ", server, " has ", count, ", expected ", n_obs)
    }
  }
  message("  Observations: ", n_obs)

  # Party 0 generates CRP
  first_server <- server_list[1]
  conn_idx <- which(server_names == first_server)

  call_expr <- call("mheInitDS",
                    party_id = 0L,
                    crp = NULL,
                    num_obs = as.integer(n_obs),
                    log_n = as.integer(log_n),
                    log_scale = as.integer(log_scale))

  result0 <- DSI::datashield.aggregate(conns = datasources[conn_idx], expr = call_expr)
  if (is.list(result0)) result0 <- result0[[1]]
  crp <- result0$crp
  message("  ", first_server, ": Party 0 initialized (CRP generated)")

  # Other parties initialize with CRP
  pk_shares <- list()
  pk_shares[[first_server]] <- result0$public_key_share

  for (i in 2:length(server_list)) {
    server <- server_list[i]
    conn_idx <- which(server_names == server)

    call_expr <- call("mheInitDS",
                      party_id = as.integer(i - 1),
                      crp = crp,
                      num_obs = as.integer(n_obs),
                      log_n = as.integer(log_n),
                      log_scale = as.integer(log_scale))

    result <- DSI::datashield.aggregate(conns = datasources[conn_idx], expr = call_expr)
    if (is.list(result)) result <- result[[1]]
    pk_shares[[server]] <- result$public_key_share
    message("  ", server, ": Party ", i-1, " initialized")
  }

  # =========================================================================
  # Phase 2: Combine public key shares into CPK
  # =========================================================================
  message("\n[Phase 2] Combining public key shares...")

  conn_idx <- which(server_names == first_server)
  call_expr <- call("mheCombineDS",
                    public_key_shares = unlist(pk_shares, use.names = FALSE),
                    crp = crp,
                    num_obs = as.integer(n_obs),
                    log_n = as.integer(log_n),
                    log_scale = as.integer(log_scale))

  combined <- DSI::datashield.aggregate(conns = datasources[conn_idx], expr = call_expr)
  if (is.list(combined)) combined <- combined[[1]]
  cpk <- combined$collective_public_key
  message("  CPK created on ", first_server)

  # Distribute CPK to other servers
  for (server in server_list[-1]) {
    conn_idx <- which(server_names == server)
    call_expr <- call("mheStoreCPKDS", cpk = cpk)
    DSI::datashield.aggregate(conns = datasources[conn_idx], expr = call_expr)
    message("  CPK stored on ", server)
  }

  # =========================================================================
  # Phase 3: Encrypt data on each server
  # =========================================================================
  message("\n[Phase 3] Encrypting data on each server...")

  server_encrypted <- list()
  for (server in server_list) {
    conn_idx <- which(server_names == server)

    call_expr <- call("mheEncryptLocalDS",
                      data_name = data_name,
                      variables = variables[[server]])

    result <- DSI::datashield.aggregate(conns = datasources[conn_idx], expr = call_expr)
    if (is.list(result)) result <- result[[1]]
    server_encrypted[[server]] <- result
    message("  ", server, ": ", length(result$encrypted_columns), " variables encrypted")
  }

  # =========================================================================
  # Phase 4: Compute local correlations (no MHE needed)
  # =========================================================================
  message("\n[Phase 4] Computing local correlations...")

  local_cors <- list()
  for (server in server_list) {
    conn_idx <- which(server_names == server)
    call_expr <- call("localCorDS",
                      data_name = data_name,
                      variables = variables[[server]])

    result <- DSI::datashield.aggregate(conns = datasources[conn_idx], expr = call_expr)
    if (is.list(result)) result <- result[[1]]
    local_cors[[server]] <- result$correlation
    message("  ", server, ": ", nrow(result$correlation), "x", ncol(result$correlation))
  }

  # =========================================================================
  # Phase 5: Cross-server correlations with THRESHOLD DECRYPTION
  # =========================================================================
  message("\n[Phase 5] Computing cross-server correlations (threshold decryption)...")

  cross_cors <- list()

  # The server that has evaluation keys (first_server) computes cross-products
  # Results are ENCRYPTED - each server partial-decrypts, then client fuses

  for (i in 1:(length(server_list) - 1)) {
    for (j in (i + 1):length(server_list)) {
      server_A <- server_list[i]
      server_B <- server_list[j]

      message("  ", server_A, " x ", server_B, ":")

      # Transfer encrypted columns from server_B to server_A using chunking
      # (cross-product runs on server_A which has the plaintext data)
      conn_idx_A <- which(server_names == server_A)
      enc_cols_B <- server_encrypted[[server_B]]$encrypted_columns
      chunk_size <- 10000  # ~10KB per chunk, safe for R parser

      for (col_k in seq_along(enc_cols_B)) {
        col_str <- enc_cols_B[[col_k]]
        n_chars <- nchar(col_str)
        chunks <- ceiling(n_chars / chunk_size)

        for (ch in seq_len(chunks)) {
          start <- (ch - 1) * chunk_size + 1
          end <- min(ch * chunk_size, n_chars)
          chunk_str <- substr(col_str, start, end)

          store_expr <- call("mheStoreEncChunkDS",
                             col_index = as.integer(col_k),
                             chunk_index = as.integer(ch),
                             chunk = chunk_str)
          DSI::datashield.aggregate(conns = datasources[conn_idx_A], expr = store_expr)
        }

        # Assemble the column from chunks
        assemble_expr <- call("mheAssembleEncColumnDS",
                              col_index = as.integer(col_k),
                              n_chunks = as.integer(chunks))
        DSI::datashield.aggregate(conns = datasources[conn_idx_A], expr = assemble_expr)
      }

      # Compute encrypted cross-product on server_A (has plaintext data, no eval keys needed)
      call_expr <- call("mheCrossProductEncDS",
                        data_name = data_name,
                        variables = variables[[server_A]],
                        n_enc_cols = as.integer(length(enc_cols_B)),
                        n_obs = as.integer(n_obs))

      enc_result <- DSI::datashield.aggregate(conns = datasources[conn_idx_A], expr = call_expr)
      if (is.list(enc_result)) enc_result <- enc_result[[1]]

      n_A <- length(variables[[server_A]])
      n_B <- length(variables[[server_B]])
      message("    Encrypted cross-product: ", n_A, "x", n_B)

      # For each encrypted result, do threshold decryption
      result_matrix <- matrix(NA, nrow = n_A, ncol = n_B)

      for (row_i in 1:n_A) {
        for (col_j in 1:n_B) {
          # Handle both matrix and list-of-lists return formats
          er <- enc_result$encrypted_results
          if (is.matrix(er)) {
            ct <- er[row_i, col_j]
          } else if (is.list(er) && length(er) >= row_i) {
            ct <- er[[row_i]][col_j]
          } else {
            # Flat vector fallback: index = (row-1)*n_B + col
            ct <- er[[(row_i - 1) * n_B + col_j]]
          }

          # Collect partial decryptions from ALL servers (chunked transfer)
          partial_shares <- list()
          for (server in server_list) {
            conn_idx <- which(server_names == server)

            # Send ciphertext in chunks
            ct_str <- ct
            n_chars_ct <- nchar(ct_str)
            n_ct_chunks <- ceiling(n_chars_ct / chunk_size)

            for (ch in seq_len(n_ct_chunks)) {
              start <- (ch - 1) * chunk_size + 1
              end <- min(ch * chunk_size, n_chars_ct)
              chunk_str <- substr(ct_str, start, end)
              store_expr <- call("mheStoreCTChunkDS",
                                 chunk_index = as.integer(ch),
                                 chunk = chunk_str)
              DSI::datashield.aggregate(conns = datasources[conn_idx], expr = store_expr)
            }

            # Compute partial decryption
            call_expr <- call("mhePartialDecryptDS",
                              n_chunks = as.integer(n_ct_chunks))
            pd <- DSI::datashield.aggregate(conns = datasources[conn_idx], expr = call_expr)
            if (is.list(pd)) pd <- pd[[1]]
            partial_shares[[server]] <- pd$decryption_share
          }

          # Client fuses partial decryptions (sum of n_obs slots = inner product)
          value <- .mheFuseLocal(ct, unlist(partial_shares, use.names = FALSE),
                                 log_n = log_n, log_scale = log_scale,
                                 num_slots = n_obs)
          result_matrix[row_i, col_j] <- value
        }
      }

      # Convert to correlation (divide by n-1)
      cross_cor <- result_matrix / (n_obs - 1)
      rownames(cross_cor) <- variables[[server_A]]
      colnames(cross_cor) <- variables[[server_B]]
      cross_cors[[paste(server_A, server_B, sep = "_")]] <- cross_cor

      message("    Threshold decryption complete")
    }
  }

  # =========================================================================
  # Phase 6: Assemble correlation matrix
  # =========================================================================
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

  list(
    correlation = full_cor,
    local_correlations = local_cors,
    cross_correlations = cross_cors,
    n_obs = n_obs,
    servers = server_list,
    variables = variables,
    mhe_params = list(log_n = log_n, log_scale = log_scale)
  )
}
