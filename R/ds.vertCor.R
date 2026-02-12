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
                       datasources = NULL) {

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

  # Helpers -------------------------------------------------------------------

  # Convert base64url to standard base64 (with padding)
  .b64url_to_b64 <- function(x) {
    x <- gsub("-", "+", x, fixed = TRUE)
    x <- gsub("_", "/", x, fixed = TRUE)
    pad <- (4 - nchar(x) %% 4) %% 4
    if (pad > 0) x <- paste0(x, paste(rep("=", pad), collapse = ""))
    x
  }

  # Find mhe-tool binary (platform-aware)
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

    # Check dsVertClient package
    bin_path <- system.file("bin", subdir, binary_name, package = "dsVertClient")
    if (bin_path != "" && file.exists(bin_path)) return(bin_path)

    # Check dsVert package
    bin_path <- system.file("bin", subdir, binary_name, package = "dsVert")
    if (bin_path != "" && file.exists(bin_path)) return(bin_path)

    # Check environment variable
    env_path <- Sys.getenv("DSVERT_MHE_TOOL")
    if (env_path != "" && file.exists(env_path)) return(env_path)

    stop("mhe-tool binary not found for platform ", subdir, ". ",
         "Install dsVert or set DSVERT_MHE_TOOL env var.", call. = FALSE)
  }

  # Fuse partial decryption shares on the client (no server involvement)
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
      stop("mhe-fuse error: ", output$error, call. = FALSE)
    }

    output$value
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
  message("\n[Phase 1] Key generation on each server...")

  # Get observation count
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
      stop("Observation count mismatch: ", server, " has ", count,
           ", expected ", n_obs, ". Ensure data is aligned.", call. = FALSE)
    }
  }
  message("  Observations: ", n_obs)

  # Party 0 generates CRP and GKG seed
  first_server <- server_list[1]
  conn_idx <- which(server_names == first_server)

  call_expr <- call("mheInitDS",
                    party_id = 0L,
                    crp = NULL,
                    gkg_seed = NULL,
                    num_obs = as.integer(n_obs),
                    log_n = as.integer(log_n),
                    log_scale = as.integer(log_scale))

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

  for (i in 2:length(server_list)) {
    server <- server_list[i]
    conn_idx <- which(server_names == server)

    call_expr <- call("mheInitDS",
                      party_id = as.integer(i - 1),
                      crp = crp,
                      gkg_seed = gkg_seed,
                      num_obs = as.integer(n_obs),
                      log_n = as.integer(log_n),
                      log_scale = as.integer(log_scale))

    result <- DSI::datashield.aggregate(conns = datasources[conn_idx], expr = call_expr)
    if (is.list(result)) result <- result[[1]]
    pk_shares[[server]] <- result$public_key_share
    gkg_shares[[server]] <- result$galois_key_shares
    message("  ", server, ": Party ", i - 1, " initialized")
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

  # Organize GKG shares as [party][galEl] for the combine step
  gkg_shares_ordered <- lapply(server_list, function(s) gkg_shares[[s]])

  call_expr <- call("mheCombineDS",
                    public_key_shares = unlist(pk_shares, use.names = FALSE),
                    crp = crp,
                    galois_key_shares = gkg_shares_ordered,
                    gkg_seed = gkg_seed,
                    num_obs = as.integer(n_obs),
                    log_n = as.integer(log_n),
                    log_scale = as.integer(log_scale))

  combined <- DSI::datashield.aggregate(conns = datasources[conn_idx], expr = call_expr)
  if (is.list(combined)) combined <- combined[[1]]
  cpk <- combined$collective_public_key
  galois_keys <- combined$galois_keys
  message("  CPK created on ", first_server)

  # Distribute CPK and Galois keys to other servers
  for (server in server_list[-1]) {
    conn_idx <- which(server_names == server)
    call_expr <- call("mheStoreCPKDS",
                      cpk = cpk,
                      galois_keys = galois_keys)
    DSI::datashield.aggregate(conns = datasources[conn_idx], expr = call_expr)
    message("  CPK stored on ", server)
  }

  # ===========================================================================
  # Phase 3: Encrypt data on each server
  # ===========================================================================
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

  # ===========================================================================
  # Phase 4: Compute local correlations (no MHE needed)
  # ===========================================================================
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

  # ===========================================================================
  # Phase 5: Cross-server correlations with THRESHOLD DECRYPTION
  # ===========================================================================
  # This is the core privacy-preserving computation. For each server pair
  # (A, B), server A computes z_A * Enc(z_B) homomorphically -- the result
  # is still encrypted. To recover the inner product (correlation), we need
  # threshold decryption: each of the K servers computes a partial decryption
  # share using its secret key, and the client fuses all K shares to recover
  # the plaintext. No single server or the client alone can decrypt.
  #
  # The loop below processes each cross-correlation element individually:
  # for each (i,j) pair, it collects K partial decryption shares and fuses
  # them. This is the dominant cost of the protocol.
  message("\n[Phase 5] Computing cross-server correlations (threshold decryption)...")

  cross_cors <- list()
  # CKKS ciphertexts are 50-200KB as base64url strings. DataSHIELD passes
  # function arguments through R's parser, which has limits on string length.
  # Chunking at 10KB per chunk stays well within parser limits on all Opal versions.
  chunk_size <- 10000

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

      # Compute encrypted cross-product on server_A
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
            ct <- er[[(row_i - 1) * n_B + col_j]]
          }

          # Threshold decryption: collect partial shares from ALL K servers.
          # Each server computes share_k = PartialDecrypt(ct, sk_k).
          # Individually, each share is indistinguishable from random noise.
          # Only when ALL K shares are fused can the plaintext be recovered.
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

          # Client-side fusion: combine all K partial shares to recover the
          # plaintext inner product. This calls the mhe-tool binary's mhe-fuse
          # command, which performs the CKKS-specific threshold decryption algebra.
          value <- .mheFuseLocal(ct, unlist(partial_shares, use.names = FALSE),
                                 log_n = log_n, log_scale = log_scale,
                                 num_slots = n_obs)
          result_matrix[row_i, col_j] <- value
        }
      }

      # The homomorphic computation yields sum(z_A * z_B) -- the raw inner
      # product of standardized columns. Dividing by (n-1) converts to
      # Pearson correlation (same as cor() on standardized data).
      cross_cor <- result_matrix / (n_obs - 1)
      rownames(cross_cor) <- variables[[server_A]]
      colnames(cross_cor) <- variables[[server_B]]
      cross_cors[[paste(server_A, server_B, sep = "_")]] <- cross_cor

      message("    Threshold decryption complete")
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
