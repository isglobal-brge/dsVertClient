#' @title Privacy-Preserving Correlation for Vertically Partitioned Data
#' @description Computes the Pearson correlation matrix for vertically partitioned
#'   data using Multiparty Homomorphic Encryption (MHE). This method correctly
#'   computes cross-server correlations while guaranteeing that:
#'   \itemize{
#'     \item Individual observations are never exposed
#'     \item The client cannot decrypt data without server cooperation
#'     \item Only aggregate statistics (correlations) are revealed
#'   }
#'
#' @param data_name Character string. Name of the (aligned) data frame on
#'   each server.
#' @param variables A named list where each name corresponds to a server name
#'   and each element is a character vector of variable names from that server.
#' @param log_n Integer. Ring dimension parameter (12, 13, or 14 recommended).
#'   Higher values support more observations but are slower. Default is 13.
#' @param log_scale Integer. Precision parameter (40 recommended).
#'   Higher values give more precision but consume more levels. Default is 40.
#' @param datasources DataSHIELD connection object or list of connections.
#'   If NULL, uses all available connections.
#'
#' @return A list with class "ds.cor" containing:
#'   \itemize{
#'     \item \code{correlation}: The correlation matrix (p x p)
#'     \item \code{var_names}: Variable names in order
#'     \item \code{n_obs}: Number of observations
#'     \item \code{method}: "MHE-CKKS" indicating the method used
#'     \item \code{servers}: Names of servers involved
#'   }
#'
#' @details
#' \subsection{Algorithm Overview}{
#' The algorithm computes correlation using homomorphic encryption:
#'
#' \enumerate{
#'   \item \strong{Key Generation}: Each server generates key shares.
#'   \item \strong{Key Combination}: Public keys are combined into a collective key.
#'   \item \strong{Standardization}: Each server standardizes its data locally
#'     (Z = (X - mean) / sd).
#'   \item \strong{Encryption}: Each server encrypts its standardized columns.
#'   \item \strong{Homomorphic Computation}: Cross-products are computed on
#'     encrypted data.
#'   \item \strong{Assembly}: Client assembles the full correlation matrix from
#'     local correlations (diagonal blocks) and cross-correlations (off-diagonal).
#' }
#' }
#'
#' \subsection{Security Guarantees}{
#' \describe{
#'   \item{Confidentiality}{Individual observations are protected by encryption.
#'     Server A sees only Enc(Z_B), not Z_B.}
#'   \item{Non-reconstruction}{The client cannot reconstruct individual data
#'     from the correlation matrix alone.}
#'   \item{Collusion resistance}{Even client + one server cannot decrypt.
#'     Decryption requires ALL servers.}
#' }
#' }
#'
#' @section Performance Notes:
#' \itemize{
#'   \item Computation scales as O(p_A * p_B * n * log(n))
#'   \item Recommended: n < 2000 for log_n=12, n < 4000 for log_n=13
#'   \item Precision: approximately 10^-3 to 10^-4 error due to CKKS approximation
#' }
#'
#' @seealso \code{\link{ds.vertPCA}} for PCA analysis
#'
#' @examples
#' \dontrun{
#' # Connect to DataSHIELD servers
#' connections <- DSI::datashield.login(builder$build())
#'
#' # Align records first
#' ref_hashes <- ds.hashId("data", "patient_id", datasource = connections["server1"])
#' ds.alignRecords("data", "patient_id", ref_hashes$hashes,
#'                 newobj = "data_aligned", datasources = connections)
#'
#' # Define variables per server
#' vars <- list(
#'   hospital_A = c("age", "bmi"),
#'   hospital_B = c("glucose", "systolic_bp")
#' )
#'
#' # Compute privacy-preserving correlation
#' result <- ds.vertCor("data_aligned", vars)
#' print(result$correlation)
#' }
#'
#' @references
#' Cheon, J.H. et al. (2017). "Homomorphic Encryption for Arithmetic of
#' Approximate Numbers". ASIACRYPT 2017.
#'
#' @importFrom DSI datashield.aggregate datashield.connections_find
#' @export
ds.vertCor <- function(data_name, variables, log_n = 13, log_scale = 40,
                       datasources = NULL) {

  # =========================================================================
  # Input Validation
  # =========================================================================
  if (!is.character(data_name) || length(data_name) != 1) {
    stop("data_name must be a single character string", call. = FALSE)
  }
  if (!is.list(variables)) {
    stop("variables must be a named list mapping server names to variable vectors",
         call. = FALSE)
  }
  if (length(variables) < 2) {
    stop("At least 2 servers required for cross-server correlation. ",
         "For single-server correlation, use standard cor() via datashield.aggregate.",
         call. = FALSE)
  }

  # Get datasources
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  if (length(datasources) == 0) {
    stop("No DataSHIELD connections found", call. = FALSE)
  }

  server_names <- names(datasources)

  # Validate all specified servers exist
  missing_servers <- setdiff(names(variables), server_names)
  if (length(missing_servers) > 0) {
    stop("Servers not found in connections: ",
         paste(missing_servers, collapse = ", "), call. = FALSE)
  }

  server_list <- names(variables)

  # =========================================================================
  # Step 1: Key Generation on each server
  # =========================================================================
  message("Initializing MHE keys on servers...")

  server_keys <- list()
  for (server in server_list) {
    conn_idx <- which(server_names == server)

    call_expr <- call("mheKeyGenDS",
                      party_id = 0L,
                      num_parties = 1L,
                      log_n = as.integer(log_n),
                      log_scale = as.integer(log_scale))

    result <- DSI::datashield.aggregate(conns = datasources[conn_idx], expr = call_expr)
    if (is.list(result) && length(result) == 1) result <- result[[1]]

    server_keys[[server]] <- result
    message("  ", server, ": keys generated")
  }

  # =========================================================================
  # Step 2: Combine keys on each server
  # =========================================================================
  message("Combining keys...")

  server_cpk <- list()
  for (server in server_list) {
    conn_idx <- which(server_names == server)

    call_expr <- call("mheCombineKeysDS",
                      public_key_shares = server_keys[[server]]$public_key_share,
                      log_n = as.integer(log_n),
                      log_scale = as.integer(log_scale))

    result <- DSI::datashield.aggregate(conns = datasources[conn_idx], expr = call_expr)
    if (is.list(result) && length(result) == 1) result <- result[[1]]

    server_cpk[[server]] <- result
    message("  ", server, ": keys combined")
  }

  # =========================================================================
  # Step 3: Get observation count and validate alignment
  # =========================================================================
  message("Verifying data alignment...")

  n_obs <- NULL
  for (server in server_list) {
    conn_idx <- which(server_names == server)
    call_expr <- call("getObsCountDS", data_name)
    result <- DSI::datashield.aggregate(conns = datasources[conn_idx], expr = call_expr)
    if (is.list(result) && length(result) == 1) result <- result[[1]]

    server_n <- if (is.list(result)) result$n_obs else result
    if (is.null(n_obs)) {
      n_obs <- server_n
    } else if (n_obs != server_n) {
      stop("Inconsistent observation counts across servers (", n_obs, " vs ", server_n, "). ",
           "Ensure data is properly aligned with ds.alignRecords()", call. = FALSE)
    }
  }
  message("  ", n_obs, " aligned observations across ", length(variables), " servers.")

  # =========================================================================
  # Step 4: Encrypt columns on each server
  # =========================================================================
  message("Encrypting data on each server...")

  all_var_names <- character()
  server_encrypted <- list()

  for (server in server_list) {
    conn_idx <- which(server_names == server)
    vars <- variables[[server]]

    call_expr <- call("mheEncryptColumnsDS",
                      data_name = data_name,
                      variables = vars,
                      collective_public_key = server_cpk[[server]]$collective_public_key,
                      log_n = as.integer(log_n),
                      log_scale = as.integer(log_scale))

    result <- DSI::datashield.aggregate(conns = datasources[conn_idx], expr = call_expr)
    if (is.list(result) && length(result) == 1) result <- result[[1]]

    server_encrypted[[server]] <- result
    all_var_names <- c(all_var_names, vars)
    message("  ", server, ": encrypted ", length(vars), " variables")
  }

  # =========================================================================
  # Step 5: Compute local correlations (diagonal blocks)
  # =========================================================================
  message("Computing local correlations...")

  local_cors <- list()
  for (server in server_list) {
    conn_idx <- which(server_names == server)
    vars <- variables[[server]]

    call_expr <- call("localCorDS",
                      data_name = data_name,
                      variables = vars)

    result <- DSI::datashield.aggregate(conns = datasources[conn_idx], expr = call_expr)
    if (is.list(result) && length(result) == 1) result <- result[[1]]

    local_cors[[server]] <- result$correlation
    message("  ", server, ": ", length(vars), "x", length(vars), " local correlation")
  }

  # =========================================================================
  # Step 6: Compute cross-server correlations (off-diagonal blocks)
  # =========================================================================
  message("Computing cross-server correlations via HE...")

  cross_cors <- list()

  for (i in 1:(length(server_list) - 1)) {
    for (j in (i + 1):length(server_list)) {
      server_A <- server_list[i]
      server_B <- server_list[j]

      conn_idx_A <- which(server_names == server_A)

      # Server A computes cross-product with Server B's encrypted data
      call_expr <- call("mheCrossProductDS",
                        plaintext_data_name = data_name,
                        plaintext_variables = variables[[server_A]],
                        encrypted_columns = server_encrypted[[server_B]]$encrypted_columns,
                        secret_key = server_keys[[server_A]]$secret_key_share,
                        n_obs = as.integer(n_obs),
                        log_n = as.integer(log_n),
                        log_scale = as.integer(log_scale))

      result <- DSI::datashield.aggregate(conns = datasources[conn_idx_A], expr = call_expr)
      if (is.list(result) && length(result) == 1) result <- result[[1]]

      cross_cors[[paste(server_A, server_B, sep = "_")]] <- result$cross_correlation

      message("  ", server_A, " x ", server_B, ": ",
              length(variables[[server_A]]), "x", length(variables[[server_B]]),
              " cross-correlation")
    }
  }

  # =========================================================================
  # Step 7: Assemble full correlation matrix
  # =========================================================================
  message("Assembling full correlation matrix...")

  p <- length(all_var_names)
  R <- matrix(0, nrow = p, ncol = p)
  rownames(R) <- all_var_names
  colnames(R) <- all_var_names

  # Fill diagonal blocks
  start_idx <- 1
  for (server in server_list) {
    vars <- variables[[server]]
    end_idx <- start_idx + length(vars) - 1

    R[start_idx:end_idx, start_idx:end_idx] <- local_cors[[server]]

    start_idx <- end_idx + 1
  }

  # Fill off-diagonal blocks
  start_i <- 1
  for (i in 1:(length(server_list) - 1)) {
    server_A <- server_list[i]
    end_i <- start_i + length(variables[[server_A]]) - 1

    start_j <- end_i + 1
    for (j in (i + 1):length(server_list)) {
      server_B <- server_list[j]
      end_j <- start_j + length(variables[[server_B]]) - 1

      key <- paste(server_A, server_B, sep = "_")
      G_AB <- cross_cors[[key]]

      R[start_i:end_i, start_j:end_j] <- G_AB
      R[start_j:end_j, start_i:end_i] <- t(G_AB)

      start_j <- end_j + 1
    }
    start_i <- end_i + 1
  }

  message("  Complete: ", p, " x ", p, " correlation matrix")

  # =========================================================================
  # Return result
  # =========================================================================
  result <- list(
    correlation = R,
    var_names = all_var_names,
    n_obs = n_obs,
    method = "MHE-CKKS",
    servers = server_list
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
