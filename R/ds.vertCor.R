#' @title Correlation Matrix for Vertically Partitioned Data
#' @description Client-side function that computes a correlation matrix
#'   across vertically partitioned data using distributed Block SVD.
#'
#' @param data_name Character string. Name of the (aligned) data frame on
#'   each server.
#' @param variables A named list where each name corresponds to a server name
#'   and each element is a character vector of variable names from that server.
#' @param datasources DataSHIELD connection object or list of connections.
#'   If NULL, uses all available connections.
#'
#' @return A correlation matrix with dimensions equal to the total number
#'   of variables across all servers. Row and column names are the
#'   variable names.
#'
#' @details
#' This function implements distributed correlation analysis using
#' Block Singular Value Decomposition:
#'
#' \enumerate{
#'   \item Each server computes U*D from SVD of its standardized variables
#'   \item Client combines [U1*D1 | U2*D2 | ...] into a single matrix
#'   \item Client computes final SVD to get V and D
#'   \item Correlation matrix = V * D^2 * V' (normalized)
#' }
#'
#' @section WARNING - Known Limitation:
#' \strong{THIS FUNCTION REQUIRES A CORRECT IMPLEMENTATION.}
#'
#' The current Block SVD algorithm does not correctly compute cross-server
#' correlations. The algorithm stacks U*D matrices and performs a final SVD,
#' but this produces results in a transformed space, NOT in the original
#' variable space. To correctly map correlations back to original variables,
#' the V matrix from each server would be needed, but sharing V along with
#' U*D would allow reconstruction of the original data, breaking privacy.
#'
#' \strong{Issue:} Cross-server correlations (e.g., cor(var_A, var_B) where
#' var_A is on server 1 and var_B is on server 2) are NOT correctly computed.
#'
#' \strong{Status:} Awaiting a privacy-preserving method for cross-covariance
#' computation. Consider using differential privacy or secure multi-party
#' computation approaches.
#'
#' @references
#' Iwen, M. & Ong, B.W. (2016). A distributed and incremental SVD algorithm
#' for agglomerative data analysis on large networks. SIAM Journal on Matrix
#' Analysis and Applications.
#'
#' @seealso \code{\link{ds.vertPCA}} for principal component analysis
#'
#' @examples
#' \dontrun{
#' # Define which variables are on which server
#' vars <- list(
#'   server1 = c("age", "weight"),
#'   server2 = c("height", "bmi"),
#'   server3 = c("glucose", "cholesterol")
#' )
#'
#' # Compute correlation matrix
#' cor_matrix <- ds.vertCor("D_aligned", vars, datasources = conns)
#' print(cor_matrix)
#' }
#'
#' @importFrom DSI datashield.aggregate datashield.connections_find
#' @importFrom stats cov2cor
#' @export
ds.vertCor <- function(data_name, variables, datasources = NULL) {
  # WARNING: Known limitation

  warning(
    "ds.vertCor: THIS FUNCTION REQUIRES A CORRECT IMPLEMENTATION.\n",
    "The current algorithm does not correctly compute cross-server correlations.\n",
    "Correlations WITHIN each server are correct, but correlations BETWEEN\n",
    "servers (variables on different servers) may be incorrect.\n",
    "See ?ds.vertCor for details.",
    call. = FALSE
  )

  # Validate inputs
  if (!is.character(data_name) || length(data_name) != 1) {
    stop("data_name must be a single character string", call. = FALSE)
  }
  if (!is.list(variables)) {
    stop("variables must be a named list mapping server names to variable vectors",
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

  # Validate that we have variables for each server
  if (!all(names(variables) %in% server_names)) {
    missing <- setdiff(names(variables), server_names)
    stop("Unknown server(s) in variables: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  # Collect U*D from each server
  UD_list <- list()
  all_var_names <- character()

  for (server in names(variables)) {
    conn_idx <- which(server_names == server)
    if (length(conn_idx) == 0) {
      stop("Server '", server, "' not found in connections", call. = FALSE)
    }

    vars <- variables[[server]]

    # Build call expression using call()
    call_expr <- call("blockSvdDS", data_name, vars, TRUE)

    # Execute on server
    result <- DSI::datashield.aggregate(
      conns = datasources[conn_idx],
      expr = call_expr
    )

    if (is.list(result) && length(result) == 1) {
      result <- result[[1]]
    }

    UD_list[[server]] <- result$UD
    all_var_names <- c(all_var_names, result$var_names)
  }

  # Combine U*D matrices column-wise
  UD_combined <- do.call(cbind, UD_list)

  # Compute final SVD
  svd_final <- svd(UD_combined)

  # Correlation matrix = V * D^2 * V'
  # Since data is standardized, this gives the correlation matrix
  n_vars <- length(svd_final$d)
  if (n_vars == 1) {
    D_squared <- matrix(svd_final$d^2, 1, 1)
  } else {
    D_squared <- diag(svd_final$d^2)
  }
  cor_matrix <- svd_final$v %*% D_squared %*% t(svd_final$v)

  # Normalize to get proper correlation matrix (diagonal = 1)
  cor_matrix <- stats::cov2cor(cor_matrix)

  # Set names
  dimnames(cor_matrix) <- list(all_var_names, all_var_names)

  return(cor_matrix)
}
