#' @title Principal Component Analysis for Vertically Partitioned Data
#' @description Performs PCA on vertically partitioned data using the
#'   privacy-preserving correlation matrix computed via Homomorphic Encryption.
#'
#' @param data_name Character string. Name of the (aligned) data frame on
#'   each server. Ignored if \code{cor_result} is provided.
#' @param variables A named list where each name corresponds to a server name
#'   and each element is a character vector of variable names from that server.
#'   Ignored if \code{cor_result} is provided.
#' @param n_components Integer. Number of principal components to return.
#'   Default is NULL (returns all).
#' @param cor_result An existing \code{ds.cor} object from \code{\link{ds.vertCor}}.
#'   If provided, the MHE protocol is NOT re-run; the PCA is computed directly
#'   from this correlation matrix. This avoids running the expensive MHE protocol
#'   twice when you already have the correlation.
#' @param log_n Integer. CKKS ring dimension parameter for MHE. Default is 12.
#'   Ignored if \code{cor_result} is provided.
#' @param log_scale Integer. CKKS precision parameter for MHE. Default is 40.
#'   Ignored if \code{cor_result} is provided.
#' @param datasources DataSHIELD connection object or list of connections.
#'   If NULL, uses all available connections. Ignored if \code{cor_result} is
#'   provided.
#'
#' @return A list with class "ds.pca" containing:
#'   \itemize{
#'     \item \code{loadings}: Matrix of variable loadings (n_vars x n_components)
#'     \item \code{eigenvalues}: Eigenvalues for each component
#'     \item \code{variance_pct}: Percentage of variance explained
#'     \item \code{cumulative_pct}: Cumulative percentage explained
#'     \item \code{var_names}: Variable names
#'     \item \code{n_obs}: Number of observations
#'     \item \code{correlation}: The correlation matrix used for PCA
#'   }
#'
#' @details
#' This function performs PCA using the correlation matrix obtained via
#' Multiparty Homomorphic Encryption (MHE). The approach is:
#'
#' \enumerate{
#'   \item Compute the privacy-preserving correlation matrix using \code{\link{ds.vertCor}}
#'     (or reuse an existing one via the \code{cor_result} parameter)
#'   \item Perform eigen decomposition on the correlation matrix
#'   \item Extract loadings (eigenvectors) and eigenvalues
#' }
#'
#' Since PCA on standardized data is equivalent to eigen decomposition of the
#' correlation matrix, this gives correct loadings and variance explained.
#'
#' \strong{Note on scores:} This function does NOT return principal component
#' scores because computing scores would require access to the raw data.
#' If you need scores, you would need to compute them on each server using
#' the loadings and aggregate the results (which is a separate operation).
#'
#' \subsection{Interpreting Loadings}{
#' Each column of the loadings matrix represents a principal component.
#' The values show how much each variable contributes to that component:
#' \itemize{
#'   \item Values close to 1 or -1 indicate strong contribution
#'   \item Values close to 0 indicate weak contribution
#'   \item Sign indicates direction of relationship
#' }
#' }
#'
#' @section Security:
#' This function inherits all security properties from \code{\link{ds.vertCor}}:
#' \itemize{
#'   \item Individual observations are never exposed
#'   \item The client cannot decrypt without all servers cooperating
#'   \item Only aggregate statistics (correlation matrix) are revealed
#' }
#'
#' @seealso \code{\link{ds.vertCor}} for correlation analysis
#'
#' @examples
#' \dontrun{
#' vars <- list(
#'   server1 = c("age", "weight"),
#'   server2 = c("height", "bmi")
#' )
#'
#' pca_result <- ds.vertPCA("D_aligned", vars, n_components = 3)
#'
#' # Or reuse an existing correlation result (avoids running MHE again):
#' cor_result <- ds.vertCor("D_aligned", vars)
#' pca_result <- ds.vertPCA(cor_result = cor_result, n_components = 3)
#'
#' # View variance explained
#' print(pca_result)
#'
#' # Biplot of loadings for first two PCs
#' plot(pca_result$loadings[, 1], pca_result$loadings[, 2],
#'      xlim = c(-1, 1), ylim = c(-1, 1),
#'      xlab = paste0("PC1 (", round(pca_result$variance_pct[1], 1), "%)"),
#'      ylab = paste0("PC2 (", round(pca_result$variance_pct[2], 1), "%)"))
#' text(pca_result$loadings[, 1], pca_result$loadings[, 2],
#'      labels = pca_result$var_names, pos = 3)
#' abline(h = 0, v = 0, lty = 2)
#' }
#'
#' @importFrom DSI datashield.aggregate datashield.connections_find
#' @export
ds.vertPCA <- function(data_name = NULL, variables = NULL, n_components = NULL,
                       cor_result = NULL,
                       log_n = 12, log_scale = 40, datasources = NULL) {

  # If an existing correlation result is provided, use it directly
  if (!is.null(cor_result)) {
    if (!inherits(cor_result, "ds.cor")) {
      stop("cor_result must be a ds.cor object from ds.vertCor()", call. = FALSE)
    }
    message("Using provided correlation matrix (skipping MHE protocol)...")
    R <- cor_result$correlation
    n_obs <- cor_result$n_obs
    var_names <- cor_result$var_names
  } else {
    # Validate inputs
    if (is.null(data_name) || !is.character(data_name) || length(data_name) != 1) {
      stop("data_name must be a single character string", call. = FALSE)
    }
    if (is.null(variables) || !is.list(variables)) {
      stop("variables must be a named list mapping server names to variable vectors",
           call. = FALSE)
    }

    # Step 1: Compute privacy-preserving correlation matrix
    message("Computing privacy-preserving correlation matrix...")
    cor_result <- ds.vertCor(data_name, variables,
                             log_n = log_n, log_scale = log_scale,
                             datasources = datasources)

    R <- cor_result$correlation
    n_obs <- cor_result$n_obs
    var_names <- cor_result$var_names
  }

  # Step 2: Eigen decomposition
  message("Performing PCA via eigen decomposition...")
  eigen_result <- eigen(R, symmetric = TRUE)

  # Extract components
  eigenvalues <- eigen_result$values
  loadings <- eigen_result$vectors

  # Handle numerical issues: set small negative eigenvalues to 0
  eigenvalues[eigenvalues < 0] <- 0

  # Variance explained
  total_var <- sum(eigenvalues)
  variance_pct <- 100 * eigenvalues / total_var
  cumulative_pct <- cumsum(variance_pct)

  # Determine number of components
  n_total <- length(eigenvalues)
  if (is.null(n_components)) {
    n_components <- n_total
  } else {
    n_components <- min(n_components, n_total)
  }

  # Set names
  pc_names <- paste0("PC", seq_len(n_components))
  colnames(loadings) <- paste0("PC", seq_len(ncol(loadings)))
  rownames(loadings) <- var_names
  names(eigenvalues) <- paste0("PC", seq_len(length(eigenvalues)))
  names(variance_pct) <- paste0("PC", seq_len(length(variance_pct)))
  names(cumulative_pct) <- paste0("PC", seq_len(length(cumulative_pct)))

  # Subset to requested components
  result <- list(
    loadings = loadings[, seq_len(n_components), drop = FALSE],
    eigenvalues = eigenvalues[seq_len(n_components)],
    variance_pct = variance_pct[seq_len(n_components)],
    cumulative_pct = cumulative_pct[seq_len(n_components)],
    var_names = var_names,
    n_obs = n_obs,
    correlation = R
  )

  class(result) <- c("ds.pca", "list")
  message("PCA complete: ", n_components, " components extracted.")
  return(result)
}

#' @title Print Method for ds.pca Objects
#' @description Prints a summary of PCA results.
#' @param x A ds.pca object
#' @param ... Additional arguments (ignored)
#' @export
print.ds.pca <- function(x, ...) {
  cat("Principal Component Analysis (Privacy-Preserving)\n")
  cat("==================================================\n\n")
  cat("Observations:", x$n_obs, "\n")
  cat("Variables:", length(x$var_names), "\n")
  cat("Components:", length(x$eigenvalues), "\n\n")

  cat("Variance Explained:\n")
  df <- data.frame(
    Component = names(x$variance_pct),
    Eigenvalue = round(x$eigenvalues, 4),
    Percent = round(x$variance_pct, 2),
    Cumulative = round(x$cumulative_pct, 2)
  )
  print(df, row.names = FALSE)

  cat("\nLoadings (top variables for first 2 PCs):\n")
  if (ncol(x$loadings) >= 2) {
    top_loadings <- x$loadings[, 1:min(2, ncol(x$loadings)), drop = FALSE]
    print(round(top_loadings, 3))
  } else {
    print(round(x$loadings, 3))
  }

  invisible(x)
}
