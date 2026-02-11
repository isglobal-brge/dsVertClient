#' @title Principal Component Analysis for Vertically Partitioned Data
#' @description Client-side function that performs PCA on vertically
#'   partitioned data using distributed Block SVD.
#'
#' @param data_name Character string. Name of the (aligned) data frame on
#'   each server.
#' @param variables A named list where each name corresponds to a server name
#'   and each element is a character vector of variable names from that server.
#' @param n_components Integer. Number of principal components to return.
#'   Default is NULL (returns all).
#' @param datasources DataSHIELD connection object or list of connections.
#'   If NULL, uses all available connections.
#'
#' @return A list with class "ds.pca" containing:
#'   \itemize{
#'     \item \code{scores}: Matrix of principal component scores (n_obs x n_components)
#'     \item \code{loadings}: Matrix of variable loadings (n_vars x n_components)
#'     \item \code{variance}: Variance explained by each component
#'     \item \code{variance_pct}: Percentage of variance explained
#'     \item \code{cumulative_pct}: Cumulative percentage explained
#'     \item \code{var_names}: Variable names
#'   }
#'
#' @details
#' This function extends the Block SVD approach for correlation to perform
#' full PCA:
#'
#' \enumerate{
#'   \item Each server computes U*D from SVD of its standardized variables
#'   \item Client combines and performs final SVD
#'   \item Principal components = U * D (combined)
#'   \item Loadings = V (right singular vectors)
#'   \item Variance = D^2 / (n-1)
#' }
#'
#' @section WARNING - Known Limitation:
#' \strong{THIS FUNCTION REQUIRES A CORRECT IMPLEMENTATION.}
#'
#' The current Block SVD algorithm has a fundamental limitation: the loadings
#' (V matrix from the final SVD) do NOT correctly map to the original variables
#' across servers. This is because the algorithm loses the per-server V matrices
#' that would be needed to transform back to variable space.
#'
#' \strong{What works:}
#' \itemize{
#'   \item Variance explained (eigenvalues) - approximately correct
#'   \item PC scores - correct for the transformed space
#' }
#'
#' \strong{What does NOT work:}
#' \itemize{
#'   \item Loadings - do NOT correctly map to original variable names
#'   \item Interpretation of which variables contribute to each PC
#' }
#'
#' \strong{Status:} Awaiting a privacy-preserving method that correctly computes
#' loadings without requiring disclosure of the V matrices from each server.
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
#' # Plot first two components
#' plot(pca_result$scores[, 1], pca_result$scores[, 2],
#'      xlab = paste0("PC1 (", round(pca_result$variance_pct[1], 1), "%)"),
#'      ylab = paste0("PC2 (", round(pca_result$variance_pct[2], 1), "%)"))
#' }
#'
#' @importFrom DSI datashield.aggregate datashield.connections_find
#' @export
ds.vertPCA <- function(data_name, variables, n_components = NULL,
                       datasources = NULL) {
  # WARNING: Known limitation
  warning(
    "ds.vertPCA: THIS FUNCTION REQUIRES A CORRECT IMPLEMENTATION.\n",
    "The current algorithm does not correctly map loadings to original variables.\n",
    "Variance explained is approximately correct, but loadings do NOT correspond\n",
    "to the original variable names when data spans multiple servers.\n",
    "See ?ds.vertPCA for details.",
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

  # Collect U*D from each server
  UD_list <- list()
  all_var_names <- character()
  n_obs <- NULL

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

    # Verify consistent observation counts
    if (is.null(n_obs)) {
      n_obs <- result$n_obs
    } else if (n_obs != result$n_obs) {
      stop("Inconsistent observation counts across servers. ",
           "Ensure data is properly aligned with ds.alignRecords()",
           call. = FALSE)
    }

    UD_list[[server]] <- result$UD
    all_var_names <- c(all_var_names, result$var_names)
  }

  # Combine U*D matrices column-wise
  UD_combined <- do.call(cbind, UD_list)

  # Compute final SVD
  svd_result <- svd(UD_combined)

  # Principal component scores
  if (length(svd_result$d) == 1) {
    scores <- matrix(svd_result$u * svd_result$d, ncol = 1)
  } else {
    scores <- svd_result$u %*% diag(svd_result$d)
  }

  # Loadings (V matrix)
  loadings <- svd_result$v

  # Variance explained
  singular_values <- svd_result$d
  total_var <- sum(singular_values^2)
  variance <- singular_values^2
  variance_pct <- 100 * variance / total_var
  cumulative_pct <- cumsum(variance_pct)

  # Limit components if requested
  n_total <- length(singular_values)
  if (is.null(n_components)) {
    n_components <- n_total
  } else {
    n_components <- min(n_components, n_total)
  }

  # Set names
  pc_names <- paste0("PC", seq_len(n_components))
  colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
  rownames(loadings) <- all_var_names
  colnames(loadings) <- paste0("PC", seq_len(ncol(loadings)))
  names(variance_pct) <- paste0("PC", seq_len(length(variance_pct)))
  names(cumulative_pct) <- paste0("PC", seq_len(length(cumulative_pct)))

  # Subset to requested components
  result <- list(
    scores = scores[, seq_len(n_components), drop = FALSE],
    loadings = loadings[, seq_len(n_components), drop = FALSE],
    variance = variance[seq_len(n_components)],
    variance_pct = variance_pct[seq_len(n_components)],
    cumulative_pct = cumulative_pct[seq_len(n_components)],
    var_names = all_var_names,
    n_obs = n_obs
  )

  class(result) <- c("ds.pca", "list")
  return(result)
}

#' @title Print Method for ds.pca Objects
#' @description Prints a summary of PCA results.
#' @param x A ds.pca object
#' @param ... Additional arguments (ignored)
#' @export
print.ds.pca <- function(x, ...) {
  cat("Principal Component Analysis (Vertically Partitioned)\n")
  cat("======================================================\n\n")
  cat("Observations:", x$n_obs, "\n")
  cat("Variables:", length(x$var_names), "\n")
  cat("Components:", length(x$variance), "\n\n")

  cat("Variance Explained:\n")
  df <- data.frame(
    Component = names(x$variance_pct),
    Variance = round(x$variance, 4),
    Percent = round(x$variance_pct, 2),
    Cumulative = round(x$cumulative_pct, 2)
  )
  print(df, row.names = FALSE)

  invisible(x)
}
