#' @title Generalized Linear Model for Vertically Partitioned Data
#' @description Client-side function that fits a Generalized Linear Model
#'   across vertically partitioned data using Block Coordinate Descent.
#'
#' @param data_name Character string. Name of the (aligned) data frame on
#'   each server.
#' @param y_var Character string. Name of the response variable (must exist
#'   on ALL servers).
#' @param x_vars A named list where each name corresponds to a server name
#'   and each element is a character vector of predictor variable names
#'   from that server.
#' @param family Character string. GLM family: "gaussian", "binomial",
#'   "poisson", "Gamma", or "inverse.gaussian". Default is "gaussian".
#' @param max_iter Integer. Maximum number of iterations. Default is 100.
#' @param tol Numeric. Convergence tolerance. Default is 1e-6.
#' @param lambda Numeric. L2 regularization parameter. Default is 1e-4.
#' @param verbose Logical. Print iteration progress. Default is TRUE.
#' @param datasources DataSHIELD connection object or list of connections.
#'   If NULL, uses all available connections.
#'
#' @return A list with class "ds.glm" containing:
#'   \itemize{
#'     \item \code{coefficients}: Named vector of coefficient estimates
#'     \item \code{iterations}: Number of iterations until convergence
#'     \item \code{converged}: Logical indicating convergence
#'     \item \code{family}: Family used
#'     \item \code{n_obs}: Number of observations
#'     \item \code{n_vars}: Number of predictor variables
#'     \item \code{deviance}: Residual deviance of the fitted model
#'     \item \code{null_deviance}: Null deviance (intercept-only model)
#'     \item \code{pseudo_r2}: McFadden's pseudo R-squared
#'     \item \code{aic}: Akaike Information Criterion
#'     \item \code{call}: The matched call
#'   }
#'
#' @details
#' This function implements the Block Coordinate Descent (BCD) algorithm
#' for privacy-preserving GLM fitting on vertically partitioned data.
#'
#' The algorithm iteratively:
#' \enumerate{
#'   \item For each partition i:
#'     \itemize{
#'       \item Compute eta_other = sum of X_j * beta_j for all j != i
#'       \item Update beta_i using IRLS with eta_other
#'       \item Share new eta_i = X_i * beta_i (NOT raw data or coefficients)
#'     }
#'   \item Check convergence (sum of |beta_new - beta_old| < tol)
#'   \item Repeat until converged or max_iter reached
#' }
#'
#' Privacy is preserved because:
#' \itemize{
#'   \item Only linear predictor contributions (eta) are shared
#'   \item Raw data never leaves servers
#'   \item Final coefficients are computed locally on each server
#' }
#'
#' @references
#' van Kesteren, E.J. et al. (2019). Privacy-preserving generalized linear
#' models using distributed block coordinate descent. arXiv:1911.03183.
#'
#' @seealso \code{\link[stats]{glm}} for standard GLM fitting
#'
#' @examples
#' \dontrun{
#' # Define predictor variables per server
#' x_vars <- list(
#'   server1 = c("age", "weight"),
#'   server2 = c("height", "bmi"),
#'   server3 = c("glucose", "cholesterol")
#' )
#'
#' # Fit Gaussian GLM (linear regression)
#' model <- ds.vertGLM("D_aligned", "outcome", x_vars, family = "gaussian")
#' print(model)
#'
#' # Fit logistic regression
#' model_logit <- ds.vertGLM("D_aligned", "binary_outcome", x_vars,
#'                           family = "binomial")
#' }
#'
#' @importFrom DSI datashield.aggregate datashield.connections_find
#' @export
ds.vertGLM <- function(data_name, y_var, x_vars, family = "gaussian",
                       max_iter = 100, tol = 1e-6, lambda = 1e-4,
                       verbose = TRUE, datasources = NULL) {
  # Capture call
  call_matched <- match.call()

  # Validate inputs
  if (!is.character(data_name) || length(data_name) != 1) {
    stop("data_name must be a single character string", call. = FALSE)
  }
  if (!is.character(y_var) || length(y_var) != 1) {
    stop("y_var must be a single character string", call. = FALSE)
  }
  if (!is.list(x_vars)) {
    stop("x_vars must be a named list mapping server names to variable vectors",
         call. = FALSE)
  }
  if (!family %in% c("gaussian", "binomial", "poisson", "Gamma", "inverse.gaussian")) {
    stop("family must be 'gaussian', 'binomial', 'poisson', 'Gamma', or 'inverse.gaussian'",
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
  if (!all(names(x_vars) %in% server_names)) {
    missing <- setdiff(names(x_vars), server_names)
    stop("Unknown server(s) in x_vars: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  n_partitions <- length(x_vars)

  # Get observation count from first server
  first_server <- names(x_vars)[1]
  conn_idx <- which(server_names == first_server)
  count_call <- call("getObsCountDS", data_name)
  count_result <- DSI::datashield.aggregate(
    conns = datasources[conn_idx],
    expr = count_call
  )
  if (is.list(count_result) && length(count_result) == 1) {
    count_result <- count_result[[1]]
  }
  n_obs <- count_result$n_obs

  # Initialize coefficients and eta for each partition
  betas <- list()
  etas <- list()
  n_vars_per_partition <- list()

  for (server in names(x_vars)) {
    p <- length(x_vars[[server]])
    betas[[server]] <- rep(0, p)
    etas[[server]] <- rep(0, n_obs)
    n_vars_per_partition[[server]] <- p
  }

  # Total number of variables
  n_vars_total <- sum(unlist(n_vars_per_partition))

  # BCD iterations
  converged <- FALSE
  iter <- 0

  if (verbose) {
    message(sprintf("Starting BCD for %s GLM with %d observations, %d variables, %d partitions",
                    family, n_obs, n_vars_total, n_partitions))
  }

  for (iter in seq_len(max_iter)) {
    betas_old <- betas
    max_diff <- 0

    # Update each partition
    for (server in names(x_vars)) {
      conn_idx <- which(server_names == server)
      vars <- x_vars[[server]]

      # Compute eta from other partitions
      other_servers <- setdiff(names(x_vars), server)
      if (length(other_servers) > 0) {
        eta_other <- Reduce(`+`, etas[other_servers])
      } else {
        eta_other <- rep(0, n_obs)
      }

      # Build call using call()
      call_expr <- call("glmPartialFitDS",
                        data_name, y_var, vars,
                        eta_other, betas[[server]],
                        family, lambda)

      # Execute on server
      result <- DSI::datashield.aggregate(
        conns = datasources[conn_idx],
        expr = call_expr
      )

      if (is.list(result) && length(result) == 1) {
        result <- result[[1]]
      }

      # Update betas and etas
      betas[[server]] <- result$beta
      etas[[server]] <- result$eta

      # Track convergence
      diff <- sum(abs(betas[[server]] - betas_old[[server]]))
      max_diff <- max(max_diff, diff)
    }

    # Check convergence
    if (max_diff < tol) {
      converged <- TRUE
      if (verbose) {
        message(sprintf("Converged after %d iterations (diff = %.2e)", iter, max_diff))
      }
      break
    }

    if (verbose && iter %% 10 == 0) {
      message(sprintf("Iteration %d: max diff = %.2e", iter, max_diff))
    }
  }

  if (!converged && verbose) {
    warning(sprintf("Did not converge after %d iterations (diff = %.2e)",
                    max_iter, max_diff))
  }

  # Compute total eta for deviance calculation
  eta_total <- Reduce(`+`, etas)

  # Calculate deviance from first server (response is same on all)
  deviance_call <- call("glmDevianceDS", data_name, y_var, eta_total, family)
  deviance_result <- DSI::datashield.aggregate(
    conns = datasources[conn_idx],
    expr = deviance_call
  )
  if (is.list(deviance_result) && length(deviance_result) == 1) {
    deviance_result <- deviance_result[[1]]
  }

  deviance <- deviance_result$deviance
  null_deviance <- deviance_result$null_deviance

  # Calculate pseudo R-squared (McFadden's)
  pseudo_r2 <- 1 - (deviance / null_deviance)

  # Calculate AIC: deviance + 2*k (k = number of parameters)
  aic <- deviance + 2 * n_vars_total

  if (verbose) {
    message(sprintf("Deviance: %.4f, Null deviance: %.4f, Pseudo R2: %.4f",
                    deviance, null_deviance, pseudo_r2))
  }

  # Combine coefficients with names
  all_coefs <- numeric()
  all_names <- character()
  for (server in names(x_vars)) {
    all_coefs <- c(all_coefs, betas[[server]])
    all_names <- c(all_names, x_vars[[server]])
  }
  names(all_coefs) <- all_names

  # Build result
  result <- list(
    coefficients = all_coefs,
    iterations = iter,
    converged = converged,
    family = family,
    n_obs = n_obs,
    n_vars = n_vars_total,
    lambda = lambda,
    deviance = deviance,
    null_deviance = null_deviance,
    pseudo_r2 = pseudo_r2,
    aic = aic,
    call = call_matched
  )

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
  cat("Regularization (lambda):", object$lambda, "\n\n")

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
  coef_df <- data.frame(
    Estimate = object$coefficients
  )
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
