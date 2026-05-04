#' @title Federated joint-softmax multinomial logistic regression
#' @description Compatibility wrapper for
#'   \code{\link{ds.vertMultinomJointNewton}}. The historical
#'   one-vs-rest/covariance-rescale implementation was removed from the
#'   product package because it did not close the softmax MLE accuracy gap.
#'   Calls now always dispatch to the non-disclosive joint Newton route for
#'   K >= 3 classes.
#'
#' @param formula R formula with the categorical outcome on the LHS.
#' @param data Aligned data-frame name.
#' @param levels Optional character vector of outcome levels (first is
#'   reference). If NULL, inferred from the outcome server.
#' @param max_iter Outer Newton iterations.
#' @param tol Convergence tolerance on max |delta beta|.
#' @param verbose Print progress.
#' @param datasources DataSHIELD connections.
#' @return A \code{ds.vertMultinomJointNewton} object.
#' @export
ds.vertMultinomJoint <- function(formula, data = NULL, levels = NULL,
                                 max_iter = 30L, tol = 1e-4,
                                 verbose = TRUE, datasources = NULL) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  y_var <- .ds_gee_extract_lhs(formula)
  server_names <- names(datasources)
  y_srv <- .ds_gee_find_server_holding(datasources, server_names, data, y_var)
  if (is.null(y_srv)) stop("y '", y_var, "' not found", call. = FALSE)

  if (is.null(levels)) {
    lv <- tryCatch({
      r <- DSI::datashield.aggregate(
        datasources[which(server_names == y_srv)],
        call(name = "dsvertOutcomeLevelsDS", data_name = data, y_var = y_var))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      r$levels
    }, error = function(e) NULL)
    if (is.null(lv)) {
      stop("dsvertOutcomeLevelsDS not available; pass levels= explicitly.",
           call. = FALSE)
    }
    levels <- as.character(lv)
  }
  if (length(levels) < 3L) {
    stop("Multinomial joint route requires >=3 outcome levels; use ",
         "ds.vertGLM for binary logistic models.",
         call. = FALSE)
  }
  if (verbose) {
    message("[ds.vertMultinomJoint] dispatching to ",
            "ds.vertMultinomJointNewton")
  }
  ds.vertMultinomJointNewton(
    formula = formula,
    data = data,
    levels = levels,
    max_outer = max_iter,
    tol = tol,
    verbose = verbose,
    datasources = datasources)
}

#' @export
print.ds.vertMultinomJoint <- function(x, ...) {
  cat(sprintf("dsVert joint-softmax multinomial (%d classes, ref='%s')\n",
              length(x$levels), x$reference))
  cat(sprintf("  N = %d\n", x$n_obs))
  cat("Coefficient matrix (rows = predictors, cols = non-reference classes):\n")
  print(round(x$coefficients, 4L))
  invisible(x)
}
