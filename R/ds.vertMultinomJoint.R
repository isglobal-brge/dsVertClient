#' @title Quarantined joint-multinomial compatibility frontdoor
#' @description This exported name is retained for API compatibility. It
#'   raises a typed \code{dsvert_route_unavailable} condition before any DSI
#'   call and returns no joint-softmax fit. Retained code after the gate is
#'   unreachable through this public frontdoor and carries no disclosure, DP,
#'   accuracy, or availability claim.
#' @param formula,data,levels,max_iter,tol,verbose,datasources,design_analysis_id
#'   Retained compatibility arguments. They are not evaluated because the
#'   public frontdoor fails locally.
#' @return No fitted object. The function raises
#'   \code{dsvert_route_unavailable} before DSI.
#' @seealso \code{\link{ds.vertMethodStatus}}
#' @export
ds.vertMultinomJoint <- function(formula, data = NULL, levels = NULL,
                                 max_iter = 30L, tol = 1e-4,
                                 verbose = TRUE, datasources = NULL,
                                 design_analysis_id = NULL) {
  .dsvert_block_retired_remote_route("multinomial")
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  y_var <- .ds_gee_extract_lhs(formula)
  server_names <- names(datasources)
  y_srv <- .ds_gee_find_server_holding(datasources, server_names, data, y_var)
  if (is.null(y_srv)) stop("y '", y_var, "' not found", call. = FALSE)

  if (is.null(levels)) {
    lv <- tryCatch({
      r <- .dsvert_aggregate_strict(
        datasources[which(server_names == y_srv)],
        call(name = "dsvertOutcomeLevelsDS", data_name = data, y_var = y_var),
        operation = "multinomial outcome-level discovery")[[1L]]
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
    design_analysis_id = design_analysis_id,
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
