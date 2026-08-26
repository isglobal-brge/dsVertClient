#' @title Sticky-DP multinomial compatibility frontdoor
#' @description The historical joint-multinomial name delegates either to
#'   intercept-only \code{ds.vertDPFrequency} post-processing or to the
#'   signed finite softmax grid accepted by \code{ds.vertMultinom()}.
#' @details This does not re-enable the joint-softmax Newton estimator.
#'   With \code{analysis_id}, an additive formula selects only from the
#'   custodian-signed finite grid. Covariance, standard errors and sampling
#'   inference remain unavailable.
#' @param formula,data,levels,max_iter,tol,verbose,datasources,design_analysis_id
#'   With \code{frequency}, only \code{formula}, \code{data}, \code{levels}
#'   and \code{verbose} are accepted. The optional levels must be exactly the
#'   signed category order; its first category is the reference.
#' @param frequency A released, validated \code{ds.vertDPFrequency} object for
#'   the outcome. It enables only the intercept-only post-processing route.
#' @param analysis_id Custodian-configured signed multinomial likelihood-grid
#'   id for an additive covariate formula.
#' @return A coefficient-only \code{ds.vertMultinom} result without sampling
#'   inference.
#' @seealso \code{\link{ds.vertMethodStatus}}
#' @export
ds.vertMultinomJoint <- function(formula, data = NULL, levels = NULL,
                                 max_iter = 30L, tol = 1e-4,
                                 verbose = TRUE, datasources = NULL,
                                 design_analysis_id = NULL,
                                 frequency = NULL, analysis_id = NULL) {
  if (!is.null(analysis_id)) {
    if (!is.null(frequency) || !is.null(levels) ||
        !identical(max_iter, 30L) || !identical(tol, 1e-4) ||
        !is.null(design_analysis_id)) {
      stop("The signed multinomial grid does not accept legacy joint controls",
           call. = FALSE)
    }
    return(ds.vertMultinom(
      formula, data = data, datasources = datasources, analysis_id = analysis_id))
  }
  if (!is.null(frequency)) {
    return(.dsvert_formal_multinom_joint_frequency_adapter(
      method = "ds.vertMultinomJoint",
      explicit_arguments = names(match.call())[-1L],
      formula = if (missing(formula)) NULL else formula,
      data = data, levels = levels, frequency = frequency))
  }
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
