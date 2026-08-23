#' @title Sticky-DP intercept-only multinomial compatibility frontdoor
#' @description The historical joint-multinomial name supports only \code{y ~ 1}
#'   when given a released, validated \code{ds.vertDPFrequency} object. It
#'   delegates to the same deterministic Jeffreys-smoothed log-odds
#'   post-processing as \code{ds.vertMultinom()} and makes no DSI request.
#' @details This does not re-enable the joint-softmax Newton estimator.
#'   Covariates, iterative controls, covariance, standard errors and inference
#'   remain unavailable until a purpose-bound protected score/design artifact
#'   exists. Calls without \code{frequency} fail locally before DSI.
#' @param formula,data,levels,max_iter,tol,verbose,datasources,design_analysis_id
#'   With \code{frequency}, only \code{formula}, \code{data}, \code{levels}
#'   and \code{verbose} are accepted. The optional levels must be exactly the
#'   signed category order; its first category is the reference.
#' @param frequency A released, validated \code{ds.vertDPFrequency} object for
#'   the outcome. It enables only the intercept-only post-processing route.
#' @return With \code{frequency}, a coefficient-only \code{ds.vertMultinom}
#'   result. Otherwise the function raises \code{dsvert_route_unavailable}
#'   before DSI.
#' @seealso \code{\link{ds.vertMethodStatus}}
#' @export
ds.vertMultinomJoint <- function(formula, data = NULL, levels = NULL,
                                 max_iter = 30L, tol = 1e-4,
                                 verbose = TRUE, datasources = NULL,
                                 design_analysis_id = NULL,
                                 frequency = NULL) {
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
