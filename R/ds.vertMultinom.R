#' @title Federated multinomial logistic regression via K-1 one-vs-rest fits
#' @description Fit a K-category multinomial logistic model by training
#'   K-1 binary logistic regressions in parallel, each contrasting one
#'   class against the pooled rest. This is NOT the full softmax-based
#'   multinomial (which jointly optimises K-1 coefficient vectors under
#'   a common denominator); it is the widely-used one-vs-rest
#'   approximation that delivers very similar coefficients in practice
#'   and avoids a new server-side softmax MPC protocol.
#'
#' @param formula R formula with the class indicator on the LHS. The
#'   class column must be a factor with K levels OR a pre-existing set
#'   of binary indicator columns named \code{paste0(class_col, "_is_",
#'   level_name)} (one per non-reference level) on a single server.
#' @param data Name of the aligned data frame on all servers.
#' @param classes Optional character vector specifying which levels to
#'   fit (default: all non-reference). The reference level is excluded
#'   and its probability is computed as \eqn{1 - \sum p_k} client-side
#'   for any subsequent prediction.
#' @param reference Optional name of the reference level.
#' @param indicator_template String format with "%s" replaced by each
#'   class name to construct indicator column names on the server.
#'   Default "\%s_ind" (e.g., stage_ind_I, stage_ind_II, ...). The
#'   indicator columns must already exist server-side.
#' @param ... passed through to each underlying \code{ds.vertGLM} call.
#'
#' @return ds.vertMultinom object: a list with per-class \code{ds.glm}
#'   fits, the level vector, the reference, and a consolidated
#'   coefficient matrix (rows = coefficients, columns = non-reference
#'   classes).
#'
#' @details Each binary fit uses the existing Ring63+Beaver+DCF pipeline,
#'   so no new MPC primitives are required. The client sees only the
#'   K-1 aggregate coefficient vectors; patient-level class indicators
#'   stay at the server that hosts the outcome.
#'
#' @export
ds.vertMultinom <- function(formula, data = NULL, classes = NULL,
                            reference = NULL, indicator_template = "%s_ind",
                            verbose = TRUE, datasources = NULL, ...) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be an R formula with class indicator on LHS",
         call. = FALSE)
  }
  class_col <- as.character(attr(terms(formula), "variables")[[2]])
  rhs <- attr(terms(formula), "term.labels")

  # Argument validation FIRST (before any DSI call) so that misuse
  # fails fast and the test suite can exercise these paths without an
  # Opal connection.
  if (is.null(classes)) {
    stop("Please pass classes = c('A','B','C') indicating the class
          names whose indicator columns should be used", call. = FALSE)
  }
  if (!is.null(reference)) {
    classes <- setdiff(classes, reference)
  }
  if (length(classes) < 2L) {
    stop("Need at least 2 non-reference classes for a multinomial fit",
         call. = FALSE)
  }

  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()

  fits <- list()
  for (k in classes) {
    ind_col <- sprintf(indicator_template, k)
    if (verbose) {
      message(sprintf("[ds.vertMultinom] fitting class '%s' vs rest (indicator '%s')",
                       k, ind_col))
    }
    fm <- as.formula(paste(ind_col, "~", paste(rhs, collapse = " + ")))
    fit <- ds.vertGLM(fm, data = data, family = "binomial",
                       verbose = verbose,
                       datasources = datasources, ...)
    fits[[k]] <- fit
  }

  # Consolidate coefficients into a matrix
  cnames <- names(fits[[1]]$coefficients)
  coef_mat <- sapply(fits, function(f) f$coefficients)
  rownames(coef_mat) <- cnames
  se_mat <- sapply(fits, function(f) f$std_errors)
  rownames(se_mat) <- cnames

  out <- list(
    fits = fits,
    classes = classes,
    reference = reference,
    coefficients = coef_mat,
    std_errors = se_mat,
    n_obs = fits[[1]]$n_obs,
    family = "multinomial (one-vs-rest)",
    call = match.call())
  class(out) <- c("ds.vertMultinom", "list")
  out
}

#' @export
print.ds.vertMultinom <- function(x, ...) {
  cat(sprintf("dsVert multinomial logistic regression (%d-class one-vs-rest)\n",
              length(x$classes) + as.integer(!is.null(x$reference))))
  if (!is.null(x$reference)) {
    cat(sprintf("  Reference class: %s\n", x$reference))
  }
  cat(sprintf("  N = %d\n\n", x$n_obs))
  cat("Coefficients (rows = predictors, columns = non-reference classes):\n")
  print(round(x$coefficients, 4L))
  invisible(x)
}
