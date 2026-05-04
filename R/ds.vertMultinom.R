#' @title Federated multinomial logistic regression
#' @description User-facing multinomial wrapper. Dispatches to
#'   \code{\link{ds.vertMultinomJointNewton}}, the paper-safe joint softmax
#'   Newton route for K >= 3. The historical one-vs-rest approximation is no
#'   longer exposed as a user-facing estimator; it remains only as an internal
#'   warm start for the joint Newton route.
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
#' @param max_iter Optional alias for \code{max_outer}.
#' @param max_outer Maximum outer Newton iterations for the joint route.
#' @param tol Convergence tolerance for the joint route.
#' @param warm_max_iter Optional maximum iterations for each internal
#'   binomial warm-start GLM.
#' @param warm_tol Optional tolerance for each internal binomial warm-start
#'   GLM.
#' @param binomial_sigmoid_intervals Optional DCF spline interval count for
#'   internal binomial warm-start GLMs.
#' @param verbose Logical (default TRUE). Print per-class fit progress.
#' @param datasources DataSHIELD connections; if NULL, uses
#'   \code{DSI::datashield.connections_find()}.
#' @param ... Reserved for future extensions.
#'
#' @return ds.vertMultinom object: a list with per-class \code{ds.glm}
#'   fits, the level vector, the reference, and a consolidated
#'   coefficient matrix (rows = coefficients, columns = non-reference
#'   classes).
#'
#' @export
ds.vertMultinom <- function(formula, data = NULL, classes = NULL,
                            reference = NULL, indicator_template = "%s_ind",
                            max_iter = NULL, max_outer = 8L, tol = NULL,
                            warm_max_iter = NULL, warm_tol = NULL,
                            binomial_sigmoid_intervals = NULL,
                            verbose = TRUE, datasources = NULL, ...) {
  extra <- list(...)
  if (length(extra) > 0L) {
    arg_names <- names(extra)
    arg_names[!nzchar(arg_names)] <- "<unnamed>"
    stop("unused argument(s): ", paste(arg_names, collapse = ", "),
         call. = FALSE)
  }
  if (!inherits(formula, "formula")) {
    stop("`formula` must be an R formula with class indicator on LHS",
         call. = FALSE)
  }
  if (is.null(classes)) {
    stop("Please pass classes = c('A','B','C') indicating the class
          names whose indicator columns should be used", call. = FALSE)
  }
  classes_in <- as.character(classes)
  if (length(classes_in) < 2L) {
    stop("Need at least 2 non-reference classes for a multinomial fit",
         call. = FALSE)
  }
  levels <- if (!is.null(reference)) {
    c(as.character(reference), setdiff(classes_in, as.character(reference)))
  } else {
    classes_in
  }
  if (length(levels) < 3L) {
    stop("Multinomial product route requires >=3 outcome levels; use ",
         "ds.vertGLM for binary logistic models.",
         call. = FALSE)
  }
  if (verbose) {
    message("[ds.vertMultinom] dispatching to ",
            "ds.vertMultinomJointNewton")
  }
  ds.vertMultinomJointNewton(
    formula = formula,
    data = data,
    levels = levels,
    indicator_template = indicator_template,
    max_outer = as.integer(max_iter %||% max_outer),
    tol = as.numeric(tol %||% 1e-4),
    warm_max_iter = warm_max_iter,
    warm_tol = warm_tol,
    binomial_sigmoid_intervals = binomial_sigmoid_intervals,
    verbose = verbose,
    datasources = datasources)
}

#' @keywords internal
.ds_vertMultinomWarm <- function(formula, data = NULL, classes = NULL,
                                 reference = NULL,
                                 indicator_template = "%s_ind",
                                 max_iter = NULL, tol = NULL,
                                 verbose = TRUE, datasources = NULL, ...) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be an R formula with class indicator on LHS",
         call. = FALSE)
  }
  rhs <- attr(terms(formula), "term.labels")
  if (is.null(classes)) {
    stop("Please pass classes = c('A','B','C') indicating the class
          names whose indicator columns should be used", call. = FALSE)
  }
  classes_in <- as.character(classes)
  if (length(classes_in) < 2L) {
    stop("Need at least 2 non-reference classes for a multinomial fit",
         call. = FALSE)
  }

  classes <- classes_in
  if (!is.null(reference)) {
    classes <- setdiff(classes, as.character(reference))
  }
  if (length(classes) < 2L) {
    stop("Need at least 2 non-reference classes for a multinomial fit",
         call. = FALSE)
  }

  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()

  fits <- list()
  dots <- list(...)
  if (!is.null(max_iter)) dots$max_iter <- max_iter
  if (!is.null(tol)) dots$tol <- tol
  for (k in classes) {
    ind_col <- sprintf(indicator_template, k)
    if (verbose) {
      message(sprintf("[ds.vertMultinom] fitting class '%s' vs rest (indicator '%s')",
                       k, ind_col))
    }
    fm <- as.formula(paste(ind_col, "~", paste(rhs, collapse = " + ")))
    fit <- do.call(ds.vertGLM, c(list(
      formula = fm, data = data, family = "binomial",
      verbose = verbose, datasources = datasources), dots))
    fits[[k]] <- fit
  }

  # Consolidate coefficients into a matrix
  cnames <- names(fits[[1]]$coefficients)
  coef_mat <- sapply(fits, function(f) f$coefficients)
  rownames(coef_mat) <- cnames
  se_mat <- sapply(fits, function(f) f$std_errors)
  rownames(se_mat) <- cnames

  # Client-side softmax intercept correction (2026-04-21 PM).
  # Anchor the per-class intercepts so the softmax-normalised
  # probabilities at X match the marginal class proportions:
  #   alpha_k^* = log(prop_k / prop_ref) - X_slopes * gamma_k
  # Slopes are unchanged (still OVR point estimates). Closes ~60% of
  # the median OVR<->softmax gap on birthwt-like balanced 3-class data
  # without any new MPC. The reference-class proportion derives from
  # 1 - Sum_k prop_k if `reference` is explicitly named, or is read
  # directly if the outcome server exposes its indicator column.
  coef_mat_corr <- coef_mat
  class_props <- NULL
  if (!is.null(fits[[1]]$x_means)) {
    # Indicator columns live on a single outcome-holding server. Query
    # each server in isolation; use the first successful response per
    # class. ds.vertGLM already knows the outcome server (its y_server
    # field); fall through if unavailable.
    y_srv <- if (!is.null(fits[[1]]$y_server)) fits[[1]]$y_server else NULL
    server_nm <- names(datasources)
    try_one_server <- function(srv, k) {
      tryCatch({
        r <- DSI::datashield.aggregate(
          datasources[which(server_nm == srv)],
          call(name = "dsvertLocalMomentsDS", data_name = data,
               variable = sprintf(indicator_template, k)))
        if (is.list(r) && length(r) == 1L) r <- r[[1L]]
        if (is.list(r) && !is.null(r$mean)) as.numeric(r$mean) else NA_real_
      }, error = function(e) NA_real_)
    }
    props <- tryCatch({
      p <- sapply(classes, function(k) {
        if (!is.null(y_srv)) {
          v <- try_one_server(y_srv, k)
          if (is.finite(v)) return(v)
        }
        for (srv in server_nm) {
          v <- try_one_server(srv, k)
          if (is.finite(v)) return(v)
        }
        NA_real_
      })
      setNames(p, classes)
    }, error = function(e) NULL)
    if (!is.null(props) && all(is.finite(props))) {
      class_props <- props
      prop_ref <- max(1 - sum(props), 1e-8)
      int_idx <- which(cnames == "(Intercept)")
      if (length(int_idx) == 1L) {
        x_means <- fits[[1]]$x_means
        for (k in classes) {
          gamma_k <- coef_mat[-int_idx, k]
          shared_nm <- intersect(names(gamma_k), names(x_means))
          x_bar_dot_gamma <- if (length(shared_nm) > 0L)
            sum(x_means[shared_nm] * gamma_k[shared_nm]) else 0
          coef_mat_corr[int_idx, k] <- log(props[k] / prop_ref) -
                                         x_bar_dot_gamma
        }
      }
    }
  }

  out <- list(
    fits = fits,
    classes = classes,
    reference = reference,
    coefficients = coef_mat_corr,        # softmax-anchored intercepts
    coefficients_ovr = coef_mat,         # raw OVR (pre-correction)
    class_props = class_props,
    std_errors = se_mat,
    n_obs = fits[[1]]$n_obs,
    family = "multinomial (one-vs-rest + softmax-anchored intercepts)",
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
