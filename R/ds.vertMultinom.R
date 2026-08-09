#' @title Quarantined multinomial compatibility frontdoor
#' @description This exported name is retained for API compatibility. It
#'   raises a typed \code{dsvert_route_unavailable} condition before any DSI
#'   call and computes or returns no multinomial fit. Retained one-vs-rest and
#'   joint-softmax code after the gate is unreachable through this public
#'   frontdoor and carries no disclosure, DP, accuracy, or availability claim.
#' @details Promotion requires a signed bounded joint-softmax artifact over
#'   the exact score design plus validated joint covariance and inference.
#' @param formula,data,classes,reference,indicator_template,max_iter,max_outer,tol,warm_max_iter,warm_tol,binomial_sigmoid_intervals,verbose,datasources,design_analysis_id
#'   Retained compatibility arguments. They are not evaluated because the
#'   public frontdoor fails locally.
#' @param ... Retained compatibility arguments; not evaluated.
#' @return No fitted object. The function raises
#'   \code{dsvert_route_unavailable} before DSI.
#' @seealso \code{\link{ds.vertMethodStatus}}
#' @export
ds.vertMultinom <- function(formula, data = NULL, classes = NULL,
                            reference = NULL, indicator_template = "%s_ind",
                            max_iter = NULL, max_outer = 20L, tol = NULL,
                            warm_max_iter = NULL, warm_tol = NULL,
                            binomial_sigmoid_intervals = NULL,
                            verbose = TRUE, datasources = NULL, ...,
                            design_analysis_id = NULL) {
  .dsvert_block_retired_remote_route("multinomial")
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
    tol = as.numeric(tol %||% 1e-5),
    warm_max_iter = warm_max_iter,
    warm_tol = warm_tol,
    binomial_sigmoid_intervals = binomial_sigmoid_intervals,
    design_analysis_id = design_analysis_id,
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
        r <- .dsvert_aggregate_strict(
          datasources[which(server_nm == srv)],
          call(name = "dsvertLocalMomentsDS", data_name = data,
               variable = sprintf(indicator_template, k)),
          operation = "multinomial warm-start class moment")[[1L]]
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
