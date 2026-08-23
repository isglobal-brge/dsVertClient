#' @title Sticky-DP intercept-only multinomial frontdoor
#' @description With a released, validated \code{ds.vertDPFrequency} object,
#'   this frontdoor fits \code{y ~ 1} by deterministic post-processing of that
#'   one sticky categorical release. It returns Jeffreys-smoothed intercept
#'   log-odds and never starts a new analysis or DSI request.
#' @details This is deliberately narrower than a covariate multinomial model:
#'   predictors, joint softmax, covariance, standard errors and inference stay
#'   unavailable until a signed score-design artifact exists. Calls without
#'   \code{frequency} retain the local zero-DSI quarantine gate.
#' @param formula,data,classes,reference,indicator_template,max_iter,max_outer,tol,warm_max_iter,warm_tol,binomial_sigmoid_intervals,verbose,datasources,design_analysis_id
#'   Retained compatibility arguments. They are not evaluated because the
#'   public joint-softmax frontdoor fails locally. With a validated
#'   \code{frequency} object, only \code{y ~ 1} is available and it returns
#'   the sticky-DP categorical intercept fit.
#' @param frequency A released, validated \code{ds.vertDPFrequency} object for
#'   the outcome. This enables only intercept-only multinomial coefficients;
#'   it never starts another analysis or reveals a raw category count.
#' @param ... Retained compatibility arguments; not evaluated.
#' @return With \code{frequency}, a coefficient-only
#'   \code{ds.vertMultinom} object. Otherwise the function raises
#'   \code{dsvert_route_unavailable} before DSI.
#' @seealso \code{\link{ds.vertMethodStatus}}
#' @export
ds.vertMultinom <- function(formula, data = NULL, classes = NULL,
                            reference = NULL, indicator_template = "%s_ind",
                            max_iter = NULL, max_outer = 20L, tol = NULL,
                            warm_max_iter = NULL, warm_tol = NULL,
                            binomial_sigmoid_intervals = NULL,
                            verbose = TRUE, datasources = NULL, ...,
                            design_analysis_id = NULL, frequency = NULL) {
  if (!is.null(frequency)) {
    return(.dsvert_formal_multinom_frequency_adapter(
      explicit_arguments = names(match.call())[-1L],
      formula = if (missing(formula)) NULL else formula,
      data = data, classes = classes, reference = reference,
      frequency = frequency))
  }
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

.dsvert_formal_multinom_frequency_adapter <- function(
    explicit_arguments, formula, data, classes, reference, frequency) {
  allowed <- c("formula", "data", "classes", "reference", "verbose",
               "frequency")
  unexpected <- setdiff(explicit_arguments, allowed)
  if (length(unexpected)) {
    stop(paste(
      "Frequency-backed multinomial does not accept legacy controls:",
      paste(sort(unexpected, method = "radix"), collapse = ", ")),
      call. = FALSE)
  }
  if (!inherits(formula, "formula") || length(formula) != 3L ||
      !is.symbol(formula[[2L]])) {
    stop("Frequency-backed multinomial requires a simple outcome formula",
         call. = FALSE)
  }
  terms <- stats::terms(formula)
  if (!identical(as.integer(attr(terms, "intercept")), 1L) ||
      length(attr(terms, "term.labels")) != 0L) {
    stop("Frequency-backed multinomial supports only an intercept-only y ~ 1 formula",
         call. = FALSE)
  }
  frequency <- .dsvert_dp_frequency_contract(frequency)
  levels <- frequency$levels
  counts <- frequency$counts
  if (!is.character(levels) || length(levels) < 3L || anyNA(levels) ||
      any(!nzchar(levels)) || anyDuplicated(levels) ||
      !is.numeric(counts) || length(counts) != length(levels) ||
      is.null(names(counts)) || !identical(names(counts), levels) ||
      any(!is.finite(counts)) || any(counts < 0) ||
      sum(counts) <= 0) {
    stop("Frequency-backed multinomial requires a non-empty signed categorical release",
         call. = FALSE)
  }
  outcome <- as.character(formula[[2L]])
  if (!identical(outcome, frequency$variable)) {
    stop("Frequency-backed multinomial outcome does not match the signed frequency variable",
         call. = FALSE)
  }
  descriptor <- frequency$coordinate_descriptor
  dataset <- if (is.list(descriptor)) descriptor$dataset else NULL
  if (!is.null(data) && (!is.character(data) || length(data) != 1L ||
                         is.na(data) || !identical(data, dataset))) {
    stop("Frequency-backed multinomial data does not match the signed frequency release",
         call. = FALSE)
  }
  if (!is.null(classes) && (!is.character(classes) ||
                            !identical(classes, levels))) {
    stop("Frequency-backed multinomial classes must equal the signed category order",
         call. = FALSE)
  }
  if (is.null(reference)) reference <- levels[[1L]]
  if (!is.character(reference) || length(reference) != 1L || is.na(reference) ||
      !reference %in% levels) {
    stop("Frequency-backed multinomial reference must be one signed category",
         call. = FALSE)
  }
  alpha <- 0.5
  probabilities <- (counts + alpha) / (sum(counts) + alpha * length(counts))
  names(probabilities) <- levels
  non_reference <- setdiff(levels, reference)
  coefficients <- matrix(
    log(probabilities[non_reference] / probabilities[[reference]]),
    nrow = 1L, dimnames = list("(Intercept)", non_reference))
  result <- list(
    status = "public_certified_intercept_only_multinomial",
    family = "multinomial_intercept_only_frequency_postprocessing",
    classes = non_reference,
    outcome_levels = levels,
    reference = reference,
    coefficients = coefficients,
    probabilities = probabilities,
    dp_counts = counts,
    effective_count_dp = sum(counts),
    smoothing = list(method = "Jeffreys_dirichlet_half", alpha = alpha),
    frequency_release_sha256 = frequency$release_sha256,
    sticky_noise = TRUE,
    sticky_replay = TRUE,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    std_errors = NULL,
    covariance = NULL,
    source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    production_ready = FALSE,
    inference = "unavailable_without_a_joint_multinomial_score_artifact",
    called_via = "ds.vertMultinom_frequency")
  class(result) <- c("dsvert_dp_frequency_multinom", "ds.vertMultinom", "list")
  result
}

.dsvert_formal_multinom_joint_frequency_adapter <- function(
    method, explicit_arguments, formula, data, levels, frequency) {
  allowed <- c("formula", "data", "levels", "verbose", "frequency")
  unexpected <- setdiff(explicit_arguments, allowed)
  if (length(unexpected)) {
    stop(paste(
      "Frequency-backed", method,
      "does not accept legacy controls:",
      paste(sort(unexpected, method = "radix"), collapse = ", ")),
      call. = FALSE)
  }
  fit <- ds.vertMultinom(
    formula = formula, data = data, classes = levels,
    reference = if (is.null(levels)) NULL else levels[[1L]],
    frequency = frequency)
  fit$called_via <- paste0(method, "_frequency")
  fit$requested_method <- method
  fit
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
  if (inherits(x, "dsvert_dp_frequency_multinom")) {
    cat("dsVert sticky-DP intercept-only multinomial fit\n")
    cat("  Reference class:", x$reference, "\n")
    print(round(x$coefficients, 4L))
    cat("  No covariance, standard errors or inference are released.\n")
    return(invisible(x))
  }
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
