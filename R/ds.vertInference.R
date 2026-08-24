#' @title Likelihood-ratio test or signed Gaussian RSS-reduction region
#' @description `type = "sampling"` computes the standard LR chi-square
#'   statistic between two nested attested ds.vertGLM fits. `type =
#'   "mechanism"` computes a bounded Gaussian residual-sum-of-squares
#'   reduction for a declared submodel from one signed Gaussian Synopsis. It
#'   is a simultaneous DP-mechanism region, not a likelihood-ratio statistic,
#'   p-value, or population-sampling test.
#'
#' @param reduced A ds.glm object with fewer coefficients for `type =
#'   "sampling"`, or a character vector of retained reduced-model terms for
#'   `type = "mechanism"`.
#' @param full A ds.glm object with more coefficients for `type = "sampling"`,
#'   or an unpenalized signed ds.vertDPGaussian result for `type =
#'   "mechanism"`.
#' @param type Either `"sampling"` (default) or `"mechanism"`.
#' @param level Simultaneous mechanism confidence. The current signed Gaussian
#'   Synopsis supports its certified level only (0.95).
#'
#' @return For `type = "sampling"`, a list with:
#'   \itemize{
#'     \item \code{statistic}: reduced$deviance - full$deviance
#'     \item \code{df}: full$n_vars - reduced$n_vars
#'     \item \code{p_value}: upper-tail chi-square p-value
#'     \item \code{deviance_reduced}, \code{deviance_full}
#'   }
#'   For `type = "mechanism"`, a list with the bounded normalized RSS
#'   reduction and its simultaneous outer interval, but no p-value or
#'   likelihood-ratio reference distribution.
#'
#' @details The sampling route is a pure client-side computation over the scalar
#'   deviances already returned by ds.vertGLM; no additional MPC round
#'   is performed and no observation-level information is exposed. It is
#'   available only for converged, unpenalized binomial or Poisson fits on the
#'   same analysis cohort. Gaussian nested-model inference requires an F test
#'   with a valid dispersion estimate and is therefore not represented as an
#'   LR chi-square test here. The currently promoted public DP GLM releases
#'   deliberately contain neither a sampling covariance nor an attested
#'   deviance, so they fail closed here until a joint inference artifact is
#'   available. The mechanism route instead resolves both nested Gaussian
#'   least-squares reductions from exactly one signed sufficient-statistic vector.
#'   It proves an outer region for that bounded, normalized RSS reduction only;
#'   it deliberately returns no chi-square reference or p-value.
#' @export
ds.vertLR <- function(reduced, full,
                      type = c("sampling", "mechanism"), level = 0.95) {
  type <- match.arg(type)
  if (identical(type, "mechanism")) {
    return(.dsvert_dp_gaussian_mechanism_lr(
      reduced = reduced, full = full, level = level))
  }
  .dsvert_reject_unavailable_public_inference(reduced)
  .dsvert_reject_unavailable_public_inference(full)
  if (!inherits(reduced, "ds.glm")) {
    stop("`reduced` must be a ds.glm object (from ds.vertGLM)", call. = FALSE)
  }
  if (!inherits(full, "ds.glm")) {
    stop("`full` must be a ds.glm object (from ds.vertGLM)", call. = FALSE)
  }
  .dsvert_validate_inference_fit(reduced, require_se = FALSE,
                                 require_covariance = FALSE)
  .dsvert_validate_inference_fit(full, require_se = FALSE,
                                 require_covariance = FALSE)
  if (!identical(reduced$family, full$family)) {
    stop("LR test requires both fits to share the same family", call. = FALSE)
  }
  if (!isTRUE(reduced$n_obs == full$n_obs)) {
    stop("LR test requires both fits on the same cohort (n_obs differs)",
         call. = FALSE)
  }
  if (!full$family %in% c("binomial", "poisson")) {
    stop("LR chi-square is supported only for binomial or Poisson fits",
         call. = FALSE)
  }
  required_deviance_type <- paste0("canonical_", full$family)
  if (!identical(reduced$deviance_type, required_deviance_type) ||
      !identical(full$deviance_type, required_deviance_type)) {
    stop("LR test requires canonical unweighted GLM deviances",
         call. = FALSE)
  }
  cohort_reduced <- reduced$cohort_id %||% reduced$alignment_manifest_hash
  cohort_full <- full$cohort_id %||% full$alignment_manifest_hash
  if (is.null(cohort_reduced) || is.null(cohort_full) ||
      !identical(cohort_reduced, cohort_full)) {
    stop("LR test requires a matching same cohort identifier on both fits",
         call. = FALSE)
  }
  if (!identical(reduced$weights %||% NULL, full$weights %||% NULL) ||
      !is.null(reduced$weights) || !is.null(full$weights)) {
    stop("LR test currently requires unweighted fits", call. = FALSE)
  }
  if (!identical(reduced$offset %||% NULL, full$offset %||% NULL)) {
    stop("LR test requires identical offsets", call. = FALSE)
  }
  same_analysis_contract <- vapply(
    c("data_name", "y_var", "missing_policy"),
    function(field) identical(reduced[[field]] %||% NULL,
                              full[[field]] %||% NULL),
    logical(1L)
  )
  if (!all(same_analysis_contract)) {
    stop("LR test requires identical data, response, and missing-data policy",
         call. = FALSE)
  }
  if (!all(names(reduced$coefficients) %in% names(full$coefficients))) {
    stop("Coefficients of `reduced` must be a subset of `full` (nested)",
         call. = FALSE)
  }

  df <- full$n_vars - reduced$n_vars
  if (df <= 0L) {
    stop("Full model must have strictly more coefficients than reduced",
         call. = FALSE)
  }

  deviances <- as.numeric(c(reduced$deviance, full$deviance))
  if (length(deviances) != 2L || any(!is.finite(deviances))) {
    stop("Both fits must expose finite deviances", call. = FALSE)
  }
  stat <- deviances[[1L]] - deviances[[2L]]
  numerical_tol <- 1e-10 * max(1, abs(deviances))
  if (stat < -numerical_tol) {
    stop("The full model has larger deviance; check convergence and nesting",
         call. = FALSE)
  }
  if (stat < 0) stat <- 0
  p <- pchisq(stat, df = df, lower.tail = FALSE)

  out <- list(
    statistic = stat,
    df = as.integer(df),
    p_value = p,
    deviance_reduced = as.numeric(reduced$deviance),
    deviance_full = as.numeric(full$deviance),
    family = full$family,
    n_obs = full$n_obs)
  class(out) <- c("ds.vertLR", "list")
  out
}

.dsvert_dp_gaussian_mechanism_lr <- function(reduced, full, level) {
  if (!inherits(full, "ds.vertDPGaussian")) {
    stop("type='mechanism' requires full to be ds.vertDPGaussian",
         call. = FALSE)
  }
  if (!is.character(reduced) || !length(reduced) || anyNA(reduced) ||
      any(!nzchar(reduced)) || anyDuplicated(reduced)) {
    stop(paste0("type='mechanism' requires reduced to be a non-empty ",
                "unique character vector of retained terms"), call. = FALSE)
  }
  ridge <- suppressWarnings(as.numeric(full$ridge))
  if (length(ridge) != 1L || !is.finite(ridge) || ridge != 0) {
    stop("type='mechanism' requires an unpenalized Gaussian fit (ridge = 0)",
         call. = FALSE)
  }

  region <- .dsvert_dp_gaussian_mechanism_region(full, parm = NULL,
                                                   level = level)
  verified <- ds.validateDPGaussianCertificate(full)
  artifact <- verified$artifact
  moment <- verified$validated_moment
  terms <- artifact$design_terms
  missing_terms <- setdiff(reduced, terms)
  if (length(missing_terms)) {
    stop("Unknown reduced-model term(s): ",
         paste(missing_terms, collapse = ", "), call. = FALSE)
  }
  if (length(reduced) >= length(terms)) {
    stop(paste0("The mechanism LR comparison requires a strict reduced ",
                "submodel"), call. = FALSE)
  }
  reduced_index <- match(reduced, terms)

  scale <- suppressWarnings(as.numeric(verified$output_lattice_scale))
  capacity <- suppressWarnings(as.numeric(verified$coordinate_capacity))
  coordinates <- verified$coordinates
  radius <- suppressWarnings(as.numeric(
    verified$accuracy_simultaneous_95$radius))
  if (length(scale) != 1L || !is.finite(scale) || scale <= 0 ||
      length(capacity) != 1L || !is.finite(capacity) || capacity <= 0 ||
      !is.numeric(coordinates) || !length(coordinates) ||
      any(!is.finite(coordinates)) || length(radius) != 1L ||
      !is.finite(radius) || radius < 0) {
    stop("The Gaussian certificate has invalid mechanism geometry",
         call. = FALSE)
  }
  count_upper <- min(capacity, max(0, coordinates[[1L]] + radius))
  quantization <- if (identical(
      artifact$version, .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_ARTIFACT_VERSION)) {
    suppressWarnings(as.numeric(
      artifact$quantization_contract$per_sum_max_abs_error_lattice_steps)) /
      scale
  } else {
    0.5 * count_upper / scale
  }
  projection <- suppressWarnings(as.numeric(
    moment$projection$frobenius_distance))
  if (length(quantization) != 1L || !is.finite(quantization) ||
      quantization < 0 || length(projection) != 1L ||
      !is.finite(projection) || projection < 0) {
    stop("The Gaussian certificate has invalid deterministic error bounds",
         call. = FALSE)
  }
  coordinate_error <- radius + quantization

  quadratic <- function(gram, cross, label) {
    dimension <- nrow(gram)
    eigenvalues <- eigen(gram, symmetric = TRUE, only.values = TRUE)$values
    tolerance <- 256 * .Machine$double.eps *
      max(1, max(abs(eigenvalues))) * dimension
    minimum <- min(eigenvalues)
    gram_error <- projection + dimension * coordinate_error
    cross_error <- projection + sqrt(dimension) * coordinate_error
    if (!is.finite(minimum) || minimum <= tolerance + gram_error) {
      .dsvert_stop_non_identifiable(
        paste0("The signed Gaussian ", label,
               " design is not stably invertible under its simultaneous ",
               "DP-mechanism region."),
        reason = "dp_mechanism_lr_not_identifiable")
    }
    inverse_norm <- 1 / minimum
    true_inverse_norm <- 1 / (minimum - gram_error)
    cross_norm <- sqrt(sum(cross^2))
    true_cross_norm <- cross_norm + cross_error
    coefficients <- as.numeric(solve(gram, cross))
    list(
      value = sum(coefficients * cross),
      error = cross_error * true_inverse_norm * true_cross_norm +
        cross_norm * true_inverse_norm * inverse_norm * gram_error *
        true_cross_norm + cross_norm * inverse_norm * cross_error)
  }
  full_quadratic <- quadratic(moment$gram, moment$cross, "full-model")
  reduced_quadratic <- quadratic(
    moment$gram[reduced_index, reduced_index, drop = FALSE],
    moment$cross[reduced_index], "reduced-model")
  statistic <- full_quadratic$value - reduced_quadratic$value
  numerical_tolerance <- 1024 * .Machine$double.eps * max(
    1, abs(full_quadratic$value), abs(reduced_quadratic$value))
  if (!is.finite(statistic) || statistic < -numerical_tolerance) {
    .dsvert_stop_non_identifiable(
      "The signed Gaussian nested-model reduction is not numerically ordered.",
      reason = "dp_mechanism_lr_not_identifiable")
  }
  statistic <- max(0, statistic)
  mechanism_radius <- full_quadratic$error + reduced_quadratic$error
  if (!is.finite(mechanism_radius) || mechanism_radius < 0) {
    .dsvert_stop_non_identifiable(
      "The signed Gaussian nested-model mechanism radius is not finite.",
      reason = "dp_mechanism_lr_not_identifiable")
  }
  lower <- max(0, statistic - mechanism_radius)
  upper <- statistic + mechanism_radius
  if (!is.finite(lower) || !is.finite(upper)) {
    .dsvert_stop_non_identifiable(
      "The signed Gaussian nested-model mechanism region is not finite.",
      reason = "dp_mechanism_lr_not_identifiable")
  }
  out <- list(
    statistic = statistic,
    lower = lower,
    upper = upper,
    mechanism_radius = mechanism_radius,
    df = as.integer(length(terms) - length(reduced)),
    reduced_terms = reduced,
    full_terms = terms,
    distribution = "simultaneous_dp_mechanism_rss_reduction",
    level = attr(region, "level"),
    sampling_inference = FALSE,
    additional_server_calls = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    estimand = paste(
      "bounded normalized Gaussian nested-model residual-sum-of-squares",
      "reduction from one signed sufficient-statistic release"))
  class(out) <- c("ds.vertLR", "list")
  out
}

.dsvert_reject_unavailable_public_inference <- function(fit) {
  if (!inherits(fit, "ds.vertDPGaussian") &&
      !inherits(fit, "dsvert_formal_dp_glm")) {
    return(invisible(NULL))
  }
  stop(structure(
    list(
      message = paste(
        "sampling inference is unavailable for this public DP GLM release;",
        "it has no attested sampling covariance or deviance. A future joint",
        "inference artifact is required."),
      call = NULL,
      reason = "joint_glm_inference_artifact_required"),
    class = c("dsvert_inference_unavailable", "error", "condition")))
}

.dsvert_validate_inference_fit <- function(fit, require_se = TRUE,
                                            require_covariance = FALSE) {
  .dsvert_reject_unavailable_public_inference(fit)
  if (!inherits(fit, "ds.glm")) {
    stop("`fit` must be a ds.glm object", call. = FALSE)
  }
  if (!isTRUE(fit$converged)) {
    stop("Inference requires a converged fit", call. = FALSE)
  }
  lambda <- fit$lambda %||% 0
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) ||
      lambda != 0) {
    stop("Wald/LR inference requires an unpenalized fit (lambda = 0)",
         call. = FALSE)
  }
  coefficients <- fit$coefficients
  if (!is.numeric(coefficients) || !length(coefficients) ||
      is.null(names(coefficients)) || any(!nzchar(names(coefficients))) ||
      anyDuplicated(names(coefficients)) || any(!is.finite(coefficients))) {
    stop("fit must expose uniquely named finite coefficients", call. = FALSE)
  }
  if (isTRUE(require_se)) {
    se <- fit$std_errors
    if (!is.numeric(se) || is.null(names(se)) ||
        !all(names(coefficients) %in% names(se)) ||
        any(!is.finite(se[names(coefficients)])) ||
        any(se[names(coefficients)] <= 0)) {
      stop("fit standard errors must be named, finite, positive and complete",
           call. = FALSE)
    }
  }
  if (isTRUE(require_covariance)) {
    covariance <- fit$covariance
    if (is.null(covariance)) {
      stop("fit does not expose the full covariance matrix", call. = FALSE)
    }
    covariance <- as.matrix(covariance)
    p <- length(coefficients)
    if (!is.numeric(covariance) || !identical(dim(covariance), c(p, p)) ||
        any(!is.finite(covariance)) ||
        max(abs(covariance - t(covariance))) >
          1e-8 * max(1, max(abs(covariance)))) {
      stop("fit covariance must be a finite symmetric coefficient matrix",
           call. = FALSE)
    }
    eigenvalues <- eigen((covariance + t(covariance)) / 2,
                         symmetric = TRUE, only.values = TRUE)$values
    covariance_scale <- max(abs(eigenvalues))
    if (!is.finite(covariance_scale) || covariance_scale <= 0 ||
        min(eigenvalues) <= sqrt(.Machine$double.eps) * covariance_scale) {
      .dsvert_stop_non_identifiable(
        paste0("The fitted coefficient covariance is not positive definite; ",
               "coefficient inference is not identifiable."),
        reason = "singular_coefficient_covariance")
    }
  }
  invisible(TRUE)
}

.dsvert_inference_reference <- function(fit) {
  if (!identical(fit$family, "gaussian")) {
    return(list(
      coefficient = "normal", contrast = "chi-square",
      df_residual = NULL
    ))
  }
  n_obs <- fit$n_obs
  n_parameters <- fit$n_parameters %||% fit$n_vars
  df_residual <- fit$df_residual %||% (n_obs - n_parameters)
  if (!is.numeric(df_residual) || length(df_residual) != 1L ||
      !is.finite(df_residual) || df_residual <= 0 ||
      df_residual != floor(df_residual)) {
    stop("Gaussian inference requires a positive residual degrees of freedom",
         call. = FALSE)
  }
  list(
    coefficient = "t", contrast = "F",
    df_residual = as.integer(df_residual)
  )
}

#' @export
print.ds.vertLR <- function(x, ...) {
  if (identical(x$distribution,
                "simultaneous_dp_mechanism_rss_reduction")) {
    cat("dsVert signed DP-mechanism Gaussian RSS-reduction region\n")
    cat(sprintf("  reduction = %.6f, certified region = [%.6f, %.6f]\n",
                x$statistic, x$lower, x$upper))
    cat("  no likelihood-ratio reference or sampling p-value\n")
    return(invisible(x))
  }
  cat("dsVert likelihood-ratio test\n")
  cat(sprintf("  Family : %s\n", x$family))
  cat(sprintf("  N      : %d\n", x$n_obs))
  cat(sprintf("  Deviance: reduced=%.4f  full=%.4f\n",
              x$deviance_reduced, x$deviance_full))
  cat(sprintf("  LR chi-sq = %.4f on %d df,  p-value = %s\n",
              x$statistic, x$df, format.pval(x$p_value, digits = 4L)))
  invisible(x)
}


#' @title Wald confidence intervals and DP-mechanism regions for coefficients
#' @description `type = "sampling"` returns Wald-type confidence intervals
#'   from a fit carrying an attested sampling covariance. For a current signed
#'   Gaussian Synopsis, `type = "mechanism"` derives a simultaneous outer
#'   region from the certificate's DP coordinate radius, quantisation bound,
#'   projection distance and the released system's inverse-norm margin. It is
#'   deterministic client post-processing: it performs no new server call,
#'   exposes no source statistic and must not be interpreted as sampling
#'   inference. Binomial and Poisson sampling inference remain unavailable
#'   until a joint signed inference artifact exists.
#'
#' @param fit A fitted model object; `type = "mechanism"` requires a signed
#'   `ds.vertDPGaussian` result.
#' @param parm Optional character vector of coefficient names to report;
#'   default all.
#' @param level Confidence level, default 0.95.
#' @param type Either `"sampling"` (the default Wald interval for an
#'   attested fit with sampling covariance) or `"mechanism"`. The latter is
#'   available only for a signed `ds.vertDPGaussian` result and returns a
#'   simultaneous DP-mechanism region, not a population-sampling interval.
#'
#' @return For `type = "sampling"`, a data frame with columns
#'   \code{estimate}, \code{std_error}, \code{lower}, \code{upper}; row
#'   names are coefficient names. For `type = "mechanism"`, a data frame
#'   with \code{estimate}, \code{lower}, \code{upper} and
#'   \code{mechanism_radius}, carrying explicit attributes that rule out a
#'   sampling-inference interpretation.
#' @export
ds.vertConfint <- function(fit, parm = NULL, level = 0.95,
                           type = c("sampling", "mechanism")) {
  type <- match.arg(type)
  if (identical(type, "mechanism")) {
    return(.dsvert_dp_gaussian_mechanism_region(fit, parm, level))
  }
  .dsvert_validate_inference_fit(fit, require_se = TRUE)
  reference <- .dsvert_inference_reference(fit)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be in (0, 1)", call. = FALSE)
  }
  coefs <- fit$coefficients
  ses <- fit$std_errors
  if (is.null(names(coefs)) || is.null(names(ses))) {
    stop("fit must expose named coefficients and std_errors", call. = FALSE)
  }
  if (!is.null(parm)) {
    missing_parm <- setdiff(parm, names(coefs))
    if (length(missing_parm) > 0L) {
      stop("Unknown coefficient(s): ", paste(missing_parm, collapse = ", "),
           call. = FALSE)
    }
    idx <- match(parm, names(coefs))
  } else {
    idx <- seq_along(coefs)
  }
  critical_value <- if (identical(reference$coefficient, "t")) {
    stats::qt((1 + level) / 2, df = reference$df_residual)
  } else {
    stats::qnorm((1 + level) / 2)
  }
  est <- as.numeric(coefs[idx])
  se <- as.numeric(ses[idx])
  out <- data.frame(
    estimate = est,
    std_error = se,
    lower = est - critical_value * se,
    upper = est + critical_value * se,
    stringsAsFactors = FALSE)
  rownames(out) <- names(coefs)[idx]
  attr(out, "level") <- level
  attr(out, "distribution") <- reference$coefficient
  attr(out, "df") <- reference$df_residual
  out
}

.dsvert_dp_gaussian_mechanism_region <- function(fit, parm, level) {
  if (!inherits(fit, "ds.vertDPGaussian")) {
    stop("type='mechanism' is available only for ds.vertDPGaussian",
         call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be in (0, 1)", call. = FALSE)
  }
  verified <- ds.validateDPGaussianCertificate(fit)
  if (!isTRUE(verified$integrity_valid) ||
      !identical(verified$authenticity, "session_transport_anchored")) {
    stop("The Gaussian mechanism region requires a transport-anchored certificate",
         call. = FALSE)
  }
  artifact <- verified$artifact
  moment <- verified$validated_moment
  accuracy <- verified$accuracy_simultaneous_95
  confidence <- suppressWarnings(as.numeric(accuracy$confidence))
  if (length(confidence) != 1L || !is.finite(confidence) ||
      abs(level - confidence) > sqrt(.Machine$double.eps) *
        max(1, abs(level), abs(confidence))) {
    stop(paste0(
      "type='mechanism' requires a current signed Synopsis certificate at ",
      "its simultaneous mechanism confidence"), call. = FALSE)
  }
  radius <- suppressWarnings(as.numeric(accuracy$radius))
  if (length(radius) != 1L || !is.finite(radius) || radius < 0) {
    stop("The Gaussian certificate has an invalid simultaneous mechanism radius",
         call. = FALSE)
  }
  ridge <- suppressWarnings(as.numeric(fit$ridge))
  if (length(ridge) != 1L || !is.finite(ridge) || ridge < 0) {
    stop("The Gaussian result has an invalid ridge value", call. = FALSE)
  }
  recomputed <- .dsvert_dp_gaussian_solve(moment, artifact, ridge)
  original <- .dsvert_dp_gaussian_original_coefficients(
    recomputed$coefficients, artifact)
  if (!isTRUE(all.equal(fit$coefficients_normalized,
                        recomputed$coefficients, tolerance = 0)) ||
      !isTRUE(all.equal(fit$coefficients_original_scale, original,
                        tolerance = 0)) ||
      !isTRUE(all.equal(fit$coefficients, original, tolerance = 0))) {
    stop("The Gaussian result does not match its signed sufficient statistics",
         call. = FALSE)
  }

  scale <- suppressWarnings(as.numeric(verified$output_lattice_scale))
  capacity <- suppressWarnings(as.numeric(verified$coordinate_capacity))
  coordinates <- verified$coordinates
  if (length(scale) != 1L || !is.finite(scale) || scale <= 0 ||
      length(capacity) != 1L || !is.finite(capacity) || capacity <= 0 ||
      !is.numeric(coordinates) || !length(coordinates) ||
      any(!is.finite(coordinates))) {
    stop("The Gaussian certificate has invalid mechanism geometry",
         call. = FALSE)
  }
  count_upper <- min(capacity, max(0, coordinates[[1L]] + radius))
  quantization <- if (identical(
      artifact$version, .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_ARTIFACT_VERSION)) {
    suppressWarnings(as.numeric(
      artifact$quantization_contract$per_sum_max_abs_error_lattice_steps)) /
      scale
  } else {
    0.5 * count_upper / scale
  }
  projection <- suppressWarnings(as.numeric(
    moment$projection$frobenius_distance))
  if (length(quantization) != 1L || !is.finite(quantization) ||
      quantization < 0 || length(projection) != 1L ||
      !is.finite(projection) || projection < 0) {
    stop("The Gaussian certificate has invalid deterministic error bounds",
         call. = FALSE)
  }

  dimension <- length(recomputed$coefficients)
  coordinate_error <- radius + quantization
  gram_error <- projection + dimension * coordinate_error
  cross_error <- projection + sqrt(dimension) * coordinate_error
  penalty <- rep(1, dimension)
  if (isTRUE(artifact$intercept)) penalty[[1L]] <- 0
  system <- moment$gram + diag(ridge * penalty, dimension)
  eigenvalues <- eigen(system, symmetric = TRUE, only.values = TRUE)$values
  minimum <- min(eigenvalues)
  tolerance <- 256 * .Machine$double.eps * max(1, max(abs(eigenvalues))) *
    dimension
  inverse_norm <- if (minimum > tolerance) 1 / minimum else Inf
  denominator <- 1 - inverse_norm * gram_error
  if (!is.finite(inverse_norm) || !is.finite(denominator) ||
      denominator <= sqrt(.Machine$double.eps)) {
    .dsvert_stop_non_identifiable(
      paste(
        "The certified simultaneous DP perturbation is too large to identify",
        "a finite Gaussian mechanism region for this fitted system."),
      reason = "dp_mechanism_region_not_identifiable")
  }
  normalized_radius <- inverse_norm * (
    cross_error + gram_error * sqrt(sum(recomputed$coefficients^2))) /
    denominator

  outcome_range <- artifact$outcome$upper - artifact$outcome$lower
  predictor_order <- artifact$predictor_order
  transform <- matrix(
    0, nrow = length(original), ncol = dimension,
    dimnames = list(names(original), names(recomputed$coefficients)))
  if (isTRUE(artifact$intercept)) {
    transform["(Intercept)", "(Intercept)"] <- outcome_range
  }
  for (predictor in predictor_order) {
    bounds <- artifact$predictors[[predictor]]
    predictor_range <- bounds$upper - bounds$lower
    if (!is.finite(predictor_range) || predictor_range <= 0) {
      stop("The Gaussian certificate has invalid predictor bounds",
           call. = FALSE)
    }
    slope_scale <- outcome_range / predictor_range
    transform[predictor, predictor] <- slope_scale
    transform["(Intercept)", predictor] <-
      transform["(Intercept)", predictor] - slope_scale * bounds$lower
  }
  mechanism_radius <- sqrt(rowSums(transform^2)) * normalized_radius
  if (any(!is.finite(mechanism_radius)) || any(mechanism_radius < 0)) {
    stop("The Gaussian mechanism region is not finite", call. = FALSE)
  }
  if (is.null(parm)) {
    selected <- names(original)
  } else {
    if (!is.character(parm) || anyNA(parm) || any(!nzchar(parm))) {
      stop("`parm` must contain non-empty coefficient names", call. = FALSE)
    }
    selected <- parm
    unknown <- setdiff(selected, names(original))
    if (length(unknown)) {
      stop("Unknown coefficient(s): ", paste(unknown, collapse = ", "),
           call. = FALSE)
    }
  }
  estimate <- original[selected]
  selected_radius <- mechanism_radius[selected]
  out <- data.frame(
    estimate = as.numeric(estimate),
    lower = as.numeric(estimate - selected_radius),
    upper = as.numeric(estimate + selected_radius),
    mechanism_radius = as.numeric(selected_radius),
    stringsAsFactors = FALSE,
    row.names = selected)
  attr(out, "level") <- confidence
  attr(out, "interval_scope") <- "simultaneous_dp_mechanism"
  attr(out, "sampling_inference") <- FALSE
  attr(out, "estimand") <- paste(
    "bounded complete-case Gaussian ridge estimator from the signed",
    "sufficient-statistic artifact")
  attr(out, "normalized_l2_radius") <- normalized_radius
  attr(out, "additional_server_calls") <- 0L
  attr(out, "additional_privacy_cost") <- c(epsilon = 0, delta = 0)
  out
}


#' @title Univariate Wald test or signed DP-mechanism null region
#' @description `type = "sampling"` tests H0: beta_j = null from an attested
#'   sampling standard error. `type = "mechanism"` instead uses the current
#'   signed Gaussian Synopsis to determine whether its simultaneous DP
#'   mechanism region excludes the null. It is deterministic post-processing,
#'   not a Wald statistic, p-value, or population-sampling test.
#'
#' @param fit A ds.glm object for sampling inference, or a signed
#'   ds.vertDPGaussian result for a mechanism region.
#' @param parm Single coefficient name.
#' @param null Null value (default 0).
#' @param type Either `"sampling"` (default) or `"mechanism"`.
#' @param level Simultaneous mechanism confidence. The current signed
#'   Gaussian Synopsis supports its certified level only (0.95).
#'
#' @return For `type = "sampling"`, a list with estimate, SE, statistic,
#'   distribution, p_value and null. For `type = "mechanism"`, a list with
#'   the certified outer interval and whether it excludes the null; it has no
#'   sampling standard error, statistic or p-value.
#' @export
ds.vertWald <- function(fit, parm, null = 0,
                        type = c("sampling", "mechanism"), level = 0.95) {
  type <- match.arg(type)
  if (!is.character(parm) || length(parm) != 1L || is.na(parm) ||
      !nzchar(parm)) {
    stop("`parm` must be a single coefficient name", call. = FALSE)
  }
  if (!is.numeric(null) || length(null) != 1L || !is.finite(null)) {
    stop("`null` must be one finite number", call. = FALSE)
  }
  if (identical(type, "mechanism")) {
    region <- ds.vertConfint(fit, parm = parm, level = level,
                             type = "mechanism")
    lower <- unname(region[parm, "lower"])
    upper <- unname(region[parm, "upper"])
    out <- list(
      parm = parm,
      estimate = unname(region[parm, "estimate"]),
      null = null,
      lower = lower,
      upper = upper,
      mechanism_radius = unname(region[parm, "mechanism_radius"]),
      null_excluded = isTRUE(null < lower || null > upper),
      distribution = "simultaneous_dp_mechanism_region",
      level = attr(region, "level"),
      sampling_inference = FALSE,
      additional_server_calls = 0L,
      additional_privacy_cost = c(epsilon = 0, delta = 0))
    class(out) <- c("ds.vertWald", "list")
    return(out)
  }
  .dsvert_validate_inference_fit(fit, require_se = TRUE)
  reference <- .dsvert_inference_reference(fit)
  if (!parm %in% names(fit$coefficients)) {
    stop("Coefficient '", parm, "' not in fit", call. = FALSE)
  }
  est <- as.numeric(fit$coefficients[parm])
  se <- as.numeric(fit$std_errors[parm])
  statistic <- (est - null) / se
  p <- if (identical(reference$coefficient, "t")) {
    2 * stats::pt(abs(statistic), df = reference$df_residual,
                  lower.tail = FALSE)
  } else {
    2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
  }
  out <- list(
    parm = parm,
    estimate = est,
    std_error = se,
    null = null,
    statistic = statistic,
    distribution = reference$coefficient,
    df = reference$df_residual,
    p_value = p)
  if (identical(reference$coefficient, "t")) {
    out$t <- statistic
  } else {
    out$z <- statistic
  }
  class(out) <- c("ds.vertWald", "list")
  out
}

#' @export
print.ds.vertWald <- function(x, ...) {
  if (identical(x$distribution, "simultaneous_dp_mechanism_region")) {
    cat(sprintf("dsVert signed DP-mechanism null region: %s = %g\n",
                x$parm, x$null))
    cat(sprintf("  certified region = [%.6f, %.6f]\n", x$lower, x$upper))
    cat(sprintf("  null excluded = %s; no sampling p-value\n",
                if (isTRUE(x$null_excluded)) "yes" else "no"))
    return(invisible(x))
  }
  cat(sprintf("dsVert Wald test: H0: %s = %g  vs  two-sided\n",
              x$parm, x$null))
  statistic_label <- if (identical(x$distribution, "t")) "t" else "z"
  cat(sprintf("  estimate = %.6f  SE = %.6f  %s = %.4f\n",
              x$estimate, x$std_error, statistic_label, x$statistic))
  if (identical(x$distribution, "t")) {
    cat(sprintf("  residual df = %d\n", x$df))
  }
  cat(sprintf("  p-value  = %s\n", format.pval(x$p_value, digits = 4L)))
  invisible(x)
}


#' @title Multi-coefficient Wald test or signed DP-mechanism contrast region
#' @description `type = "sampling"` tests H0: K * beta = m with an attested
#'   covariance. `type = "mechanism"` propagates the one signed simultaneous
#'   Gaussian coefficient region through K. It can prove that individual
#'   contrast components exclude their nulls, but it is not a Wald statistic,
#'   p-value, or population-sampling test.
#'
#' @param fit A ds.glm object with a non-NULL `covariance` slot for sampling
#'   inference, or a signed ds.vertDPGaussian result for a mechanism region.
#' @param K   Contrast matrix: numeric matrix with ncol equal to the
#'   number of coefficients. Rows define the contrasts under test.
#'   Alternatively a named-coef character vector (treated as indicator
#'   rows) or a character RHS parsed against the coefficient names.
#' @param m   Null vector (length nrow(K)); default zero.
#' @param type Either `"sampling"` (default) or `"mechanism"`.
#' @param level Simultaneous mechanism confidence. The current signed
#'   Gaussian Synopsis supports its certified level only (0.95).
#'
#' @return For `type = "sampling"`, a list of class ds.vertContrast with
#'   covariance-derived statistics and a p-value. For `type = "mechanism"`,
#'   a list with componentwise certified outer bounds and null-exclusion flags;
#'   it contains no sampling covariance, statistic or p-value.
#' @export
ds.vertContrast <- function(fit, K, m = NULL,
                            type = c("sampling", "mechanism"),
                            level = 0.95) {
  type <- match.arg(type)
  if (identical(type, "mechanism")) {
    return(.dsvert_dp_gaussian_mechanism_contrast(
      fit = fit, K = K, m = m, level = level))
  }
  .dsvert_validate_inference_fit(fit, require_se = TRUE,
                                 require_covariance = TRUE)
  reference <- .dsvert_inference_reference(fit)
  cov <- as.matrix(fit$covariance)
  coef <- as.numeric(fit$coefficients)
  names(coef) <- names(fit$coefficients)

  # Coerce K into a numeric matrix with one row per contrast and one
  # column per coefficient.
  if (is.character(K)) {
    # Character vector of coefficient names -> identity-like contrast.
    miss <- setdiff(K, names(coef))
    if (length(miss)) {
      stop("Unknown coefficient(s) in K: ", paste(miss, collapse = ", "),
           call. = FALSE)
    }
    Kmat <- matrix(0, nrow = length(K), ncol = length(coef),
                   dimnames = list(K, names(coef)))
    for (i in seq_along(K)) Kmat[i, K[i]] <- 1
    K <- Kmat
  } else if (is.vector(K)) {
    K <- matrix(K, nrow = 1L, dimnames = list(NULL, names(coef)))
  }
  K <- as.matrix(K)

  if (!is.numeric(K) || any(!is.finite(K)) || nrow(K) < 1L) {
    stop("K must be a non-empty finite numeric contrast matrix",
         call. = FALSE)
  }

  if (ncol(K) != length(coef)) {
    stop("ncol(K) = ", ncol(K), " must equal number of coefficients (",
         length(coef), ")", call. = FALSE)
  }
  if (is.null(m)) m <- rep(0, nrow(K))
  if (!is.numeric(m) || length(m) != nrow(K) || any(!is.finite(m))) {
    stop("length(m) must equal nrow(K)", call. = FALSE)
  }

  estimate <- drop(K %*% coef) - m
  var_mat <- K %*% cov %*% t(K)
  var_mat <- (var_mat + t(var_mat)) / 2  # enforce symmetry
  inv_var <- .dsvert_solve_identifiable(
    var_mat,
    context = "The requested K * Cov * K^T contrast",
    reason = "singular_contrast_covariance",
    symmetric = TRUE)
  wald_stat <- drop(t(estimate) %*% inv_var %*% estimate)
  df <- nrow(K)
  if (identical(reference$contrast, "F")) {
    stat <- wald_stat / df
    p <- stats::pf(stat, df1 = df, df2 = reference$df_residual,
                   lower.tail = FALSE)
  } else {
    stat <- wald_stat
    p <- stats::pchisq(stat, df = df, lower.tail = FALSE)
  }

  out <- list(
    estimate = estimate,
    variance = var_mat,
    statistic = as.numeric(stat),
    wald_statistic = as.numeric(wald_stat),
    distribution = reference$contrast,
    df = as.integer(df),
    df_residual = reference$df_residual,
    p_value = as.numeric(p),
    K = K,
    null = m)
  class(out) <- c("ds.vertContrast", "list")
  out
}

.dsvert_dp_gaussian_mechanism_contrast <- function(fit, K, m, level) {
  region <- ds.vertConfint(fit, level = level, type = "mechanism")
  coefficients <- rownames(region)
  if (is.character(K)) {
    missing_coefficients <- setdiff(K, coefficients)
    if (length(missing_coefficients)) {
      stop("Unknown coefficient(s) in K: ",
           paste(missing_coefficients, collapse = ", "), call. = FALSE)
    }
    Kmat <- matrix(0, nrow = length(K), ncol = length(coefficients),
                   dimnames = list(K, coefficients))
    for (index in seq_along(K)) Kmat[index, K[index]] <- 1
    K <- Kmat
  } else if (is.vector(K)) {
    K <- matrix(K, nrow = 1L, dimnames = list(NULL, coefficients))
  }
  K <- as.matrix(K)
  if (!is.numeric(K) || any(!is.finite(K)) || nrow(K) < 1L) {
    stop("K must be a non-empty finite numeric contrast matrix",
         call. = FALSE)
  }
  if (ncol(K) != length(coefficients)) {
    stop("ncol(K) = ", ncol(K), " must equal number of coefficients (",
         length(coefficients), ")", call. = FALSE)
  }
  if (is.null(m)) m <- rep(0, nrow(K))
  if (!is.numeric(m) || length(m) != nrow(K) || any(!is.finite(m))) {
    stop("length(m) must equal nrow(K)", call. = FALSE)
  }
  estimate <- drop(K %*% as.numeric(region$estimate)) - m
  radius <- drop(abs(K) %*% as.numeric(region$mechanism_radius))
  lower <- estimate - radius
  upper <- estimate + radius
  if (any(!is.finite(c(estimate, radius, lower, upper))) || any(radius < 0)) {
    stop("The signed Gaussian mechanism contrast is not finite", call. = FALSE)
  }
  out <- list(
    estimate = estimate,
    lower = lower,
    upper = upper,
    mechanism_radius = radius,
    null = m,
    null_excluded = lower > 0 | upper < 0,
    all_null_components_excluded = all(lower > 0 | upper < 0),
    K = K,
    distribution = "simultaneous_dp_mechanism_region",
    level = attr(region, "level"),
    sampling_inference = FALSE,
    additional_server_calls = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0))
  class(out) <- c("ds.vertContrast", "list")
  out
}

#' @export
print.ds.vertContrast <- function(x, ...) {
  if (identical(x$distribution, "simultaneous_dp_mechanism_region")) {
    cat("dsVert signed DP-mechanism linear contrast region\n")
    cat("  componentwise null exclusion (not a sampling test):\n")
    print(data.frame(estimate = x$estimate, lower = x$lower, upper = x$upper,
                     null_excluded = x$null_excluded,
                     check.names = FALSE))
    return(invisible(x))
  }
  cat(sprintf("dsVert multi-coefficient Wald / linear contrast test\n"))
  if (identical(x$distribution, "F")) {
    cat(sprintf("  F = %.4f on %d and %d df,  p-value = %s\n",
                x$statistic, x$df, x$df_residual,
                format.pval(x$p_value, digits = 4L)))
  } else {
    cat(sprintf("  chi-sq = %.4f on %d df,  p-value = %s\n",
                x$statistic, x$df,
                format.pval(x$p_value, digits = 4L)))
  }
  cat("  K * beta - m estimates:\n")
  print(x$estimate)
  invisible(x)
}
