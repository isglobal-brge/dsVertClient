# Deterministic L1 post-processing for a validated Gaussian DP capsule.
# Source authentication and scale conversion are kept separate from this
# numerical core so the solver never needs a connection or a server call.

.dsvert_lasso_dp_moments <- function(
    gram_projected, cross_projected, outcome_square_projected, n_obs) {
  gram <- as.matrix(gram_projected)
  cross_names <- names(cross_projected)
  cross <- as.numeric(cross_projected)
  names(cross) <- cross_names
  terms <- rownames(gram)
  p <- nrow(gram)
  valid <- is.numeric(gram) && length(dim(gram)) == 2L && p > 0L &&
    ncol(gram) == p && all(is.finite(gram)) &&
    !is.null(terms) && !anyNA(terms) && !any(!nzchar(terms)) &&
    !anyDuplicated(terms) && identical(colnames(gram), terms) &&
    length(cross) == p && all(is.finite(cross)) &&
    identical(names(cross), terms) &&
    is.numeric(outcome_square_projected) &&
    length(outcome_square_projected) == 1L &&
    is.finite(outcome_square_projected) &&
    is.numeric(n_obs) && length(n_obs) == 1L && is.finite(n_obs) &&
    n_obs > 0
  if (!isTRUE(valid)) {
    stop("The certified DP Gaussian projected moments are malformed",
         call. = FALSE)
  }
  augmented <- rbind(cbind(gram, cross), c(cross, outcome_square_projected))
  scale <- max(1, max(abs(augmented)))
  symmetry_tolerance <- 256 * .Machine$double.eps * scale * (p + 1L)
  if (max(abs(augmented - t(augmented))) > symmetry_tolerance) {
    stop("The certified DP Gaussian augmented moment matrix is not symmetric",
         call. = FALSE)
  }
  eigenvalues <- eigen(
    augmented, symmetric = TRUE, only.values = TRUE)$values
  psd_tolerance <- 256 * .Machine$double.eps *
    max(1, max(abs(eigenvalues))) * (p + 1L)
  if (min(eigenvalues) < -psd_tolerance) {
    .dsvert_stop_non_identifiable(
      paste0(
        "The certified DP Gaussian augmented moment matrix is not ",
        "positive semidefinite; no client-side reprojection was applied."),
      reason = "non_psd_dp_lasso_augmented_moments")
  }
  list(
    n_obs = as.numeric(n_obs),
    gram_sum = gram,
    cross_sum = cross,
    outcome_square_sum = as.numeric(outcome_square_projected),
    gram = gram / n_obs,
    cross = cross / n_obs,
    outcome_square = as.numeric(outcome_square_projected) / n_obs,
    normalization = paste0(
      "projected_DP_moment_sums_divided_by_the_separate_positive_",
      "DP_count_coordinate_without_count_Gram11_averaging"),
    augmented_certificate = list(
      positive_semidefinite = TRUE,
      minimum_eigenvalue = min(eigenvalues),
      maximum_eigenvalue = max(eigenvalues),
      psd_tolerance = psd_tolerance,
      client_reprojection_applied = FALSE))
}

.dsvert_lasso_dp_artifact_scales <- function(artifact) {
  predictor_order <- artifact$predictor_order
  predictors <- artifact$predictors
  outcome <- artifact$outcome
  valid <- is.list(artifact) && is.character(predictor_order) &&
    length(predictor_order) > 0L && !anyNA(predictor_order) &&
    !any(!nzchar(predictor_order)) && !anyDuplicated(predictor_order) &&
    is.list(predictors) && identical(names(predictors), predictor_order) &&
    is.list(outcome) && is.numeric(outcome$lower) &&
    length(outcome$lower) == 1L && is.finite(outcome$lower) &&
    is.numeric(outcome$upper) && length(outcome$upper) == 1L &&
    is.finite(outcome$upper) && outcome$upper > outcome$lower &&
    is.logical(artifact$intercept) && length(artifact$intercept) == 1L &&
    !is.na(artifact$intercept)
  if (!isTRUE(valid)) {
    stop("The validated DP Gaussian scale contract is malformed",
         call. = FALSE)
  }
  lowers <- vapply(predictor_order, function(variable) {
    value <- predictors[[variable]]
    if (!is.list(value) || !is.numeric(value$lower) ||
        length(value$lower) != 1L || !is.finite(value$lower) ||
        !is.numeric(value$upper) || length(value$upper) != 1L ||
        !is.finite(value$upper) || value$upper <= value$lower) {
      stop("The validated DP Gaussian predictor bounds are malformed",
           call. = FALSE)
    }
    as.numeric(value$lower)
  }, numeric(1L))
  uppers <- vapply(
    predictor_order,
    function(variable) as.numeric(predictors[[variable]]$upper),
    numeric(1L))
  predictor_ranges <- uppers - lowers
  outcome_range <- as.numeric(outcome$upper - outcome$lower)
  list(
    predictor_order = predictor_order,
    predictor_lowers = lowers,
    predictor_ranges = predictor_ranges,
    outcome_lower = as.numeric(outcome$lower),
    outcome_range = outcome_range,
    intercept = isTRUE(artifact$intercept),
    normalized_terms = c(
      if (isTRUE(artifact$intercept)) "(Intercept)" else character(),
      predictor_order),
    original_penalty_weights = predictor_ranges / outcome_range)
}

.dsvert_lasso_dp_original_scale <- function(theta, artifact) {
  scales <- .dsvert_lasso_dp_artifact_scales(artifact)
  if (!is.numeric(theta) || length(theta) != length(scales$normalized_terms) ||
      any(!is.finite(theta)) || is.null(names(theta)) ||
      !setequal(names(theta), scales$normalized_terms)) {
    stop("Normalized coefficients do not match the signed DP design terms",
         call. = FALSE)
  }
  theta <- theta[scales$normalized_terms]
  slopes <- scales$outcome_range *
    theta[scales$predictor_order] / scales$predictor_ranges
  normalized_intercept <- if (scales$intercept) {
    theta[["(Intercept)"]]
  } else {
    0
  }
  intercept <- scales$outcome_lower +
    scales$outcome_range * normalized_intercept -
    sum(slopes * scales$predictor_lowers)
  c("(Intercept)" = intercept, slopes)
}

.dsvert_lasso_dp_normalized_scale <- function(coefficients, artifact) {
  scales <- .dsvert_lasso_dp_artifact_scales(artifact)
  expected <- c("(Intercept)", scales$predictor_order)
  if (!is.numeric(coefficients) || any(!is.finite(coefficients)) ||
      is.null(names(coefficients)) || anyNA(names(coefficients)) ||
      any(!nzchar(names(coefficients))) || anyDuplicated(names(coefficients)) ||
      !setequal(names(coefficients), expected)) {
    stop("Original-scale coefficients do not match the signed DP model",
         call. = FALSE)
  }
  coefficients <- coefficients[expected]
  theta_slopes <- coefficients[scales$predictor_order] *
    scales$predictor_ranges / scales$outcome_range
  if (scales$intercept) {
    theta_intercept <- (
      coefficients[["(Intercept)"]] - scales$outcome_lower +
        sum(coefficients[scales$predictor_order] *
              scales$predictor_lowers)) / scales$outcome_range
    return(c("(Intercept)" = theta_intercept, theta_slopes))
  }
  implied <- .dsvert_lasso_dp_original_scale(theta_slopes, artifact)
  tolerance <- 1024 * .Machine$double.eps *
    max(1, abs(implied[["(Intercept)"]]))
  if (abs(coefficients[["(Intercept)"]] -
          implied[["(Intercept)"]]) > tolerance) {
    stop(paste0(
      "For a no-intercept normalized model, the original-scale intercept ",
      "must equal the deterministic bound-transform offset."),
    call. = FALSE)
  }
  theta_slopes
}

.dsvert_lasso_dp_lambda_max <- function(
    gram, cross, keep_intercept = TRUE) {
  gram <- as.matrix(gram)
  cross_names <- names(cross)
  cross <- as.numeric(cross)
  names(cross) <- cross_names
  terms <- rownames(gram)
  p <- nrow(gram)
  if (!is.numeric(gram) || length(dim(gram)) != 2L || p < 1L ||
      ncol(gram) != p || any(!is.finite(gram)) ||
      is.null(terms) || !identical(colnames(gram), terms) ||
      length(cross) != p || any(!is.finite(cross)) ||
      !identical(names(cross), terms) ||
      !is.logical(keep_intercept) || length(keep_intercept) != 1L ||
      is.na(keep_intercept)) {
    stop("Invalid DP normal equations for lambda_max", call. = FALSE)
  }
  intercept <- match("(Intercept)", terms, nomatch = 0L)
  unpenalised <- if (keep_intercept && intercept > 0L) {
    intercept
  } else {
    integer()
  }
  penalised <- setdiff(seq_len(p), unpenalised)
  if (!length(penalised)) return(0)
  beta_unpenalised <- numeric(length(unpenalised))
  if (length(unpenalised)) {
    curvature <- gram[unpenalised, unpenalised]
    scale <- max(1, abs(curvature), abs(cross[unpenalised]))
    tolerance <- 256 * .Machine$double.eps * scale
    if (curvature <= tolerance) {
      if (abs(cross[unpenalised]) > tolerance) {
        .dsvert_stop_non_identifiable(
          "The unpenalized DP intercept objective is unbounded.",
          reason = "unbounded_dp_lasso_intercept")
      }
      beta_unpenalised[] <- 0
    } else {
      beta_unpenalised <- cross[unpenalised] / curvature
    }
  }
  residual_score <- cross[penalised]
  if (length(unpenalised)) {
    residual_score <- residual_score -
      gram[penalised, unpenalised, drop = FALSE] %*% beta_unpenalised
  }
  as.numeric(max(abs(residual_score)))
}

.dsvert_lasso_dp_default_lambda <- function(
    gram, cross, keep_intercept = TRUE, length_out = 50L,
    minimum_ratio = 1e-3) {
  if (!is.numeric(length_out) || length(length_out) != 1L ||
      !is.finite(length_out) || length_out < 1 ||
      length_out != floor(length_out) ||
      !is.numeric(minimum_ratio) || length(minimum_ratio) != 1L ||
      !is.finite(minimum_ratio) || minimum_ratio <= 0 ||
      minimum_ratio > 1) {
    stop("Invalid default lambda-grid contract", call. = FALSE)
  }
  maximum <- .dsvert_lasso_dp_lambda_max(
    gram, cross, keep_intercept = keep_intercept)
  if (!is.finite(maximum) || maximum < 0) {
    stop("The signed DP lambda_max is invalid", call. = FALSE)
  }
  if (maximum == 0 || length_out == 1L) return(maximum)
  exp(seq(log(maximum), log(maximum * minimum_ratio),
          length.out = as.integer(length_out)))
}

.dsvert_lasso_dp_solver <- function(
    gram, cross, lambda, max_iter, tol, keep_intercept,
    warm_start = NULL) {
  gram <- as.matrix(gram)
  cross_names <- names(cross)
  cross <- as.numeric(cross)
  names(cross) <- cross_names
  term_names <- rownames(gram)
  p <- nrow(gram)
  valid_matrix <- is.numeric(gram) && length(dim(gram)) == 2L &&
    p > 0L && ncol(gram) == p && all(is.finite(gram)) &&
    !is.null(term_names) && !is.null(colnames(gram)) &&
    !anyNA(term_names) && !any(!nzchar(term_names)) &&
    !anyDuplicated(term_names) &&
    identical(colnames(gram), term_names)
  valid_cross <- length(cross) == p && all(is.finite(cross)) &&
    !is.null(cross_names) && identical(cross_names, term_names)
  if (!isTRUE(valid_matrix) || !isTRUE(valid_cross)) {
    stop("The validated DP LASSO normal equations are malformed",
         call. = FALSE)
  }
  if (!is.numeric(lambda) || length(lambda) != 1L ||
      !is.finite(lambda) || lambda < 0) {
    stop("lambda must be one finite non-negative number", call. = FALSE)
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L ||
      !is.finite(max_iter) || max_iter < 1 || max_iter != floor(max_iter)) {
    stop("max_iter must be one positive integer", call. = FALSE)
  }
  max_iter <- as.integer(max_iter)
  if (!is.numeric(tol) || length(tol) != 1L ||
      !is.finite(tol) || tol <= 0) {
    stop("tol must be one finite positive number", call. = FALSE)
  }
  if (!is.logical(keep_intercept) || length(keep_intercept) != 1L ||
      is.na(keep_intercept)) {
    stop("keep_intercept must be TRUE or FALSE", call. = FALSE)
  }

  matrix_scale <- max(1, max(abs(gram)))
  symmetry_tolerance <- 256 * .Machine$double.eps * matrix_scale * p
  if (max(abs(gram - t(gram))) > symmetry_tolerance) {
    stop("The certified DP Gaussian Gram matrix is not symmetric",
         call. = FALSE)
  }
  # The Gaussian producer already performs its declared augmented-moment PSD
  # projection. Do not project, ridge or otherwise alter that matrix here.
  eigenvalues <- eigen(gram, symmetric = TRUE, only.values = TRUE)$values
  psd_tolerance <- 256 * .Machine$double.eps *
    max(1, max(abs(eigenvalues))) * p
  if (min(eigenvalues) < -psd_tolerance) {
    .dsvert_stop_non_identifiable(
      paste0(
        "The certified DP Gaussian Gram matrix is not positive ",
        "semidefinite; no client-side PSD or ridge repair was applied."),
      reason = "non_psd_dp_lasso_information")
  }
  strong_convexity <- min(eigenvalues) > psd_tolerance
  intercept_index <- match("(Intercept)", term_names, nomatch = 0L)
  unpenalised <- if (isTRUE(keep_intercept) && intercept_index > 0L) {
    intercept_index
  } else {
    integer()
  }
  penalised <- setdiff(seq_len(p), unpenalised)

  if (is.null(warm_start)) {
    beta <- stats::setNames(numeric(p), term_names)
  } else {
    if (!is.numeric(warm_start) || length(warm_start) != p ||
        any(!is.finite(warm_start))) {
      stop("warm_start must have one finite value per DP design term",
           call. = FALSE)
    }
    if (!is.null(names(warm_start))) {
      if (anyNA(names(warm_start)) || any(!nzchar(names(warm_start))) ||
          anyDuplicated(names(warm_start)) ||
          !setequal(names(warm_start), term_names)) {
        stop("named warm_start values must match the DP design terms",
             call. = FALSE)
      }
      warm_start <- warm_start[term_names]
    }
    beta <- stats::setNames(as.numeric(warm_start), term_names)
  }

  soft_threshold <- function(value, threshold) {
    sign(value) * pmax(abs(value) - threshold, 0)
  }
  diagonal <- diag(gram)
  coordinate_tolerance <- psd_tolerance
  converged_by_step <- FALSE
  iterations <- max_iter
  for (iteration in seq_len(max_iter)) {
    previous <- beta
    for (index in seq_len(p)) {
      curvature <- diagonal[[index]]
      if (curvature <= coordinate_tolerance) {
        # A zero diagonal in a PSD matrix implies a zero row. The deterministic
        # zero representative is KKT-valid whenever the objective is bounded.
        beta[[index]] <- 0
        next
      }
      partial <- cross[[index]] -
        sum(gram[index, -index, drop = TRUE] * beta[-index])
      beta[[index]] <- if (index %in% unpenalised) {
        partial / curvature
      } else {
        soft_threshold(partial, lambda) / curvature
      }
    }
    step <- max(abs(beta - previous))
    if (step <= tol * max(1, max(abs(beta)))) {
      converged_by_step <- TRUE
      iterations <- iteration
      break
    }
  }

  gradient <- as.numeric(gram %*% beta - cross)
  names(gradient) <- term_names
  support_tolerance <- max(tol, sqrt(.Machine$double.eps))
  violations <- numeric(p)
  if (length(unpenalised)) {
    violations[unpenalised] <- abs(gradient[unpenalised])
  }
  active <- penalised[abs(beta[penalised]) > support_tolerance]
  inactive <- setdiff(penalised, active)
  if (length(active)) {
    violations[active] <- abs(
      gradient[active] + lambda * sign(beta[active]))
  }
  if (length(inactive)) {
    violations[inactive] <- pmax(abs(gradient[inactive]) - lambda, 0)
  }
  names(violations) <- term_names
  kkt_scale <- max(1, max(abs(cross)), max(abs(gram %*% beta)), lambda)
  kkt_tolerance <- max(100 * tol, 1024 * .Machine$double.eps) * kkt_scale
  max_violation <- max(violations)
  kkt_satisfied <- is.finite(max_violation) &&
    max_violation <= kkt_tolerance
  if (!isTRUE(kkt_satisfied)) {
    zero_curvature_violation <- which(
      diagonal <= coordinate_tolerance & violations > kkt_tolerance)
    if (length(zero_curvature_violation)) {
      .dsvert_stop_non_identifiable(
        paste0(
          "The signed DP LASSO objective has effectively zero curvature ",
          "with an incompatible score for: ",
          paste(term_names[zero_curvature_violation], collapse = ", "),
          ". A bounded identifiable solution cannot be certified at the ",
          "declared numeric tolerance."),
        reason = "zero_curvature_dp_lasso_score")
    }
    .dsvert_stop_numeric(
      "lasso_non_convergence",
      paste0(
        "The deterministic DP LASSO solver did not satisfy its KKT ",
        "certificate after ", max_iter, " iterations."),
      reason = "dp_lasso_kkt_not_satisfied")
  }

  quadratic <- 0.5 * sum(beta * as.numeric(gram %*% beta)) -
    sum(cross * beta)
  l1 <- lambda * sum(abs(beta[penalised]))
  list(
    coefficients = beta,
    objective_without_constant = quadratic + l1,
    quadratic_without_constant = quadratic,
    l1_penalty = l1,
    iterations = as.integer(iterations),
    converged = isTRUE(converged_by_step) || isTRUE(kkt_satisfied),
    kkt = list(
      satisfied = TRUE, max_violation = max_violation,
      tolerance = kkt_tolerance, gradient = gradient,
      coordinate_violations = violations),
    curvature = list(
      positive_semidefinite = TRUE,
      minimum_eigenvalue = min(eigenvalues),
      maximum_eigenvalue = max(eigenvalues),
      psd_tolerance = psd_tolerance,
      strong_convexity = strong_convexity,
      solution_unique = if (strong_convexity) TRUE else NA,
      solution_uniqueness_certified = strong_convexity,
      warning = if (strong_convexity) NULL else paste(
        "The DP Gram matrix is PSD but not strongly convex; this",
        "deterministic solution satisfies KKT, but neither uniqueness nor",
        "non-uniqueness is certified.")),
    implicit_ridge = FALSE,
    client_psd_projection_applied = FALSE)
}

.dsvert_lasso_dp_lambda_path <- function(
    gram, cross, lambda, max_iter, tol, keep_intercept,
    warm_start = NULL) {
  if (!is.numeric(lambda) || !length(lambda) || any(!is.finite(lambda)) ||
      any(lambda < 0)) {
    stop("lambda must contain finite non-negative values", call. = FALSE)
  }
  solutions <- vector("list", length(lambda))
  names(solutions) <- make.unique(vapply(
    lambda, function(value) sprintf("%.17g", value), character(1L)))
  current <- warm_start
  for (index in seq_along(lambda)) {
    solutions[[index]] <- .dsvert_lasso_dp_solver(
      gram = gram, cross = cross, lambda = lambda[[index]],
      max_iter = max_iter, tol = tol, keep_intercept = keep_intercept,
      warm_start = current)
    current <- solutions[[index]]$coefficients
  }
  solutions
}

.dsvert_lasso_dp_pseudo_ic <- function(
    moments, solutions, lambda, criterion = c("BIC", "AIC", "EBIC"),
    ebic_gamma = 0.5, keep_intercept = TRUE,
    support_tolerance = sqrt(.Machine$double.eps),
    parsimonious_delta = 0.02) {
  criterion <- match.arg(criterion)
  if (!is.list(moments) || !is.numeric(moments$n_obs) ||
      length(moments$n_obs) != 1L || !is.finite(moments$n_obs) ||
      moments$n_obs <= 0 || !is.matrix(moments$gram) ||
      !is.numeric(moments$cross) ||
      !is.numeric(moments$outcome_square) ||
      length(moments$outcome_square) != 1L ||
      !is.finite(moments$outcome_square)) {
    stop("Invalid normalized DP moments for pseudo-IC", call. = FALSE)
  }
  if (!is.list(solutions) || !length(solutions) ||
      !is.numeric(lambda) || length(lambda) != length(solutions) ||
      any(!is.finite(lambda)) || any(lambda < 0)) {
    stop("Invalid DP LASSO path for pseudo-IC", call. = FALSE)
  }
  if (!is.numeric(ebic_gamma) || length(ebic_gamma) != 1L ||
      !is.finite(ebic_gamma) || ebic_gamma < 0 || ebic_gamma > 1) {
    stop("ebic_gamma must be one number in [0, 1]", call. = FALSE)
  }
  if (!is.logical(keep_intercept) || length(keep_intercept) != 1L ||
      is.na(keep_intercept)) {
    stop("keep_intercept must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(support_tolerance) || length(support_tolerance) != 1L ||
      !is.finite(support_tolerance) || support_tolerance < 0 ||
      !is.numeric(parsimonious_delta) || length(parsimonious_delta) != 1L ||
      !is.finite(parsimonious_delta) || parsimonious_delta < 0) {
    stop("Invalid DP pseudo-IC tolerance", call. = FALSE)
  }
  terms <- rownames(moments$gram)
  intercept <- match("(Intercept)", terms, nomatch = 0L)
  unpenalised <- if (keep_intercept && intercept > 0L) {
    intercept
  } else {
    integer()
  }
  penalised <- setdiff(seq_along(terms), unpenalised)
  coefficients <- lapply(solutions, function(solution) {
    value <- if (is.list(solution)) solution$coefficients else solution
    if (!is.numeric(value) || length(value) != length(terms) ||
        any(!is.finite(value)) || is.null(names(value)) ||
        !setequal(names(value), terms)) {
      stop("A DP LASSO path coefficient vector is malformed",
           call. = FALSE)
    }
    value[terms]
  })
  rss_mean <- vapply(coefficients, function(beta) {
    as.numeric(moments$outcome_square -
      2 * sum(moments$cross * beta) +
      sum(beta * as.numeric(moments$gram %*% beta)))
  }, numeric(1L))
  active_penalised <- vapply(coefficients, function(beta) {
    sum(abs(beta[penalised]) > support_tolerance)
  }, integer(1L))
  degrees_freedom <- active_penalised + length(unpenalised)
  unavailable <- function(reason) {
    list(
      available = FALSE, reason = reason,
      criterion = paste0("DP_projected_pseudo_", criterion),
      classical_information_criterion = FALSE,
      cross_validation = FALSE, one_standard_error_rule = FALSE,
      ic = rep(NA_real_, length(solutions)),
      rss_mean_dp_projected = rss_mean,
      rss_sum_dp_projected = moments$n_obs * rss_mean,
      df = degrees_freedom,
      active_penalised = active_penalised,
      lambda.min = NULL, lambda.parsimonious = NULL,
      lambda.1se = NULL, parsimonious_delta = parsimonious_delta)
  }
  if (moments$n_obs <= 1) {
    return(unavailable("dp_effective_count_not_greater_than_one"))
  }
  rss_tolerance <- 1024 * .Machine$double.eps * max(
    1, abs(moments$outcome_square),
    max(abs(vapply(coefficients, function(beta) {
      sum(beta * as.numeric(moments$gram %*% beta))
    }, numeric(1L)))))
  if (any(!is.finite(rss_mean)) || any(rss_mean <= rss_tolerance)) {
    return(unavailable("non_positive_dp_projected_residual_mean_square"))
  }
  base <- moments$n_obs * log(rss_mean)
  ic <- switch(
    criterion,
    AIC = base + 2 * degrees_freedom,
    BIC = base + log(moments$n_obs) * degrees_freedom,
    EBIC = base + log(moments$n_obs) * degrees_freedom +
      2 * ebic_gamma * vapply(active_penalised, function(active) {
        lchoose(length(penalised), active)
      }, numeric(1L)))
  minimum_index <- which.min(ic)
  viable <- which(ic <= ic[[minimum_index]] + parsimonious_delta)
  parsimonious_index <- viable[[which.max(lambda[viable])]]
  list(
    available = TRUE, reason = NULL,
    criterion = paste0("DP_projected_pseudo_", criterion),
    classical_information_criterion = FALSE,
    cross_validation = FALSE, one_standard_error_rule = FALSE,
    ic = ic,
    rss_mean_dp_projected = rss_mean,
    rss_sum_dp_projected = moments$n_obs * rss_mean,
    df = degrees_freedom,
    active_penalised = active_penalised,
    lambda.min = lambda[[minimum_index]],
    lambda.parsimonious = lambda[[parsimonious_index]],
    lambda.1se = lambda[[parsimonious_index]],
    minimum_index = as.integer(minimum_index),
    parsimonious_index = as.integer(parsimonious_index),
    parsimonious_delta = parsimonious_delta)
}

.dsvert_lasso_dp_source <- function(fit, trusted_pinset = NULL) {
  if (!inherits(fit, "ds.vertDPGaussian") || !is.list(fit) ||
      !identical(fit$status, "ok") || !identical(fit$family, "gaussian") ||
      !identical(fit$source_values_exposed, FALSE) ||
      !identical(fit$intermediate_values_exposed, FALSE) ||
      !is.list(fit$signed_artifact) ||
      !is.list(fit$sufficient_statistics_dp)) {
    stop("fit is not a valid ds.vertDPGaussian release", call. = FALSE)
  }
  verification <- tryCatch(
    ds.validateDPGaussianCertificate(fit, trusted_pinset = trusted_pinset),
    error = function(error) {
      stop(
        "The Gaussian DP provenance certificate is invalid: ",
        conditionMessage(error), call. = FALSE)
    })
  if (!is.list(verification) ||
      !identical(verification$integrity_valid, TRUE)) {
    stop("The Gaussian DP provenance certificate failed integrity validation",
         call. = FALSE)
  }
  if (!identical(verification$authenticity, "caller_anchored") &&
      !identical(verification$authenticity, "session_transport_anchored")) {
    stop("The Gaussian DP provenance certificate is not anchored to trusted peers",
         call. = FALSE)
  }
  statistics <- fit$sufficient_statistics_dp
  required <- c(
    "gram_projected", "cross_projected", "outcome_square_projected")
  if (is.null(names(statistics)) ||
      !all(required %in% names(statistics))) {
    stop("The Gaussian DP sufficient statistics are incomplete",
         call. = FALSE)
  }
  moments <- .dsvert_lasso_dp_moments(
    gram_projected = statistics$gram_projected,
    cross_projected = statistics$cross_projected,
    outcome_square_projected = statistics$outcome_square_projected,
    n_obs = fit$n_obs)
  scales <- .dsvert_lasso_dp_artifact_scales(fit$signed_artifact)
  if (!identical(rownames(moments$gram), scales$normalized_terms)) {
    stop("The Gaussian DP moments do not match the signed design terms",
         call. = FALSE)
  }
  list(
    moments = moments, artifact = fit$signed_artifact, scales = scales,
    verification = verification)
}

.dsvert_lasso_dp_warm_start <- function(warm_start, source) {
  if (is.null(warm_start)) return(NULL)
  expected <- c("(Intercept)", source$scales$predictor_order)
  if (!is.numeric(warm_start) || length(warm_start) != length(expected) ||
      any(!is.finite(warm_start))) {
    stop(paste(
      "For a ds.vertDPGaussian fit, warm_start must have one finite",
      "original-scale value for the intercept and every predictor"),
    call. = FALSE)
  }
  if (is.null(names(warm_start))) {
    names(warm_start) <- expected
  } else if (anyNA(names(warm_start)) || any(!nzchar(names(warm_start))) ||
             anyDuplicated(names(warm_start)) ||
             !setequal(names(warm_start), expected)) {
    stop(paste(
      "For a ds.vertDPGaussian fit, named warm_start values must match",
      "the original-scale coefficient names"), call. = FALSE)
  }
  .dsvert_lasso_dp_normalized_scale(
    warm_start[expected], source$artifact)
}

.dsvert_lasso_dp_proximal <- function(
    fit, lambda, max_iter, tol, keep_intercept, warm_start, accelerate,
    trusted_pinset = NULL) {
  if (!is.numeric(lambda) || length(lambda) != 1L ||
      !is.finite(lambda) || lambda < 0) {
    stop("lambda must be one finite non-negative number", call. = FALSE)
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L ||
      !is.finite(max_iter) || max_iter < 1 || max_iter != floor(max_iter)) {
    stop("max_iter must be one positive integer", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L ||
      !is.finite(tol) || tol <= 0) {
    stop("tol must be one finite positive number", call. = FALSE)
  }
  if (!is.logical(keep_intercept) || length(keep_intercept) != 1L ||
      is.na(keep_intercept)) {
    stop("keep_intercept must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(accelerate) || length(accelerate) != 1L ||
      is.na(accelerate)) {
    stop("accelerate must be TRUE or FALSE", call. = FALSE)
  }
  source <- .dsvert_lasso_dp_source(fit, trusted_pinset = trusted_pinset)
  warm_normalized <- .dsvert_lasso_dp_warm_start(warm_start, source)
  solution <- .dsvert_lasso_dp_solver(
    gram = source$moments$gram,
    cross = source$moments$cross,
    lambda = as.numeric(lambda), max_iter = as.integer(max_iter),
    tol = as.numeric(tol), keep_intercept = keep_intercept,
    warm_start = warm_normalized)
  coefficients <- .dsvert_lasso_dp_original_scale(
    solution$coefficients, source$artifact)
  eigenvalues <- eigen(
    source$moments$gram, symmetric = TRUE, only.values = TRUE)$values
  support <- which(abs(solution$coefficients) > tol)
  source_penalty <- if (is.numeric(fit$ridge) && length(fit$ridge) == 1L &&
                        is.finite(fit$ridge)) {
    as.numeric(fit$ridge)
  } else {
    NA_real_
  }
  out <- list(
    coefficients = coefficients,
    coefficients_normalized = solution$coefficients,
    support = support,
    support_terms = names(solution$coefficients)[support],
    lambda = as.numeric(lambda),
    lambda_scale = "signed_bound_normalized_gaussian_objective",
    converged = solution$converged,
    iterations = solution$iterations,
    L = max(eigenvalues), step_size = NA_real_,
    XtX_over_n = source$moments$gram,
    cross_over_n = source$moments$cross,
    outcome_square_over_n = source$moments$outcome_square,
    sigma2_hat = NA_real_,
    objective = solution$objective_without_constant,
    n_obs = source$moments$n_obs,
    family = "gaussian",
    estimand = "bounded_normalized_DP_gaussian_lasso",
    solver = "deterministic_coordinate_descent_on_signed_DP_moments",
    source_fit_penalty = source_penalty,
    source_fit_penalty_affects_released_moments = FALSE,
    input_provenance = "signed_dp_gaussian_capsule",
    source_certificate_validation = source$verification,
    source_query_contract_sha256 = fit$query_contract_sha256,
    source_release_contract_hash = fit$release_contract_hash,
    source_accuracy = fit$accuracy,
    kkt = solution$kkt,
    curvature = solution$curvature,
    augmented_moment_certificate =
      source$moments$augmented_certificate,
    coefficient_scale = list(
      normalized_terms = source$scales$normalized_terms,
      original_predictor_penalty_weights =
        source$scales$original_penalty_weights,
      normalized_intercept_penalized =
        !isTRUE(keep_intercept) && isTRUE(source$scales$intercept)),
    inference = list(
      classical_standard_errors = NULL, p_values = NULL,
      confidence_intervals = NULL, sampling_inference_available = FALSE),
    coefficient_regions_available = FALSE,
    additional_server_calls_after_capsule = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    implicit_ridge = FALSE,
    client_psd_projection_applied = FALSE,
    compatibility = list(
      accelerate_argument = accelerate, accelerate_effective = FALSE,
      legacy_ds_glm_route_used = FALSE),
    comparison = list(
      coefficients_soft = NULL,
      note = "Naive soft-thresholding is not reported for DP moments"),
    disclosure_guard = list(
      satisfied = TRUE,
      basis = "deterministic_postprocessing_of_one_validated_DP_capsule"))
  class(out) <- c("ds.vertLASSOProximal", "ds.vertDPLASSO", "list")
  out
}

.dsvert_lasso_dp_select <- function(
    fit, lambda_grid, criterion, ebic_gamma, keep_intercept,
    se_threshold, trusted_pinset = NULL) {
  if (!is.numeric(ebic_gamma) || length(ebic_gamma) != 1L ||
      !is.finite(ebic_gamma) || ebic_gamma < 0 || ebic_gamma > 1) {
    stop("ebic_gamma must be one number in [0, 1]", call. = FALSE)
  }
  if (!is.logical(keep_intercept) || length(keep_intercept) != 1L ||
      is.na(keep_intercept)) {
    stop("keep_intercept must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(se_threshold) || length(se_threshold) != 1L ||
      !is.finite(se_threshold) || se_threshold < 0) {
    stop("se_threshold must be one finite non-negative number",
         call. = FALSE)
  }
  source <- .dsvert_lasso_dp_source(fit, trusted_pinset = trusted_pinset)
  if (is.null(lambda_grid)) {
    lambda_grid <- .dsvert_lasso_dp_default_lambda(
      source$moments$gram, source$moments$cross,
      keep_intercept = keep_intercept)
  }
  if (!is.numeric(lambda_grid) || !length(lambda_grid) ||
      any(!is.finite(lambda_grid)) || any(lambda_grid < 0) ||
      anyDuplicated(lambda_grid)) {
    stop("lambda_grid must contain unique finite non-negative values",
         call. = FALSE)
  }
  lambda_grid <- as.numeric(lambda_grid)
  solutions <- .dsvert_lasso_dp_lambda_path(
    source$moments$gram, source$moments$cross, lambda_grid,
    max_iter = 5000L, tol = 1e-10,
    keep_intercept = keep_intercept)
  selected <- .dsvert_lasso_dp_pseudo_ic(
    source$moments, solutions, lambda_grid, criterion = criterion,
    ebic_gamma = ebic_gamma, keep_intercept = keep_intercept,
    parsimonious_delta = 0)
  if (isTRUE(selected$available)) {
    minimum_index <- which.min(selected$ic)
    relative_delta <- abs(selected$ic[[minimum_index]]) * se_threshold
    viable <- which(
      selected$ic <= selected$ic[[minimum_index]] + relative_delta)
    parsimonious_index <- viable[[which.max(lambda_grid[viable])]]
    selected$minimum_index <- as.integer(minimum_index)
    selected$parsimonious_index <- as.integer(parsimonious_index)
    selected$lambda.min <- lambda_grid[[minimum_index]]
    selected$lambda.parsimonious <- lambda_grid[[parsimonious_index]]
    selected$lambda.1se <- selected$lambda.parsimonious
    selected$parsimonious_delta <- relative_delta
  }
  normalized_paths <- lapply(solutions, `[[`, "coefficients")
  original_paths <- lapply(normalized_paths, function(coefficients) {
    .dsvert_lasso_dp_original_scale(coefficients, source$artifact)
  })
  names(original_paths) <- names(normalized_paths) <- names(solutions)
  minimum_index <- selected$minimum_index
  parsimonious_index <- selected$parsimonious_index
  out <- list(
    lambda = lambda_grid,
    ic = selected$ic,
    df = selected$df,
    criterion = criterion,
    lambda.min = selected$lambda.min,
    lambda.parsimonious = selected$lambda.parsimonious,
    lambda.1se = selected$lambda.1se,
    beta.min = if (is.null(minimum_index)) NULL else
      original_paths[[minimum_index]],
    beta.parsimonious = if (is.null(parsimonious_index)) NULL else
      original_paths[[parsimonious_index]],
    beta.1se = if (is.null(parsimonious_index)) NULL else
      original_paths[[parsimonious_index]],
    beta.normalized.min = if (is.null(minimum_index)) NULL else
      normalized_paths[[minimum_index]],
    beta.normalized.parsimonious =
      if (is.null(parsimonious_index)) NULL else
        normalized_paths[[parsimonious_index]],
    selection_available = selected$available,
    selection_unavailable_reason = selected$reason,
    selection_method = "DP_projected_pseudo_information_criterion",
    classical_information_criterion = FALSE,
    cross_validation = FALSE,
    one_standard_error_rule = FALSE,
    relative_ic_tolerance = se_threshold,
    absolute_ic_tolerance = selected$parsimonious_delta,
    paths = original_paths,
    paths_normalized = normalized_paths,
    path_certificates = lapply(solutions, function(solution) {
      list(kkt = solution$kkt, curvature = solution$curvature)
    }),
    fit = fit,
    input_provenance = "signed_dp_gaussian_capsule",
    source_certificate_validation = source$verification,
    source_query_contract_sha256 = fit$query_contract_sha256,
    source_release_contract_hash = fit$release_contract_hash,
    inference = list(
      classical_standard_errors = NULL, p_values = NULL,
      confidence_intervals = NULL, sampling_inference_available = FALSE),
    additional_server_calls_after_capsule = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    disclosure_guard = list(
      satisfied = TRUE,
      basis = "deterministic_postprocessing_of_one_validated_DP_capsule"))
  class(out) <- c("ds.vertLASSOCV", "ds.vertDPLASSOSelect", "list")
  out
}
