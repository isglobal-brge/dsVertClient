#' @title Signed independent-working GEE point estimate
#' @description With a custodian-configured \code{formal_analysis_id}, this
#'   frontdoor returns the completed, two-authority-certified binomial or
#'   Poisson GLM point estimate under an independence working correlation.
#'   \code{fresh_formal_analysis_id} first completes that same configured
#'   durable GLM analysis, then consumes its public result.
#'   With \code{dp_analysis_id}, it reads the signed Gaussian Synopsis fit
#'   under that same score equation. Neither route chooses privacy controls or
#'   exposes a cluster statistic.
#' @details This is deliberately narrower than a general GEE: cluster ids,
#'   exchangeable and AR(1) correlations, sandwich covariance, standard errors
#'   and inference remain unavailable until their contribution-bounded,
#'   protected artifacts exist. Calls without a configured selector keep
#'   failing locally before DSI.
#' @param formula,data,family,corstr,verbose,datasources Formula, registered
#'   data name, Gaussian/binomial/Poisson family, working correlation, progress flag and
#'   DataSHIELD connections. Only \code{corstr = "independence"} is available.
#' @param formal_analysis_id Custodian-configured completed formal GLM
#'   certificate selector.
#' @param fresh_formal_analysis_id Custodian-configured binomial/Poisson
#'   durable GLM selector. It is mutually exclusive with the completed
#'   certificate and Gaussian Synopsis selectors.
#' @param dp_analysis_id Custodian-configured signed Gaussian Synopsis
#'   artifact selector. It is mutually exclusive with both formal selectors.
#' @param id_col,order_col,max_iter,tol,lambda,working_max_iter,ring,binomial_sigmoid_intervals
#'   Retained clustered-GEE controls. They are unavailable with a formal point
#'   release and are never silently ignored.
#' @return A \code{ds.vertGEE} point-estimate object without covariance,
#'   standard errors, working-correlation estimate or inference.
#' @seealso \code{\link{ds.vertMethodStatus}}
#' @export
ds.vertGEE <- function(formula, data = NULL,
                       family = c("gaussian", "binomial", "poisson"),
                       id_col = NULL,
                       order_col = NULL,
                       corstr = c("independence", "exchangeable", "ar1"),
                       max_iter = 100L, tol = 1e-4, lambda = 0,
                       working_max_iter = NULL,
                       ring = 63L,
                       binomial_sigmoid_intervals = NULL,
                       verbose = TRUE, datasources = NULL,
                       formal_analysis_id = NULL,
                       fresh_formal_analysis_id = NULL,
                       dp_analysis_id = NULL) {
  selected_analysis_ids <- sum(!vapply(
    list(formal_analysis_id, fresh_formal_analysis_id, dp_analysis_id),
    is.null, logical(1L)))
  if (selected_analysis_ids > 1L) {
    stop(paste(
      "formal_analysis_id, fresh_formal_analysis_id and dp_analysis_id are",
      "mutually exclusive"),
         call. = FALSE)
  }
  if (!is.null(dp_analysis_id)) {
    family <- match.arg(family)
    corstr <- match.arg(corstr)
    return(.dsvert_gaussian_gee_independence_adapter(
      explicit_arguments = names(match.call())[-1L],
      formula = if (missing(formula)) NULL else formula,
      data = data, family = family, id_col = id_col, order_col = order_col,
      corstr = corstr, verbose = verbose, datasources = datasources,
      dp_analysis_id = dp_analysis_id))
  }
  if (!is.null(formal_analysis_id)) {
    corstr <- match.arg(corstr)
    return(.dsvert_formal_gee_independence_adapter(
      explicit_arguments = names(match.call())[-1L],
      formula = if (missing(formula)) NULL else formula,
      data = data, family = family, id_col = id_col, order_col = order_col,
      corstr = corstr, verbose = verbose, datasources = datasources,
      analysis_id = formal_analysis_id))
  }
  if (!is.null(fresh_formal_analysis_id)) {
    corstr <- match.arg(corstr)
    return(.dsvert_formal_gee_independence_adapter(
      explicit_arguments = names(match.call())[-1L],
      formula = if (missing(formula)) NULL else formula,
      data = data, family = family, id_col = id_col, order_col = order_col,
      corstr = corstr, verbose = verbose, datasources = datasources,
      analysis_id = fresh_formal_analysis_id,
      selector_name = "fresh_formal_analysis_id"))
  }
  return(.dsvert_block_retired_remote_route("gee", .allow_test = FALSE))
}

.dsvert_gaussian_gee_independence_adapter <- function(
    explicit_arguments, formula, data, family, id_col, order_col, corstr,
    verbose, datasources, dp_analysis_id) {
  if (!identical(family, "gaussian")) {
    stop("dp_analysis_id GEE supports only family='gaussian'", call. = FALSE)
  }
  if (!identical(corstr, "independence")) {
    stop(paste(
      "dp_analysis_id GEE supports only corstr='independence';",
      "cluster working correlations remain unavailable"), call. = FALSE)
  }
  if (!is.null(id_col) || !is.null(order_col)) {
    stop(paste(
      "dp_analysis_id GEE does not accept cluster id_col or order_col;",
      "robust clustered covariance remains unavailable"), call. = FALSE)
  }
  allowed <- c("formula", "data", "family", "corstr", "verbose",
               "datasources", "dp_analysis_id")
  unexpected <- setdiff(explicit_arguments, allowed)
  if (length(unexpected)) {
    stop(paste(
      "dp_analysis_id GEE does not accept legacy controls:",
      paste(sort(unexpected, method = "radix"), collapse = ", ")),
      call. = FALSE)
  }
  fit <- ds.vertGLM(
    formula = formula, data = data, family = "gaussian", verbose = verbose,
    datasources = datasources, dp_analysis_id = dp_analysis_id)
  if (!inherits(fit, "ds.vertDPGaussian") ||
      !identical(fit$family, "gaussian") ||
      !identical(fit$analysis_id, dp_analysis_id) ||
      !is.character(fit$certificate_sha256) ||
      length(fit$certificate_sha256) != 1L || is.na(fit$certificate_sha256) ||
      !grepl("^[0-9a-f]{64}$", fit$certificate_sha256) ||
      !is.numeric(fit$coefficients) || !length(fit$coefficients) ||
      is.null(names(fit$coefficients)) || any(!is.finite(fit$coefficients)) ||
      !isTRUE(fit$source_values_exposed == FALSE) ||
      !isTRUE(fit$intermediate_values_exposed == FALSE) ||
      !identical(fit$additional_server_calls_after_synopsis, 0L) ||
      !identical(fit$additional_privacy_cost, c(epsilon = 0, delta = 0))) {
    stop("signed Gaussian release cannot support independent GEE", call. = FALSE)
  }
  result <- list(
    status = "public_certified_independence_gee_gaussian",
    family = "gaussian",
    corstr = "independence",
    coefficients = fit$coefficients,
    analysis_id = fit$analysis_id,
    dp_analysis_id = dp_analysis_id,
    certificate_sha256 = fit$certificate_sha256,
    robust_covariance = NULL,
    std_errors = NULL,
    cluster_correlation_estimated = FALSE,
    cluster_columns = NULL,
    source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    production_ready = FALSE,
    additional_server_calls_after_synopsis = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    inference = "unavailable_without_protected_cluster_score_and_meat",
    called_via = "ds.vertGEE_dp_analysis_id")
  class(result) <- c("dsvert_dp_gaussian_gee", "ds.vertGEE", "list")
  result
}

.dsvert_formal_gee_independence_adapter <- function(
    explicit_arguments, formula, data, family, id_col, order_col, corstr,
    verbose, datasources, analysis_id,
    selector_name = "formal_analysis_id") {
  selector_name <- match.arg(
    selector_name, c("formal_analysis_id", "fresh_formal_analysis_id"))
  if (!is.character(family) || length(family) != 1L || is.na(family) ||
      !family %in% c("binomial", "poisson")) {
    stop(paste(
      selector_name, "GEE supports only family='binomial' or",
      "family='poisson'"), call. = FALSE)
  }
  if (!is.character(corstr) || length(corstr) != 1L || is.na(corstr) ||
      !identical(corstr, "independence")) {
    stop(paste(
      selector_name, "GEE supports only corstr='independence';",
      "cluster working correlations remain unavailable"), call. = FALSE)
  }
  if (!is.null(id_col) || !is.null(order_col)) {
    stop(paste(
      selector_name, "GEE does not accept cluster id_col or order_col;",
      "robust clustered covariance remains unavailable"), call. = FALSE)
  }
  allowed <- c("formula", "data", "family", "corstr", "verbose",
               "datasources", selector_name)
  unexpected <- setdiff(explicit_arguments, allowed)
  if (length(unexpected)) {
    stop(paste(
      selector_name, "GEE does not accept legacy controls:",
      paste(sort(unexpected, method = "radix"), collapse = ", ")),
      call. = FALSE)
  }
  glm_arguments <- list(
    formula = formula, data = data, family = family, verbose = verbose,
    datasources = datasources)
  glm_arguments[[selector_name]] <- analysis_id
  fit <- do.call(ds.vertGLM, glm_arguments)
  expected_called_via <- if (identical(
    selector_name, "fresh_formal_analysis_id")) {
    "ds.vertGLM_fresh_formal_analysis_id"
  } else "ds.vertGLM_formal_analysis_id"
  if (!inherits(fit, "dsvert_formal_dp_glm") ||
      !identical(fit$family, family) ||
      !is.numeric(fit$coefficients) || !length(fit$coefficients) ||
      is.null(names(fit$coefficients)) || any(!is.finite(fit$coefficients)) ||
      !is.null(fit$covariance) || !is.null(fit$std_errors) ||
      !isTRUE(fit$source_values_exposed == FALSE) ||
      !isTRUE(fit$intermediate_values_exposed == FALSE) ||
      !identical(fit$called_via, expected_called_via)) {
    stop("formal GLM point release cannot support formal independent GEE",
         call. = FALSE)
  }
  result <- list(
    status = "public_certified_independence_gee_point_estimate",
    family = family,
    corstr = "independence",
    coefficients = fit$coefficients,
    artifact_id = fit$artifact_id,
    certificate_sha256 = fit$certificate_sha256,
    formal_analysis_id = if (identical(selector_name, "formal_analysis_id")) {
      fit$formal_analysis_id
    } else NULL,
    fresh_formal_analysis_id = if (identical(
      selector_name, "fresh_formal_analysis_id")) analysis_id else NULL,
    formula_sha256 = fit$formula_sha256,
    robust_covariance = NULL,
    std_errors = NULL,
    cluster_correlation_estimated = FALSE,
    cluster_columns = NULL,
    source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    production_ready = FALSE,
    inference = "unavailable_without_protected_cluster_score_and_meat",
    called_via = paste0("ds.vertGEE_", selector_name))
  class(result) <- c("dsvert_formal_dp_gee", "ds.vertGEE", "list")
  result
}

#' @keywords internal
.ds_gee_extract_lhs <- function(formula) {
  if (inherits(formula, "formula")) {
    as.character(attr(terms(formula), "variables")[[2]])
  } else if (is.character(formula) && grepl("~", formula)) {
    as.character(attr(terms(as.formula(formula)), "variables")[[2]])
  } else stop("cannot extract y from formula", call. = FALSE)
}

#' @keywords internal
.ds_gee_find_server_holding <- function(datasources, server_names,
                                        data_name, var) {
  results <- .dsvert_aggregate_strict(
    datasources,
    call(name = "dsvertColNamesDS", data_name = data_name),
    operation = "vertical column discovery")
  present <- server_names[vapply(
    results, function(value) is.list(value) && var %in% value$columns,
    logical(1L))]
  if (length(present) > 1L) {
    stop("variable '", var,
         "' is present on more than one server; choose an unambiguous vertical partition",
         call. = FALSE)
  }
  if (length(present) == 1L) present[[1L]] else NULL
}

#' @keywords internal
.ds_gee_standardized_parameters <- function(fit, features) {
  beta_orig <- as.numeric(fit$coefficients[features])
  names(beta_orig) <- features
  x_sds <- as.numeric(fit$x_sds[features])
  x_means <- as.numeric(fit$x_means[features])
  names(x_sds) <- names(x_means) <- features

  if (fit$family == "gaussian" && isTRUE(fit$standardize_y) &&
      !is.null(fit$y_sd) && is.finite(fit$y_sd) && fit$y_sd > 0) {
    beta_std <- beta_orig * x_sds / fit$y_sd
    intercept_std <- (as.numeric(fit$coefficients["(Intercept)"]) +
                        sum(beta_orig * x_means) - fit$y_mean) / fit$y_sd
  } else {
    beta_std <- beta_orig * x_sds
    intercept_std <- as.numeric(fit$coefficients["(Intercept)"]) +
                       sum(beta_orig * x_means)
  }
  names(beta_std) <- features
  list(intercept = as.numeric(intercept_std), beta = beta_std)
}

#' @keywords internal
.ds_gee_standardization_jacobian <- function(fit, features) {
  p <- length(features)
  out_names <- c("(Intercept)", features)
  J <- diag(p + 1L)
  dimnames(J) <- list(out_names, out_names)
  x_sds <- as.numeric(fit$x_sds[features])
  x_means <- as.numeric(fit$x_means[features])

  if (fit$family == "gaussian" && isTRUE(fit$standardize_y) &&
      !is.null(fit$y_sd) && is.finite(fit$y_sd) && fit$y_sd > 0) {
    for (jj in seq_len(p)) {
      J[jj + 1L, jj + 1L] <- fit$y_sd / x_sds[jj]
      J[1L, jj + 1L] <- -fit$y_sd * x_means[jj] / x_sds[jj]
    }
    J[1L, 1L] <- fit$y_sd
  } else {
    for (jj in seq_len(p)) {
      J[jj + 1L, jj + 1L] <- 1 / x_sds[jj]
      J[1L, jj + 1L] <- -x_means[jj] / x_sds[jj]
    }
  }
  J
}

#' @keywords internal
.ds_gee_unstandardized_parameters <- function(fit, beta_std, features) {
  beta_std <- as.numeric(beta_std[c("(Intercept)", features)])
  names(beta_std) <- c("(Intercept)", features)
  x_sds <- as.numeric(fit$x_sds[features])
  x_means <- as.numeric(fit$x_means[features])
  names(x_sds) <- names(x_means) <- features

  if (fit$family == "gaussian" && isTRUE(fit$standardize_y) &&
      !is.null(fit$y_sd) && is.finite(fit$y_sd) && fit$y_sd > 0) {
    beta_orig <- beta_std[features] * fit$y_sd / x_sds
    intercept <- beta_std["(Intercept)"] * fit$y_sd + fit$y_mean -
      sum(beta_orig * x_means)
  } else {
    beta_orig <- beta_std[features] / x_sds
    intercept <- beta_std["(Intercept)"] -
      sum(beta_orig * x_means)
  }
  out <- c("(Intercept)" = as.numeric(intercept), beta_orig)
  names(out) <- c("(Intercept)", features)
  out
}

#' @keywords internal
.ds_gee_solve_exchangeable_from_cluster_stats <- function(
    sx, sy, sxy, xx, syy, beta_start, max_iter = 100L, tol = 1e-6) {
  q <- ncol(sx)
  n_clusters <- nrow(sx)
  cluster_sizes <- as.numeric(sx[, 1L])
  if (q < 1L || n_clusters < 1L) {
    stop("empty cluster statistics", call. = FALSE)
  }
  if (any(!is.finite(cluster_sizes)) || any(cluster_sizes < 2)) {
    stop("invalid exchangeable cluster sizes", call. = FALSE)
  }
  n_obs <- sum(cluster_sizes)
  beta <- as.numeric(beta_start)
  if (length(beta) != q || any(!is.finite(beta))) {
    beta <- rep(0, q)
  }
  alpha <- 0
  lower_alpha <- -1 / (max(cluster_sizes) - 1) + 1e-6
  upper_alpha <- 0.95

  safe_solve <- function(A, b = NULL) {
    .dsvert_solve_identifiable(
      A, b,
      context = "The exchangeable Gaussian GEE/GLS system",
      reason = "singular_exchangeable_gee_information",
      symmetric = TRUE)
  }

  residual_summaries <- function(beta) {
    rss <- numeric(n_clusters)
    rsum <- numeric(n_clusters)
    for (cc in seq_len(n_clusters)) {
      xxy <- xx[cc, , , drop = FALSE][1L, , ]
      rss[cc] <- syy[cc] - 2 * sum(beta * sxy[cc, ]) +
        drop(t(beta) %*% xxy %*% beta)
      rsum[cc] <- sy[cc] - sum(sx[cc, ] * beta)
    }
    rss <- pmax(rss, 0)
    list(rss = rss, rsum = rsum)
  }

  update_alpha <- function(beta) {
    rs <- residual_summaries(beta)
    phi <- sum(rs$rss) / max(n_obs - q, 1)
    if (!is.finite(phi) || phi <= 0) return(0)
    pair_num <- sum((rs$rsum^2 - rs$rss) / 2)
    pair_den <- sum(cluster_sizes * (cluster_sizes - 1) / 2) * phi
    if (!is.finite(pair_den) || pair_den <= 0) return(0)
    alpha_new <- pair_num / pair_den
    min(max(alpha_new, lower_alpha), upper_alpha)
  }

  fixed_alpha_fit <- function(alpha) {
    A <- matrix(0, q, q)
    rhs <- numeric(q)
    for (cc in seq_len(n_clusters)) {
      m <- cluster_sizes[cc]
      a <- 1 / (1 - alpha)
      b <- alpha / ((1 - alpha) * (1 - alpha + alpha * m))
      A <- A + a * xx[cc, , ] - b * tcrossprod(sx[cc, ])
      rhs <- rhs + a * sxy[cc, ] - b * sx[cc, ] * sy[cc]
    }
    list(beta = as.numeric(safe_solve(A, rhs)), A = A)
  }

  converged <- FALSE
  iter <- 0L
  for (iter in seq_len(as.integer(max_iter))) {
    beta_old <- beta
    alpha_old <- alpha
    alpha <- update_alpha(beta)
    fa <- fixed_alpha_fit(alpha)
    beta <- fa$beta
    if (max(abs(beta - beta_old), abs(alpha - alpha_old)) <= tol) {
      converged <- TRUE
      break
    }
  }
  fa <- fixed_alpha_fit(alpha)
  beta <- fa$beta
  A <- (fa$A + t(fa$A)) / 2
  bread <- safe_solve(A)

  rs <- residual_summaries(beta)
  phi_num <- 0
  scores <- matrix(0, nrow = n_clusters, ncol = q)
  for (cc in seq_len(n_clusters)) {
    m <- cluster_sizes[cc]
    a <- 1 / (1 - alpha)
    b <- alpha / ((1 - alpha) * (1 - alpha + alpha * m))
    rsum <- rs$rsum[cc]
    phi_num <- phi_num + a * rs$rss[cc] - b * rsum^2
    scores[cc, ] <- a * (sxy[cc, ] - xx[cc, , ] %*% beta) -
      b * sx[cc, ] * rsum
  }
  phi <- phi_num / max(n_obs - q, 1)
  if (!is.finite(phi) || phi <= 0) {
    phi <- sum(rs$rss) / max(n_obs - q, 1)
  }
  model_cov <- as.numeric(phi) * bread
  robust_cov <- bread %*% crossprod(scores) %*% bread
  model_cov <- (model_cov + t(model_cov)) / 2
  robust_cov <- (robust_cov + t(robust_cov)) / 2

  list(beta = beta, alpha = alpha, phi = as.numeric(phi),
       model_cov = model_cov, robust_cov = robust_cov,
       iterations = as.integer(iter), converged = converged,
       scores = scores)
}

#' @keywords internal
.ds_gee_residual_quad_from_stats <- function(stats, beta) {
  stats$yy - 2 * sum(beta * stats$xy) +
    drop(t(beta) %*% stats$xx %*% beta)
}

#' @keywords internal
.ds_gee_adjacent_residual_cross <- function(adj, beta) {
  adj$yy - sum(beta * adj$xy_forward) - sum(beta * adj$xy_backward) +
    drop(t(beta) %*% adj$xx %*% beta)
}

#' @keywords internal
.ds_gee_solve_ar1_from_stats <- function(
    total, interior, nonlast, adj, beta_start, max_iter = 100L,
    tol = 1e-6, lag_stats = NULL) {
  q <- ncol(total$xx)
  n_obs <- as.numeric(total$xx[1L, 1L])
  if (q < 1L || !is.finite(n_obs) || n_obs < q + 1L) {
    stop("invalid AR1 sufficient statistics", call. = FALSE)
  }
  beta <- as.numeric(beta_start)
  if (length(beta) != q || any(!is.finite(beta))) beta <- rep(0, q)
  rho <- 0
  lower_rho <- -0.95
  upper_rho <- 0.95

  safe_solve <- function(A, b = NULL) {
    .dsvert_solve_identifiable(
      A, b,
      context = "The Gaussian AR1 GEE/GLS system",
      reason = "singular_ar1_gee_information",
      symmetric = TRUE)
  }

  if (is.null(lag_stats) || length(lag_stats) == 0L) {
    lag_stats <- list(adj)
  }

  update_rho <- function(beta_current) {
    if (length(lag_stats) == 1L) {
      num <- .ds_gee_adjacent_residual_cross(adj, beta_current)
      den <- .ds_gee_residual_quad_from_stats(nonlast, beta_current)
      if (!is.finite(num) || !is.finite(den) || den <= 0) return(0)
      return(min(max(num / den, lower_rho), upper_rho))
    }
    phi <- .ds_gee_residual_quad_from_stats(total, beta_current) /
      max(n_obs - q, 1)
    if (!is.finite(phi) || phi <= 0) return(0)
    counts <- vapply(lag_stats, function(x) as.numeric(x$count %||%
                                                        x$xx[1L, 1L]),
                     numeric(1L))
    rho_targets <- vapply(lag_stats, function(x) {
      .ds_gee_adjacent_residual_cross(x, beta_current) /
        max(as.numeric(x$count %||% x$xx[1L, 1L]) * phi, .Machine$double.eps)
    }, numeric(1L))
    keep <- is.finite(counts) & counts > 0 & is.finite(rho_targets)
    if (!any(keep)) return(0)
    counts <- counts[keep]
    rho_targets <- rho_targets[keep]
    lags <- seq_along(lag_stats)[keep]
    objective <- function(rho_current) {
      sum(counts * (rho_targets - rho_current^lags)^2)
    }
    out <- tryCatch(stats::optimize(objective,
                                    interval = c(lower_rho, upper_rho))$minimum,
                    error = function(e) 0)
    if (!is.finite(out)) out <- 0
    min(max(out, lower_rho), upper_rho)
  }

  fixed_rho_fit <- function(rho_current) {
    cfac <- 1 / (1 - rho_current^2)
    A <- cfac * (total$xx + rho_current^2 * interior$xx -
                   rho_current * (adj$xx + t(adj$xx)))
    rhs <- cfac * (total$xy + rho_current^2 * interior$xy -
                     rho_current * (adj$xy_forward + adj$xy_backward))
    list(beta = as.numeric(safe_solve(A, rhs)), A = A)
  }

  converged <- FALSE
  iter <- 0L
  for (iter in seq_len(as.integer(max_iter))) {
    beta_old <- beta
    rho_old <- rho
    rho <- update_rho(beta)
    fr <- fixed_rho_fit(rho)
    beta <- fr$beta
    if (max(abs(beta - beta_old), abs(rho - rho_old)) <= tol) {
      converged <- TRUE
      break
    }
  }
  fr <- fixed_rho_fit(rho)
  beta <- fr$beta
  A <- (fr$A + t(fr$A)) / 2
  bread <- safe_solve(A)

  rss_total <- .ds_gee_residual_quad_from_stats(total, beta)
  phi <- rss_total / max(n_obs - q, 1)
  if (!is.finite(phi) || phi <= 0) phi <- 1
  model_cov <- as.numeric(phi) * bread
  model_cov <- (model_cov + t(model_cov)) / 2
  list(beta = beta, rho = rho, phi = as.numeric(phi), A = A,
       bread = bread, model_cov = model_cov,
       iterations = as.integer(iter), converged = converged)
}

#' @keywords internal
print.ds.vertGEE <- function(x, ...) {
  if (inherits(x, "dsvert_formal_dp_gee")) {
    cat("dsVert formal independent-working GEE point estimate\n")
    cat("  Family :", x$family, "\n")
    print(x$coefficients)
    cat("  No robust covariance, standard errors or inference are released.\n")
    return(invisible(x))
  }
  cat("dsVert GEE (", x$corstr, " working correlation)\n", sep = "")
  cat(sprintf("  Family: %s   N = %d\n", x$family, x$n_obs))
  if (!is.null(x$estimand)) {
    cat(sprintf("  Estimand: %s (lambda = %.4g)\n",
                x$estimand, x$lambda %||% 0))
  }
  cat("\nCoefficients:\n")
  df <- data.frame(
    Estimate = x$coefficients,
    `SE (model)` = x$model_se,
    `SE (robust)` = x$robust_se,
    z = x$coefficients / x$robust_se,
    check.names = FALSE)
  df$p <- 2 * stats::pnorm(-abs(df$z))
  print(round(df, 5L))
  invisible(x)
}
