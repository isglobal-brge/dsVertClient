# DP-aware categorical inference is deliberately client-only.  The statistic
# is calibrated against the already released sticky biomedical vector; it
# never asks a peer for another table, margin, seed or mechanism draw.
.DSVERT_DP_CHISQ_BOOTSTRAP_VERSION <-
  "dsvert-dp-categorical-parametric-bootstrap-v1"
.DSVERT_DP_CHISQ_EXACT_BOOTSTRAP_VERSION <-
  "dsvert-dp-categorical-parametric-bootstrap-one-draw-v2"
.DSVERT_DP_FISHER_BOOTSTRAP_VERSION <-
  "dsvert-dp-2x2-conditional-hypergeometric-bootstrap-v1"
.DSVERT_DP_FISHER_EXACT_BOOTSTRAP_VERSION <-
  "dsvert-dp-2x2-conditional-hypergeometric-bootstrap-one-draw-v2"

#' @title DP-aware independence test on a 2-way contingency table
#' @description Calibrate a categorical independence statistic against the
#'   sampling distribution and the privacy mechanism. The input is either an
#'   already released `ds.vertDPContingency` object (zero DSI calls) or the
#'   original data/variable identifiers, in which case the one immutable
#'   signed Synopsis projection is obtained through `ds.vertDPContingency`. Ordinary
#'   chi-square or Fisher reference laws are never applied to the noisy cells.
#'
#' @param data_name A released `ds.vertDPContingency` object, the name of the
#'   protected data frame, or a reusable `ds.vertFederation`.
#' @param var1,var2 Row and column variables when `data_name` is a character
#'   string. They may be omitted for an already released object.
#' @param server Optional owner-peer assertion forwarded to
#'   \code{ds.vertDPContingency}. No separate column-discovery request is
#'   performed. It must be omitted for an existing release.
#' @param correct Logical. If TRUE (default for 2x2 tables), apply Yates'
#'   continuity correction to the statistic. Its null law is still obtained
#'   by the same DP-aware bootstrap rather than a chi-square approximation.
#' @param datasources DataSHIELD connections.
#' @param simulations Number of parametric-bootstrap replicates. The default
#'   gives p-value resolution `1 / 10000`.
#' @param mc_confidence Confidence level for the reported exact binomial
#'   interval around the Monte Carlo p-value.
#'
#' @return An object of class \code{ds.vertChisq} with elements
#'   \code{statistic}, \code{df}, \code{p_value}, \code{observed},
#'   \code{expected}, \code{residuals}, \code{n}, \code{correct}, and the
#'   bootstrap/mechanism uncertainty contract.
#'
#' @details The test follows the Monte Carlo independence-testing principle of
#'   Gaboardi et al. (2016), adapted to dsVert's signed mechanism. It estimates
#'   the latent contributing sample size and null margins from the released
#'   table, simulates multinomial tables under the fitted independence model,
#'   adds the signed exact-GC plan's dyadic one-draw reference or two
#'   independent peer discrete-Laplace draws per cell, applies the same
#'   common-lattice clamp,
#'   and repeats the nuisance fit for every
#'   replicate. This is a parametric-bootstrap test: it is asymptotically
#'   calibrated under positive cell probabilities, but is not a finite-sample
#'   exact conditional test. Tables whose fitted expected count is below five
#'   return a structured non-tested result rather than an unreliable p-value.
#'
#'   The production sampler is a finite binary-geometric approximation to its
#'   signed-plan dyadic one-draw or two-peer-convolution reference law. Its
#'   total-variation certificate is propagated into a separate calibration
#'   interval; the exact-GC branch also charges an outward-rounded coupling
#'   bound for R's numeric stop probability. Bootstrap randomness
#'   is deterministically derived from public release commitments solely for
#'   reproducible post-processing; it is not DP mechanism randomness and adds
#'   no privacy cost.
#'
#' @importFrom stats pchisq
#' @export
ds.vertChisq <- function(data_name, var1 = NULL, var2 = NULL, server = NULL,
                         correct = TRUE, datasources = NULL,
                         simulations = 9999L, mc_confidence = 0.95) {
  if (!is.logical(correct) || length(correct) != 1L || is.na(correct)) {
    stop("correct must be one non-missing logical value", call. = FALSE)
  }
  .dsvert_dp_chisq_validate_simulation_inputs(simulations, mc_confidence)

  if (inherits(data_name, "ds.vertDPContingency")) {
    release <- data_name
    if (!is.null(datasources) || !is.null(server)) {
      stop("server and datasources must be omitted for an existing DP release",
           call. = FALSE)
    }
    assertions <- list(row = var1, column = var2)
    expected <- list(row = release$row_var, column = release$col_var)
    for (name in names(assertions)) {
      value <- assertions[[name]]
      if (!is.null(value) && (!is.character(value) || length(value) != 1L ||
          is.na(value) || !identical(value, expected[[name]]))) {
        stop(name, " variable assertion does not match the DP release",
             call. = FALSE)
      }
    }
  } else {
    if (inherits(data_name, "ds.vertFederation")) {
      resolved <- .dsvert_federation_argument(data_name, datasources)
      data_name <- resolved$value
      datasources <- resolved$datasources
    }
    for (value in list(data_name = data_name, var1 = var1, var2 = var2)) {
      if (!is.character(value) || length(value) != 1L || is.na(value) ||
          !nzchar(value)) {
        stop("data_name, var1 and var2 must be non-empty strings",
             call. = FALSE)
      }
    }
    release <- ds.vertDPContingency(
      data_name = data_name, row_var = var1, col_var = var2,
      server = server, datasources = datasources)
  }
  .dsvert_dp_chisq_from_release(
    release, correct = correct, simulations = simulations,
    mc_confidence = mc_confidence)
}

.dsvert_dp_chisq_validate_simulation_inputs <- function(
    simulations, mc_confidence) {
  if (!is.numeric(simulations) || length(simulations) != 1L ||
      is.na(simulations) || !is.finite(simulations) || simulations < 1L ||
      simulations > .Machine$integer.max ||
      simulations != floor(simulations)) {
    stop("simulations must be one positive integer", call. = FALSE)
  }
  if (!is.numeric(mc_confidence) || length(mc_confidence) != 1L ||
      is.na(mc_confidence) || !is.finite(mc_confidence) ||
      mc_confidence <= 0 || mc_confidence >= 1) {
    stop("mc_confidence must lie strictly between zero and one",
         call. = FALSE)
  }
  invisible(TRUE)
}

.dsvert_dp_chisq_seed <- function(x, correct, calibration) {
  release_identity <- if (.dsvert_vector_hex(x$artifact_key)) {
    x$artifact_key
  } else x$capsule_id
  fields <- c(
    calibration, release_identity,
    x$final_vector_root, x$coordinate_order_sha256,
    x$row_var, x$col_var, if (isTRUE(correct)) "yates" else "pearson")
  if (anyNA(fields) || any(!nzchar(fields))) {
    stop("The DP release lacks a public reproducibility commitment",
         call. = FALSE)
  }
  hash <- digest::digest(
    paste(enc2utf8(fields), collapse = "|"), algo = "sha256",
    serialize = FALSE)
  as.integer(strtoi(substr(hash, 1L, 7L), base = 16L)) + 1L
}

.dsvert_dp_table_source_binding <- function(x) {
  identity <- if (.dsvert_vector_hex(x$artifact_key)) {
    list(artifact_key = x$artifact_key)
  } else {
    list(capsule_id = x$capsule_id)
  }
  c(identity, list(
    final_vector_root = x$final_vector_root,
    manifest_sha256 = x$manifest_sha256,
    coordinate_order_sha256 = x$coordinate_order_sha256))
}

.dsvert_dp_chisq_with_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv,
                                inherits = FALSE) else NULL
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv,
                      inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  force(code)
}

.dsvert_dp_chisq_project_simplex <- function(values, total) {
  values <- as.numeric(values)
  if (!length(values) || anyNA(values) || any(!is.finite(values)) ||
      any(values < 0) || !is.numeric(total) || length(total) != 1L ||
      is.na(total) || !is.finite(total) || total <= 0) return(NULL)
  ordered <- sort(values, decreasing = TRUE, method = "radix")
  candidate <- ordered - (cumsum(ordered) - total) / seq_along(ordered)
  positive <- which(candidate > 0)
  if (!length(positive)) return(NULL)
  rho <- max(positive)
  threshold <- (sum(ordered[seq_len(rho)]) - total) / rho
  projected <- pmax(values - threshold, 0)
  drift <- total - sum(projected)
  if (abs(drift) > 0) projected[[which.max(projected)]] <-
    projected[[which.max(projected)]] + drift
  if (any(!is.finite(projected)) || any(projected < 0)) return(NULL)
  projected
}

.dsvert_dp_chisq_fit <- function(table, capacity, minimum_expected = 5) {
  dims <- dim(table)
  total_dp <- sum(table)
  if (length(dims) != 2L || any(dims < 2L) || !is.finite(total_dp) ||
      total_dp <= 0) {
    return(list(ok = FALSE, status = "not_tested_degenerate_dp_table"))
  }
  n_estimate <- min(capacity, floor(total_dp + 0.5))
  if (!is.finite(n_estimate) || n_estimate < 1) {
    return(list(ok = FALSE, status = "not_tested_degenerate_dp_table"))
  }
  latent <- .dsvert_dp_chisq_project_simplex(table, n_estimate)
  if (is.null(latent)) {
    return(list(ok = FALSE, status = "not_tested_degenerate_dp_table"))
  }
  latent <- matrix(latent, nrow = dims[[1L]], ncol = dims[[2L]])
  row_probability <- rowSums(latent) / n_estimate
  col_probability <- colSums(latent) / n_estimate
  if (any(row_probability <= 0) || any(col_probability <= 0)) {
    return(list(ok = FALSE, status = "not_tested_degenerate_dp_table"))
  }
  expected <- n_estimate * outer(row_probability, col_probability)
  if (any(!is.finite(expected)) || any(expected < minimum_expected)) {
    return(list(
      ok = FALSE, status = "not_tested_minimum_expected_count",
      n = n_estimate, latent = latent, expected = expected,
      row_probability = row_probability,
      col_probability = col_probability))
  }
  list(
    ok = TRUE, status = "ok", n = n_estimate, latent = latent,
    expected = expected, row_probability = row_probability,
    col_probability = col_probability)
}

.dsvert_dp_chisq_laplace_mode <- function(x) {
  profile <- .dsvert_dp_table_vector_profile(x)
  if (!is.list(profile) || isTRUE(profile$gaussian) ||
      !identical(x$mechanism, profile$release_mechanism)) return(NULL)
  if (isTRUE(profile$exact_gc)) {
    "signed_one_draw_discrete_laplace"
  } else {
    "signed_two_peer_discrete_laplace"
  }
}

.dsvert_dp_chisq_exact_stop_probability <- function(plan) {
  if (!is.list(plan)) {
    stop("The exact-GC release lacks its signed finite-sampler plan",
         call. = FALSE)
  }
  if (!identical(suppressWarnings(as.numeric(plan$stop_bits)), 128)) {
    stop("The signed exact-GC stop probability is not representable in R",
         call. = FALSE)
  }
  interval <- .dsvert_dp_vector_dyadic_tail_context(plan)
  value <- interval$q
  if (!is.finite(value) || value <= 0 || value >= 1 ||
      value < interval$q_lower || value > interval$q_upper) {
    stop("The signed exact-GC stop probability is not representable in R",
         call. = FALSE)
  }
  round_up <- function(x) {
    for (index in seq_len(8L)) x <- .dsvert_dp_vector_next_up(x)
    x
  }
  error_upper <- max(
    round_up(value - interval$q_lower),
    round_up(interval$q_upper - value))
  # Shared Bernoulli uniforms couple Geom(q) and Geom(value) with mismatch
  # probability at most |q-value|/max(q,value) <= error_upper/value. A signed
  # two-sided draw uses two geometricals, followed by a product over every
  # released coordinate.
  one_geometric_tv <- round_up(error_upper / value)
  numeric_vector_tv <- round_up(
    2 * as.numeric(plan$total_coordinate_count) * one_geometric_tv)
  if (!is.finite(numeric_vector_tv) || numeric_vector_tv < 0 ||
      numeric_vector_tv >= 1) {
    stop("The signed exact-GC stop probability cannot certify R simulation",
         call. = FALSE)
  }
  list(
    value = value,
    interval = c(lower = interval$q_lower, upper = interval$q_upper),
    source = "signed_exact_gc_dyadic_plan",
    numeric_vector_tv_upper_bound = numeric_vector_tv)
}

.dsvert_dp_chisq_noise_contract <- function(x) {
  mode <- .dsvert_dp_chisq_laplace_mode(x)
  draw_count <- if (identical(mode, "signed_one_draw_discrete_laplace")) {
    1L
  } else 2L
  scale <- as.numeric(x$output_lattice_scale)
  sensitivity_steps <- as.numeric(x$l1_sensitivity) * scale
  rate <- as.numeric(x$epsilon) / sensitivity_steps
  stop_contract <- if (draw_count == 1L) {
    .dsvert_dp_chisq_exact_stop_probability(x$mechanism_plan)
  } else {
    value <- -expm1(-rate)
    list(
      value = value, interval = c(lower = value, upper = value),
      source = "declared_epsilon_exponential_law",
      numeric_vector_tv_upper_bound = 0)
  }
  stop_probability <- stop_contract$value
  continuation_probability <- 1 - stop_probability
  if (!is.finite(scale) || scale <= 0 || !is.finite(sensitivity_steps) ||
      sensitivity_steps <= 0 || !is.finite(stop_probability) ||
      stop_probability <= 0 || stop_probability >= 1) {
    stop("The signed DP noise calibration is not representable in R",
         call. = FALSE)
  }
  # Each complete draw is the difference of two iid geometricals. Exact GC
  # emits one joint draw; the convolution backend emits two peer draws.
  variance_stop_probability <- if (draw_count == 1L) {
    stop_contract$interval[["lower"]]
  } else stop_probability
  variance <- if (draw_count == 1L) {
    outward <- function(value, upward) {
      for (index in seq_len(8L)) {
        value <- if (isTRUE(upward)) {
          .dsvert_dp_vector_next_up(value)
        } else {
          .dsvert_dp_vector_next_down(value)
        }
      }
      value
    }
    continuation_upper <- outward(1 - variance_stop_probability, TRUE)
    numerator_upper <- outward(2 * continuation_upper, TRUE)
    probability_squared_lower <- outward(
      variance_stop_probability * variance_stop_probability, FALSE)
    scale_squared_lower <- outward(scale * scale, FALSE)
    if (!is.finite(probability_squared_lower) ||
        probability_squared_lower <= 0 ||
        !is.finite(scale_squared_lower) || scale_squared_lower <= 0) {
      stop("The signed DP noise variance is not representable in R",
           call. = FALSE)
    }
    steps_upper <- outward(
      numerator_upper / probability_squared_lower, TRUE)
    outward(steps_upper / scale_squared_lower, TRUE)
  } else {
    variance_steps <- 2 * draw_count *
      (1 - variance_stop_probability) / variance_stop_probability^2
    variance_steps / scale^2
  }
  vector_tv <- if (draw_count == 1L) {
    .dsvert_dp_vector_sampler_tv_upper(x$mechanism_plan, TRUE)
  } else {
    min(1, 2 * .dsvert_dp_vector_fraction(x$implementation_delta))
  }
  numeric_tv <- stop_contract$numeric_vector_tv_upper_bound
  reference_tv <- if (draw_count == 1L) {
    value <- vector_tv + numeric_tv
    for (index in seq_len(8L)) {
      value <- .dsvert_dp_vector_next_up(value)
    }
    min(1, value)
  } else vector_tv
  calibration_tv <- if (draw_count == 1L) {
    value <- 2 * reference_tv
    for (index in seq_len(8L)) {
      value <- .dsvert_dp_vector_next_up(value)
    }
    min(1, value)
  } else min(1, 2 * vector_tv)
  list(
    scale = scale, sensitivity_steps = sensitivity_steps,
    stop_probability = stop_probability,
    stop_probability_interval = stop_contract$interval,
    stop_probability_source = stop_contract$source,
    continuation_probability = continuation_probability,
    variance_stop_probability = variance_stop_probability,
    variance_is_upper_bound = draw_count == 1L,
    mode = mode, draw_count = draw_count,
    chisq_calibration = if (draw_count == 1L) {
      .DSVERT_DP_CHISQ_EXACT_BOOTSTRAP_VERSION
    } else .DSVERT_DP_CHISQ_BOOTSTRAP_VERSION,
    fisher_calibration = if (draw_count == 1L) {
      .DSVERT_DP_FISHER_EXACT_BOOTSTRAP_VERSION
    } else .DSVERT_DP_FISHER_BOOTSTRAP_VERSION,
    variance = variance,
    vector_tv_upper_bound = vector_tv,
    numeric_reference_tv_upper_bound = numeric_tv,
    one_transfer_tv_upper_bound = reference_tv,
    calibration_tv_upper_bound = calibration_tv)
}

.dsvert_dp_chisq_sample_one_laplace <- function(count, noise) {
  zero_probability <- noise$stop_probability /
    (2 - noise$stop_probability)
  nonzero <- stats::runif(count) >= zero_probability
  result <- numeric(count)
  number <- sum(nonzero)
  if (number > 0L) {
    magnitude <- 1 + stats::rgeom(
      number, prob = noise$stop_probability)
    sign <- ifelse(stats::runif(number) < 0.5, -1, 1)
    result[nonzero] <- sign * magnitude
  }
  if (any(!is.finite(result))) {
    stop("The DP bootstrap noise is outside R's numeric simulation domain",
         call. = FALSE)
  }
  result
}

.dsvert_dp_chisq_sample_noise <- function(count, noise) {
  draw_count <- if (is.null(noise$draw_count)) 2L else noise$draw_count
  if (!identical(draw_count, 1L) && !identical(draw_count, 2L)) {
    stop("Invalid signed DP noise draw count", call. = FALSE)
  }
  result <- .dsvert_dp_chisq_sample_one_laplace(count, noise)
  if (draw_count == 2L) {
    result <- result + .dsvert_dp_chisq_sample_one_laplace(count, noise)
  }
  result / noise$scale
}

.dsvert_dp_chisq_apply_release_clamp <- function(values, noise, capacity) {
  pmin(pmax(values + noise, 0), capacity)
}

.dsvert_dp_chisq_multinomial <- function(replicates, size, probability) {
  cells <- length(probability)
  result <- matrix(0, nrow = cells, ncol = replicates)
  remaining_size <- rep(size, replicates)
  remaining_probability <- 1
  if (cells > 1L) {
    for (cell in seq_len(cells - 1L)) {
      conditional <- min(1, max(0,
        probability[[cell]] / remaining_probability))
      draw <- stats::rbinom(replicates, remaining_size, conditional)
      result[cell, ] <- draw
      remaining_size <- remaining_size - draw
      remaining_probability <- remaining_probability - probability[[cell]]
    }
  }
  result[cells, ] <- remaining_size
  result
}

.dsvert_dp_chisq_statistic <- function(table, fit, correct) {
  expected <- fit$expected
  difference <- abs(table - expected)
  corrected <- isTRUE(correct) && identical(dim(table), c(2L, 2L))
  if (corrected) difference <- pmax(0, difference - 0.5)
  # This is Q_D from the Monte Carlo independence test: the privacy noise is
  # represented in the bootstrap reference distribution, not hidden in an
  # ordinary chi-square tail approximation.
  denominator <- expected
  list(
    statistic = sum(difference^2 / denominator),
    residuals = (table - expected) / sqrt(denominator),
    correct = corrected)
}

.dsvert_dp_chisq_mc_interval <- function(exceedances, simulations,
                                         confidence) {
  alpha <- 1 - confidence
  lower <- if (exceedances == 0L) 0 else stats::qbeta(
    alpha / 2, exceedances, simulations - exceedances + 1)
  upper <- if (exceedances == simulations) 1 else stats::qbeta(
    1 - alpha / 2, exceedances + 1, simulations - exceedances)
  c(lower = (1 + simulations * lower) / (simulations + 1),
    upper = (1 + simulations * upper) / (simulations + 1))
}

.dsvert_dp_chisq_unavailable <- function(
    x, fit, correct, simulations, mc_confidence, noise) {
  expected <- if (is.null(fit$expected)) {
    matrix(NA_real_, nrow = nrow(x$table), ncol = ncol(x$table))
  } else {
    fit$expected
  }
  dimnames(expected) <- dimnames(x$table)
  result <- list(
    status = fit$status, decision_available = FALSE,
    statistic = NA_real_, df = as.integer((nrow(x$table) - 1L) *
                                            (ncol(x$table) - 1L)),
    p_value = NA_real_, observed = x$table, expected = expected,
    residuals = matrix(NA_real_, nrow = nrow(x$table),
                       ncol = ncol(x$table), dimnames = dimnames(x$table)),
    n = if (is.null(fit$n)) NA_real_ else fit$n,
    n_na = NA_integer_,
    correct = isTRUE(correct) && identical(dim(x$table), c(2L, 2L)),
    server = x$server, var1 = x$row_var, var2 = x$col_var,
    simulations = as.integer(simulations), mc_confidence = mc_confidence,
    minimum_expected_count = 5,
    calibration = noise$chisq_calibration,
    reference_mechanism = if (noise$draw_count == 1L) paste(
      "signed-plan dyadic one joint complete discrete-Laplace draw per cell;",
      "finite-sampler and R-parameter TV accounted separately;",
      "exact public clamp")
    else paste(
      "ideal two-independent complete discrete-Laplace draws per cell;",
      "signed finite-sampler TV accounted separately; exact public clamp"),
    finite_sampler_bit_exact_reference = FALSE,
    mechanism_noise_variance_per_cell = noise$variance,
    mechanism_noise_variance_per_cell_is_upper_bound =
      noise$variance_is_upper_bound,
    mechanism_vector_tv_upper_bound = noise$vector_tv_upper_bound,
    mechanism_numeric_reference_tv_upper_bound =
      noise$numeric_reference_tv_upper_bound,
    mechanism_one_transfer_tv_upper_bound =
      noise$one_transfer_tv_upper_bound,
    mechanism_reference_tv_upper_bound =
      noise$calibration_tv_upper_bound,
    additional_server_calls = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    source_dp_release = x,
    inferential_contract = paste(
      "No p-value is returned because the fitted independence model does",
      "not meet its predeclared identifiability/expected-count contract"))
  class(result) <- c("ds.vertChisq", "list")
  result
}

.dsvert_dp_chisq_mechanism_dispatch <- function(x) {
  profile <- .dsvert_dp_table_vector_profile(x)
  mode <- .dsvert_dp_chisq_laplace_mode(x)
  if (!is.null(mode)) return(mode)
  gaussian <- is.list(profile) && isTRUE(profile$gaussian)
  detail <- if (gaussian) {
    paste(
      "The TV-bounded dyadic discrete-Gaussian contingency release is valid for bounded",
      "estimation, but ds.vertChisq does not yet carry a bit-exact public",
      "reference sampler and calibration-error certificate for that mechanism.")
  } else {
    paste(
      "ds.vertChisq currently supports only signed Ring128 one-draw or",
      "two-peer discrete-Laplace contingency mechanisms.")
  }
  .dsvert_stop_numeric(
    "dsvert_dp_chisq_mechanism_not_certified", detail,
    reason = if (gaussian) {
      "gaussian_parametric_reference_not_certified"
    } else {
      "unsupported_dp_chisq_mechanism"
    })
}

.dsvert_dp_chisq_from_release <- function(
    x, correct = TRUE, simulations = 9999L, mc_confidence = 0.95) {
  .dsvert_dp_chisq_validate_simulation_inputs(simulations, mc_confidence)
  x <- .dsvert_dp_table_contract(x)
  .dsvert_dp_chisq_mechanism_dispatch(x)
  noise <- .dsvert_dp_chisq_noise_contract(x)
  fit <- .dsvert_dp_chisq_fit(x$table, x$coordinate_maximum)
  if (!isTRUE(fit$ok)) {
    return(.dsvert_dp_chisq_unavailable(
      x, fit, correct, simulations, mc_confidence, noise))
  }
  observed_stat <- .dsvert_dp_chisq_statistic(
    x$table, fit, correct)
  probability <- as.numeric(outer(
    fit$row_probability, fit$col_probability))
  seed <- .dsvert_dp_chisq_seed(
    x, observed_stat$correct, noise$chisq_calibration)
  cells <- length(x$table)
  batch_limit <- max(1L, floor(1000000 / cells))

  simulation <- .dsvert_dp_chisq_with_seed(seed, {
    exceedances <- 0L
    invalid <- 0L
    completed <- 0L
    while (completed < simulations) {
      batch <- as.integer(min(batch_limit, simulations - completed))
      latent <- .dsvert_dp_chisq_multinomial(
        batch, fit$n, probability)
      mechanism_noise <- matrix(
        .dsvert_dp_chisq_sample_noise(cells * batch, noise),
        nrow = cells, ncol = batch)
      noisy <- .dsvert_dp_chisq_apply_release_clamp(
        latent, mechanism_noise, x$coordinate_maximum)
      for (index in seq_len(batch)) {
        table <- matrix(noisy[, index], nrow = nrow(x$table),
                        ncol = ncol(x$table))
        replicate_fit <- .dsvert_dp_chisq_fit(
          table, x$coordinate_maximum)
        if (!isTRUE(replicate_fit$ok)) {
          # MCIndep fails to reject if a bootstrap nuisance fit is invalid.
          # Treating it as +Inf is the equivalent conservative tail action.
          invalid <- invalid + 1L
          exceedances <- exceedances + 1L
        } else {
          statistic <- .dsvert_dp_chisq_statistic(
            table, replicate_fit, observed_stat$correct)$statistic
          tolerance <- 64 * .Machine$double.eps *
            max(1, abs(statistic), abs(observed_stat$statistic))
          if (statistic >= observed_stat$statistic - tolerance) {
            exceedances <- exceedances + 1L
          }
        }
      }
      completed <- completed + batch
    }
    list(exceedances = exceedances, invalid = invalid)
  })
  simulations <- as.integer(simulations)
  p_value <- (simulation$exceedances + 1) / (simulations + 1)
  mc_interval <- .dsvert_dp_chisq_mc_interval(
    simulation$exceedances, simulations, mc_confidence)
  calibration_interval <- c(
    lower = max(0, mc_interval[["lower"]] -
                  noise$calibration_tv_upper_bound),
    upper = min(1, mc_interval[["upper"]] +
                  noise$calibration_tv_upper_bound))
  dimnames(fit$expected) <- dimnames(x$table)
  dimnames(observed_stat$residuals) <- dimnames(x$table)
  result <- list(
    status = "ok", decision_available = TRUE,
    statistic = observed_stat$statistic,
    statistic_name = paste(
      "Pearson distance with refitted DP-aware parametric-bootstrap null"),
    df = as.integer((nrow(x$table) - 1L) * (ncol(x$table) - 1L)),
    p_value = p_value,
    p_value_mc_se = sqrt(p_value * (1 - p_value) /
                           (simulations + 1)),
    p_value_mc_interval = mc_interval,
    p_value_calibration_interval = calibration_interval,
    observed = x$table, expected = fit$expected,
    residuals = observed_stat$residuals,
    n = fit$n, n_na = NA_integer_, correct = observed_stat$correct,
    server = x$server, var1 = x$row_var, var2 = x$col_var,
    simulations = simulations, exceedances = simulation$exceedances,
    invalid_bootstrap_fits = simulation$invalid,
    mc_confidence = mc_confidence, minimum_expected_count = 5,
    n_estimation_method = paste(
      "nearest feasible integer to the noisy table total followed by",
      "Euclidean simplex projection"),
    calibration = noise$chisq_calibration,
    reference_mechanism = if (noise$draw_count == 1L) paste(
      "multinomial null plus one joint complete signed-plan dyadic",
      "discrete-Laplace draw per",
      "cell on the common lattice and exact [0, capacity] clamp") else paste(
      "multinomial null plus two independent complete discrete-Laplace",
      "draws per cell on the common lattice and exact [0, capacity] clamp"),
    finite_sampler_handling = if (noise$draw_count == 1L) paste(
      "signed-plan dyadic target; exact-GC finite-sampler and outward-bounded",
      "R-parameter conversion TV propagated twice for the observed-release",
      "and bootstrap-reference transfers to p_value_calibration_interval")
    else paste(
      "ideal geometric reference; signed two-peer finite-sampler total",
      "variation propagated to p_value_calibration_interval"),
    finite_sampler_bit_exact_reference = FALSE,
    mechanism_noise_variance_per_cell = noise$variance,
    mechanism_noise_variance_per_cell_is_upper_bound =
      noise$variance_is_upper_bound,
    mechanism_vector_tv_upper_bound = noise$vector_tv_upper_bound,
    mechanism_numeric_reference_tv_upper_bound =
      noise$numeric_reference_tv_upper_bound,
    mechanism_one_transfer_tv_upper_bound =
      noise$one_transfer_tv_upper_bound,
    mechanism_reference_tv_upper_bound =
      noise$calibration_tv_upper_bound,
    bootstrap_seed_source = paste(
      "SHA-256 of public artifact/vector/orientation commitments;",
      "not analyst-controlled and not privacy randomness"),
    bootstrap_seed_id = digest::digest(
      paste(noise$chisq_calibration, seed, sep = "|"),
      algo = "sha256", serialize = FALSE),
    additional_server_calls = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    inferential_contract = paste(
      "DP-aware parametric bootstrap with plug-in nuisance parameters;",
      "asymptotically calibrated for positive cell probabilities; not a",
      "finite-sample exact conditional test; MC and sampler-TV uncertainty",
      "reported separately"),
    source_release = .dsvert_dp_table_source_binding(x),
    source_dp_release = x)
  class(result) <- c("ds.vertChisq", "list")
  result
}

#' @title Compute chi-square statistics from observed counts
#' @description Pure helper exposed internally so unit tests can exercise the
#'   math path without a DataSHIELD round trip.
#'
#' @param observed Integer matrix of cell counts.
#' @param row_margins Integer vector of row sums.
#' @param col_margins Integer vector of column sums.
#' @param n Total observation count.
#' @param correct Apply Yates continuity correction for 2x2 tables.
#' @return list with statistic, df, p_value, observed, expected,
#'   residuals, n, correct.
#' @keywords internal
.dsvert_chisq_compute <- function(observed, row_margins, col_margins, n,
                                  correct = TRUE) {
  observed <- as.matrix(observed)
  n <- as.numeric(n)
  if (n < 1 || sum(observed) < 1) {
    stop("Contingency table has no observations", call. = FALSE)
  }
  row_m <- as.numeric(row_margins)
  col_m <- as.numeric(col_margins)
  if (length(row_m) != nrow(observed) ||
      length(col_m) != ncol(observed) ||
      any(!is.finite(row_m)) || any(!is.finite(col_m)) ||
      any(row_m <= 0) || any(col_m <= 0)) {
    .dsvert_stop_non_identifiable(
      paste0("The contingency-table independence model has an empty or ",
             "invalid row/column margin; its requested degrees of freedom ",
             "are not identifiable."),
      reason = "degenerate_contingency_margins")
  }
  expected <- outer(row_m, col_m) / n

  dims <- dim(observed)
  is_2x2 <- identical(dims, c(2L, 2L))
  correct <- isTRUE(correct) && is_2x2

  eps <- .Machine$double.eps
  if (correct) {
    diff <- abs(observed - expected) - 0.5
    diff[diff < 0] <- 0
    stat <- sum(diff^2 / pmax(expected, eps))
  } else {
    stat <- sum((observed - expected)^2 / pmax(expected, eps))
  }
  df <- (dims[1] - 1L) * (dims[2] - 1L)
  p <- pchisq(stat, df = df, lower.tail = FALSE)
  residuals <- (observed - expected) / sqrt(pmax(expected, eps))

  list(
    statistic = stat,
    df = as.integer(df),
    p_value = p,
    observed = observed,
    expected = expected,
    residuals = residuals,
    n = n,
    correct = correct)
}

#' @title DP-aware conditional test for a 2-by-2 contingency release
#' @description Test association by conditioning a plug-in latent table on
#'   deterministic integer margins, then reproducing the signed Synopsis noise
#'   and clamp in a Monte Carlo reference law. The input is either one already
#'   released `ds.vertDPContingency` object (zero DSI calls) or data and variable
#'   identifiers, in which case `ds.vertDPContingency` is called exactly once.
#'
#' @param data_name A released `ds.vertDPContingency` object, the protected
#'   data-frame name, or a reusable `ds.vertFederation`.
#' @param var1,var2 Row and column variables for a character `data_name`. They
#'   may be omitted for an existing release.
#' @param server Optional owner-peer assertion forwarded to
#'   \code{ds.vertDPContingency}. No separate column-discovery request is
#'   performed. It must be omitted for an existing release.
#' @param alternative One of `"two.sided"`, `"greater"`, or `"less"`. The
#'   directional alternatives refer to an odds ratio above or below one in the
#'   displayed row/column orientation.
#' @param conf.int Logical compatibility argument. A classical conditional
#'   odds-ratio interval is deliberately not returned from a DP release.
#' @param conf.level Compatibility confidence level recorded in the result.
#' @param datasources DataSHIELD connections. Omit for an existing release.
#' @param simulations Number of Monte Carlo replicates.
#' @param mc_confidence Confidence level of the exact binomial Monte Carlo
#'   interval.
#'
#' @return An object of class \code{ds.vertFisher}. It retains the historical
#'   `p_value`, `odds_ratio`, and `conf_int` fields, but the p-value is from the
#'   explicitly labelled DP-aware plug-in conditional bootstrap and
#'   `conf_int` is always `NULL`.
#'
#' @details The projected DP table supplies a fitted contributing total and
#'   margins. Those margins are rounded deterministically to a feasible
#'   positive integer 2-by-2 configuration. Under the odds-ratio-one null, the
#'   upper-left latent cell is drawn from its hypergeometric law. Each simulated
#'   latent table then receives the signed plan's dyadic one-draw reference or
#'   the two-peer discrete-Laplace law and public clamp, after which
#'   nuisance
#'   margins and the signed root-Pearson statistic are refitted.
#'
#'   This is not Fisher's finite-sample exact test for the confidential table:
#'   noisy projected margins are neither the confidential margins nor exact
#'   nuisance-sufficient statistics. The plug-in calibration is asymptotic.
#'   Monte Carlo error, signed finite-sampler total variation, and the
#'   exact-GC numeric-reference transfer are reported separately. Degenerate
#'   fitted margins return a structured
#'   non-tested result. Only the signed Ring128 discrete-Laplace Synopsis artifact is
#'   currently certified; other mechanisms fail with a typed condition rather
#'   than being silently approximated.
#'
#' @export
ds.vertFisher <- function(data_name, var1 = NULL, var2 = NULL, server = NULL,
                          alternative = c("two.sided", "greater", "less"),
                          conf.int = TRUE, conf.level = 0.95,
                          datasources = NULL, simulations = 9999L,
                          mc_confidence = 0.95) {
  alternative <- match.arg(alternative)
  if (!is.logical(conf.int) || length(conf.int) != 1L || is.na(conf.int)) {
    stop("conf.int must be one non-missing logical value", call. = FALSE)
  }
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      is.na(conf.level) || !is.finite(conf.level) || conf.level <= 0 ||
      conf.level >= 1) {
    stop("conf.level must lie strictly between zero and one", call. = FALSE)
  }
  .dsvert_dp_chisq_validate_simulation_inputs(simulations, mc_confidence)

  if (inherits(data_name, "ds.vertDPContingency")) {
    release <- data_name
    if (!is.null(datasources) || !is.null(server)) {
      stop("server and datasources must be omitted for an existing DP release",
           call. = FALSE)
    }
    assertions <- list(row = var1, column = var2)
    expected <- list(row = release$row_var, column = release$col_var)
    for (name in names(assertions)) {
      value <- assertions[[name]]
      if (!is.null(value) && (!is.character(value) || length(value) != 1L ||
          is.na(value) || !identical(value, expected[[name]]))) {
        stop(name, " variable assertion does not match the DP release",
             call. = FALSE)
      }
    }
  } else {
    if (inherits(data_name, "ds.vertFederation")) {
      resolved <- .dsvert_federation_argument(data_name, datasources)
      data_name <- resolved$value
      datasources <- resolved$datasources
    }
    for (value in list(data_name = data_name, var1 = var1, var2 = var2)) {
      if (!is.character(value) || length(value) != 1L || is.na(value) ||
          !nzchar(value)) {
        stop("data_name, var1 and var2 must be non-empty strings",
             call. = FALSE)
      }
    }
    release <- ds.vertDPContingency(
      data_name = data_name, row_var = var1, col_var = var2,
      server = server, datasources = datasources)
  }

  .dsvert_dp_fisher_from_release(
    release, alternative = alternative, conf.int = conf.int,
    conf.level = conf.level, simulations = simulations,
    mc_confidence = mc_confidence)
}

.dsvert_dp_fisher_mechanism_dispatch <- function(x) {
  profile <- .dsvert_dp_table_vector_profile(x)
  mode <- .dsvert_dp_chisq_laplace_mode(x)
  if (!is.null(mode)) return(mode)
  gaussian <- is.list(profile) && isTRUE(profile$gaussian)
  detail <- if (gaussian) {
    paste(
      "The Gaussian contingency mechanism does not yet carry the conditional",
      "reference-law and total-variation certificate required by ds.vertFisher.")
  } else {
    paste(
      "ds.vertFisher currently supports only signed Ring128 one-draw or",
      "two-peer discrete-Laplace contingency mechanisms.")
  }
  .dsvert_stop_numeric(
    "dsvert_dp_fisher_mechanism_not_certified", detail,
    reason = if (gaussian) {
      "gaussian_conditional_reference_not_certified"
    } else {
      "unsupported_dp_fisher_mechanism"
    })
}

.dsvert_dp_fisher_seed <- function(x, alternative, calibration) {
  release_identity <- if (.dsvert_vector_hex(x$artifact_key)) {
    x$artifact_key
  } else x$capsule_id
  fields <- list(
    calibration, release_identity,
    x$final_vector_root, x$coordinate_order_sha256,
    x$row_var, x$col_var, alternative)
  valid <- vapply(fields, function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) &&
      nzchar(value)
  }, logical(1L))
  if (!all(valid)) {
    stop("The DP release lacks a public reproducibility commitment",
         call. = FALSE)
  }
  fields <- vapply(fields, identity, character(1L))
  hash <- digest::digest(
    paste(enc2utf8(fields), collapse = "|"), algo = "sha256",
    serialize = FALSE)
  as.integer(strtoi(substr(hash, 1L, 7L), base = 16L)) + 1L
}

.dsvert_dp_fisher_fit <- function(table, capacity) {
  fit <- .dsvert_dp_chisq_fit(table, capacity, minimum_expected = 0)
  if (!isTRUE(fit$ok)) return(fit)
  if (!identical(dim(table), c(2L, 2L)) || fit$n < 2 ||
      fit$n > 2^53) {
    fit$ok <- FALSE
    fit$status <- "not_tested_degenerate_conditional_margins"
    return(fit)
  }
  nearest_margin <- function(value, total) {
    min(total - 1, max(1, floor(value + 0.5)))
  }
  row_one <- nearest_margin(sum(fit$latent[1L, ]), fit$n)
  col_one <- nearest_margin(sum(fit$latent[, 1L]), fit$n)
  lower <- max(0, row_one + col_one - fit$n)
  upper <- min(row_one, col_one)
  if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
    fit$ok <- FALSE
    fit$status <- "not_tested_degenerate_conditional_support"
    return(fit)
  }
  variance <- row_one * (fit$n - row_one) * col_one *
    (fit$n - col_one) / (fit$n^2 * (fit$n - 1))
  if (!is.finite(variance) || variance <= 0) {
    fit$ok <- FALSE
    fit$status <- "not_tested_degenerate_conditional_support"
    return(fit)
  }
  fit$integer_row_margins <- c(row_one, fit$n - row_one)
  fit$integer_col_margins <- c(col_one, fit$n - col_one)
  names(fit$integer_row_margins) <- rownames(table)
  names(fit$integer_col_margins) <- colnames(table)
  fit$hypergeometric_support <- c(lower = lower, upper = upper)
  fit$hypergeometric_variance <- variance
  fit
}

.dsvert_dp_fisher_direction <- function(table) {
  left <- table[1L, 1L] * table[2L, 2L]
  right <- table[1L, 2L] * table[2L, 1L]
  if (is.finite(left) && is.finite(right)) return(sign(left - right))
  log_product <- function(first, second) {
    if (first <= 0 || second <= 0) -Inf else log(first) + log(second)
  }
  sign(log_product(table[1L, 1L], table[2L, 2L]) -
         log_product(table[1L, 2L], table[2L, 1L]))
}

.dsvert_dp_fisher_statistic <- function(table, fit, alternative) {
  pearson <- sum((table - fit$expected)^2 / fit$expected)
  if (!is.finite(pearson) || pearson < 0) return(NULL)
  signed_root <- .dsvert_dp_fisher_direction(table) * sqrt(pearson)
  tail_statistic <- switch(
    alternative,
    two.sided = abs(signed_root),
    greater = signed_root,
    less = -signed_root)
  list(
    statistic = tail_statistic,
    signed_root_pearson = signed_root,
    pearson = pearson,
    direction = .dsvert_dp_fisher_direction(table))
}

.dsvert_dp_fisher_latent_tables <- function(replicates, fit) {
  row_one <- fit$integer_row_margins[[1L]]
  col_one <- fit$integer_col_margins[[1L]]
  upper_left <- stats::rhyper(
    replicates, m = col_one, n = fit$n - col_one, k = row_one)
  rbind(
    upper_left,
    row_one - upper_left,
    col_one - upper_left,
    fit$n - row_one - col_one + upper_left)[c(1L, 3L, 2L, 4L), ,
                                               drop = FALSE]
}

.dsvert_dp_fisher_odds_ratio <- function(table) {
  values <- as.numeric(table)
  if (length(values) != 4L || anyNA(values) || any(!is.finite(values)) ||
      any(values < 0)) return(NA_real_)
  numerator_zero <- values[[1L]] <= 0 || values[[4L]] <= 0
  denominator_zero <- values[[2L]] <= 0 || values[[3L]] <= 0
  if (numerator_zero && denominator_zero) return(NA_real_)
  if (denominator_zero) return(Inf)
  if (numerator_zero) return(0)
  log_ratio <- log(values[[1L]]) + log(values[[4L]]) -
    log(values[[2L]]) - log(values[[3L]])
  if (log_ratio > log(.Machine$double.xmax)) return(Inf)
  if (log_ratio < log(.Machine$double.xmin)) return(0)
  exp(log_ratio)
}

.dsvert_dp_fisher_unavailable <- function(
    x, fit, alternative, conf.int, conf.level, simulations, mc_confidence,
    noise) {
  latent <- if (is.null(fit$latent)) {
    matrix(NA_real_, 2L, 2L, dimnames = dimnames(x$table))
  } else {
    fit$latent
  }
  dimnames(latent) <- dimnames(x$table)
  result <- list(
    status = fit$status, decision_available = FALSE,
    method = "DP-aware conditional hypergeometric plug-in Monte Carlo test",
    statistic = NA_real_, statistic_name = "signed root Pearson",
    signed_root_pearson = NA_real_, p_value = NA_real_,
    p_value_mc_se = NA_real_,
    p_value_mc_interval = c(lower = NA_real_, upper = NA_real_),
    p_value_calibration_interval = c(lower = NA_real_, upper = NA_real_),
    odds_ratio = .dsvert_dp_fisher_odds_ratio(latent), conf_int = NULL,
    conf_int_requested = conf.int, conf_level = conf.level,
    confidence_interval_available = FALSE,
    confidence_interval_reason = paste(
      "No classical conditional odds-ratio interval is valid for noisy",
      "projected margins."),
    alternative = alternative, null_odds_ratio = 1,
    observed = x$table, latent_projected_table = latent,
    expected = if (is.null(fit$expected)) {
      matrix(NA_real_, 2L, 2L, dimnames = dimnames(x$table))
    } else fit$expected,
    n = if (is.null(fit$n)) NA_real_ else fit$n, n_na = NA_integer_,
    server = x$server, var1 = x$row_var, var2 = x$col_var,
    simulations = as.integer(simulations), mc_confidence = mc_confidence,
    calibration = noise$fisher_calibration,
    finite_sample_exact = FALSE, raw_table_fisher_exact = FALSE,
    finite_sample_conservative_available = FALSE,
    mechanism_noise_variance_per_cell = noise$variance,
    mechanism_noise_variance_per_cell_is_upper_bound =
      noise$variance_is_upper_bound,
    mechanism_vector_tv_upper_bound = noise$vector_tv_upper_bound,
    mechanism_numeric_reference_tv_upper_bound =
      noise$numeric_reference_tv_upper_bound,
    mechanism_one_transfer_tv_upper_bound =
      noise$one_transfer_tv_upper_bound,
    mechanism_reference_tv_upper_bound =
      noise$calibration_tv_upper_bound,
    additional_server_calls = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    source_dp_release = x,
    inferential_contract = paste(
      "No p-value is returned because the fitted conditional null is",
      "degenerate; no exact-raw or unproved fallback was applied."))
  class(result) <- c("ds.vertFisher", "list")
  result
}

.dsvert_dp_fisher_from_release <- function(
    x, alternative = c("two.sided", "greater", "less"), conf.int = TRUE,
    conf.level = 0.95, simulations = 9999L, mc_confidence = 0.95) {
  alternative <- match.arg(alternative)
  .dsvert_dp_chisq_validate_simulation_inputs(simulations, mc_confidence)
  x <- .dsvert_dp_table_contract(x)
  if (!identical(dim(x$table), c(2L, 2L))) {
    .dsvert_stop_numeric(
      "dsvert_dp_fisher_unsupported_dimension",
      "ds.vertFisher requires a fixed-domain 2-by-2 DP contingency release.",
      reason = "dp_fisher_requires_2x2")
  }
  .dsvert_dp_fisher_mechanism_dispatch(x)
  noise <- .dsvert_dp_chisq_noise_contract(x)
  fit <- .dsvert_dp_fisher_fit(x$table, x$coordinate_maximum)
  if (!isTRUE(fit$ok)) {
    return(.dsvert_dp_fisher_unavailable(
      x, fit, alternative, conf.int, conf.level, simulations,
      mc_confidence, noise))
  }
  observed_stat <- .dsvert_dp_fisher_statistic(
    x$table, fit, alternative)
  if (is.null(observed_stat)) {
    fit$ok <- FALSE
    fit$status <- "not_tested_non_finite_conditional_statistic"
    return(.dsvert_dp_fisher_unavailable(
      x, fit, alternative, conf.int, conf.level, simulations,
      mc_confidence, noise))
  }
  seed <- .dsvert_dp_fisher_seed(
    x, alternative, noise$fisher_calibration)
  batch_limit <- 250000L

  simulation <- .dsvert_dp_chisq_with_seed(seed, {
    exceedances <- 0L
    invalid <- 0L
    completed <- 0L
    while (completed < simulations) {
      batch <- as.integer(min(batch_limit, simulations - completed))
      latent <- .dsvert_dp_fisher_latent_tables(batch, fit)
      mechanism_noise <- matrix(
        .dsvert_dp_chisq_sample_noise(4L * batch, noise),
        nrow = 4L, ncol = batch)
      noisy <- .dsvert_dp_chisq_apply_release_clamp(
        latent, mechanism_noise, x$coordinate_maximum)
      for (index in seq_len(batch)) {
        table <- matrix(noisy[, index], nrow = 2L, ncol = 2L)
        replicate_fit <- .dsvert_dp_fisher_fit(
          table, x$coordinate_maximum)
        statistic <- if (isTRUE(replicate_fit$ok)) {
          .dsvert_dp_fisher_statistic(
            table, replicate_fit, alternative)
        } else NULL
        if (is.null(statistic)) {
          invalid <- invalid + 1L
          exceedances <- exceedances + 1L
        } else {
          tolerance <- 64 * .Machine$double.eps * max(
            1, abs(statistic$statistic), abs(observed_stat$statistic))
          if (statistic$statistic >= observed_stat$statistic - tolerance) {
            exceedances <- exceedances + 1L
          }
        }
      }
      completed <- completed + batch
    }
    list(exceedances = exceedances, invalid = invalid)
  })
  simulations <- as.integer(simulations)
  p_value <- (simulation$exceedances + 1) / (simulations + 1)
  mc_interval <- .dsvert_dp_chisq_mc_interval(
    simulation$exceedances, simulations, mc_confidence)
  calibration_interval <- c(
    lower = max(0, mc_interval[["lower"]] -
                  noise$calibration_tv_upper_bound),
    upper = min(1, mc_interval[["upper"]] +
                  noise$calibration_tv_upper_bound))
  dimnames(fit$latent) <- dimnames(x$table)
  dimnames(fit$expected) <- dimnames(x$table)

  result <- list(
    status = "ok", decision_available = TRUE,
    method = "DP-aware conditional hypergeometric plug-in Monte Carlo test",
    statistic = observed_stat$statistic,
    statistic_name = paste(
      "tail-oriented signed root Pearson with refitted noisy margins"),
    signed_root_pearson = observed_stat$signed_root_pearson,
    pearson_statistic = observed_stat$pearson,
    p_value = p_value,
    p_value_mc_se = sqrt(p_value * (1 - p_value) / (simulations + 1)),
    p_value_mc_interval = mc_interval,
    p_value_calibration_interval = calibration_interval,
    odds_ratio = .dsvert_dp_fisher_odds_ratio(fit$latent),
    odds_ratio_method =
      "cross-product ratio of the continuous simplex-projected DP table",
    conf_int = NULL, conf_int_requested = conf.int,
    conf_level = conf.level, confidence_interval_available = FALSE,
    confidence_interval_reason = paste(
      "A classical conditional odds-ratio interval is not valid for noisy",
      "projected margins; no unproved substitute is returned."),
    alternative = alternative, null_odds_ratio = 1,
    observed = x$table, latent_projected_table = fit$latent,
    expected = fit$expected,
    n = fit$n, n_na = NA_integer_,
    integer_row_margins = fit$integer_row_margins,
    integer_col_margins = fit$integer_col_margins,
    hypergeometric_support = fit$hypergeometric_support,
    hypergeometric_variance = fit$hypergeometric_variance,
    server = x$server, var1 = x$row_var, var2 = x$col_var,
    simulations = simulations, exceedances = simulation$exceedances,
    invalid_bootstrap_fits = simulation$invalid,
    mc_confidence = mc_confidence,
    n_estimation_method = paste(
      "nearest feasible integer to noisy total plus simplex projection;",
      "projected margins rounded deterministically to positive integers"),
    calibration = noise$fisher_calibration,
    conditional_null = paste(
      "upper-left latent cell is hypergeometric under odds ratio one,",
      "conditional on deterministic plug-in integer margins"),
    reference_mechanism = if (noise$draw_count == 1L) paste(
      "conditional hypergeometric latent table plus one joint complete",
      "signed-plan dyadic",
      "discrete-Laplace draw per cell and exact public clamp") else paste(
      "conditional hypergeometric latent table plus two independent complete",
      "discrete-Laplace draws per cell and exact public clamp"),
    finite_sampler_handling = if (noise$draw_count == 1L) paste(
      "signed-plan dyadic target; exact-GC finite-sampler and outward-bounded",
      "R-parameter conversion TV propagated twice for the observed-release",
      "and bootstrap-reference transfers to p_value_calibration_interval")
    else paste(
      "ideal geometric reference; signed two-peer finite-sampler total",
      "variation propagated to p_value_calibration_interval"),
    finite_sampler_bit_exact_reference = FALSE,
    finite_sample_exact = FALSE, raw_table_fisher_exact = FALSE,
    finite_sample_conservative_available = FALSE,
    finite_sample_conservative_reason = paste(
      "No nuisance-uniform finite-sample inversion has been proved for the",
      "noisy clamped plug-in margins."),
    mechanism_noise_variance_per_cell = noise$variance,
    mechanism_noise_variance_per_cell_is_upper_bound =
      noise$variance_is_upper_bound,
    mechanism_vector_tv_upper_bound = noise$vector_tv_upper_bound,
    mechanism_numeric_reference_tv_upper_bound =
      noise$numeric_reference_tv_upper_bound,
    mechanism_one_transfer_tv_upper_bound =
      noise$one_transfer_tv_upper_bound,
    mechanism_reference_tv_upper_bound =
      noise$calibration_tv_upper_bound,
    bootstrap_seed_source = paste(
      "SHA-256 of public artifact/vector/orientation/alternative commitments;",
      "not analyst-controlled and not privacy randomness"),
    bootstrap_seed_id = digest::digest(
      paste(noise$fisher_calibration, seed, sep = "|"),
      algo = "sha256", serialize = FALSE),
    additional_server_calls = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    disclosure_guard = list(
      satisfied = TRUE,
      basis = "validated formal sticky DP artifact post-processing"),
    inferential_contract = paste(
      "DP-aware conditional hypergeometric plug-in bootstrap; asymptotically",
      "calibrated under positive margins, not Fisher-exact for confidential",
      "data; Monte Carlo and signed sampler-TV uncertainty reported"),
    source_release = .dsvert_dp_table_source_binding(x),
    source_dp_release = x)
  class(result) <- c("ds.vertFisher", "list")
  result
}

#' @export
print.ds.vertFisher <- function(x, ...) {
  cat(sprintf(
    "dsVert DP-aware conditional 2x2 test on %s x %s (server: %s)\n",
    x$var1, x$var2, x$server))
  cat("  fitted contributing n =", format(x$n), "\n")
  if (!isTRUE(x$decision_available)) {
    cat("  no p-value:", x$status, "\n")
  } else {
    cat(sprintf("  signed-root statistic = %.4f  DP-aware p-value = %s\n",
                x$signed_root_pearson,
                format.pval(x$p_value, digits = 4L)))
    cat("  Monte Carlo interval: [",
        format(x$p_value_mc_interval[["lower"]], digits = 4L), ", ",
        format(x$p_value_mc_interval[["upper"]], digits = 4L), "]\n",
        sep = "")
    cat("  mechanism+MC calibration interval: [",
        format(x$p_value_calibration_interval[["lower"]], digits = 4L),
        ", ",
        format(x$p_value_calibration_interval[["upper"]], digits = 4L),
        "]\n", sep = "")
    if (!is.na(x$odds_ratio)) {
      cat("  projected-table cross-product odds ratio =",
          format(x$odds_ratio, digits = 5L), "\n")
    }
  }
  cat("  not a finite-sample Fisher-exact test for confidential data\n")
  cat("  client-only post-processing; additional privacy cost: (0, 0)\n")
  cat("\nReleased DP counts:\n")
  print(x$observed)
  invisible(x)
}

#' @export
print.ds.vertChisq <- function(x, ...) {
  dp_aware <- is.character(x$calibration) &&
    length(x$calibration) == 1L && !is.na(x$calibration) &&
    x$calibration %in% c(
      .DSVERT_DP_CHISQ_BOOTSTRAP_VERSION,
      .DSVERT_DP_CHISQ_EXACT_BOOTSTRAP_VERSION)
  heading <- if (dp_aware) "dsVert DP-aware independence test" else
    "dsVert chi-square test"
  cat(sprintf("%s on %s x %s (server: %s)\n",
              heading, x$var1, x$var2, x$server))
  if (x$correct) cat("  (Yates continuity correction applied)\n")
  if (dp_aware) {
    cat("  fitted contributing n =", format(x$n), "\n")
    if (!isTRUE(x$decision_available)) {
      cat("  no p-value:", x$status, "\n")
    } else {
      cat(sprintf("  statistic = %.4f  nominal df = %d  DP-aware p-value = %s\n",
                  x$statistic, x$df,
                  format.pval(x$p_value, digits = 4L)))
      cat("  Monte Carlo interval: [",
          format(x$p_value_mc_interval[["lower"]], digits = 4L), ", ",
          format(x$p_value_mc_interval[["upper"]], digits = 4L), "]\n",
          sep = "")
      cat("  mechanism+MC calibration interval: [",
          format(x$p_value_calibration_interval[["lower"]], digits = 4L),
          ", ",
          format(x$p_value_calibration_interval[["upper"]], digits = 4L),
          "]\n", sep = "")
    }
    cat("  client-only post-processing; additional privacy cost: (0, 0)\n")
    cat("\nReleased DP counts:\n")
  } else {
    cat(sprintf("  n = %d (n_na = %d)\n", as.integer(x$n), x$n_na))
    cat(sprintf("  chi-sq = %.4f  df = %d  p-value = %s\n",
                x$statistic, x$df,
                format.pval(x$p_value, digits = 4L)))
    cat("\nObserved counts:\n")
  }
  print(x$observed)
  invisible(x)
}
