# Plaintext specification oracle for the formal GLM artifact.  This code is
# deliberately internal and accepts only explicitly supplied public fixtures.
# Production data must never be routed through it.

.dsvert_formal_glm_fixture_stop <- function(message, reason) {
  condition <- structure(
    list(message = message, call = NULL, reason = reason),
    class = c("dsvert_formal_glm_fixture_error",
              "dsvert_numeric_condition", "error", "condition"))
  stop(condition)
}

.dsvert_formal_glm_rat_try <- function(value) {
  tryCatch(.dsvert_glm_rat(value), error = function(error) NULL)
}

.dsvert_formal_glm_clamp <- function(value, lower, upper) {
  if (.dsvert_glm_rat_cmp(value, lower) < 0L) return(.dsvert_glm_rat(lower))
  if (.dsvert_glm_rat_cmp(value, upper) > 0L) return(.dsvert_glm_rat(upper))
  .dsvert_glm_rat(value)
}

.dsvert_formal_glm_fixture_rows <- function(compilation, data) {
  .dsvert_formal_glm_validate_compilation(compilation)
  if (!is.data.frame(data)) {
    .dsvert_formal_glm_fixture_stop(
      "The formal GLM oracle accepts only an explicit data-frame fixture",
      "invalid_reference_fixture")
  }
  schema <- compilation$unsigned_schema
  estimand <- schema$estimand
  capacity <- as.integer(estimand$capacity_normalizer)
  if (nrow(data) > capacity) {
    .dsvert_formal_glm_fixture_stop(
      "The reference fixture exceeds the signed patient capacity",
      "reference_fixture_exceeds_capacity")
  }
  required <- unique(c(
    estimand$response,
    unlist(estimand$predictors, use.names = FALSE),
    estimand$weights$column,
    estimand$offset$column))
  required <- required[!vapply(required, is.null, logical(1L))]
  required <- as.character(required)
  if (any(!required %in% names(data))) {
    .dsvert_formal_glm_fixture_stop(
      "The reference fixture omits a registered model column",
      "invalid_reference_fixture")
  }
  if (nrow(data) < capacity) {
    padding <- as.data.frame(
      setNames(replicate(ncol(data), rep(NA, capacity - nrow(data)),
                         simplify = FALSE), names(data)),
      stringsAsFactors = FALSE)
    data <- rbind(data, padding)
  }
  numeric <- schema$numeric
  coefficients <- schema$estimand$coefficients
  columns <- estimand$column_registry
  family <- estimand$family
  rows <- vector("list", capacity)

  for (row_index in seq_len(capacity)) {
    admitted <- TRUE
    weight <- .dsvert_glm_rat("1")
    if (identical(estimand$weights$mode, "bounded_column")) {
      raw <- .dsvert_formal_glm_rat_try(
        data[[estimand$weights$column]][[row_index]])
      if (is.null(raw)) {
        admitted <- FALSE
      } else {
        weight <- .dsvert_formal_glm_clamp(
          raw, "0", estimand$weights$source_maximum_patient_weight)
        weight <- .dsvert_glm_rat_round_dyadic(
          weight, numeric$x_fraction_bits)
      }
    }

    response <- .dsvert_formal_glm_rat_try(
      data[[estimand$response]][[row_index]])
    if (is.null(response)) {
      admitted <- FALSE
      response <- .dsvert_glm_rat("0")
    } else if (identical(family, "binomial")) {
      if (.dsvert_glm_rat_cmp(response, "0") != 0L &&
          .dsvert_glm_rat_cmp(response, "1") != 0L) {
        admitted <- FALSE
        response <- .dsvert_glm_rat("0")
      }
    } else {
      rounded <- .dsvert_glm_rat_round_dyadic(response, 0L)
      if (.dsvert_glm_rat_cmp(response, rounded) != 0L) {
        admitted <- FALSE
        response <- .dsvert_glm_rat("0")
      } else {
        response <- .dsvert_formal_glm_clamp(
          response, "0", columns[[estimand$response]]$upper)
      }
    }

    offset <- .dsvert_glm_rat("0")
    if (!identical(estimand$offset$mode, "none")) {
      raw <- .dsvert_formal_glm_rat_try(
        data[[estimand$offset$column]][[row_index]])
      if (is.null(raw)) {
        admitted <- FALSE
      } else if (identical(estimand$offset$mode, "bounded_offset")) {
        offset <- .dsvert_formal_glm_clamp(
          raw, estimand$offset$source_lower, estimand$offset$source_upper)
        offset <- .dsvert_glm_rat_round_dyadic(
          offset, numeric$offset_fraction_bits)
      } else {
        source <- columns[[estimand$offset$column]]
        if (.dsvert_glm_rat_cmp(raw, "0") <= 0L) {
          admitted <- FALSE
        } else {
          raw <- .dsvert_formal_glm_clamp(
            raw, source$lower, source$upper)
          interval <- .dsvert_glm_rat_log_interval(
            raw, numeric$reference_precision_bits)
          offset <- .dsvert_glm_rat_round_dyadic(
            .dsvert_glm_rat_interval_midpoint(interval),
            numeric$offset_fraction_bits)
          offset <- .dsvert_formal_glm_clamp(
            offset, estimand$offset$lower, estimand$offset$upper)
        }
      }
    }

    design <- vector("list", length(coefficients))
    names(design) <- vapply(coefficients, `[[`, character(1L), "name")
    for (index in seq_along(coefficients)) {
      term <- coefficients[[index]]$term
      if (identical(term$kind, "intercept")) {
        design[[index]] <- .dsvert_glm_rat("1")
      } else if (identical(term$kind, "numeric")) {
        raw <- .dsvert_formal_glm_rat_try(
          data[[term$source_column]][[row_index]])
        if (is.null(raw)) {
          admitted <- FALSE
          design[[index]] <- .dsvert_glm_rat("0")
        } else {
          design[[index]] <- .dsvert_glm_rat_round_dyadic(
            .dsvert_formal_glm_clamp(
              raw, term$clipping_lower, term$clipping_upper),
            numeric$x_fraction_bits)
        }
      } else {
        raw <- data[[term$source_column]][[row_index]]
        registered <- unlist(columns[[term$source_column]]$levels,
                             use.names = FALSE)
        if (length(raw) != 1L || is.na(raw) ||
            !as.character(raw) %in% registered) {
          admitted <- FALSE
          design[[index]] <- .dsvert_glm_rat("0")
        } else {
          design[[index]] <- .dsvert_glm_rat(
            if (identical(as.character(raw), term$source_level)) "1" else "0")
        }
      }
    }
    if (!admitted || .dsvert_glm_rat_cmp(weight, "0") == 0L) {
      weight <- .dsvert_glm_rat("0")
      response <- .dsvert_glm_rat("0")
      offset <- .dsvert_glm_rat("0")
      design <- lapply(design, function(value) .dsvert_glm_rat("0"))
      names(design) <- vapply(coefficients, `[[`, character(1L), "name")
    }
    rows[[row_index]] <- list(
      weight = weight, design = design, response = response, offset = offset)
  }
  rows
}

.dsvert_formal_glm_pwl_mean <- function(eta, table) {
  knots <- lapply(table$knots, .dsvert_glm_rat)
  values <- lapply(table$values, .dsvert_glm_rat)
  slopes <- lapply(table$slopes, .dsvert_glm_rat)
  if (.dsvert_glm_rat_cmp(eta, knots[[1L]]) < 0L ||
      .dsvert_glm_rat_cmp(eta, knots[[length(knots)]]) > 0L) {
    .dsvert_formal_glm_fixture_stop(
      "The reference linear predictor escaped its public certified domain",
      "internal_reference_bound_violation")
  }
  if (.dsvert_glm_rat_cmp(eta, knots[[1L]]) == 0L) return(values[[1L]])
  if (.dsvert_glm_rat_cmp(eta, knots[[length(knots)]]) == 0L) {
    return(values[[length(values)]])
  }
  segment <- which(vapply(knots[-1L], function(knot) {
    .dsvert_glm_rat_cmp(eta, knot) <= 0L
  }, logical(1L)))[[1L]]
  .dsvert_glm_rat_add(
    values[[segment]],
    .dsvert_glm_rat_mul(
      slopes[[segment]], .dsvert_glm_rat_sub(eta, knots[[segment]])))
}

.dsvert_formal_glm_rational_oracle <- function(compilation, data) {
  rows <- .dsvert_formal_glm_fixture_rows(compilation, data)
  schema <- compilation$unsigned_schema
  coefficients <- schema$estimand$coefficients
  coefficient_names <- vapply(coefficients, `[[`, character(1L), "name")
  capacity <- as.integer(schema$estimand$capacity_normalizer)
  iterations <- as.integer(schema$optimizer$iterations)
  alpha <- .dsvert_glm_rat(schema$optimizer$alpha)
  boxes <- lapply(coefficients, function(value) {
    .dsvert_glm_rat(value$box_abs)
  })
  ridge <- lapply(coefficients, function(value) {
    .dsvert_glm_rat(value$ridge)
  })
  beta <- rep(list(.dsvert_glm_rat("0")), length(coefficients))
  names(beta) <- coefficient_names

  for (iteration in seq_len(iterations)) {
    gradient_sum <- rep(list(.dsvert_glm_rat("0")), length(coefficients))
    for (row in rows) {
      eta <- row$offset
      for (index in seq_along(beta)) {
        eta <- .dsvert_glm_rat_add(
          eta, .dsvert_glm_rat_mul(row$design[[index]], beta[[index]]))
      }
      mean <- .dsvert_formal_glm_pwl_mean(
        eta, schema$link_surrogate)
      residual <- .dsvert_glm_rat_sub(mean, row$response)
      for (index in seq_along(beta)) {
        contribution <- .dsvert_glm_rat_mul(
          .dsvert_glm_rat_mul(row$weight, row$design[[index]]), residual)
        gradient_sum[[index]] <- .dsvert_glm_rat_add(
          gradient_sum[[index]], contribution)
      }
    }
    for (index in seq_along(beta)) {
      gradient <- .dsvert_glm_rat_add(
        .dsvert_glm_rat_div(gradient_sum[[index]], as.character(capacity)),
        .dsvert_glm_rat_mul(ridge[[index]], beta[[index]]))
      candidate <- .dsvert_glm_rat_sub(
        beta[[index]], .dsvert_glm_rat_mul(alpha, gradient))
      beta[[index]] <- .dsvert_formal_glm_clamp(
        candidate, .dsvert_glm_rat_neg(boxes[[index]]), boxes[[index]])
    }
  }

  rational <- lapply(beta, .dsvert_glm_rat_json)
  decoded <- vapply(beta, .dsvert_glm_rat_double, numeric(1L))
  names(rational) <- names(decoded) <- coefficient_names
  if (any(!is.finite(decoded))) {
    .dsvert_stop_numeric_unrepresentable(
      "The rational oracle result cannot be decoded as a finite R double",
      list(status = "numeric_backend_unrepresentable",
           required_ring_bits =
             schema$numeric$plan$required_signed_ring_bits_upper,
           required_frac_bits = schema$numeric$beta_fraction_bits,
           numeric_error_budget = "exact rational result remains available"))
  }
  structure(list(
    status = "estimable_regularized_reference",
    coefficients = decoded,
    coefficients_rational = rational,
    schema_sha256 = compilation$sha256,
    family = schema$estimand$family,
    estimand = schema$estimand$target,
    iterations = as.numeric(iterations),
    logical_work = list(
      rows_per_iteration = as.numeric(capacity),
      coefficient_coordinates = as.numeric(length(beta)),
      link_segments = schema$link_surrogate$segment_count,
      data_dependent_stopping = FALSE),
    certificate = schema$theorem_certificate),
    class = "dsvert_formal_glm_reference")
}

.dsvert_formal_glm_ordinary_reference <- function(compilation, data) {
  rows <- .dsvert_formal_glm_fixture_rows(compilation, data)
  schema <- compilation$unsigned_schema
  p <- length(schema$estimand$coefficients)
  design <- matrix(0, nrow = length(rows), ncol = p)
  outcome <- offset <- weight <- numeric(length(rows))
  for (row_index in seq_along(rows)) {
    design[row_index, ] <- vapply(
      rows[[row_index]]$design, .dsvert_glm_rat_double, numeric(1L))
    outcome[[row_index]] <- .dsvert_glm_rat_double(rows[[row_index]]$response)
    offset[[row_index]] <- .dsvert_glm_rat_double(rows[[row_index]]$offset)
    weight[[row_index]] <- .dsvert_glm_rat_double(rows[[row_index]]$weight)
  }
  active <- weight > 0
  certificate <- list(
    status = "ordinary_reference_checked",
    family = schema$estimand$family,
    target = "ordinary_unbounded_unpenalized_glm_comparator",
    schema_sha256 = compilation$sha256,
    not_the_private_release_estimand = TRUE)
  if (!any(active)) {
    .dsvert_stop_non_identifiable(
      "The ordinary GLM comparator has no active complete patient tuple",
      reason = "no_complete_patient_tuple", certificate = certificate)
  }
  rank <- qr(design[active, , drop = FALSE], tol = 1e-12)$rank
  if (rank < p) {
    .dsvert_stop_non_identifiable(
      "The ordinary GLM comparator is collinear",
      reason = "rank_deficient_design", certificate = certificate)
  }
  if (identical(schema$estimand$family, "poisson") &&
      all(outcome[active] == 0)) {
    .dsvert_stop_non_identifiable(
      "The ordinary Poisson comparator has no finite intercept MLE",
      reason = "all_zero_poisson_boundary", certificate = certificate)
  }
  warnings <- character()
  family <- if (identical(schema$estimand$family, "binomial")) {
    stats::binomial(link = "logit")
  } else {
    stats::poisson(link = "log")
  }
  fit <- tryCatch(
    withCallingHandlers(
      stats::glm.fit(
        x = design, y = outcome, weights = weight, offset = offset,
        family = family, control = stats::glm.control(maxit = 100L,
                                                      epsilon = 1e-12)),
      warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }),
    error = function(error) error)
  failed <- inherits(fit, "error") || !isTRUE(fit$converged) ||
    any(!is.finite(fit$coefficients)) ||
    any(grepl("numerically 0 or 1|numerically 0 occurred", warnings))
  if (failed) {
    reason <- if (identical(schema$estimand$family, "binomial")) {
      "separation_or_no_finite_binomial_mle"
    } else {
      "no_finite_poisson_mle"
    }
    .dsvert_stop_non_identifiable(
      "The ordinary unpenalized GLM comparator has no certified finite fit",
      reason = reason, certificate = certificate)
  }
  names(fit$coefficients) <- vapply(
    schema$estimand$coefficients, `[[`, character(1L), "name")
  list(status = "ordinary_reference_estimable",
       coefficients = fit$coefficients, warnings = warnings,
       certificate = certificate)
}

# Plaintext specification oracle for the separately released inference vector.
# It is intentionally fixture-only.  Production protected rows must never be
# passed to this function or to any client-side equivalent.
.dsvert_formal_glm_pwl_slope_at <- function(eta, table) {
  knots <- lapply(table$knots, .dsvert_glm_rat)
  slopes <- lapply(table$slopes, .dsvert_glm_rat)
  if (.dsvert_glm_rat_cmp(eta, knots[[1L]]) < 0L ||
      .dsvert_glm_rat_cmp(eta, knots[[length(knots)]]) > 0L) {
    .dsvert_formal_glm_fixture_stop(
      "The inference oracle escaped its public PWL domain",
      "internal_reference_bound_violation")
  }
  # The same convention as .dsvert_formal_glm_pwl_mean(): an internal knot
  # belongs to the segment on its left.  Endpoint choices are immaterial for
  # the primitive but are material for the released working information.
  segment <- which(vapply(knots[-1L], function(knot) {
    .dsvert_glm_rat_cmp(eta, knot) <= 0L
  }, logical(1L)))[[1L]]
  slopes[[segment]]
}

.dsvert_formal_glm_pwl_primitive_double <- function(value, table) {
  knots <- vapply(table$knots, .dsvert_glm_rat_double, numeric(1L))
  means <- vapply(table$values, .dsvert_glm_rat_double, numeric(1L))
  slopes <- vapply(table$slopes, .dsvert_glm_rat_double, numeric(1L))
  if (!is.finite(value) || value < knots[[1L]] ||
      value > knots[[length(knots)]] || 0 < knots[[1L]] ||
      0 > knots[[length(knots)]]) {
    .dsvert_formal_glm_fixture_stop(
      "The inference primitive escaped its public PWL domain",
      "internal_reference_bound_violation")
  }
  from_lower <- function(target) {
    total <- 0
    for (segment in seq_along(slopes)) {
      left <- knots[[segment]]
      right <- min(target, knots[[segment + 1L]])
      if (right > left) {
        width <- right - left
        total <- total + means[[segment]] * width +
          slopes[[segment]] * width^2 / 2
      }
      if (target <= knots[[segment + 1L]]) break
    }
    total
  }
  from_lower(value) - from_lower(0)
}

.dsvert_formal_glm_joint_inference_oracle <- function(
    point_result, compilation, data) {
  if (!inherits(point_result, "dsvert_formal_dp_glm")) {
    .dsvert_formal_glm_fixture_stop(
      "The inference oracle requires a validated formal DP beta",
      "invalid_reference_fixture")
  }
  .dsvert_formal_glm_validate_compilation(compilation)
  if (!identical(point_result$schema_sha256, compilation$sha256)) {
    .dsvert_formal_glm_fixture_stop(
      "The inference oracle beta and compilation do not match",
      "invalid_reference_fixture")
  }
  rows <- .dsvert_formal_glm_fixture_rows(compilation, data)
  schema <- compilation$unsigned_schema
  coefficients <- unlist(
    schema$estimand$coefficient_order, use.names = FALSE)
  p <- length(coefficients)
  capacity <- as.integer(schema$estimand$capacity_normalizer)
  scale <- .dsvert_glm_rat(point_result$output_lattice_scale)
  beta <- lapply(point_result$coefficient_lattice_steps, function(value) {
    .dsvert_glm_rat_div(value, scale)
  })
  names(beta) <- coefficients
  information <- matrix(0, p, p,
                        dimnames = list(coefficients, coefficients))
  meat <- information
  score <- stats::setNames(numeric(p), coefficients)
  canonical_loglik <- 0
  surrogate_loss <- 0
  admitted_n <- 0
  for (row in rows) {
    eta <- row$offset
    for (index in seq_len(p)) {
      eta <- .dsvert_glm_rat_add(
        eta, .dsvert_glm_rat_mul(row$design[[index]], beta[[index]]))
    }
    mean <- .dsvert_formal_glm_pwl_mean(eta, schema$link_surrogate)
    slope <- .dsvert_formal_glm_pwl_slope_at(
      eta, schema$link_surrogate)
    x <- vapply(row$design, .dsvert_glm_rat_double, numeric(1L))
    weight <- .dsvert_glm_rat_double(row$weight)
    y <- .dsvert_glm_rat_double(row$response)
    eta_double <- .dsvert_glm_rat_double(eta)
    mean_double <- .dsvert_glm_rat_double(mean)
    slope_double <- .dsvert_glm_rat_double(slope)
    information <- information + weight * slope_double * tcrossprod(x)
    patient_score <- weight * x * (y - mean_double)
    score <- score + patient_score
    meat <- meat + tcrossprod(patient_score)
    if (weight > 0) {
      admitted_n <- admitted_n + 1
      canonical_loglik <- canonical_loglik + weight * if (
        identical(schema$estimand$family, "binomial")) {
        y * eta_double - if (eta_double > 0) {
          eta_double + log1p(exp(-eta_double))
        } else log1p(exp(eta_double))
      } else {
        y * eta_double - exp(eta_double) - lgamma(y + 1)
      }
      surrogate_loss <- surrogate_loss + weight * (
        .dsvert_formal_glm_pwl_primitive_double(
          eta_double, schema$link_surrogate) - y * eta_double)
    }
  }
  ridge <- vapply(schema$estimand$coefficients, function(coefficient) {
    .dsvert_glm_rat_double(coefficient$ridge)
  }, numeric(1L))
  beta_double <- vapply(beta, .dsvert_glm_rat_double, numeric(1L))
  information <- information + capacity * diag(ridge, p)
  surrogate_loss <- surrogate_loss +
    capacity * sum(ridge * beta_double^2) / 2
  list(
    version = "dsvert-formal-glm-joint-inference-fixture-oracle-v1",
    schema_sha256 = compilation$sha256,
    beta_public_release_sha256 = point_result$public_release_sha256,
    coefficient_order = coefficients,
    information = information,
    score_meat = meat,
    score = score,
    canonical_bounded_log_likelihood_at_dp_beta = canonical_loglik,
    integrated_pwl_surrogate_loss_at_dp_beta = surrogate_loss,
    admitted_n = as.numeric(admitted_n),
    information_scale = "patient_sum_scale",
    protected_data_used_in_production = FALSE)
}

.dsvert_formal_glm_artifact_requirements <- function(compilation) {
  .dsvert_formal_glm_validate_compilation(compilation)
  schema <- compilation$unsigned_schema
  p <- length(schema$estimand$coefficients)
  list(
    version = "formal-glm-phase0-artifact-requirements-v1",
    schema_sha256 = compilation$sha256,
    reusable_capsule_components = list(
      "immutable capsule identity and cross-signed manifest",
      "two-recipient typed encrypted source streams",
      "sticky joint vector ledger and one-global-draw Gaussian planner",
      "dynamic-ring canonical integer encoding"),
    protected_input_artifact = list(
      kind = "fixed_capacity_patient_tuple_shares",
      capacity = schema$estimand$capacity_normalizer,
      design_width = as.numeric(p),
      row_coordinates = as.numeric(p + 3L),
      semantics = as.list(c("weight", "design", "outcome", "offset")),
      finite_sufficient_statistic = FALSE),
    public_dp_output_artifact = list(
      kind = "one_sticky_joint_coefficient_vector",
      coordinate_count = as.numeric(p),
      intermediate_release_count = 0,
      score_hessian_deviance_residual_release = FALSE),
    exact_gc_required = as.list(c(
      "signed_segment_comparison_and_mux",
      "checked_multiply_and_exact_truncate",
      "fixed_order_checked_ring_accumulation",
      "coefficient_box_projection",
      "one_global_draw_two_peer_noise_and_final_projection")),
    cannot_reuse_as_exact_glm_sufficient_statistic = as.list(c(
      "gaussian_second_moments", "categorical_contingency_counts",
      "ordinary_low_order_moment_vector")),
    production_phase = "phase1_not_connected")
}
