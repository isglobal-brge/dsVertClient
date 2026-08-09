.dsvert_numeric_condition <- function(class, message, certificate = NULL,
                                      reason = NULL) {
  structure(
    list(message = message, call = NULL, certificate = certificate,
         reason = reason),
    class = c(class, "dsvert_numeric_condition", "error", "condition")
  )
}

.dsvert_stop_numeric <- function(class, message, certificate = NULL,
                                 reason = NULL) {
  stop(.dsvert_numeric_condition(class, message, certificate, reason))
}

.dsvert_stop_numeric_unrepresentable <- function(message, certificate) {
  condition <- .dsvert_numeric_condition(
    "numeric_backend_unrepresentable", message,
    certificate = certificate,
    reason = "certified_ring_or_precision_exhausted")
  class(condition) <- c(
    "dsvert_numeric_backend_unrepresentable", class(condition))
  condition$required_ring_bits <- certificate$required_ring_bits
  condition$required_frac_bits <- certificate$required_frac_bits
  condition$numeric_error_budget <- certificate$numeric_error_budget
  stop(condition)
}

.dsvert_non_identifiable <- function(message, reason = "singular_information",
                                     certificate = NULL) {
  if (!is.null(certificate)) {
    certificate$status <- "non_identifiable"
    certificate$estimability_status <- "non_identifiable"
    certificate$estimability_reason <- reason
    certificate$estimability_evaluated <- TRUE
    certificate$certificate_phase <- "identifiability_failed"
    certificate$numerically_certified <- FALSE
  }
  .dsvert_numeric_condition(
    "non_identifiable", message, certificate = certificate, reason = reason)
}

.dsvert_stop_non_identifiable <- function(message,
                                           reason = "singular_information",
                                           certificate = NULL) {
  stop(.dsvert_non_identifiable(
    message, reason = reason, certificate = certificate))
}

.dsvert_solve_identifiable <- function(a, b = NULL, context = "model",
                                       reason = "singular_information_matrix",
                                       tol = sqrt(.Machine$double.eps),
                                       symmetric = FALSE,
                                       certificate = NULL) {
  a <- as.matrix(a)
  if (!is.numeric(a) || length(dim(a)) != 2L || nrow(a) != ncol(a) ||
      nrow(a) < 1L || any(!is.finite(a))) {
    .dsvert_stop_non_identifiable(
      paste0(context, " has a malformed or non-finite information matrix."),
      reason = "non_finite_information_matrix", certificate = certificate)
  }
  if (!is.null(b) && (!is.numeric(b) || any(!is.finite(b)) ||
                      NROW(b) != nrow(a))) {
    .dsvert_stop_non_identifiable(
      paste0(context, " has a malformed or non-finite score vector."),
      reason = "non_finite_score", certificate = certificate)
  }
  if (isTRUE(symmetric)) {
    a <- (a + t(a)) / 2
    eigenvalues <- tryCatch(
      eigen(a, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) NULL)
    scale <- if (is.null(eigenvalues)) NA_real_ else max(abs(eigenvalues))
    if (is.null(eigenvalues) || any(!is.finite(eigenvalues)) ||
        !is.finite(scale) || scale <= 0 || min(eigenvalues) <= tol * scale) {
      .dsvert_stop_non_identifiable(
        paste0(context, " is not identifiable: its information matrix is ",
               "not positive definite. No ridge, Firth, or pseudoinverse ",
               "fallback was applied."),
        reason = reason, certificate = certificate)
    }
  }
  rank <- qr(a, tol = tol)$rank
  if (rank < ncol(a)) {
    .dsvert_stop_non_identifiable(
      paste0(context, " is not identifiable: its information matrix has ",
             "rank ", rank, " of ", ncol(a),
             ". No ridge, Firth, or pseudoinverse fallback was applied."),
      reason = reason, certificate = certificate)
  }
  out <- tryCatch(
    if (is.null(b)) solve(a) else solve(a, b),
    error = function(e) NULL)
  if (is.null(out) || any(!is.finite(out))) {
    .dsvert_stop_non_identifiable(
      paste0(context, " is not identifiable: its information system could ",
             "not be solved. No ridge, Firth, or pseudoinverse fallback ",
             "was applied."),
      reason = reason, certificate = certificate)
  }
  out
}

.dsvert_numeric_mark_estimable <- function(
    certificate, reason = "positive_definite_full_rank_information") {
  if (!inherits(certificate, "dsvert_numeric_certificate") ||
      !is.character(reason) || length(reason) != 1L || is.na(reason) ||
      !nzchar(reason)) {
    stop("invalid numeric estimability certificate", call. = FALSE)
  }
  certificate$estimability_status <- "estimable"
  certificate$estimability_reason <- reason
  certificate$estimability_evaluated <- TRUE
  certificate
}

.dsvert_numeric_scalar <- function(value, name, lower = 0,
                                   integer = FALSE, allow_zero = FALSE) {
  valid_shape <- is.numeric(value) && length(value) == 1L &&
    !is.na(value) && is.finite(value)
  bad_lower <- if (!valid_shape) TRUE else if (isTRUE(allow_zero)) {
    value < lower
  } else {
    value <= lower
  }
  if (!valid_shape || bad_lower ||
      (isTRUE(integer) &&
       (value != floor(value) || value > .Machine$integer.max))) {
    qualifier <- if (isTRUE(integer)) "integer" else "number"
    relation <- if (isTRUE(allow_zero)) "at least" else "greater than"
    stop(name, " must be one finite ", qualifier, " ", relation, " ", lower,
         call. = FALSE)
  }
  if (isTRUE(integer)) as.integer(value) else as.numeric(value)
}

.dsvert_numeric_safe_product <- function(...) {
  values <- unlist(list(...), use.names = FALSE)
  if (!length(values) || any(!is.finite(values)) || any(values < 0)) {
    return(Inf)
  }
  positive <- values[values > 0]
  if (!length(positive)) return(0)
  log_value <- sum(log(positive))
  if (!is.finite(log_value) || log_value > log(.Machine$double.xmax)) {
    return(Inf)
  }
  exp(log_value)
}

.dsvert_numeric_log2_sum <- function(...) {
  values <- unlist(list(...), use.names = FALSE)
  if (!length(values) || anyNA(values) || any(values == Inf)) return(Inf)
  values <- values[is.finite(values)]
  if (!length(values)) return(-Inf)
  largest <- max(values)
  largest + log2(sum(2^(values - largest)))
}

.dsvert_numeric_from_log2_bound <- function(value) {
  if (identical(value, -Inf)) return(0)
  if (!is.finite(value) ||
      value >= log2(.Machine$double.xmax)) return(Inf)
  2^value
}

.dsvert_numeric_canonicalize <- function(value) {
  if (is.list(value)) {
    value_names <- names(value)
    if (!is.null(value_names)) {
      if (anyNA(value_names) || anyDuplicated(value_names)) {
        stop("numeric contract objects require unique non-missing names",
             call. = FALSE)
      }
      value <- value[order(value_names)]
    }
    return(lapply(value, .dsvert_numeric_canonicalize))
  }
  value_names <- names(value)
  if (!is.null(value_names)) {
    if (anyNA(value_names) || anyDuplicated(value_names)) {
      stop("numeric contract vectors require unique non-missing names",
           call. = FALSE)
    }
    value <- value[order(value_names)]
  }
  value
}

.dsvert_numeric_hash <- function(value) {
  canonical <- jsonlite::toJSON(
    .dsvert_numeric_canonicalize(value),
    auto_unbox = TRUE, null = "null", na = "null", digits = 17,
    pretty = FALSE)
  digest::digest(canonical, algo = "sha256", serialize = FALSE)
}

.dsvert_numeric_policy_id <- function(policy) {
  policy$policy_id <- NULL
  .dsvert_numeric_hash(policy)
}

.dsvert_numeric_glm_workload <- function(n_obs, n_predictors, family,
                                         max_iter, compute_se = TRUE,
                                         compute_deviance = TRUE,
                                         weights_active = FALSE,
                                         offset_active = FALSE) {
  n_obs <- .dsvert_numeric_scalar(
    n_obs, "n_obs", integer = TRUE)
  n_predictors <- .dsvert_numeric_scalar(
    n_predictors, "n_predictors", integer = TRUE, allow_zero = TRUE)
  max_iter <- .dsvert_numeric_scalar(
    max_iter, "max_iter", integer = TRUE, allow_zero = TRUE)
  if (n_predictors >= .Machine$integer.max) {
    stop("n_predictors is too large to add the intercept safely",
         call. = FALSE)
  }
  if (!is.character(family) || length(family) != 1L || is.na(family) ||
      !family %in% c("gaussian", "binomial", "poisson")) {
    stop("family must be gaussian, binomial, or poisson", call. = FALSE)
  }
  logicals <- list(compute_se = compute_se,
                   compute_deviance = compute_deviance,
                   weights_active = weights_active,
                   offset_active = offset_active)
  if (any(!vapply(logicals, function(x) {
    is.logical(x) && length(x) == 1L && !is.na(x)
  }, logical(1L)))) {
    stop("workload flags must be non-missing logical scalars", call. = FALSE)
  }

  # Keep the workload arithmetic in exactly represented doubles.  R integer
  # addition would otherwise produce NA near .Machine$integer.max before the
  # custodian's public workload limit can reject the request cleanly.
  n_parameters <- as.numeric(n_predictors) + 1
  approximation_degree <- switch(
    family, gaussian = 0L, binomial = 29L, poisson = 30L)
  evaluations <- as.numeric(max_iter) +
    if (isTRUE(compute_se)) n_parameters + 1 else 0
  evaluations <- evaluations + if (isTRUE(compute_deviance)) 1 else 0
  evaluations <- max(1, evaluations)
  inner_terms <- .dsvert_numeric_safe_product(
    n_obs, n_parameters, evaluations)
  operation_count <- .dsvert_numeric_safe_product(
    inner_terms, 2L + approximation_degree)
  truncation_count <- .dsvert_numeric_safe_product(
    inner_terms, 1L + approximation_degree)
  # The current GLM inverse links are fixed-domain Chebyshev evaluations in
  # the share domain.  Domain admission is proved from custodian-owned input
  # bounds plus the public beta/intercept/offset guard before each evaluation;
  # this route performs no secret comparison.
  comparison_count <- 0

  list(
    workload = "glm",
    family = family,
    n_obs = n_obs,
    n_predictors = n_predictors,
    n_parameters = n_parameters,
    max_iter = max_iter,
    operation_count = operation_count,
    truncation_count = truncation_count,
    comparison_count = comparison_count,
    max_sequential_truncations = 1 + as.numeric(approximation_degree),
    quantized_terms_per_path = n_parameters + 4,
    accumulator_terms = n_obs,
    weights_active = isTRUE(weights_active),
    offset_active = isTRUE(offset_active),
    approximation_domain_guard =
      "custodian_input_bounds_plus_public_coefficient_bound",
    requires_exact_truncation = truncation_count > 0,
    requires_exact_comparison = FALSE
  )
}

.dsvert_numeric_validate_policy <- function(policy, server) {
  fail <- function(message) {
    .dsvert_stop_numeric(
      "numeric_backend_unavailable",
      paste0("Invalid numeric policy from ", server, ": ", message),
      reason = "malformed_custodian_policy")
  }
  if (!is.list(policy) || !identical(policy$schema_version, 1L) ||
      !identical(policy$workload, "glm")) {
    fail("unsupported or malformed policy header")
  }
  if (!is.character(policy$policy_version) ||
      length(policy$policy_version) != 1L || is.na(policy$policy_version) ||
      !nzchar(policy$policy_version)) {
    fail("invalid policy version")
  }
  if (!isTRUE(policy$enabled)) {
    .dsvert_stop_numeric(
      "numeric_backend_unavailable",
      paste0("Numeric releases are disabled by custodian ", server),
      reason = "custodian_numeric_policy_disabled")
  }
  if (!is.character(policy$policy_id) || length(policy$policy_id) != 1L ||
      is.na(policy$policy_id) ||
      !grepl("^[0-9a-f]{64}$", policy$policy_id)) {
    fail("invalid policy identifier")
  }
  computed_policy_id <- tryCatch(
    .dsvert_numeric_policy_id(policy), error = function(e) NA_character_)
  if (is.na(computed_policy_id) ||
      !identical(computed_policy_id, policy$policy_id)) {
    fail("policy identifier does not authenticate the returned policy body")
  }

  bounds <- policy$bounds
  required_bounds <- c(
    "max_abs_predictor_input", "max_abs_predictor",
    "max_abs_response_input", "max_abs_response",
    "max_abs_linear_predictor", "max_abs_approximation_intermediate",
    "max_abs_offset", "max_abs_weight",
    "max_observations",
    "max_predictors", "max_iterations", "max_numeric_error")
  if (!is.list(bounds) || is.null(names(bounds)) || anyNA(names(bounds)) ||
      anyDuplicated(names(bounds)) ||
      !all(required_bounds %in% names(bounds))) {
    fail("missing GLM bounds")
  }
  positive_scalar <- function(x) {
    is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) && x > 0
  }
  if (any(!vapply(bounds[c(
    "max_abs_predictor_input", "max_abs_predictor", "max_abs_offset",
    "max_abs_weight", "max_observations", "max_predictors",
    "max_iterations", "max_numeric_error")],
    positive_scalar, logical(1L)))) {
    fail("non-positive or non-finite GLM bound")
  }
  families <- c("gaussian", "binomial", "poisson")
  for (field in c("max_abs_response_input", "max_abs_response",
                  "max_abs_linear_predictor",
                  "max_abs_approximation_intermediate")) {
    value <- bounds[[field]]
    if (!(is.list(value) || is.numeric(value)) ||
        is.null(names(value)) || anyNA(names(value)) ||
        anyDuplicated(names(value)) ||
        !identical(sort(names(value)), sort(families)) ||
        any(!vapply(value[families], positive_scalar, logical(1L)))) {
      fail(paste0("invalid family-specific bound: ", field))
    }
  }

  approximation <- policy$approximation
  if (!is.list(approximation) || !all(families %in% names(approximation))) {
    fail("missing approximation contracts")
  }
  for (family in families) {
    contract <- approximation[[family]]
    if (!is.list(contract) || !is.numeric(contract$max_abs_error) ||
        length(contract$max_abs_error) != 1L ||
        !is.finite(contract$max_abs_error) || contract$max_abs_error < 0) {
      fail(paste0("invalid approximation error for ", family))
    }
    if (!is.null(contract$domain) &&
        (!is.numeric(contract$domain) || length(contract$domain) != 2L ||
         any(!is.finite(contract$domain)) ||
         contract$domain[[1L]] >= contract$domain[[2L]])) {
      fail(paste0("invalid approximation domain for ", family))
    }
  }

  backend_names <- c("ring63", "ring127", "exact_gc", "multiprecision")
  if (!is.list(policy$capabilities) ||
      is.null(names(policy$capabilities)) ||
      anyNA(names(policy$capabilities)) ||
      anyDuplicated(names(policy$capabilities)) ||
      !identical(sort(names(policy$capabilities)), sort(backend_names))) {
    fail("missing backend capabilities")
  }
  logical_fields <- c(
    "available", "allowed", "e2e_verified", "canonical_encoding",
    "fail_closed_overflow", "runtime_bounds_enforced",
    "workload_adapter_e2e_verified",
    "public_scalar_mul_truncate_e2e_verified",
    "full_iteration_e2e_verified", "exact_truncation", "exact_comparison")
  for (backend in backend_names) {
    capability <- policy$capabilities[[backend]]
    if (!is.list(capability) ||
        any(!vapply(capability[logical_fields], function(x) {
          is.logical(x) && length(x) == 1L && !is.na(x)
        }, logical(1L)))) {
      fail(paste0("invalid capability flags for ", backend))
    }
    if (!is.character(capability$truncation_semantics) ||
        length(capability$truncation_semantics) != 1L ||
        is.na(capability$truncation_semantics)) {
      fail(paste0("invalid truncation contract for ", backend))
    }
    integer_vector <- function(x, allow_empty = FALSE) {
      is.numeric(x) && !anyNA(x) && all(is.finite(x)) &&
        all(x == floor(x)) && all(x >= 2) &&
        (isTRUE(allow_empty) || length(x) > 0L) && !anyDuplicated(x)
    }
    optional_integer <- function(x, lower = 1L) {
      is.numeric(x) && length(x) == 1L &&
        (is.na(x) || (is.finite(x) && x == floor(x) && x >= lower))
    }
    if (!integer_vector(
      capability$supported_ring_bits,
      allow_empty = identical(backend, "multiprecision")) ||
        !optional_integer(capability$ring_bits, 2L) ||
        !optional_integer(capability$frac_bits, 1L) ||
        !optional_integer(capability$max_ring_bits, 2L) ||
        !optional_integer(capability$max_frac_bits, 1L)) {
      fail(paste0("invalid ring/precision metadata for ", backend))
    }
    if (backend %in% c("ring63", "ring127")) {
      expected_ring <- if (identical(backend, "ring63")) 63L else 127L
      expected_frac <- if (identical(backend, "ring63")) 20L else 50L
      if (!identical(as.integer(capability$ring_bits), expected_ring) ||
          !identical(as.integer(capability$frac_bits), expected_frac) ||
          !identical(as.integer(capability$supported_ring_bits),
                     expected_ring) ||
          !identical(as.integer(capability$max_ring_bits), expected_ring) ||
          !identical(as.integer(capability$max_frac_bits), expected_frac)) {
        fail(paste0("fixed wire metadata disagrees with ", backend))
      }
    }
    if (identical(backend, "exact_gc") &&
        (any(capability$supported_ring_bits < 63L) ||
         any(capability$supported_ring_bits > 4096L) ||
         (is.finite(capability$max_ring_bits) &&
          capability$max_ring_bits > 4096L) ||
         (is.finite(capability$max_frac_bits) &&
          capability$max_frac_bits > 4095L))) {
      fail("exact_gc metadata exceeds the Ring63-through-Ring4096 protocol")
    }
    if (isTRUE(capability$e2e_verified) &&
        (!isTRUE(capability$available) ||
         !isTRUE(capability$canonical_encoding) ||
         !isTRUE(capability$fail_closed_overflow) ||
         !isTRUE(capability$runtime_bounds_enforced) ||
         !isTRUE(capability$workload_adapter_e2e_verified) ||
         !isTRUE(capability$public_scalar_mul_truncate_e2e_verified) ||
         !isTRUE(capability$full_iteration_e2e_verified))) {
      fail(paste0("contradictory E2E-verification flags for ", backend))
    }
    if ((isTRUE(capability$exact_truncation) ||
         isTRUE(capability$exact_comparison)) &&
        !isTRUE(capability$e2e_verified)) {
      fail(paste0("exact primitive advertised without E2E verification for ",
                  backend))
    }
    if (isTRUE(capability$exact_truncation) &&
        !capability$truncation_semantics %in% c("floor", "toward_zero")) {
      fail(paste0("exact truncation semantics are invalid for ", backend))
    }
  }
  policy
}

.dsvert_numeric_combine_bounds <- function(policies, workload) {
  get_bound <- function(name, family = NULL, fun = max) {
    values <- vapply(policies, function(policy) {
      value <- policy$bounds[[name]]
      if (!is.null(family)) value <- value[[family]]
      as.numeric(value)
    }, numeric(1L))
    fun(values)
  }

  limits <- list(
    max_observations = get_bound("max_observations", fun = min),
    max_predictors = get_bound("max_predictors", fun = min),
    max_iterations = get_bound("max_iterations", fun = min))
  if (workload$n_obs > limits$max_observations ||
      workload$n_predictors > limits$max_predictors ||
      workload$max_iter > limits$max_iterations) {
    certificate <- list(
      schema_version = 1L, status = "numeric_backend_unavailable",
      workload = workload$workload, family = workload$family,
      requested_backend = NA_character_, effective_backend = NA_character_,
      reason = "workload_exceeds_custodian_bounds", limits = limits)
    class(certificate) <- c("dsvert_numeric_certificate", "list")
    .dsvert_stop_numeric(
      "numeric_backend_unavailable",
      "The requested GLM workload exceeds a custodian-owned public bound",
      certificate, "workload_exceeds_custodian_bounds")
  }

  family <- workload$family
  predictor <- get_bound("max_abs_predictor")
  predictor_input <- get_bound("max_abs_predictor_input")
  response <- get_bound("max_abs_response", family)
  response_input <- get_bound("max_abs_response_input", family)
  eta <- get_bound("max_abs_linear_predictor", family)
  approximation_intermediate <- get_bound(
    "max_abs_approximation_intermediate", family)
  offset <- if (workload$offset_active) get_bound("max_abs_offset") else 0
  if (offset > eta) {
    certificate <- list(
      schema_version = 1L, status = "numeric_backend_unavailable",
      workload = workload$workload, family = family,
      requested_backend = NA_character_, effective_backend = NA_character_,
      public_linear_predictor_bound = eta,
      public_offset_bound = offset,
      reason = "offset_bound_exceeds_linear_predictor_bound")
    class(certificate) <- c("dsvert_numeric_certificate", "list")
    .dsvert_stop_numeric(
      "numeric_backend_unavailable",
      "The public offset bound exceeds the total linear-predictor bound",
      certificate, "offset_bound_exceeds_linear_predictor_bound")
  }
  weight <- if (workload$weights_active) get_bound("max_abs_weight") else 1
  mu_log2 <- switch(family,
    gaussian = log2(eta),
    binomial = 0,
    # Computing exp(eta) first would erase a finite public bound as Inf for
    # eta > log(.Machine$double.xmax).  Planning in log2 space retains the
    # ring requirement and permits a typed Ring4096 boundary failure.
    poisson = eta / log(2))
  residual_log2 <- .dsvert_numeric_log2_sum(log2(response), mu_log2)
  curvature_log2 <- switch(family,
    gaussian = 0,
    binomial = -2,
    poisson = mu_log2)
  gradient_log2 <- log2(workload$n_obs) + log2(predictor) +
    residual_log2 + log2(weight)
  hessian_log2 <- log2(workload$n_obs) + 2 * log2(predictor) +
    curvature_log2 + log2(weight)
  deviance_log2 <- log2(workload$n_obs) + 2 * residual_log2 + log2(weight)
  gradient <- .dsvert_numeric_from_log2_bound(gradient_log2)
  hessian <- .dsvert_numeric_from_log2_bound(hessian_log2)
  deviance <- .dsvert_numeric_from_log2_bound(deviance_log2)
  public_log2_bound <- max(
    log2(eta), log2(approximation_intermediate), mu_log2,
    residual_log2, gradient_log2, hessian_log2, deviance_log2)
  public_bound <- .dsvert_numeric_from_log2_bound(public_log2_bound)

  product_log2_bound <- max(
    log2(predictor) + residual_log2 + log2(weight),
    2 * log2(predictor) + curvature_log2 + log2(weight),
    2 * residual_log2 + log2(weight),
    2 * log2(approximation_intermediate),
    # Matvec/deviance kernels can accumulate products while they are still at
    # scale 2^(2f), before the one exact truncation that restores scale 2^f.
    # Their raw-ring capacity must therefore cover the complete sum, not only
    # one observation's product.  Treating these only as post-truncation
    # magnitudes would understate the required ring width by roughly f bits.
    gradient_log2, hessian_log2, deviance_log2)
  product_bound <- .dsvert_numeric_from_log2_bound(product_log2_bound)

  domains <- lapply(policies, function(policy) {
    policy$approximation[[family]]$domain
  })
  non_null_domains <- Filter(Negate(is.null), domains)
  domain <- if (!length(non_null_domains)) NULL else c(
    max(vapply(non_null_domains, `[[`, numeric(1L), 1L)),
    min(vapply(non_null_domains, `[[`, numeric(1L), 2L)))
  if (!is.null(domain) &&
      (domain[[1L]] >= domain[[2L]] || eta > min(abs(domain)))) {
    certificate <- list(
      schema_version = 1L, status = "approximation_domain_failure",
      workload = workload$workload, family = family,
      requested_backend = NA_character_, effective_backend = NA_character_,
      approximation_domain = domain,
      public_linear_predictor_bound = eta,
      reason = "linear_predictor_bound_outside_approximation_domain")
    class(certificate) <- c("dsvert_numeric_certificate", "list")
    .dsvert_stop_numeric(
      "approximation_domain_failure",
      "The custodian-owned linear-predictor bound exceeds the validated approximation domain",
      certificate, "linear_predictor_bound_outside_approximation_domain")
  }

  list(
    public_predictor_input_bound = predictor_input,
    public_response_input_bound = response_input,
    public_predictor_bound = predictor,
    public_response_bound = response,
    public_offset_bound = offset,
    public_weight_bound = weight,
    public_magnitude_bound = public_bound,
    public_product_bound = product_bound,
    public_magnitude_log2_bound = public_log2_bound,
    public_product_log2_bound = product_log2_bound,
    public_linear_predictor_bound = eta,
    public_approximation_intermediate_bound = approximation_intermediate,
    approximation_domain = domain,
    approximation_error = max(vapply(policies, function(policy) {
      policy$approximation[[family]]$max_abs_error
    }, numeric(1L))),
    numeric_error_budget = min(vapply(policies, function(policy) {
      policy$bounds$max_numeric_error
    }, numeric(1L))),
    limits = limits
  )
}

.dsvert_numeric_magnitude_bits <- function(value, log2_bound = NULL) {
  if (is.null(log2_bound)) {
    if (!is.finite(value)) return(Inf)
    if (value <= 1) return(0)
    log2_bound <- log2(value)
  }
  if (!is.numeric(log2_bound) || length(log2_bound) != 1L ||
      is.na(log2_bound) || log2_bound == Inf) return(Inf)
  if (log2_bound <= 0) return(0)
  # A tiny upward guard makes a floating log-domain bound conservative at an
  # exact power-of-two boundary.  At worst it costs one extra ring bit.
  guard <- 16 * .Machine$double.eps * max(1, abs(log2_bound))
  ceiling(log2_bound + guard)
}

.dsvert_numeric_quantization_error <- function(path_terms, frac_bits) {
  if (!is.finite(path_terms) || path_terms <= 0 ||
      !is.finite(frac_bits) || frac_bits < 0) return(Inf)
  value <- 2^(log2(path_terms) - frac_bits)
  if (value == 0) .Machine$double.xmin * .Machine$double.eps else value
}

.dsvert_numeric_output_required_bits <- function(public_bound, product_bound,
                                                 frac_bits,
                                                 public_log2_bound = NULL,
                                                 product_log2_bound = NULL) {
  2 + frac_bits + .dsvert_numeric_magnitude_bits(
    max(public_bound, product_bound),
    if (is.null(public_log2_bound) || is.null(product_log2_bound)) NULL
    else max(public_log2_bound, product_log2_bound))
}

.dsvert_numeric_raw_product_required_bits <- function(product_bound,
                                                      frac_bits,
                                                      product_log2_bound = NULL) {
  2 + 2 * frac_bits + .dsvert_numeric_magnitude_bits(
    product_bound, product_log2_bound)
}

.dsvert_numeric_required_bits <- function(public_bound, product_bound,
                                          frac_bits,
                                          public_log2_bound = NULL,
                                          product_log2_bound = NULL) {
  max(
    .dsvert_numeric_output_required_bits(
      public_bound, product_bound, frac_bits,
      public_log2_bound, product_log2_bound),
    .dsvert_numeric_raw_product_required_bits(
      product_bound, frac_bits, product_log2_bound)
  )
}

.dsvert_numeric_integer_or_infinity <- function(value) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      value < 0 || (is.finite(value) && value != floor(value))) {
    return(NA_real_)
  }
  if (value == Inf || value > .Machine$integer.max) return(as.numeric(value))
  as.integer(value)
}

.dsvert_numeric_integer_or_na <- function(value) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value != floor(value) || value < -.Machine$integer.max ||
      value > .Machine$integer.max) {
    return(NA_integer_)
  }
  as.integer(value)
}

.dsvert_numeric_raw_product_limit <- function(ring_bits, frac_bits,
                                               planner_headroom = FALSE) {
  if (!is.finite(ring_bits) || !is.finite(frac_bits)) return(NA_real_)
  exponent <- ring_bits -
    (if (isTRUE(planner_headroom)) 2 else 1) - 2 * frac_bits
  if (exponent > log2(.Machine$double.xmax)) return(Inf)
  2^exponent
}

.dsvert_numeric_truncated_output_limit <- function(
    ring_bits, frac_bits, planner_headroom = FALSE) {
  if (!is.finite(ring_bits) || !is.finite(frac_bits)) return(NA_real_)
  exponent <- ring_bits -
    (if (isTRUE(planner_headroom)) 2 else 1) - frac_bits
  if (exponent > log2(.Machine$double.xmax)) return(Inf)
  2^exponent
}

.dsvert_numeric_backend_order <- function(requested_backend, requested_ring,
                                          family) {
  if (identical(requested_backend, "auto")) {
    first <- if (requested_ring == 127L || !identical(family, "gaussian")) {
      "ring127"
    } else {
      "ring63"
    }
  } else {
    first <- requested_backend
  }
  switch(first,
    ring63 = c("ring63", "ring127", "exact_gc", "multiprecision"),
    ring127 = c("ring127", "exact_gc", "multiprecision"),
    exact_gc = c("exact_gc", "multiprecision"),
    multiprecision = "multiprecision")
}

.dsvert_numeric_capability_consensus <- function(policies, backend) {
  capabilities <- lapply(policies, function(policy) {
    policy$capabilities[[backend]]
  })
  all_true <- function(field) all(vapply(
    capabilities, function(capability) isTRUE(capability[[field]]),
    logical(1L)))
  fixed_number <- function(field) {
    values <- vapply(capabilities, function(capability) {
      as.numeric(capability[[field]])[1L]
    }, numeric(1L))
    if (any(!is.finite(values)) || length(unique(values)) != 1L) NA_real_
    else values[[1L]]
  }
  supported <- Reduce(intersect, lapply(capabilities, function(capability) {
    as.integer(capability$supported_ring_bits)
  }))
  list(
    available = all_true("available"),
    allowed = all_true("allowed"),
    e2e_verified = all_true("e2e_verified"),
    canonical_encoding = all_true("canonical_encoding"),
    fail_closed_overflow = all_true("fail_closed_overflow"),
    runtime_bounds_enforced = all_true("runtime_bounds_enforced"),
    workload_adapter_e2e_verified =
      all_true("workload_adapter_e2e_verified"),
    public_scalar_mul_truncate_e2e_verified =
      all_true("public_scalar_mul_truncate_e2e_verified"),
    full_iteration_e2e_verified = all_true("full_iteration_e2e_verified"),
    exact_truncation = all_true("exact_truncation"),
    exact_comparison = all_true("exact_comparison"),
    ring_bits = fixed_number("ring_bits"),
    frac_bits = fixed_number("frac_bits"),
    supported_ring_bits = sort(unique(supported)),
    max_ring_bits = min(vapply(capabilities, function(capability) {
      value <- as.numeric(capability$max_ring_bits)[1L]
      if (is.finite(value)) value else Inf
    }, numeric(1L))),
    max_frac_bits = min(vapply(capabilities, function(capability) {
      value <- as.numeric(capability$max_frac_bits)[1L]
      if (is.finite(value)) value else Inf
    }, numeric(1L))),
    truncation_semantics = unique(vapply(capabilities, function(capability) {
      capability$truncation_semantics
    }, character(1L)))
  )
}

.dsvert_numeric_candidate <- function(backend, policies, workload, bounds,
                                      requested_ring) {
  capability <- .dsvert_numeric_capability_consensus(policies, backend)
  reasons <- character()
  required_flags <- c(
    "available", "allowed", "e2e_verified", "canonical_encoding",
    "fail_closed_overflow", "runtime_bounds_enforced",
    "workload_adapter_e2e_verified",
    "public_scalar_mul_truncate_e2e_verified",
    "full_iteration_e2e_verified")
  for (field in required_flags) {
    if (!isTRUE(capability[[field]])) reasons <- c(reasons, field)
  }
  if (isTRUE(workload$requires_exact_truncation) &&
      !isTRUE(capability$exact_truncation)) {
    reasons <- c(reasons, "exact_truncation")
  }
  if (isTRUE(workload$requires_exact_comparison) &&
      !isTRUE(capability$exact_comparison)) {
    reasons <- c(reasons, "exact_comparison")
  }
  if (length(capability$truncation_semantics) != 1L ||
      (isTRUE(workload$requires_exact_truncation) &&
       !capability$truncation_semantics %in% c("floor", "toward_zero"))) {
    reasons <- c(reasons, "truncation_semantics")
  }

  minimum_ring <- if (!identical(workload$family, "gaussian") ||
                      requested_ring == 127L) 127L else 63L
  residual_budget <- bounds$numeric_error_budget -
    bounds$approximation_error
  path_terms <- workload$quantized_terms_per_path +
    workload$max_sequential_truncations
  required_frac_bits <- if (is.finite(residual_budget) &&
                            residual_budget > 0 &&
                            is.finite(path_terms) && path_terms > 0) {
    log2_requirement <- log2(path_terms) - log2(residual_budget)
    precision <- max(0, ceiling(log2_requirement))
    while (.dsvert_numeric_quantization_error(path_terms, precision) >
           residual_budget) {
      precision <- precision + 1
    }
    while (precision > 0 &&
           .dsvert_numeric_quantization_error(path_terms, precision - 1) <=
           residual_budget) {
      precision <- precision - 1
    }
    precision
  } else {
    Inf
  }
  ring_bits <- frac_bits <- NA_real_
  if (backend %in% c("ring63", "ring127")) {
    ring_bits <- capability$ring_bits
    frac_bits <- capability$frac_bits
  } else if (identical(backend, "exact_gc")) {
    choices <- sort(capability$supported_ring_bits[
      capability$supported_ring_bits >= minimum_ring &
      capability$supported_ring_bits <= 4096L &
      capability$supported_ring_bits <= capability$max_ring_bits])
    max_frac_bits <- min(capability$max_frac_bits, 4095L)
    if (!is.finite(required_frac_bits)) {
      reasons <- c(reasons, "numeric_error_budget")
    } else if (required_frac_bits > max_frac_bits) {
      reasons <- c(reasons, "fractional_precision")
    }
    if (length(choices) && is.finite(required_frac_bits)) {
      choice_fracs <- vapply(choices, function(choice) {
        fast_frac <- if (choice == 63L) 20L else 50L
        max(fast_frac, required_frac_bits)
      }, numeric(1L))
      choice_ok <- vapply(seq_along(choices), function(index) {
        choice <- choices[[index]]
        choice_frac <- choice_fracs[[index]]
        # The direct-wide exact circuit computes the product in an internal
        # 2*w type.  Its wire ring therefore needs to contain the exactly
        # truncated result, not the untruncated 2*f product.
        choice_required <- .dsvert_numeric_output_required_bits(
          bounds$public_magnitude_bound, bounds$public_product_bound,
          choice_frac, bounds$public_magnitude_log2_bound,
          bounds$public_product_log2_bound)
        choice_error <- .dsvert_numeric_quantization_error(
          path_terms, choice_frac) + bounds$approximation_error
        choice_frac < choice &&
          choice_frac <= max_frac_bits &&
          choice_required <= choice &&
          choice_error <= bounds$numeric_error_budget
      }, logical(1L))
      selected_index <- if (any(choice_ok)) {
        which(choice_ok)[[1L]]
      } else {
        length(choices)
      }
      ring_bits <- choices[[selected_index]]
      frac_bits <- choice_fracs[[selected_index]]
      if (!any(choice_ok)) reasons <- c(reasons, "ring_width")
    }
  } else {
    if (!is.finite(required_frac_bits)) {
      reasons <- c(reasons, "numeric_error_budget")
    } else {
      frac_bits <- max(50, required_frac_bits)
      if (frac_bits > capability$max_frac_bits) {
        reasons <- c(reasons, "fractional_precision")
      }
      required <- .dsvert_numeric_required_bits(
        bounds$public_magnitude_bound, bounds$public_product_bound, frac_bits,
        bounds$public_magnitude_log2_bound,
        bounds$public_product_log2_bound)
      ring_bits <- max(minimum_ring, required)
      if (ring_bits > capability$max_ring_bits) {
        reasons <- c(reasons, "ring_width")
      }
    }
  }

  output_required_bits <- if (is.finite(frac_bits)) {
    .dsvert_numeric_output_required_bits(
      bounds$public_magnitude_bound, bounds$public_product_bound, frac_bits,
      bounds$public_magnitude_log2_bound,
      bounds$public_product_log2_bound)
  } else {
    Inf
  }
  raw_required_bits <- if (is.finite(frac_bits)) {
    .dsvert_numeric_raw_product_required_bits(
      bounds$public_product_bound, frac_bits,
      bounds$public_product_log2_bound)
  } else {
    Inf
  }
  required_bits <- if (identical(backend, "exact_gc")) {
    output_required_bits
  } else {
    max(output_required_bits, raw_required_bits)
  }
  if (!is.finite(ring_bits) || required_bits > ring_bits) {
    reasons <- c(reasons, "ring_width")
  }
  raw_product_limit <- .dsvert_numeric_raw_product_limit(
    ring_bits, frac_bits, planner_headroom = FALSE)
  raw_product_fits <- is.finite(ring_bits) && is.finite(frac_bits) &&
    !is.na(bounds$public_product_log2_bound) &&
    bounds$public_product_log2_bound != Inf &&
    bounds$public_product_log2_bound < ring_bits - 1 - 2 * frac_bits
  if (!identical(backend, "exact_gc") && !isTRUE(raw_product_fits)) {
    reasons <- c(reasons, "raw_product_capacity")
  }
  product_headroom_proof <- if (identical(backend, "exact_gc") &&
                                !isTRUE(raw_product_fits)) {
    "truncated_output_headroom_direct_wide"
  } else {
    "raw_product_headroom"
  }
  product_headroom_verified <- is.finite(ring_bits) &&
    required_bits <= ring_bits &&
    (identical(product_headroom_proof,
               "truncated_output_headroom_direct_wide") ||
     isTRUE(raw_product_fits))
  quantization_error <- .dsvert_numeric_quantization_error(
    path_terms, frac_bits)
  total_error <- quantization_error + bounds$approximation_error
  if (!is.finite(total_error) || total_error > bounds$numeric_error_budget) {
    reasons <- c(reasons, "numeric_error_budget")
  }

  list(
    backend = backend,
    eligible = !length(unique(reasons)),
    reasons = unique(reasons),
    ring_bits = .dsvert_numeric_integer_or_na(ring_bits),
    frac_bits = .dsvert_numeric_integer_or_na(frac_bits),
    required_frac_bits =
      .dsvert_numeric_integer_or_infinity(required_frac_bits),
    required_ring_bits =
      .dsvert_numeric_integer_or_infinity(required_bits),
    truncated_output_required_ring_bits =
      .dsvert_numeric_integer_or_infinity(output_required_bits),
    raw_product_required_ring_bits =
      .dsvert_numeric_integer_or_infinity(raw_required_bits),
    supported_ring_ceiling = if (identical(backend, "exact_gc")) {
      min(capability$max_ring_bits, 4096L)
    } else {
      capability$max_ring_bits
    },
    raw_ring_product_limit_exclusive = .dsvert_numeric_raw_product_limit(
      ring_bits, frac_bits, planner_headroom = FALSE),
    truncated_output_limit_exclusive =
      .dsvert_numeric_truncated_output_limit(
        ring_bits, frac_bits, planner_headroom = FALSE),
    planner_product_limit_exclusive = if (identical(backend, "exact_gc")) {
      .dsvert_numeric_truncated_output_limit(
        ring_bits, frac_bits, planner_headroom = TRUE)
    } else {
      .dsvert_numeric_raw_product_limit(
        ring_bits, frac_bits, planner_headroom = TRUE)
    },
    public_product_within_raw_ring = isTRUE(raw_product_fits),
    product_headroom_proof = product_headroom_proof,
    product_headroom_verified = isTRUE(product_headroom_verified),
    workload_adapter_e2e_verified =
      capability$workload_adapter_e2e_verified,
    public_scalar_mul_truncate_e2e_verified =
      capability$public_scalar_mul_truncate_e2e_verified,
    full_iteration_e2e_verified = capability$full_iteration_e2e_verified,
    quantization_error = quantization_error,
    approximation_error = bounds$approximation_error,
    total_error = total_error,
    truncation_semantics = if (length(capability$truncation_semantics) == 1L) {
      capability$truncation_semantics
    } else {
      "inconsistent"
    }
  )
}

.dsvert_numeric_certificate <- function(status, workload, requested_backend,
                                         requested_ring, bounds = NULL,
                                         selected = NULL,
                                         failed_candidate = NULL,
                                         policies = NULL,
                                         attempts = NULL, reason = NULL) {
  diagnostic <- if (is.null(selected)) failed_candidate else selected
  certificate <- list(
    schema_version = 1L,
    contract_version = "dsvert-numeric-certificate-v1",
    status = status,
    workload = workload$workload,
    family = workload$family,
    requested_backend = requested_backend,
    effective_backend = if (is.null(selected)) NA_character_ else selected$backend,
    failed_backend = if (is.null(failed_candidate)) {
      NA_character_
    } else {
      failed_candidate$backend
    },
    failed_ring_bits = if (is.null(failed_candidate)) {
      NA_integer_
    } else {
      failed_candidate$ring_bits
    },
    failed_frac_bits = if (is.null(failed_candidate)) {
      NA_integer_
    } else {
      failed_candidate$frac_bits
    },
    requested_ring = requested_ring,
    ring_bits = if (is.null(selected)) NA_integer_ else selected$ring_bits,
    frac_bits = if (is.null(selected)) NA_integer_ else selected$frac_bits,
    required_frac_bits = if (is.null(diagnostic)) {
      NA_integer_
    } else {
      diagnostic$required_frac_bits
    },
    required_ring_bits = if (is.null(diagnostic)) {
      NA_integer_
    } else {
      diagnostic$required_ring_bits
    },
    supported_ring_ceiling = if (is.null(diagnostic)) {
      NA_real_
    } else {
      diagnostic$supported_ring_ceiling
    },
    truncated_output_required_ring_bits = if (is.null(diagnostic)) {
      NA_integer_
    } else {
      diagnostic$truncated_output_required_ring_bits
    },
    raw_product_required_ring_bits = if (is.null(diagnostic)) {
      NA_integer_
    } else {
      diagnostic$raw_product_required_ring_bits
    },
    raw_ring_product_limit_exclusive = if (is.null(selected)) {
      NA_real_
    } else {
      selected$raw_ring_product_limit_exclusive
    },
    planner_product_limit_exclusive = if (is.null(selected)) {
      NA_real_
    } else {
      selected$planner_product_limit_exclusive
    },
    truncated_output_limit_exclusive = if (is.null(selected)) {
      NA_real_
    } else {
      selected$truncated_output_limit_exclusive
    },
    public_product_within_raw_ring = if (is.null(selected)) {
      FALSE
    } else {
      isTRUE(selected$public_product_within_raw_ring)
    },
    product_headroom_proof = if (is.null(selected)) {
      NA_character_
    } else {
      selected$product_headroom_proof
    },
    product_headroom_verified = if (is.null(selected)) {
      FALSE
    } else {
      isTRUE(selected$product_headroom_verified)
    },
    workload_adapter_e2e_verified = if (is.null(selected)) {
      FALSE
    } else {
      isTRUE(selected$workload_adapter_e2e_verified)
    },
    public_scalar_mul_truncate_e2e_verified = if (is.null(selected)) {
      FALSE
    } else {
      isTRUE(selected$public_scalar_mul_truncate_e2e_verified)
    },
    full_iteration_e2e_verified = if (is.null(selected)) {
      FALSE
    } else {
      isTRUE(selected$full_iteration_e2e_verified)
    },
    public_predictor_input_bound = if (is.null(bounds)) NA_real_ else bounds$public_predictor_input_bound,
    public_response_input_bound = if (is.null(bounds)) NA_real_ else bounds$public_response_input_bound,
    public_predictor_bound = if (is.null(bounds)) NA_real_ else bounds$public_predictor_bound,
    public_response_bound = if (is.null(bounds)) NA_real_ else bounds$public_response_bound,
    public_offset_bound = if (is.null(bounds)) NA_real_ else bounds$public_offset_bound,
    public_weight_bound = if (is.null(bounds)) NA_real_ else bounds$public_weight_bound,
    public_magnitude_bound = if (is.null(bounds)) NA_real_ else bounds$public_magnitude_bound,
    public_product_bound = if (is.null(bounds)) NA_real_ else bounds$public_product_bound,
    public_magnitude_log2_bound = if (is.null(bounds)) {
      NA_real_
    } else {
      bounds$public_magnitude_log2_bound
    },
    public_product_log2_bound = if (is.null(bounds)) {
      NA_real_
    } else {
      bounds$public_product_log2_bound
    },
    operation_count = workload$operation_count,
    operation_count_semantics = "conservative_workload_estimate",
    truncation_count = workload$truncation_count,
    comparison_count = workload$comparison_count,
    approximation_domain_guard = workload$approximation_domain_guard,
    max_sequential_truncations = workload$max_sequential_truncations,
    planned_quantization_error = if (is.null(selected)) NA_real_ else selected$quantization_error,
    planned_approximation_error = if (is.null(selected)) NA_real_ else selected$approximation_error,
    planned_total_numeric_error = if (is.null(selected)) NA_real_ else selected$total_error,
    quantization_error_max = NA_real_,
    approximation_error_max = NA_real_,
    total_numeric_error_max = NA_real_,
    error_bound_scope = "one_element_nonlinear_arithmetic_path",
    aggregate_output_error_bound = NA_real_,
    estimator_error_bound = NA_real_,
    numeric_error_budget = if (is.null(bounds)) NA_real_ else bounds$numeric_error_budget,
    approximation_domain = if (is.null(bounds)) NULL else bounds$approximation_domain,
    public_linear_predictor_bound = if (is.null(bounds)) NA_real_ else bounds$public_linear_predictor_bound,
    public_approximation_intermediate_bound = if (is.null(bounds)) {
      NA_real_
    } else {
      bounds$public_approximation_intermediate_bound
    },
    truncation_semantics = if (is.null(selected)) NA_character_ else selected$truncation_semantics,
    bounds_custodian_owned = !is.null(policies),
    policy_body_hash_verified = !is.null(policies),
    capability_claims_received = !is.null(policies),
    capability_claims_validated = !is.null(policies),
    capabilities_attested = FALSE,
    backend_policy_e2e_claim = identical(status, "preflight_eligible"),
    backend_e2e_verified = FALSE,
    numeric_error_bound_certified = FALSE,
    preflight_eligible = identical(status, "preflight_eligible"),
    certificate_phase = "preflight",
    runtime_input_bounds_attested = FALSE,
    runtime_intermediate_bounds_attested = FALSE,
    execution_binding_ids = character(),
    result_postcondition_verified = FALSE,
    numerically_certified = FALSE,
    estimator_unchanged = TRUE,
    estimability_status = "not_evaluated",
    estimability_reason = "information_matrix_not_evaluated",
    estimability_evaluated = FALSE,
    policy_ids = if (is.null(policies)) character() else vapply(
      policies, `[[`, character(1L), "policy_id"),
    attempts = attempts,
    reason = reason
  )
  class(certificate) <- c("dsvert_numeric_certificate", "list")
  certificate
}

.dsvert_numeric_attestation_binding <- function(
    kind, policy_id, session_id, data_name, variables, family, ring, n) {
  .dsvert_numeric_hash(list(
    schema_version = 1L,
    kind = kind,
    policy_id = policy_id,
    session_id = session_id,
    data_name = data_name,
    variables = as.character(variables),
    family = family,
    ring = as.integer(ring),
    n = as.integer(n)))
}

.dsvert_numeric_validate_attestation <- function(
    attestation, certificate, server, kind, session_id, data_name,
    variables, family, ring, n) {
  fail <- function(message) {
    .dsvert_stop_numeric(
      "numeric_backend_unavailable",
      paste0("Invalid numeric execution attestation from ", server,
             ": ", message),
      certificate = certificate,
      reason = "malformed_execution_attestation")
  }
  required <- c(
    "schema_version", "kind", "policy_id", "binding_id", "ring", "n",
    "checks", "runtime_input_bounds_enforced",
    "runtime_intermediate_bounds_enforced", "observed_extrema_released",
    "attestation_scope")
  if (!is.list(attestation) || is.null(names(attestation)) ||
      anyNA(names(attestation)) || anyDuplicated(names(attestation)) ||
      !identical(sort(names(attestation)), sort(required)) ||
      !identical(attestation$schema_version, 1L) ||
      !identical(attestation$kind, kind) ||
      !is.character(attestation$policy_id) ||
      length(attestation$policy_id) != 1L ||
      !grepl("^[0-9a-f]{64}$", attestation$policy_id) ||
      !is.character(attestation$binding_id) ||
      length(attestation$binding_id) != 1L ||
      !grepl("^[0-9a-f]{64}$", attestation$binding_id) ||
      !is.numeric(attestation$ring) || length(attestation$ring) != 1L ||
      is.na(attestation$ring) || !is.finite(attestation$ring) ||
      attestation$ring != floor(attestation$ring) ||
      !identical(as.integer(attestation$ring), as.integer(ring)) ||
      !is.numeric(attestation$n) || length(attestation$n) != 1L ||
      is.na(attestation$n) || !is.finite(attestation$n) ||
      attestation$n != floor(attestation$n) ||
      !identical(as.integer(attestation$n), as.integer(n)) ||
      !is.list(attestation$checks) || !length(attestation$checks) ||
      any(!vapply(attestation$checks, isTRUE, logical(1L))) ||
      !isTRUE(attestation$runtime_input_bounds_enforced) ||
      isTRUE(attestation$runtime_intermediate_bounds_enforced) ||
      isTRUE(attestation$observed_extrema_released) ||
      !is.character(attestation$attestation_scope) ||
      length(attestation$attestation_scope) != 1L ||
      is.na(attestation$attestation_scope) ||
      !nzchar(attestation$attestation_scope)) {
    fail("malformed fields")
  }
  if (is.null(names(certificate$policy_ids)) ||
      !server %in% names(certificate$policy_ids)) {
    fail("preflight certificate has no policy for this server")
  }
  expected_policy_id <- unname(certificate$policy_ids[server])
  if (length(expected_policy_id) != 1L || is.na(expected_policy_id) ||
      !identical(attestation$policy_id, expected_policy_id)) {
    fail("policy_id does not match this server's preflight policy")
  }
  expected_binding <- .dsvert_numeric_attestation_binding(
    kind, expected_policy_id, session_id, data_name, variables, family, ring, n)
  if (!identical(attestation$binding_id, expected_binding)) {
    fail("binding does not match the session/dataset/variables/workload")
  }
  attestation$binding_id
}

.dsvert_numeric_attach_attestations <- function(certificate, binding_ids,
                                                 all_inputs = FALSE) {
  if (!inherits(certificate, "dsvert_numeric_certificate") ||
      !is.character(binding_ids) || anyNA(binding_ids) ||
      any(!grepl("^[0-9a-f]{64}$", binding_ids))) {
    stop("invalid numeric certificate/attestation binding", call. = FALSE)
  }
  certificate$execution_binding_ids <- unique(c(
    certificate$execution_binding_ids, unname(binding_ids)))
  if (isTRUE(all_inputs)) {
    certificate$runtime_input_bounds_attested <- TRUE
    certificate$certificate_phase <- "input_attested"
  }
  certificate$runtime_intermediate_bounds_attested <- FALSE
  certificate$execution_attestation_trust <-
    "custodian assertion bound to this execution; observed extrema withheld"
  certificate
}

.dsvert_numeric_assert_eta_bound <- function(beta, intercept, certificate) {
  if (is.null(certificate)) return(invisible(TRUE))
  if (!is.numeric(beta) || any(!is.finite(beta)) ||
      !is.numeric(intercept) || length(intercept) != 1L ||
      !is.finite(intercept) ||
      !is.numeric(certificate$public_predictor_bound) ||
      length(certificate$public_predictor_bound) != 1L ||
      !is.finite(certificate$public_predictor_bound) ||
      !is.numeric(certificate$public_linear_predictor_bound) ||
      length(certificate$public_linear_predictor_bound) != 1L ||
      !is.finite(certificate$public_linear_predictor_bound) ||
      !is.numeric(certificate$public_offset_bound) ||
      length(certificate$public_offset_bound) != 1L ||
      !is.finite(certificate$public_offset_bound)) {
    .dsvert_stop_numeric(
      "numeric_backend_unavailable",
      "Cannot establish a finite public linear-predictor bound",
      certificate, "malformed_runtime_eta_bound")
  }
  eta_bound <- abs(intercept) +
    certificate$public_predictor_bound * sum(abs(beta)) +
    certificate$public_offset_bound
  if (!is.finite(eta_bound) ||
      eta_bound > certificate$public_linear_predictor_bound) {
    .dsvert_stop_numeric(
      "approximation_domain_failure",
      paste0(
        "The public coefficient/input bound permits |eta| beyond the ",
        "custodian-owned linear-predictor/approximation domain"),
      certificate, "runtime_public_eta_bound_exceeded")
  }
  invisible(TRUE)
}

.dsvert_numeric_finalize_certificate <- function(
    certificate, coefficients, converged, returned_numeric = list()) {
  vectors <- c(list(coefficients = coefficients), returned_numeric)
  invalid <- names(vectors)[vapply(vectors, function(value) {
    !is.null(value) && (!is.numeric(value) || any(is.infinite(value)) ||
      any(is.nan(value)))
  }, logical(1L))]
  if (length(invalid) || !is.numeric(coefficients) ||
      !length(coefficients) || any(!is.finite(coefficients))) {
    failed <- certificate
    failed$status <- "numeric_backend_unavailable"
    failed$reason <- "non_finite_result_postcondition"
    .dsvert_stop_numeric(
      "numeric_backend_unavailable",
      paste0("The numeric result postcondition failed",
             if (length(invalid)) paste0(": ", paste(invalid, collapse = ", "))
             else ""),
      failed, "non_finite_result_postcondition")
  }
  certificate$result_postcondition_verified <- TRUE
  certificate$result_fields_ieee_valid <- TRUE
  certificate$estimator_converged <- isTRUE(converged)
  certificate$certificate_phase <- "result_validated"
  # Input attestations and finite output checks do not prove that every hidden
  # intermediate stayed in range or establish an end-to-end estimator error
  # bound.  Those requirements remain fail-closed until the executing backend
  # supplies an operation-complete runtime attestation and certified aggregate
  # and estimator bounds for this execution.
  finite_nonnegative_scalar <- function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value >= 0
  }
  certificate$numerically_certified <-
    isTRUE(certificate$preflight_eligible) &&
    isTRUE(certificate$backend_policy_e2e_claim) &&
    isTRUE(certificate$workload_adapter_e2e_verified) &&
    isTRUE(certificate$public_scalar_mul_truncate_e2e_verified) &&
    isTRUE(certificate$full_iteration_e2e_verified) &&
    isTRUE(certificate$product_headroom_verified) &&
    isTRUE(certificate$runtime_input_bounds_attested) &&
    isTRUE(certificate$runtime_intermediate_bounds_attested) &&
    isTRUE(certificate$result_postcondition_verified) &&
    isTRUE(certificate$estimator_converged) &&
    isTRUE(certificate$estimability_evaluated) &&
    identical(certificate$estimability_status, "estimable") &&
    finite_nonnegative_scalar(certificate$planned_quantization_error) &&
    finite_nonnegative_scalar(certificate$planned_approximation_error) &&
    finite_nonnegative_scalar(certificate$planned_total_numeric_error) &&
    finite_nonnegative_scalar(certificate$numeric_error_budget) &&
    certificate$planned_total_numeric_error <=
      certificate$numeric_error_budget &&
    finite_nonnegative_scalar(certificate$aggregate_output_error_bound) &&
    finite_nonnegative_scalar(certificate$estimator_error_bound) &&
    certificate$aggregate_output_error_bound <=
      certificate$numeric_error_budget &&
    certificate$estimator_error_bound <= certificate$numeric_error_budget
  if (isTRUE(certificate$numerically_certified)) {
    certificate$status <- "certified"
    certificate$capabilities_attested <- TRUE
    certificate$backend_e2e_verified <- TRUE
    certificate$numeric_error_bound_certified <- TRUE
    certificate$quantization_error_max <-
      certificate$planned_quantization_error
    certificate$approximation_error_max <-
      certificate$planned_approximation_error
    certificate$total_numeric_error_max <-
      certificate$planned_total_numeric_error
  }
  certificate
}

.dsvert_numeric_preflight_from_policies <- function(
    policies, workload, requested_backend = "auto", requested_ring = 63L) {
  if (!is.list(policies) || !length(policies)) {
    .dsvert_stop_numeric(
      "numeric_backend_unavailable", "No custodian numeric policies returned",
      reason = "missing_custodian_policy")
  }
  if (is.null(names(policies))) names(policies) <- paste0("server", seq_along(policies))
  policies <- Map(.dsvert_numeric_validate_policy, policies, names(policies))
  bounds <- .dsvert_numeric_combine_bounds(policies, workload)
  order <- .dsvert_numeric_backend_order(
    requested_backend, requested_ring, workload$family)
  attempts <- lapply(order, .dsvert_numeric_candidate,
    policies = policies, workload = workload, bounds = bounds,
    requested_ring = requested_ring)
  names(attempts) <- order
  eligible <- which(vapply(attempts, `[[`, logical(1L), "eligible"))
  if (!length(eligible)) {
    representability_reasons <- c(
      "fractional_precision", "ring_width", "raw_product_capacity",
      "numeric_error_budget")
    unrepresentable <- which(vapply(attempts, function(attempt) {
      attempt$backend %in% c("exact_gc", "multiprecision") &&
        length(attempt$reasons) > 0L &&
        all(attempt$reasons %in% representability_reasons)
    }, logical(1L)))
    if (length(unrepresentable)) {
      failed_candidate <- attempts[[unrepresentable[[1L]]]]
      certificate <- .dsvert_numeric_certificate(
        "numeric_backend_unrepresentable", workload, requested_backend,
        requested_ring, bounds = bounds,
        failed_candidate = failed_candidate, policies = policies,
        attempts = attempts,
        reason = "certified_ring_or_precision_exhausted")
      bit_label <- function(value) {
        if (is.finite(value)) format(value, scientific = FALSE, trim = TRUE)
        else "unbounded"
      }
      .dsvert_stop_numeric_unrepresentable(
        paste0(
          "The public numeric contract requires Ring",
          bit_label(certificate$required_ring_bits), " at planned f=",
          bit_label(certificate$failed_frac_bits), " (minimum f=",
          bit_label(certificate$required_frac_bits), ")",
          ", beyond the jointly certified backend ceiling Ring",
          bit_label(certificate$supported_ring_ceiling), "."),
        certificate)
    }
    certificate <- .dsvert_numeric_certificate(
      "numeric_backend_unavailable", workload, requested_backend,
      requested_ring, bounds = bounds, policies = policies,
      attempts = attempts, reason = "no_certified_backend")
    .dsvert_stop_numeric(
      "numeric_backend_unavailable",
      "No backend jointly attested by all custodians satisfies the numeric contract",
      certificate, "no_certified_backend")
  }
  .dsvert_numeric_certificate(
    "preflight_eligible", workload, requested_backend, requested_ring,
    bounds = bounds, selected = attempts[[eligible[[1L]]]],
    policies = policies, attempts = attempts)
}

.dsvert_fetch_numeric_policies <- function(datasources,
                                           .aggregate = DSI::datashield.aggregate) {
  expected <- names(datasources)
  if (is.null(expected) || anyNA(expected) || any(!nzchar(expected)) ||
      anyDuplicated(expected)) {
    .dsvert_stop_numeric(
      "numeric_backend_unavailable",
      "Numeric preflight requires uniquely named DataSHIELD connections",
      reason = "invalid_datasource_set")
  }
  result <- .dsvert_aggregate_strict(
    conns = datasources,
    expr = call(name = "dsvertNumericPolicyDS"),
    operation = "numeric policy handshake",
    .aggregate = .aggregate)
  if (!is.list(result) || is.null(names(result)) ||
      anyNA(names(result)) || any(!nzchar(names(result))) ||
      anyDuplicated(names(result)) ||
      !identical(sort(names(result)), sort(expected))) {
    .dsvert_stop_numeric(
      "numeric_backend_unavailable",
      "The numeric-policy handshake returned an invalid server set",
      reason = "invalid_policy_server_set")
  }
  result[expected]
}

.dsvert_require_numeric_policies <- function(datasources) {
  policies <- tryCatch(
    .dsvert_fetch_numeric_policies(datasources),
    error = function(error) {
      if (inherits(error, "dsvert_numeric_condition")) stop(error)
      .dsvert_stop_numeric(
        "numeric_backend_unavailable",
        "The custodian numeric-policy endpoint is unavailable",
        reason = "custodian_policy_endpoint_unavailable")
    })
  Map(.dsvert_numeric_validate_policy, policies, names(policies))
}

#' Preflight a vertically partitioned GLM numeric backend
#'
#' Negotiates only server-published bounds and runtime capabilities.  The
#' analyst may request a preferred backend/ring, but cannot supply or relax a
#' bound.  A stronger backend is selected automatically when required.
#'
#' @param n_obs Public observation count.
#' @param n_predictors Number of predictor columns.
#' @param family GLM family.
#' @param max_iter Maximum optimizer iterations.
#' @param compute_se Whether the workload includes Hessian evaluations.
#' @param compute_deviance Whether the workload includes deviance evaluation.
#' @param weights_active Whether protected weights are active.
#' @param offset_active Whether a protected offset is active.
#' @param backend Preferred backend: auto, ring63, ring127, exact_gc, or
#'   multiprecision.
#' @param ring Preferred fast ring, 63 or 127.
#' @param datasources DataSHIELD connections.
#' @return A numeric certificate.
#' @export
ds.vertNumericPreflight <- function(
    n_obs, n_predictors, family = "gaussian", max_iter = 100L,
    compute_se = TRUE, compute_deviance = TRUE,
    weights_active = FALSE, offset_active = FALSE,
    backend = "auto", ring = 63L, datasources = NULL) {
  backend <- tolower(as.character(backend))
  if (length(backend) != 1L || is.na(backend) ||
      !backend %in% c("auto", "ring63", "ring127", "exact_gc",
                     "multiprecision")) {
    stop("backend must be auto, ring63, ring127, exact_gc, or multiprecision",
         call. = FALSE)
  }
  if (!is.numeric(ring) || length(ring) != 1L || is.na(ring) ||
      !is.finite(ring) || ring != floor(ring) ||
      !as.integer(ring) %in% c(63L, 127L)) {
    stop("ring must be 63 or 127", call. = FALSE)
  }
  ring <- as.integer(ring)
  workload <- .dsvert_numeric_glm_workload(
    n_obs = n_obs, n_predictors = n_predictors, family = family,
    max_iter = max_iter, compute_se = compute_se,
    compute_deviance = compute_deviance,
    weights_active = weights_active, offset_active = offset_active)
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  if (!length(datasources)) {
    .dsvert_stop_numeric(
      "numeric_backend_unavailable", "No DataSHIELD connections found",
      reason = "missing_datasources")
  }

  fetched <- tryCatch(
    list(ok = TRUE, value = .dsvert_fetch_numeric_policies(datasources)),
    error = function(error) {
      if (inherits(error, "dsvert_numeric_condition")) stop(error)
      list(ok = FALSE, error = error)
    })
  if (!isTRUE(fetched$ok)) {
    .dsvert_stop_numeric(
      "numeric_backend_unavailable",
      "The custodian numeric-policy endpoint is unavailable",
      reason = "custodian_policy_endpoint_unavailable")
  }

  .dsvert_numeric_preflight_from_policies(
    fetched$value, workload, requested_backend = backend,
    requested_ring = ring)
}
