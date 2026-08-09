# Pure, signature-ready compiler for the formal binomial/Poisson GLM artifact.
# It consumes public custodian policy only.  It neither contacts a server nor
# accepts observations, fitted values, gradients or analyst sensitivity fields.

.DSVERT_FORMAL_GLM_SCHEMA_VERSION <- "dsvert-formal-glm-schema-v1"
.DSVERT_FORMAL_GLM_COMPILER_VERSION <- "dsvert-formal-glm-compiler-v1"
.DSVERT_FORMAL_GLM_THEOREM_VERSION <-
  "projected-full-gradient-contraction-rational-v1"
.DSVERT_FORMAL_GLM_SIGNATURE_DOMAIN <-
  "dsVert/dp/formal-glm-schema/v1|"

.dsvert_formal_glm_spec_condition <- function(message, reason,
                                              certificate = NULL) {
  structure(
    list(message = message, call = NULL, reason = reason,
         certificate = certificate),
    class = c("dsvert_formal_glm_specification_error",
              "dsvert_numeric_condition", "error", "condition"))
}

.dsvert_formal_glm_stop <- function(message, reason,
                                    certificate = NULL) {
  stop(.dsvert_formal_glm_spec_condition(
    message, reason = reason, certificate = certificate))
}

.dsvert_formal_glm_identifier <- function(value, what) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !grepl("^[A-Za-z][A-Za-z0-9_.]{0,127}$", value)) {
    .dsvert_formal_glm_stop(
      paste0(what, " must be one bounded ASCII identifier"),
      "invalid_public_identifier")
  }
  value
}

.dsvert_formal_glm_string <- function(value, what, maximum = 256L) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value) || nchar(value, type = "bytes") > maximum ||
      grepl("\\x00", value)) {
    .dsvert_formal_glm_stop(
      paste0(what, " must be one bounded UTF-8 string"),
      "invalid_public_string")
  }
  enc2utf8(value)
}

.dsvert_formal_glm_hash <- function(value) {
  digest::digest(.dsvert_joint_dp_client_json(value), algo = "sha256",
                 serialize = FALSE)
}

.dsvert_formal_glm_rat_field <- function(value, what,
                                         lower = NULL, upper = NULL,
                                         lower_open = FALSE,
                                         upper_open = FALSE) {
  rational <- tryCatch(.dsvert_glm_rat(value), error = function(error) NULL)
  valid <- !is.null(rational)
  if (valid && !is.null(lower)) {
    comparison <- .dsvert_glm_rat_cmp(rational, lower)
    valid <- comparison > 0L || (!lower_open && comparison == 0L)
  }
  if (valid && !is.null(upper)) {
    comparison <- .dsvert_glm_rat_cmp(rational, upper)
    valid <- comparison < 0L || (!upper_open && comparison == 0L)
  }
  if (!valid) {
    .dsvert_formal_glm_stop(
      paste0("Invalid public rational field: ", what),
      "invalid_public_rational")
  }
  rational
}

.dsvert_formal_glm_integer <- function(value, what, lower, upper) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value != floor(value) || value < lower || value > upper) {
    .dsvert_formal_glm_stop(
      paste0("Invalid public integer field: ", what),
      "invalid_public_integer")
  }
  as.integer(value)
}

.dsvert_formal_glm_exact_names <- function(value, expected, what) {
  if (!is.list(value) || !identical(sort(names(value)), sort(expected))) {
    .dsvert_formal_glm_stop(
      paste0("Invalid ", what, " fields"), "invalid_public_schema")
  }
  invisible(TRUE)
}

.dsvert_formal_glm_formula <- function(formula) {
  if (is.character(formula)) {
    if (length(formula) != 1L || is.na(formula) || !nzchar(formula) ||
        nchar(formula, type = "bytes") > 4096L) {
      .dsvert_formal_glm_stop(
        "formula must be one bounded two-sided formula",
        "formula_not_certified")
    }
    formula <- tryCatch(stats::as.formula(formula, env = baseenv()),
                        error = function(error) NULL)
  }
  if (!inherits(formula, "formula") || length(formula) != 3L ||
      !is.symbol(formula[[2L]])) {
    .dsvert_formal_glm_stop(
      "Only a two-sided formula with one plain response is certified",
      "formula_not_certified")
  }
  response <- as.character(formula[[2L]])
  .dsvert_formal_glm_identifier(response, "formula response")
  terms <- tryCatch(stats::terms(formula), error = function(error) NULL)
  if (is.null(terms) || "." %in% all.names(formula[[3L]], functions = FALSE)) {
    .dsvert_formal_glm_stop(
      "Data-dependent formula expansion is not certified",
      "formula_not_certified")
  }
  labels <- attr(terms, "term.labels")
  orders <- attr(terms, "order")
  plain <- vapply(labels, function(label) {
    expression <- tryCatch(parse(text = label, keep.source = FALSE),
                           error = function(error) NULL)
    length(expression) == 1L && is.symbol(expression[[1L]]) &&
      identical(as.character(expression[[1L]]), label) &&
      grepl("^[A-Za-z][A-Za-z0-9_.]{0,127}$", label)
  }, logical(1L))
  if (any(!plain) || any(orders != 1L) || anyDuplicated(labels)) {
    .dsvert_formal_glm_stop(
      paste0("Only unique additive registered columns are certified; ",
             "interactions, transformations and splines require a distinct ",
             "artifact"),
      "formula_not_certified")
  }
  predictors <- sort(labels, method = "radix")
  intercept <- identical(as.integer(attr(terms, "intercept")), 1L)
  right <- c(if (intercept) "1" else "0", predictors)
  list(
    response = response,
    predictors = predictors,
    intercept = intercept,
    canonical = paste(response, "~", paste(right, collapse = " + ")))
}

.dsvert_formal_glm_levels <- function(value, what) {
  if (!is.character(value) || !length(value) || anyNA(value) ||
      any(!nzchar(value)) || any(nchar(value, type = "bytes") > 256L) ||
      anyDuplicated(value)) {
    .dsvert_formal_glm_stop(
      paste0("Invalid fixed categorical domain for ", what),
      "categorical_domain_not_certified")
  }
  enc2utf8(unname(value))
}

.dsvert_formal_glm_column <- function(value, name) {
  if (!is.list(value) || !"kind" %in% names(value) ||
      !is.character(value$kind) || length(value$kind) != 1L) {
    .dsvert_formal_glm_stop(
      paste0("Invalid registered column ", name), "invalid_column_schema")
  }
  kind <- value$kind
  expected <- switch(kind,
    numeric = c("kind", "owner", "lower", "upper"),
    categorical = c("kind", "owner", "levels", "reference", "contrast"),
    binary = c("kind", "owner"),
    count = c("kind", "owner", "upper"),
    weight = c("kind", "owner", "upper"),
    offset = c("kind", "owner", "lower", "upper"),
    exposure = c("kind", "owner", "lower", "upper"),
    NULL)
  if (is.null(expected) || !identical(sort(names(value)), sort(expected))) {
    .dsvert_formal_glm_stop(
      paste0("Unsupported registered column schema for ", name),
      "invalid_column_schema")
  }
  owner <- .dsvert_formal_glm_identifier(value$owner,
                                         paste0(name, " owner"))
  if (identical(kind, "categorical")) {
    levels <- .dsvert_formal_glm_levels(value$levels, name)
    reference <- .dsvert_formal_glm_string(
      value$reference, paste0(name, " reference"))
    if (!reference %in% levels || !identical(value$contrast, "treatment")) {
      .dsvert_formal_glm_stop(
        paste0("Only a registered treatment contrast is certified for ", name),
        "categorical_contrast_not_certified")
    }
    return(list(kind = kind, owner = owner, levels = as.list(levels),
                reference = reference, contrast = "treatment"))
  }
  if (identical(kind, "binary")) {
    return(list(kind = kind, owner = owner, domain = list(0, 1)))
  }
  if (kind %in% c("count", "weight")) {
    upper <- .dsvert_formal_glm_rat_field(
      value$upper, paste0(name, " upper"), lower = "0", lower_open = TRUE)
    if (identical(kind, "count") &&
        .dsvert_glm_rat_cmp(upper, .dsvert_glm_rat_round_dyadic(upper, 0L)) !=
          0L) {
      .dsvert_formal_glm_stop(
        "Poisson patient-level count caps must be integers",
        "count_domain_not_certified")
    }
    return(list(kind = kind, owner = owner,
                upper = .dsvert_glm_rat_json(upper)))
  }
  lower <- .dsvert_formal_glm_rat_field(
    value$lower, paste0(name, " lower"))
  upper <- .dsvert_formal_glm_rat_field(
    value$upper, paste0(name, " upper"))
  if (.dsvert_glm_rat_cmp(lower, upper) >= 0L ||
      (identical(kind, "exposure") &&
       .dsvert_glm_rat_cmp(lower, "0") <= 0L)) {
    .dsvert_formal_glm_stop(
      paste0("Invalid signed bounds for ", name),
      "column_bounds_not_certified")
  }
  list(kind = kind, owner = owner,
       lower = .dsvert_glm_rat_json(lower),
       upper = .dsvert_glm_rat_json(upper))
}

.dsvert_formal_glm_expand_terms <- function(parsed, columns,
                                            x_fraction_bits) {
  result <- list()
  if (parsed$intercept) {
    result[[1L]] <- list(
      coefficient = "(Intercept)", source_column = NULL,
      source_level = NULL, kind = "intercept", owner = NULL,
      lower = .dsvert_glm_rat_json("1"),
      upper = .dsvert_glm_rat_json("1"),
      abs_bound = .dsvert_glm_rat_json("1"))
  }
  for (name in parsed$predictors) {
    column <- columns[[name]]
    if (is.null(column) || !column$kind %in% c("numeric", "categorical")) {
      .dsvert_formal_glm_stop(
        paste0("Predictor ", name, " is not a registered numeric or factor ",
               "column"),
        "predictor_not_certified")
    }
    if (identical(column$kind, "numeric")) {
      clipping_lower <- .dsvert_glm_rat(column$lower)
      clipping_upper <- .dsvert_glm_rat(column$upper)
      lower <- .dsvert_glm_rat_round_dyadic(
        clipping_lower, x_fraction_bits)
      upper <- .dsvert_glm_rat_round_dyadic(
        clipping_upper, x_fraction_bits)
      bound <- .dsvert_glm_rat_max(
        .dsvert_glm_rat_abs(lower), .dsvert_glm_rat_abs(upper))
      result[[length(result) + 1L]] <- list(
        coefficient = name, source_column = name, source_level = NULL,
        kind = "numeric", owner = column$owner,
        clipping_lower = column$lower, clipping_upper = column$upper,
        quantized_lower = .dsvert_glm_rat_json(lower),
        quantized_upper = .dsvert_glm_rat_json(upper),
        abs_bound = .dsvert_glm_rat_json(bound))
    } else {
      levels <- unlist(column$levels, use.names = FALSE)
      included <- if (parsed$intercept) {
        levels[levels != column$reference]
      } else {
        levels
      }
      for (level in included) {
        result[[length(result) + 1L]] <- list(
          coefficient = paste0(name, "[", level, "]"),
          source_column = name, source_level = level,
          kind = "categorical_indicator", owner = column$owner,
          lower = .dsvert_glm_rat_json("0"),
          upper = .dsvert_glm_rat_json("1"),
          abs_bound = .dsvert_glm_rat_json("1"))
      }
    }
  }
  names(result) <- vapply(result, `[[`, character(1L), "coefficient")
  result
}

.dsvert_formal_glm_expand_parameter <- function(value, names, what,
                                                positive = TRUE) {
  if (length(value) == 1L && is.null(names(value))) {
    value <- rep(value, length(names))
    names(value) <- names
  }
  if (is.list(value)) value <- unlist(value, use.names = TRUE)
  if (is.null(names(value)) || anyDuplicated(names(value)) ||
      !setequal(names(value), names)) {
    .dsvert_formal_glm_stop(
      paste0(what, " must be scalar or named in canonical coefficient order"),
      "coefficient_parameter_not_certified")
  }
  value <- value[names]
  result <- lapply(seq_along(value), function(index) {
    .dsvert_formal_glm_rat_field(
      value[[index]], paste0(what, " for ", names[[index]]),
      lower = "0", lower_open = positive)
  })
  names(result) <- names
  result
}

.dsvert_formal_glm_link_table <- function(family, maximum_eta, segments,
                                          fraction_bits,
                                          precision_bits) {
  maximum_eta <- .dsvert_glm_rat(maximum_eta)
  width <- .dsvert_glm_rat_div(
    .dsvert_glm_rat_mul("2", maximum_eta), as.character(segments))
  knots <- vector("list", segments + 1L)
  values <- vector("list", segments + 1L)
  knot_errors <- vector("list", segments + 1L)
  true_intervals <- vector("list", segments + 1L)
  for (index in 0:segments) {
    eta <- .dsvert_glm_rat_add(
      .dsvert_glm_rat_neg(maximum_eta),
      .dsvert_glm_rat_mul(as.character(index), width))
    exponential <- .dsvert_glm_rat_exp_interval(eta, precision_bits)
    truth <- if (identical(family, "poisson")) {
      exponential
    } else {
      list(
        lower = .dsvert_glm_rat_div(
          exponential$lower,
          .dsvert_glm_rat_add("1", exponential$lower)),
        upper = .dsvert_glm_rat_div(
          exponential$upper,
          .dsvert_glm_rat_add("1", exponential$upper)))
    }
    approximation <- .dsvert_glm_rat_round_dyadic(
      .dsvert_glm_rat_interval_midpoint(truth), fraction_bits)
    error <- .dsvert_glm_rat_max(
      .dsvert_glm_rat_abs(
        .dsvert_glm_rat_sub(approximation, truth$lower)),
      .dsvert_glm_rat_abs(
        .dsvert_glm_rat_sub(truth$upper, approximation)))
    knots[[index + 1L]] <- eta
    values[[index + 1L]] <- approximation
    knot_errors[[index + 1L]] <- error
    true_intervals[[index + 1L]] <- truth
  }
  if (any(vapply(seq_len(segments), function(index) {
    .dsvert_glm_rat_cmp(values[[index]], values[[index + 1L]]) > 0L
  }, logical(1L)))) {
    .dsvert_formal_glm_stop(
      "The certified link grid is not monotone", "link_grid_not_certified")
  }
  if (identical(family, "poisson") &&
      .dsvert_glm_rat_cmp(values[[1L]], "0") <= 0L) {
    .dsvert_formal_glm_stop(
      "The Poisson link precision rounds its lower endpoint to zero",
      "link_grid_not_certified")
  }
  slopes <- lapply(seq_len(segments), function(index) {
    .dsvert_glm_rat_div(
      .dsvert_glm_rat_sub(values[[index + 1L]], values[[index]]), width)
  })
  maximum_slope <- Reduce(.dsvert_glm_rat_max, slopes)
  maximum_knot_error <- Reduce(.dsvert_glm_rat_max, knot_errors)
  second_derivative <- if (identical(family, "binomial")) {
    .dsvert_glm_rat(list(numerator = "1", denominator = "4"))
  } else {
    true_intervals[[length(true_intervals)]]$upper
  }
  interpolation_error <- .dsvert_glm_rat_div(
    .dsvert_glm_rat_mul(
      .dsvert_glm_rat_mul(width, width), second_derivative), "8")
  error <- .dsvert_glm_rat_add(maximum_knot_error, interpolation_error)
  table <- list(
    version = "monotone-dyadic-pwl-integrated-score-v1",
    family = family,
    domain_lower = .dsvert_glm_rat_json(.dsvert_glm_rat_neg(maximum_eta)),
    domain_upper = .dsvert_glm_rat_json(maximum_eta),
    segment_count = as.numeric(segments),
    selection_topology = "oblivious_balanced_exact_gc_v1",
    fraction_bits = as.numeric(fraction_bits),
    construction_precision_bits = as.numeric(precision_bits),
    knots = lapply(knots, .dsvert_glm_rat_json),
    values = lapply(values, .dsvert_glm_rat_json),
    slopes = lapply(slopes, .dsvert_glm_rat_json),
    range_lower = .dsvert_glm_rat_json(values[[1L]]),
    range_upper = .dsvert_glm_rat_json(values[[length(values)]]),
    max_nonnegative_slope = .dsvert_glm_rat_json(maximum_slope),
    uniform_mean_error_upper = .dsvert_glm_rat_json(error),
    error_proof = list(
      endpoint_interval_method = "exact-rational-taylor-with-geometric-tail",
      maximum_knot_error = .dsvert_glm_rat_json(maximum_knot_error),
      interpolation_remainder = .dsvert_glm_rat_json(interpolation_error),
      second_derivative_upper = .dsvert_glm_rat_json(second_derivative),
      continuity = TRUE, monotone = TRUE))
  table$sha256 <- .dsvert_formal_glm_hash(table)
  table
}

.dsvert_formal_glm_sum <- function(values) {
  Reduce(.dsvert_glm_rat_add, values, init = .dsvert_glm_rat("0"))
}

.dsvert_formal_glm_max <- function(values) {
  Reduce(.dsvert_glm_rat_max, values)
}

.dsvert_formal_glm_numeric_plan <- function(bounds, numeric) {
  maximum <- .dsvert_formal_glm_max(lapply(bounds, .dsvert_glm_rat))
  scaled <- .dsvert_glm_rat_ceil_abs_scaled(
    maximum, numeric$working_fraction_bits)
  digits <- nchar(as.character(scaled), type = "bytes")
  required <- as.integer(ceiling(digits * log2(10))) + 3L
  backend <- if (required <= 63L) {
    "Ring63"
  } else if (required <= 127L) {
    "Ring127"
  } else if (required <= 4096L) {
    paste0("Ring", required)
  } else {
    certificate <- list(
      status = "numeric_backend_unrepresentable",
      required_ring_bits = required,
      required_frac_bits = numeric$working_fraction_bits,
      numeric_error_budget = "reference_rational_exact")
    .dsvert_stop_numeric_unrepresentable(
      "The public formal GLM interval plan exceeds Ring4096", certificate)
  }
  list(
    version = "formal-glm-public-interval-plan-v1",
    reference_backend = "openssl-bignum-exact-rational",
    production_backend_candidate = backend,
    required_signed_ring_bits_upper = as.numeric(required),
    bit_bound_method = "decimal-digit upper bound with three guard bits",
    working_fraction_bits = as.numeric(numeric$working_fraction_bits),
    reference_update_error_rho = .dsvert_glm_rat_json("0"),
    production_update_error_rho_status = "phase1_exact_dag_required",
    production_release_ready = FALSE,
    interval_nodes = bounds)
}

.dsvert_formal_glm_compile <- function(model, authority) {
  required_model <- c(
    "formula", "family", "capacity", "columns", "adjacency",
    "patient_collapse", "missingness", "clipping", "weights", "offset",
    "coefficient_box", "ridge", "optimizer", "numeric", "link", "privacy")
  .dsvert_formal_glm_exact_names(model, required_model,
                                 "formal GLM model")
  .dsvert_formal_glm_exact_names(
    authority,
    c("consortium_id", "capsule_id", "logical_snapshot",
      "policy_contract_sha256", "custodian_peers", "designated_peers",
      "pinset_sha256"),
    "formal GLM authority")
  if (!is.character(model$family) || length(model$family) != 1L ||
      is.na(model$family) || !model$family %in% c("binomial", "poisson")) {
    .dsvert_formal_glm_stop(
      "Only binomial/logit and Poisson/log are certified",
      "family_not_certified")
  }
  family <- model$family
  parsed <- .dsvert_formal_glm_formula(model$formula)
  capacity <- .dsvert_formal_glm_integer(
    model$capacity, "capacity", 1L, 100000000L)
  if (!is.character(model$adjacency) || length(model$adjacency) != 1L ||
      is.na(model$adjacency) ||
      !model$adjacency %in% c("add_remove", "replace_one")) {
    .dsvert_formal_glm_stop(
      "Adjacency must be add_remove or replace_one",
      "adjacency_not_certified")
  }
  adjacency <- model$adjacency
  gamma <- if (identical(adjacency, "add_remove")) 1L else 2L

  if (!is.list(model$columns) || is.null(names(model$columns)) ||
      anyNA(names(model$columns)) || anyDuplicated(names(model$columns))) {
    .dsvert_formal_glm_stop(
      "columns must be one uniquely named custodian registry",
      "invalid_column_schema")
  }
  column_names <- sort(names(model$columns), method = "radix")
  if (any(!vapply(column_names, function(name) {
    grepl("^[A-Za-z][A-Za-z0-9_.]{0,127}$", name)
  }, logical(1L)))) {
    .dsvert_formal_glm_stop(
      "Every registered column needs a bounded ASCII identifier",
      "invalid_column_schema")
  }
  columns <- lapply(column_names, function(name) {
    .dsvert_formal_glm_column(model$columns[[name]], name)
  })
  names(columns) <- column_names
  response <- columns[[parsed$response]]
  expected_response <- if (identical(family, "binomial")) "binary" else "count"
  if (is.null(response) || !identical(response$kind, expected_response)) {
    .dsvert_formal_glm_stop(
      "The registered outcome domain does not match the GLM family",
      "outcome_domain_not_certified")
  }

  .dsvert_formal_glm_exact_names(
    model$patient_collapse,
    c("unit", "repeated_records", "row_order_invariant",
      "max_records_per_unit", "conflict_policy"),
    "patient collapse")
  collapse <- model$patient_collapse
  collapse$max_records_per_unit <- .dsvert_formal_glm_integer(
    collapse$max_records_per_unit, "max_records_per_unit", 1L, 1000000L)
  if (!identical(collapse$unit, "aligned_patient") ||
      !identical(collapse$repeated_records, "reject_duplicates") ||
      !identical(collapse$row_order_invariant, TRUE) ||
      collapse$max_records_per_unit != 1L ||
      !identical(collapse$conflict_policy, "zero_weight")) {
    .dsvert_formal_glm_stop(
      paste0("Phase 0 certifies one aligned record per patient; duplicates ",
             "or cross-source conflicts map to the zero-weight tuple"),
      "patient_collapse_not_certified")
  }
  if (!identical(model$missingness, "complete_tuple_zero_weight")) {
    .dsvert_formal_glm_stop(
      "Only fixed-shape complete-tuple zero-weight missingness is certified",
      "missingness_not_certified")
  }
  required_clipping <- list(
    numeric = "clamp_then_quantize",
    categorical = "registered_level_or_zero_weight",
    binomial = "binary_or_zero_weight",
    poisson_count = "integer_then_patient_cap",
    offset = "clamp_then_quantize",
    weight = "clamp_then_quantize")
  if (!is.list(model$clipping) ||
      !identical(sort(names(model$clipping)), sort(names(required_clipping))) ||
      !identical(model$clipping[names(required_clipping)], required_clipping)) {
    .dsvert_formal_glm_stop(
      "The clipping contract is not the certified fixed policy",
      "clipping_not_certified")
  }

  .dsvert_formal_glm_exact_names(
    model$numeric,
    c("x_fraction_bits", "offset_fraction_bits", "beta_fraction_bits",
      "link_fraction_bits", "working_fraction_bits",
      "reference_precision_bits"),
    "numeric policy")
  numeric <- lapply(names(model$numeric), function(name) {
    .dsvert_formal_glm_integer(model$numeric[[name]], name, 0L, 256L)
  })
  names(numeric) <- names(model$numeric)
  if (numeric$reference_precision_bits < 64L ||
      numeric$working_fraction_bits < max(
        numeric$x_fraction_bits, numeric$offset_fraction_bits,
        numeric$beta_fraction_bits, numeric$link_fraction_bits)) {
    .dsvert_formal_glm_stop(
      "The numeric grid does not meet the reference precision contract",
      "numeric_grid_not_certified")
  }

  .dsvert_formal_glm_exact_names(model$weights, c("mode", "column"),
                                 "weight policy")
  if (!model$weights$mode %in% c("unit", "bounded_column")) {
    .dsvert_formal_glm_stop(
      "Unsupported analysis-weight policy", "weights_not_certified")
  }
  weight_column <- model$weights$column
  if (identical(model$weights$mode, "unit")) {
    if (!is.null(weight_column)) {
      .dsvert_formal_glm_stop(
        "Unit weights cannot name a source column", "weights_not_certified")
    }
    source_maximum_weight <- maximum_weight <- .dsvert_glm_rat("1")
  } else {
    weight_column <- .dsvert_formal_glm_identifier(
      weight_column, "weight column")
    if (is.null(columns[[weight_column]]) ||
        !identical(columns[[weight_column]]$kind, "weight")) {
      .dsvert_formal_glm_stop(
        "The bounded weight column is not registered",
        "weights_not_certified")
    }
    source_maximum_weight <- .dsvert_glm_rat(
      columns[[weight_column]]$upper)
    maximum_weight <- .dsvert_glm_rat_round_dyadic(
      source_maximum_weight, numeric$x_fraction_bits)
    if (.dsvert_glm_rat_cmp(maximum_weight, "0") <= 0L) {
      .dsvert_formal_glm_stop(
        "The bounded patient weight vanishes on the signed numeric grid",
        "weights_not_certified")
    }
  }

  .dsvert_formal_glm_exact_names(model$offset, c("mode", "column"),
                                 "offset policy")
  if (!model$offset$mode %in% c("none", "bounded_offset", "log_exposure")) {
    .dsvert_formal_glm_stop(
      "Unsupported offset/exposure policy", "offset_not_certified")
  }
  offset_column <- model$offset$column
  offset_error <- .dsvert_glm_rat("0")
  if (identical(model$offset$mode, "none")) {
    if (!is.null(offset_column)) {
      .dsvert_formal_glm_stop(
        "A model without an offset cannot name an offset column",
        "offset_not_certified")
    }
    offset_source_lower <- offset_source_upper <-
      offset_lower <- offset_upper <- .dsvert_glm_rat("0")
  } else {
    offset_column <- .dsvert_formal_glm_identifier(
      offset_column, "offset/exposure column")
    expected_kind <- if (identical(model$offset$mode, "bounded_offset")) {
      "offset"
    } else {
      "exposure"
    }
    column <- columns[[offset_column]]
    if (is.null(column) || !identical(column$kind, expected_kind)) {
      .dsvert_formal_glm_stop(
        "The offset/exposure column is not registered with matching semantics",
        "offset_not_certified")
    }
    if (identical(expected_kind, "offset")) {
      offset_source_lower <- .dsvert_glm_rat(column$lower)
      offset_source_upper <- .dsvert_glm_rat(column$upper)
      offset_lower <- .dsvert_glm_rat_round_dyadic(
        offset_source_lower, numeric$offset_fraction_bits)
      offset_upper <- .dsvert_glm_rat_round_dyadic(
        offset_source_upper, numeric$offset_fraction_bits)
      offset_error <- .dsvert_glm_rat_div(
        "1", .dsvert_glm_rat_pow(
          "2", numeric$offset_fraction_bits + 1L))
    } else {
      offset_source_lower <- .dsvert_glm_rat(column$lower)
      offset_source_upper <- .dsvert_glm_rat(column$upper)
      intervals <- tryCatch(list(
        lower = .dsvert_glm_rat_log_interval(
          .dsvert_glm_rat(column$lower), numeric$reference_precision_bits),
        upper = .dsvert_glm_rat_log_interval(
          .dsvert_glm_rat(column$upper), numeric$reference_precision_bits)),
        error = function(error) NULL)
      if (is.null(intervals)) {
        .dsvert_stop_numeric_unrepresentable(
          "The registered exposure range exceeds the Phase-0 log oracle",
          list(status = "numeric_backend_unrepresentable",
               required_ring_bits = NA_integer_,
               required_frac_bits = numeric$offset_fraction_bits,
               numeric_error_budget =
                 "bounded high-precision exposure transform required"))
      }
      lower_interval <- intervals$lower
      upper_interval <- intervals$upper
      offset_lower <- .dsvert_glm_rat_round_dyadic(
        lower_interval$lower, numeric$offset_fraction_bits)
      offset_upper <- .dsvert_glm_rat_round_dyadic(
        upper_interval$upper, numeric$offset_fraction_bits)
      offset_error <- .dsvert_glm_rat_add(
        .dsvert_glm_rat_div(
          "1", .dsvert_glm_rat_pow(
            "2", numeric$offset_fraction_bits + 1L)),
        .dsvert_glm_rat_div(
          "1", .dsvert_glm_rat_pow(
            "2", numeric$reference_precision_bits)))
    }
  }

  artifact_column_names <- sort(unique(c(
    parsed$response, parsed$predictors, weight_column, offset_column)),
    method = "radix")
  artifact_columns <- columns[artifact_column_names]
  terms <- .dsvert_formal_glm_expand_terms(
    parsed, artifact_columns, numeric$x_fraction_bits)
  coefficient_names <- names(terms)
  if (!length(coefficient_names)) {
    .dsvert_formal_glm_stop(
      "The compiled design has no coefficients", "formula_not_certified")
  }
  boxes <- .dsvert_formal_glm_expand_parameter(
    model$coefficient_box, coefficient_names, "coefficient_box")
  ridge <- .dsvert_formal_glm_expand_parameter(
    model$ridge, coefficient_names, "ridge")

  .dsvert_formal_glm_exact_names(
    model$optimizer, c("alpha", "iterations", "start"), "optimizer")
  if (!identical(model$optimizer$start, "zero")) {
    .dsvert_formal_glm_stop(
      "Only the public zero initial coefficient vector is certified",
      "optimizer_not_certified")
  }
  alpha <- .dsvert_formal_glm_rat_field(
    model$optimizer$alpha, "optimizer alpha", lower = "0", lower_open = TRUE)
  iterations <- .dsvert_formal_glm_integer(
    model$optimizer$iterations, "optimizer iterations", 1L, 4096L)
  .dsvert_formal_glm_exact_names(
    model$link, c("segments", "construction"), "link policy")
  if (!identical(model$link$construction,
                 "uniform_monotone_dyadic_linear")) {
    .dsvert_formal_glm_stop(
      "The link construction is not certified", "link_grid_not_certified")
  }
  segments <- .dsvert_formal_glm_integer(
    model$link$segments, "link segments", 2L, 256L)

  a_values <- lapply(terms, function(term) .dsvert_glm_rat(term$abs_bound))
  r_squared <- .dsvert_formal_glm_sum(lapply(a_values, function(value) {
    .dsvert_glm_rat_mul(value, value)
  }))
  r_upper <- .dsvert_formal_glm_sum(a_values)
  maximum_offset <- .dsvert_glm_rat_max(
    .dsvert_glm_rat_abs(offset_lower),
    .dsvert_glm_rat_abs(offset_upper))
  maximum_eta <- .dsvert_glm_rat_add(
    maximum_offset,
    .dsvert_formal_glm_sum(Map(.dsvert_glm_rat_mul, a_values, boxes)))
  if (!is.finite(.dsvert_glm_rat_double(maximum_eta)) ||
      .dsvert_glm_rat_double(maximum_eta) > 64) {
    .dsvert_stop_numeric_unrepresentable(
      "The public link domain exceeds the Phase-0 high-precision oracle",
      list(status = "numeric_backend_unrepresentable",
           required_ring_bits = NA_integer_,
           required_frac_bits = numeric$link_fraction_bits,
           numeric_error_budget = "M_eta must not exceed 64 in Phase 0"))
  }
  link <- .dsvert_formal_glm_link_table(
    family, maximum_eta, segments, numeric$link_fraction_bits,
    numeric$reference_precision_bits)

  lambda_min <- Reduce(.dsvert_glm_rat_min, ridge)
  lambda_max <- Reduce(.dsvert_glm_rat_max, ridge)
  slope <- .dsvert_glm_rat(link$max_nonnegative_slope)
  smoothness <- .dsvert_glm_rat_add(
    lambda_max,
    .dsvert_glm_rat_mul(
      .dsvert_glm_rat_mul(maximum_weight, r_squared), slope))
  step_limit <- .dsvert_glm_rat_div(
    "2", .dsvert_glm_rat_add(lambda_min, smoothness))
  if (.dsvert_glm_rat_cmp(alpha, step_limit) > 0L) {
    .dsvert_formal_glm_stop(
      "The public step exceeds 2/(m+L)", "optimizer_not_contractive")
  }
  q <- .dsvert_glm_rat_max(
    .dsvert_glm_rat_abs(.dsvert_glm_rat_sub(
      "1", .dsvert_glm_rat_mul(alpha, lambda_min))),
    .dsvert_glm_rat_abs(.dsvert_glm_rat_sub(
      "1", .dsvert_glm_rat_mul(alpha, smoothness))))
  if (.dsvert_glm_rat_cmp(q, "1") >= 0L) {
    .dsvert_formal_glm_stop(
      "The public optimizer is not strictly contractive",
      "optimizer_not_contractive")
  }
  outcome_upper <- if (identical(family, "binomial")) {
    .dsvert_glm_rat("1")
  } else {
    .dsvert_glm_rat(response$upper)
  }
  link_upper <- .dsvert_glm_rat(link$range_upper)
  residual_upper <- if (identical(family, "binomial")) {
    .dsvert_glm_rat("1")
  } else {
    .dsvert_glm_rat_max(link_upper, outcome_upper)
  }
  score_upper <- .dsvert_glm_rat_mul(
    .dsvert_glm_rat_mul(maximum_weight, r_upper), residual_upper)
  fixed_gradient_difference <- .dsvert_glm_rat_div(
    .dsvert_glm_rat_mul(as.character(gamma), score_upper),
    as.character(capacity))
  q_power <- .dsvert_glm_rat_pow(q, iterations)
  sensitivity <- .dsvert_glm_rat_mul(
    .dsvert_glm_rat_mul(alpha, fixed_gradient_difference),
    .dsvert_glm_rat_div(
      .dsvert_glm_rat_sub("1", q_power),
      .dsvert_glm_rat_sub("1", q)))
  diameter_upper <- .dsvert_glm_rat_mul(
    "2", .dsvert_formal_glm_sum(boxes))
  optimization_error <- .dsvert_glm_rat_mul(q_power, diameter_upper)
  link_error <- .dsvert_glm_rat_div(
    .dsvert_glm_rat_mul(
      .dsvert_glm_rat_mul(maximum_weight, r_upper),
      .dsvert_glm_rat(link$uniform_mean_error_upper)),
    lambda_min)

  interval_bounds <- list(
    coefficient_abs = .dsvert_glm_rat_json(Reduce(.dsvert_glm_rat_max, boxes)),
    offset_abs = .dsvert_glm_rat_json(maximum_offset),
    eta_abs = .dsvert_glm_rat_json(maximum_eta),
    link_abs = .dsvert_glm_rat_json(link_upper),
    residual_abs = .dsvert_glm_rat_json(residual_upper),
    row_score_l2_upper = .dsvert_glm_rat_json(score_upper),
    score_accumulator_l2_upper = .dsvert_glm_rat_json(
      .dsvert_glm_rat_mul(as.character(capacity), score_upper)),
    regularized_gradient_l2_upper = .dsvert_glm_rat_json(
      .dsvert_glm_rat_add(
        score_upper,
        .dsvert_glm_rat_mul(lambda_max,
                            .dsvert_formal_glm_sum(boxes)))),
    projected_update_abs = .dsvert_glm_rat_json(
      Reduce(.dsvert_glm_rat_max, boxes)))
  numeric_plan <- .dsvert_formal_glm_numeric_plan(interval_bounds, numeric)

  .dsvert_formal_glm_exact_names(
    model$privacy, c("mechanism", "epsilon", "delta", "allocation"),
    "privacy policy")
  if (!identical(model$privacy$mechanism,
                 "joint_discrete_gaussian_one_global_draw") ||
      !identical(model$privacy$allocation, "one_stacked_capsule_vector")) {
    .dsvert_formal_glm_stop(
      "The privacy mechanism is not the formal one-shot vector contract",
      "privacy_mechanism_not_certified")
  }
  epsilon <- .dsvert_formal_glm_rat_field(
    model$privacy$epsilon, "epsilon", lower = "0", lower_open = TRUE)
  delta <- .dsvert_formal_glm_rat_field(
    model$privacy$delta, "delta", lower = "0", lower_open = TRUE,
    upper = "1", upper_open = TRUE)

  designated <- authority$designated_peers
  if (!is.character(designated) || length(designated) != 2L ||
      anyNA(designated) || anyDuplicated(designated)) {
    .dsvert_formal_glm_stop(
      "Exactly two distinct designated pinned peers are required",
      "invalid_peer_authority")
  }
  designated <- unname(sort(vapply(
    designated, .dsvert_formal_glm_identifier,
    character(1L), what = "designated peer"), method = "radix"))
  custodians <- authority$custodian_peers
  if (!is.character(custodians) || length(custodians) < 2L ||
      anyNA(custodians) || anyDuplicated(custodians)) {
    .dsvert_formal_glm_stop(
      "The authority must bind the complete K>=2 custodian set",
      "invalid_peer_authority")
  }
  custodians <- unname(sort(vapply(
    custodians, .dsvert_formal_glm_identifier,
    character(1L), what = "custodian peer"), method = "radix"))
  owners <- unique(vapply(artifact_columns, `[[`, character(1L), "owner"))
  if (!all(designated %in% custodians) || !all(owners %in% custodians)) {
    .dsvert_formal_glm_stop(
      "Every designated peer and registered column owner must be in the pinset",
      "invalid_peer_authority")
  }
  hashes <- c("policy_contract_sha256", "pinset_sha256")
  if (any(!vapply(authority[hashes], function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) &&
      grepl("^[0-9a-f]{64}$", value)
  }, logical(1L)))) {
    .dsvert_formal_glm_stop(
      "Authority commitments must be lowercase SHA-256 values",
      "invalid_peer_authority")
  }

  coefficients <- lapply(seq_along(coefficient_names), function(index) {
    list(
      index = as.numeric(index), name = coefficient_names[[index]],
      term = terms[[index]],
      box_abs = .dsvert_glm_rat_json(boxes[[index]]),
      ridge = .dsvert_glm_rat_json(ridge[[index]]))
  })
  unsigned <- .dsvert_joint_dp_client_canonical(list(
    version = .DSVERT_FORMAL_GLM_SCHEMA_VERSION,
    compiler_version = .DSVERT_FORMAL_GLM_COMPILER_VERSION,
    theorem_version = .DSVERT_FORMAL_GLM_THEOREM_VERSION,
    authority = list(
      consortium_id = .dsvert_formal_glm_string(
        authority$consortium_id, "consortium_id"),
      capsule_id = .dsvert_formal_glm_string(
        authority$capsule_id, "capsule_id"),
      logical_snapshot = .dsvert_formal_glm_string(
        authority$logical_snapshot, "logical_snapshot", 1024L),
      policy_contract_sha256 = authority$policy_contract_sha256,
      custodian_peers = as.list(custodians),
      designated_peers = as.list(designated),
      pinset_sha256 = authority$pinset_sha256,
      schema_authority = "two-custodian-cross-signed-before-materialization"),
    estimand = list(
      family = family,
      link = if (identical(family, "binomial")) "logit" else "log",
      target = paste0(
        "capacity-normalized-box-constrained-fully-ridge-regularized-",
        "quantized-pwl-fixed-iteration-", family),
      formula = parsed$canonical,
      response = parsed$response,
      predictors = as.list(parsed$predictors),
      intercept = parsed$intercept,
      coefficient_order = as.list(coefficient_names),
      coefficients = coefficients,
      capacity_normalizer = as.numeric(capacity),
      missingness = model$missingness,
      clipping = model$clipping,
      patient_collapse = collapse,
      weights = list(
        mode = model$weights$mode, column = weight_column,
        source_maximum_patient_weight =
          .dsvert_glm_rat_json(source_maximum_weight),
        maximum_patient_weight = .dsvert_glm_rat_json(maximum_weight)),
      offset = list(
        mode = model$offset$mode, column = offset_column,
        source_lower = .dsvert_glm_rat_json(offset_source_lower),
        source_upper = .dsvert_glm_rat_json(offset_source_upper),
        lower = .dsvert_glm_rat_json(offset_lower),
        upper = .dsvert_glm_rat_json(offset_upper),
        log_quantization_error_upper =
          .dsvert_glm_rat_json(offset_error)),
      column_registry = artifact_columns),
    adjacency = list(
      definition = adjacency, unit = "one_aligned_patient",
      fixed_padded_capacity = as.numeric(capacity),
      triangle_factor_gamma = as.numeric(gamma)),
    optimizer = list(
      algorithm = "fixed_iteration_projected_full_gradient_v1",
      start = "zero", iterations = as.numeric(iterations),
      alpha = .dsvert_glm_rat_json(alpha),
      early_stopping = FALSE, line_search = FALSE,
      data_dependent_branching = FALSE,
      reduction_order = "capacity_slot_then_coefficient_v1"),
    link_surrogate = link,
    numeric = c(numeric, list(
      plan = numeric_plan,
      quantized_snapshot_estimand = TRUE,
      production_fixed_point_certificate_pending_phase1 = TRUE)),
    theorem_certificate = list(
      status = "reference_rational_certified",
      production_release_ready = FALSE,
      arithmetic_error_rho = .dsvert_glm_rat_json("0"),
      arithmetic_error_scope = "exact rational reference oracle only",
      r_squared = .dsvert_glm_rat_json(r_squared),
      r_l2_upper = .dsvert_glm_rat_json(r_upper),
      r_l2_upper_method = "l1 upper bound over coordinate absolute bounds",
      maximum_eta = .dsvert_glm_rat_json(maximum_eta),
      strong_convexity_m = .dsvert_glm_rat_json(lambda_min),
      smoothness_L = .dsvert_glm_rat_json(smoothness),
      contraction_q = .dsvert_glm_rat_json(q),
      patient_score_l2_upper = .dsvert_glm_rat_json(score_upper),
      fixed_beta_gradient_difference =
        .dsvert_glm_rat_json(fixed_gradient_difference),
      global_l2_sensitivity_reference =
        .dsvert_glm_rat_json(sensitivity),
      optimization_error_reference =
        .dsvert_glm_rat_json(optimization_error),
      link_coefficient_error_upper = .dsvert_glm_rat_json(link_error),
      numeric_feature_quantization_abs_error_per_coordinate =
        .dsvert_glm_rat_json(.dsvert_glm_rat_div(
          "1", .dsvert_glm_rat_pow(
            "2", numeric$x_fraction_bits + 1L))),
      weight_quantization_abs_error =
        .dsvert_glm_rat_json(.dsvert_glm_rat_div(
          "1", .dsvert_glm_rat_pow(
            "2", numeric$x_fraction_bits + 1L))),
      offset_or_log_exposure_quantization_abs_error =
        .dsvert_glm_rat_json(offset_error),
      clipping_is_estimand_change_without_universal_error_bound = TRUE,
      ridge_and_box_are_estimand_changes = TRUE,
      estimability = "unique_for_every_admissible_dataset_by_full_ridge_and_box",
      no_ordinary_mle_claim = TRUE),
    privacy = list(
      mechanism = model$privacy$mechanism,
      allocation = model$privacy$allocation,
      epsilon = .dsvert_glm_rat_json(epsilon),
      delta = .dsvert_glm_rat_json(delta),
      history_can_deny_operation = FALSE,
      request_limit = FALSE,
      operation_limit = FALSE,
      calibration_status = "phase2_pending"),
    execution_requirements = list(
      protected_input_artifact = "fixed_capacity_patient_tuple_shares_v1",
      protected_row_coordinates = as.numeric(length(coefficient_names) + 3L),
      finite_exact_glm_sufficient_statistic = FALSE,
      exact_mpc_primitives = as.list(c(
        "signed_segment_comparison_and_mux",
        "checked_multiply_and_exact_truncate",
        "fixed_order_checked_ring_accumulation",
        "coefficient_box_projection",
        "one_global_draw_two_peer_noise_and_final_projection")),
      legacy_dcf_or_local_truncation_allowed = FALSE,
      plaintext_intermediate_opening_allowed = FALSE,
      production_phase = "phase1_not_connected"),
    transcript = list(
      protected_row_count = as.numeric(capacity),
      coefficient_count = as.numeric(length(coefficient_names)),
      optimizer_iterations = as.numeric(iterations),
      link_segment_count = as.numeric(segments),
      source_shape_is_public_fixed = TRUE,
      data_dependent_stop_or_error = FALSE,
      public_release_coordinates = "dp_coefficients_only")))
  json <- .dsvert_joint_dp_client_json(unsigned)
  hash <- digest::digest(json, algo = "sha256", serialize = FALSE)
  payload <- charToRaw(paste0(.DSVERT_FORMAL_GLM_SIGNATURE_DOMAIN, json))
  structure(list(
    unsigned_schema = unsigned,
    canonical_json = json,
    sha256 = hash,
    signature_domain = .DSVERT_FORMAL_GLM_SIGNATURE_DOMAIN,
    signature_payload = payload,
    signature_status = "unsigned_phase0_requires_two_custodian_signatures",
    production_release_ready = FALSE),
    class = "dsvert_formal_glm_compilation")
}

.dsvert_formal_glm_validate_compilation <- function(compilation) {
  if (!inherits(compilation, "dsvert_formal_glm_compilation") ||
      !is.list(compilation$unsigned_schema) ||
      !identical(compilation$unsigned_schema$version,
                 .DSVERT_FORMAL_GLM_SCHEMA_VERSION)) {
    .dsvert_formal_glm_stop(
      "Invalid formal GLM compilation object", "compiled_schema_tampered")
  }
  json <- tryCatch(
    .dsvert_joint_dp_client_json(compilation$unsigned_schema),
    error = function(error) NULL)
  hash <- if (is.null(json)) NULL else digest::digest(
    json, algo = "sha256", serialize = FALSE)
  payload <- if (is.null(json)) raw() else charToRaw(paste0(
    .DSVERT_FORMAL_GLM_SIGNATURE_DOMAIN, json))
  if (is.null(json) || !identical(json, compilation$canonical_json) ||
      !identical(hash, compilation$sha256) ||
      !identical(payload, compilation$signature_payload) ||
      !identical(compilation$production_release_ready, FALSE)) {
    .dsvert_formal_glm_stop(
      "The formal GLM compilation commitment was modified",
      "compiled_schema_tampered")
  }
  invisible(TRUE)
}

.dsvert_formal_glm_verify_cross_signatures <- function(
    compilation, signatures, pinset) {
  .dsvert_formal_glm_validate_compilation(compilation)
  peers <- unlist(
    compilation$unsigned_schema$authority$designated_peers,
    use.names = FALSE)
  custodians <- unlist(
    compilation$unsigned_schema$authority$custodian_peers,
    use.names = FALSE)
  if (!is.list(signatures) || is.null(names(signatures)) ||
      anyDuplicated(names(signatures)) || !setequal(names(signatures), peers) ||
      !is.character(pinset) || is.null(names(pinset)) ||
      anyDuplicated(names(pinset)) || anyDuplicated(unname(pinset)) ||
      !setequal(names(pinset), custodians)) {
    .dsvert_formal_glm_stop(
      "The formal GLM schema requires both designated pinned signatures",
      "schema_signature_set_invalid")
  }
  valid_pinset <- vapply(custodians, function(peer) {
    length(tryCatch(.dsvert_joint_dp_client_b64url(
      unname(pinset[[peer]]), 32L, "formal GLM custodian identity key"),
      error = function(error) raw())) == 32L
  }, logical(1L))
  if (any(!valid_pinset)) {
    .dsvert_formal_glm_stop(
      "The formal GLM schema pinset contains an invalid identity key",
      "schema_signature_set_invalid")
  }
  signatures <- signatures[peers]
  for (peer in peers) {
    public <- tryCatch(.dsvert_joint_dp_client_b64url(
      unname(pinset[[peer]]), 32L, "formal GLM schema identity key"),
      error = function(error) raw())
    signature <- tryCatch(.dsvert_joint_dp_client_b64url(
      signatures[[peer]], 64L, "formal GLM schema signature"),
      error = function(error) raw())
    valid <- length(public) == 32L && length(signature) == 64L &&
      tryCatch(openssl::ed25519_verify(
        compilation$signature_payload, signature,
        openssl::read_ed25519_pubkey(public)),
        error = function(error) FALSE)
    if (!isTRUE(valid)) {
      .dsvert_formal_glm_stop(
        paste0("Invalid formal GLM schema signature from pinned peer ", peer),
        "schema_signature_invalid")
    }
  }
  value <- .dsvert_joint_dp_client_canonical(list(
    version = "dsvert-formal-glm-cross-signed-schema-v1",
    schema_sha256 = compilation$sha256,
    designated_peers = as.list(peers),
    signatures = signatures))
  structure(value, class = c("dsvert_formal_glm_signed_schema", "list"))
}
