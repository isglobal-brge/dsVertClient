#' @title Deterministic categorical multiple-imputation core
#' @description Completes a bounded categorical count vector from one signed
#'   sticky synopsis release. It is pure post-processing: neither private
#'   rows nor an analyst seed enter the calculation.
#' @noRd
.dsvert_mi_complete_counts_v1 <- function(
    observed_counts, admitted_count, m, release_sha256) {
  valid_counts <- is.numeric(observed_counts) && length(observed_counts) >= 2L &&
    !is.null(names(observed_counts)) && !anyNA(names(observed_counts)) &&
    all(nzchar(names(observed_counts))) && !anyDuplicated(names(observed_counts)) &&
    !anyNA(observed_counts) && all(is.finite(observed_counts)) &&
    all(observed_counts >= 0)
  if (!isTRUE(valid_counts)) {
    stop("observed_counts must be named non-negative DP projected counts",
         call. = FALSE)
  }
  valid_total <- is.numeric(admitted_count) && length(admitted_count) == 1L &&
    !is.na(admitted_count) && is.finite(admitted_count) &&
    admitted_count >= 0
  if (!isTRUE(valid_total)) {
    stop("admitted_count must be a non-negative DP projected count", call. = FALSE)
  }
  valid_m <- is.numeric(m) && length(m) == 1L && !is.na(m) &&
    is.finite(m) && m == floor(m) && m >= 2L && m <= 100L
  if (!isTRUE(valid_m)) stop("m must be an integer in [2, 100]", call. = FALSE)
  if (!is.character(release_sha256) || length(release_sha256) != 1L ||
      is.na(release_sha256) || !grepl("^[0-9a-f]{64}$", release_sha256)) {
    stop("release_sha256 must be a canonical signed release root", call. = FALSE)
  }

  observed_total <- sum(observed_counts)
  completed_total <- max(admitted_count, observed_total)
  missing_raw <- admitted_count - observed_total
  missing_count <- max(0, missing_raw)
  uniform <- function(draw, coordinate, purpose) {
    digest <- openssl::sha256(charToRaw(paste(
      "dsVert/mi/categorical-mcar-draw/v1", release_sha256, draw,
      coordinate, purpose, sep = "|")))
    integer <- sum(as.numeric(digest[seq_len(4L)]) *
      c(2^24, 2^16, 2^8, 1))
    (integer + 0.5) / 2^32
  }
  draws <- matrix(0, nrow = as.integer(m), ncol = length(observed_counts),
                  dimnames = list(NULL, names(observed_counts)))
  for (draw in seq_len(as.integer(m))) {
    shape <- observed_counts + 0.5
    gamma <- vapply(seq_along(shape), function(index) {
      stats::qgamma(uniform(draw, index, "dirichlet"), shape = shape[[index]])
    }, numeric(1L))
    probabilities <- gamma / sum(gamma)
    draws[draw, ] <- observed_counts + missing_count * probabilities
  }
  result <- list(
    observed_counts_dp = observed_counts,
    admitted_count_dp = admitted_count,
    missing_count_dp_raw = missing_raw,
    missing_count_dp = missing_count,
    completed_count_dp = completed_total,
    missing_count_projection = if (missing_raw < 0) {
      "observed_total_exceeds_noisy_admitted_count_use_observed_total_v1"
    } else "none",
    completed_counts = draws,
    pooled_probabilities = colMeans(draws) / completed_total)
  result$draws_sha256 <- digest::digest(
    list(version = "dsvert-mi-categorical-mcar-draws-v1",
         release_sha256 = release_sha256, completed_counts = draws),
    algo = "sha256", serialize = TRUE, serializeVersion = 3L)
  result
}

#' @keywords internal
.dsvert_mi_strict_missingness_policy_v1 <- paste(
  "missing_values_have_no_marginal_cell_and_unknown_or_conflicting",
  "nonmissing_values_reject_before_release_v1", sep = "_")

#' @keywords internal
.dsvert_mi_response_columns_v1 <- function(formula) {
  valid_formula <- inherits(formula, "formula") && length(formula) == 3L &&
    length(attr(stats::terms(formula), "term.labels")) == 0L &&
    identical(attr(stats::terms(formula), "intercept"), 1L)
  response <- if (isTRUE(valid_formula)) formula[[2L]] else NULL
  columns <- if (is.symbol(response)) {
    as.character(response)
  } else if (is.call(response) && identical(as.character(response[[1L]]), "cbind")) {
    values <- as.list(response)[-1L]
    if (length(values) < 2L || !all(vapply(values, is.symbol, logical(1L)))) {
      character()
    } else {
      vapply(values, as.character, character(1L))
    }
  } else {
    character()
  }
  if (!length(columns) || anyNA(columns) || any(!nzchar(columns)) ||
      anyDuplicated(columns)) {
    stop(paste(
      "ds.vertMI supports only outcome ~ 1 or",
      "cbind(outcome1, outcome2, ...) ~ 1 with untransformed columns"),
         call. = FALSE)
  }
  columns
}

#' @keywords internal
.dsvert_mi_variable_draw_root_v1 <- function(release_sha256, variable) {
  digest::digest(
    list(version = "dsvert-mi-independent-marginal-draw-root-v1",
         release_sha256 = release_sha256, variable = variable),
    algo = "sha256", serialize = TRUE, serializeVersion = 3L)
}

#' @keywords internal
.dsvert_mi_synopsis_result_v1 <- function(
    formula, data_name, impute_columns, m, family, datasources, .aggregate,
    .run = .dsvert_dp_synopsis_vector_run,
    .count = .dsvert_dp_count_synopsis_result_v1,
    .frequency = .dsvert_dp_frequency_synopsis_result_v1) {
  outcomes <- .dsvert_mi_response_columns_v1(formula)
  if (!is.character(data_name) || length(data_name) != 1L ||
      is.na(data_name) || !nzchar(data_name)) {
    stop("data must name one signed protected dataset", call. = FALSE)
  }
  if (is.null(impute_columns)) impute_columns <- outcomes
  if (!is.character(impute_columns) || length(impute_columns) != length(outcomes) ||
      anyNA(impute_columns) || any(!nzchar(impute_columns)) ||
      anyDuplicated(impute_columns) || !identical(impute_columns, outcomes)) {
    stop(paste(
      "impute_columns must exactly match the categorical response columns",
      "in formula order"), call. = FALSE)
  }
  if (!is.numeric(m) || length(m) != 1L || is.na(m) || !is.finite(m) ||
      m != floor(m) || m < 2L || m > 100L) {
    stop("m must be an integer in [2, 100]", call. = FALSE)
  }
  family <- match.arg(family, c("auto", "binomial", "multinomial"))
  if (length(outcomes) > 1L && !identical(family, "auto")) {
    stop(paste(
      "multivariable categorical MI requires family = 'auto' because",
      "each signed marginal determines its own binary or multinomial domain"),
         call. = FALSE)
  }
  if (!is.function(.run) || !is.function(.count) || !is.function(.frequency)) {
    stop("Invalid categorical MI Synopsis dependency", call. = FALSE)
  }
  run <- .run(datasources, .aggregate = .aggregate)
  replay <- function(...) run
  count <- .count(
    data_name, NULL, datasources, .aggregate = .aggregate, .run = replay)
  binding_fields <- c(
    "artifact_key", "execution_id", "manifest_sha256", "contract_sha256",
    "attempt_sha256", "source_contract_sha256", "result_set_sha256",
    "final_vector_root", "coordinate_order_sha256", "release_sha256")
  variable_result <- function(outcome) {
    frequency <- .frequency(
      data_name, outcome, count$source_owner, datasources,
      .aggregate = .aggregate, .run = replay)
    binding_ok <- all(vapply(binding_fields, function(field) {
      identical(count[[field]], frequency[[field]])
    }, logical(1L))) &&
      identical(count$dataset, data_name) &&
      identical(frequency$variable, outcome) &&
      identical(frequency$source_owner, count$source_owner) &&
      identical(frequency$missingness_policy,
                .dsvert_mi_strict_missingness_policy_v1) &&
      is.list(frequency$coordinate_descriptor) &&
      identical(frequency$coordinate_descriptor$missingness_policy,
                .dsvert_mi_strict_missingness_policy_v1)
    if (!isTRUE(binding_ok)) {
      stop("The signed categorical marginal is not bound to strict missingness",
           call. = FALSE)
    }
    levels <- frequency$levels
    counts <- frequency$counts
    if (!is.character(levels) || length(levels) < 2L || anyNA(levels) ||
        any(!nzchar(levels)) || anyDuplicated(levels) ||
        !is.numeric(counts) || !identical(names(counts), levels) ||
        anyNA(counts) || any(!is.finite(counts)) || any(counts < 0) ||
        !is.numeric(count$value) || length(count$value) != 1L ||
        is.na(count$value) || !is.finite(count$value) || count$value < 0) {
      stop("The signed categorical MI coordinates are invalid", call. = FALSE)
    }
    selected_family <- if (identical(family, "auto")) {
      if (length(levels) == 2L) "binomial" else "multinomial"
    } else family
    if (identical(selected_family, "binomial") && length(levels) != 2L) {
      stop("family='binomial' requires exactly two signed outcome levels",
           call. = FALSE)
    }
    draw_root <- if (length(outcomes) == 1L) frequency$release_sha256 else
      .dsvert_mi_variable_draw_root_v1(frequency$release_sha256, outcome)
    completion <- .dsvert_mi_complete_counts_v1(
      counts, count$value, as.integer(m), draw_root)
    total <- completion$completed_count_dp
    status <- if (total > 0) "ok" else "dp_effective_count_zero"
    probabilities <- completion$pooled_probabilities
    coefficients <- if (identical(status, "ok")) {
      floor_probability <- 1 / (2 * total)
      if (identical(selected_family, "binomial")) {
        stats::setNames(stats::qlogis(min(max(probabilities[[2L]],
                                                floor_probability),
                                            1 - floor_probability)),
                        "(Intercept)")
      } else {
        stats::setNames(log(pmax(probabilities[-1L], floor_probability) /
                              max(probabilities[[1L]], floor_probability)),
                        paste0("(Intercept):", levels[-1L]))
      }
    } else NULL
    list(
      status = status, family = selected_family, outcome = outcome,
      outcome_levels = levels, reference_level = levels[[1L]],
      coefficients = coefficients, probabilities = probabilities,
      observed_counts_dp = completion$observed_counts_dp,
      admitted_count_dp = completion$admitted_count_dp,
      missing_count_dp = completion$missing_count_dp,
      completed_count_dp = completion$completed_count_dp,
      missing_count_projection = completion$missing_count_projection,
      completed_draws_sha256 = completion$draws_sha256,
      release_sha256 = frequency$release_sha256)
  }
  variables <- stats::setNames(lapply(outcomes, variable_result), outcomes)
  if (length(variables) == 1L) {
    variable <- variables[[1L]]
    result <- c(variable, list(
      method = "signed_categorical_mcar_intercept_only_v1",
      formula = stats::as.formula(formula), m = as.integer(m),
      source_owner = count$source_owner,
      artifact_key = count$artifact_key, execution_id = count$execution_id,
      manifest_sha256 = count$manifest_sha256, contract_sha256 = count$contract_sha256,
      attempt_sha256 = count$attempt_sha256,
      source_contract_sha256 = count$source_contract_sha256,
      result_set_sha256 = count$result_set_sha256,
      final_vector_root = count$final_vector_root,
      coordinate_order_sha256 = count$coordinate_order_sha256))
  } else {
    result <- list(
      status = if (all(vapply(variables, `[[`, character(1L), "status") == "ok")) {
        "ok"
      } else "some_variables_dp_effective_count_zero",
      method = "signed_categorical_mcar_independent_marginals_v2",
      joint_model = "independent_marginal_completion_no_joint_imputation_v1",
      formula = stats::as.formula(formula), impute_columns = outcomes,
      variables = variables, m = as.integer(m), source_owner = count$source_owner,
      artifact_key = count$artifact_key, execution_id = count$execution_id,
      manifest_sha256 = count$manifest_sha256, contract_sha256 = count$contract_sha256,
      attempt_sha256 = count$attempt_sha256,
      source_contract_sha256 = count$source_contract_sha256,
      result_set_sha256 = count$result_set_sha256,
      final_vector_root = count$final_vector_root,
      coordinate_order_sha256 = count$coordinate_order_sha256,
      completed_draws_sha256 = stats::setNames(vapply(
        variables, `[[`, character(1L), "completed_draws_sha256"), outcomes))
  }
  result <- c(result, list(
    sticky_replay = TRUE,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    assumption = paste(
      "MCAR categorical missingness under the signed strict-missingness",
      "contract; multivariable completions are independent marginals"),
    inference_scope = "No classical or Rubin sampling inference is provided",
    standard_errors = NULL, covariance = NULL, p_values = NULL,
    source_values_exposed = FALSE, intermediate_values_exposed = FALSE,
    production_ready = FALSE))
  class(result) <- c("ds.vertMI", "list")
  result
}

#' @title Signed categorical multiple-imputation compatibility route
#' @description Returns categorical MCAR completions from one signed sticky
#'   Synopsis release. It accepts either \code{outcome ~ 1} or
#'   \code{cbind(outcome1, outcome2, ...) ~ 1}; the multivariable form
#'   completes signed marginals independently and never claims a joint model.
#'   It never mutates source tables and its deterministic completion draws are
#'   post-processing of the released vector.
#' @param formula An intercept-only categorical response formula exactly of the
#'   form \code{outcome ~ 1} or \code{cbind(outcome1, outcome2, ...) ~ 1}.
#' @param data Signed protected dataset name or federation.
#' @param impute_columns Must be omitted or exactly match the response columns
#'   in formula order.
#' @param m Number of deterministic categorical completion draws, from 2 to 100.
#' @param family One of \code{"auto"}, \code{"binomial"}, or
#'   \code{"multinomial"}. The multivariable route requires \code{"auto"}
#'   so each signed marginal determines whether it is binary or multinomial.
#' @param max_iter,tol,lambda,intercept_only,verbose,seed Compatibility controls.
#'   Only \code{lambda = 0}, \code{intercept_only = "aggregate"}, the default
#'   iteration values, and \code{seed = NULL} are supported.
#' @param datasources DataSHIELD connections.
#' @return A \code{ds.vertMI} object with DP-projected completed-category
#'   probabilities and no classical or Rubin sampling inference. Multi-response
#'   calls return one protected marginal completion per named response, not a
#'   completed joint microdata set.
#' @seealso \code{\link{ds.vertMethodStatus}}
#' @export
ds.vertMI <- function(formula, data = NULL, impute_columns = NULL,
                      m = 20L, family = c("auto", "binomial", "multinomial"),
                      max_iter = 50L, tol = 1e-4, lambda = 0,
                      intercept_only = "aggregate",
                      verbose = TRUE, datasources = NULL, seed = NULL) {
  if (!identical(max_iter, 50L) || !identical(tol, 1e-4) ||
      !is.numeric(lambda) || length(lambda) != 1L || is.na(lambda) ||
      !is.finite(lambda) || lambda != 0 ||
      !identical(intercept_only, "aggregate") || !is.logical(verbose) ||
      length(verbose) != 1L || is.na(verbose) || !is.null(seed)) {
    stop(paste(
      "ds.vertMI supports only the signed categorical MCAR route with",
      "default iteration controls, lambda=0, intercept_only='aggregate',",
      "and no analyst seed"), call. = FALSE)
  }
  resolved <- .dsvert_federation_argument(data, datasources)
  .dsvert_mi_synopsis_result_v1(
    formula, resolved$value, impute_columns, m, family,
    resolved$datasources, DSI::datashield.aggregate)
}

#' @export
print.ds.vertMI <- function(x, ...) {
  cat("dsVert signed categorical MCAR completion\n")
  if (is.list(x$variables)) {
    cat("  Independent signed marginals:",
        paste(names(x$variables), collapse = ", "), "\n")
    for (name in names(x$variables)) {
      variable <- x$variables[[name]]
      cat("   -", name, "| family:", variable$family,
          "| missing (DP):", variable$missing_count_dp, "\n")
    }
    cat("  No joint completed data, associations or Rubin inference are released.\n")
  } else {
    cat("  Outcome:", x$outcome, "| family:", x$family,
        "| missing (DP):", x$missing_count_dp, "\n")
    if (identical(x$status, "ok")) print(x$coefficients)
  }
  cat("  No classical or Rubin sampling inference is provided.\n")
  invisible(x)
}
