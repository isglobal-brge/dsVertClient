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
.dsvert_mi_strict_joint_missingness_policy_v1 <- paste(
  "missing_values_have_no_joint_cell_and_unknown_or_conflicting",
  "nonmissing_values_reject_before_release_v1", sep = "_")

#' @keywords internal
.dsvert_mi_complete_joint_counts_v1 <- function(
    observed_table, admitted_count, m, release_sha256) {
  valid_table <- is.matrix(observed_table) && is.numeric(observed_table) &&
    nrow(observed_table) >= 2L && ncol(observed_table) >= 2L &&
    !is.null(rownames(observed_table)) && !is.null(colnames(observed_table)) &&
    !anyNA(rownames(observed_table)) && !anyNA(colnames(observed_table)) &&
    all(nzchar(rownames(observed_table))) && all(nzchar(colnames(observed_table))) &&
    !anyDuplicated(rownames(observed_table)) && !anyDuplicated(colnames(observed_table)) &&
    !anyNA(observed_table) && all(is.finite(observed_table)) &&
    all(observed_table >= 0)
  if (!isTRUE(valid_table)) {
    stop("observed_table must be a named non-negative DP joint table",
         call. = FALSE)
  }
  completed <- .dsvert_mi_complete_counts_v1(
    stats::setNames(as.numeric(observed_table),
                    sprintf("cell-%08d", seq_len(length(observed_table)))),
    admitted_count, m, release_sha256)
  draws <- array(
    completed$completed_counts,
    dim = c(nrow(completed$completed_counts), nrow(observed_table),
            ncol(observed_table)),
    dimnames = list(NULL, rownames(observed_table), colnames(observed_table)))
  pooled <- matrix(
    completed$pooled_probabilities, nrow = nrow(observed_table),
    ncol = ncol(observed_table),
    dimnames = dimnames(observed_table))
  list(
    observed_table_dp = observed_table,
    admitted_count_dp = completed$admitted_count_dp,
    missing_count_dp = completed$missing_count_dp,
    completed_count_dp = completed$completed_count_dp,
    missing_count_projection = completed$missing_count_projection,
    pooled_joint_probabilities = pooled,
    draws_sha256 = digest::digest(
      list(version = "dsvert-mi-categorical-joint-pair-draws-v1",
           release_sha256 = release_sha256, completed_counts = draws),
      algo = "sha256", serialize = TRUE, serializeVersion = 3L))
}

#' @keywords internal
.dsvert_mi_joint_pair_result_v1 <- function(
    formula, data_name, outcomes, m, datasources, .aggregate, .contingency) {
  release <- .contingency(
    data_name, outcomes[[1L]], outcomes[[2L]], NULL, datasources, .aggregate)
  .dsvert_mi_joint_pair_completion_v1(formula, outcomes, m, release)
}

#' @keywords internal
.dsvert_mi_joint_pair_completion_v1 <- function(formula, outcomes, m, release) {
  binding_fields <- c(
    "artifact_key", "execution_id", "manifest_sha256", "contract_sha256",
    "attempt_sha256", "source_contract_sha256", "result_set_sha256",
    "final_vector_root", "coordinate_order_sha256", "release_sha256")
  binding_ok <- is.list(release) &&
    all(vapply(binding_fields, function(field) {
      is.character(release[[field]]) && length(release[[field]]) == 1L &&
        grepl("^[0-9a-f]{64}$", release[[field]])
    }, logical(1L))) &&
    identical(release$row_var, outcomes[[1L]]) &&
    identical(release$col_var, outcomes[[2L]]) &&
    identical(release$missingness_policy,
              .dsvert_mi_strict_joint_missingness_policy_v1) &&
    is.numeric(release$admitted_count_dp) && length(release$admitted_count_dp) == 1L &&
    !is.na(release$admitted_count_dp) && is.finite(release$admitted_count_dp) &&
    release$admitted_count_dp >= 0
  if (!isTRUE(binding_ok)) {
    stop("The signed categorical pair is not bound to strict missingness",
         call. = FALSE)
  }
  completion <- .dsvert_mi_complete_joint_counts_v1(
    release$table, release$admitted_count_dp, as.integer(m),
    release$release_sha256)
  status <- if (completion$completed_count_dp > 0) "ok" else
    "dp_effective_count_zero"
  variables <- stats::setNames(list(
    list(levels = rownames(completion$pooled_joint_probabilities),
         probabilities = rowSums(completion$pooled_joint_probabilities)),
    list(levels = colnames(completion$pooled_joint_probabilities),
         probabilities = colSums(completion$pooled_joint_probabilities))), outcomes)
  result <- c(list(
    status = status,
    method = "signed_categorical_mcar_joint_pair_v3",
    joint_model = "strict_missing_signed_joint_pair_completion_v1",
    formula = stats::as.formula(formula), impute_columns = outcomes,
    variables = variables,
    joint_probabilities = completion$pooled_joint_probabilities,
    observed_joint_table_dp = completion$observed_table_dp,
    admitted_count_dp = completion$admitted_count_dp,
    missing_count_dp = completion$missing_count_dp,
    completed_count_dp = completion$completed_count_dp,
    missing_count_projection = completion$missing_count_projection,
    completed_draws_sha256 = completion$draws_sha256,
    m = as.integer(m)), release[binding_fields])
  result <- c(result, list(
    sticky_replay = TRUE,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    assumption = paste(
      "MCAR pair missingness under the signed strict-joint-missingness",
      "contract; the released joint categorical dependence is retained"),
    inference_scope = "No classical or Rubin sampling inference is provided",
    standard_errors = NULL, covariance = NULL, p_values = NULL,
    source_values_exposed = FALSE, intermediate_values_exposed = FALSE,
    production_ready = FALSE))
  class(result) <- c("ds.vertMI", "list")
  result
}

#' @keywords internal
.dsvert_mi_star_pair_release_v1 <- function(run, data_name, row_var, col_var) {
  references <- lapply(
    list(row_var, col_var), .dsvert_dp_synopsis_pair_reference_v1,
    what = "A categorical star variable")
  if (identical(references[[1L]]$reference, references[[2L]]$reference)) {
    stop("A categorical star pair needs two distinct variables", call. = FALSE)
  }
  context <- .dsvert_dp_vector_context(run, allow_synopsis = TRUE)
  accepts <- function(reference, side) {
    identical(reference$column, side$column) &&
      (is.null(reference$owner_peer) ||
         identical(reference$owner_peer, side$owner_peer))
  }
  block <- .dsvert_dp_capsule_single_block(
    context$layout, "categorical_pairs", dataset = data_name,
    predicate = function(candidate) {
      descriptor <- candidate$descriptor
      sides <- tryCatch(list(
        descriptor$left[c("dataset", "column", "owner_peer")],
        descriptor$right[c("dataset", "column", "owner_peer")]),
        error = function(error) list())
      if (length(sides) != 2L) return(FALSE)
      direct <- accepts(references[[1L]], sides[[1L]]) &&
        accepts(references[[2L]], sides[[2L]])
      reverse <- accepts(references[[1L]], sides[[2L]]) &&
        accepts(references[[2L]], sides[[1L]])
      xor(direct, reverse)
    },
    description = paste0("signed categorical star pair for '", row_var,
                         "' and '", col_var, "'"))
  descriptor <- block$descriptor
  canonical <- .dsvert_dp_capsule_vector_values(context$release, block)
  row_levels <- .dsvert_dp_capsule_manifest_strings(
    descriptor$left$levels, "left categorical star levels", sorted = TRUE)
  col_levels <- .dsvert_dp_capsule_manifest_strings(
    descriptor$right$levels, "right categorical star levels", sorted = TRUE)
  expected <- as.double(length(row_levels)) * length(col_levels)
  count_block <- .dsvert_dp_capsule_single_block(
    context$layout, "admitted_count",
    description = "signed admitted-count capacity block")
  admitted_count <- .dsvert_dp_capsule_vector_values(
    context$release, count_block)
  capacity <- .dsvert_dp_vector_block_capacity(count_block)
  if (length(canonical) != expected || any(!is.finite(canonical)) ||
      any(canonical < 0) || any(canonical > capacity) ||
      length(admitted_count) != 1L || !is.finite(admitted_count) ||
      admitted_count < 0 || admitted_count > capacity) {
    stop("The signed categorical star release violates its bound", call. = FALSE)
  }
  table <- matrix(canonical, nrow = length(row_levels),
                  ncol = length(col_levels),
                  dimnames = list(row_levels, col_levels))
  left <- descriptor$left[c("dataset", "column", "owner_peer")]
  right <- descriptor$right[c("dataset", "column", "owner_peer")]
  if (accepts(references[[1L]], left) && accepts(references[[2L]], right)) {
    mapped <- table
  } else if (accepts(references[[1L]], right) &&
             accepts(references[[2L]], left)) {
    mapped <- t(table)
  } else {
    stop("The signed categorical star descriptor changed during mapping",
         call. = FALSE)
  }
  c(.dsvert_dp_vector_public_metadata(context), list(
    row_var = row_var, col_var = col_var, table = mapped,
    admitted_count_dp = unname(admitted_count),
    release_sha256 = context$release$final_vector_root,
    missingness_policy = descriptor$missingness_policy))
}

#' @keywords internal
.dsvert_mi_joint_star_result_v1 <- function(
    formula, data_name, outcomes, m, run, .pair_release,
    root = outcomes[[1L]], include_root = TRUE) {
  if (!is.character(root) || length(root) != 1L || is.na(root) ||
      !nzchar(root) || !is.logical(include_root) || length(include_root) != 1L ||
      is.na(include_root) || !is.function(.pair_release)) {
    stop("Multivariable categorical MI requires signed strict pair releases",
         call. = FALSE)
  }
  children <- if (isTRUE(include_root)) outcomes[-1L] else outcomes
  if (!length(children) || anyNA(children) || any(!nzchar(children)) ||
      anyDuplicated(children) || root %in% children) {
    stop("Multivariable categorical MI requires signed strict pair releases",
         call. = FALSE)
  }
  pairs <- stats::setNames(lapply(children, function(child) {
    .dsvert_mi_joint_pair_completion_v1(
      formula, c(root, child), m,
      .pair_release(run, data_name, root, child))
  }), children)
  binding_fields <- c(
    "artifact_key", "execution_id", "manifest_sha256", "contract_sha256",
    "attempt_sha256", "source_contract_sha256", "result_set_sha256",
    "final_vector_root", "coordinate_order_sha256", "release_sha256",
    "admitted_count_dp")
  first <- pairs[[1L]]
  same_release <- vapply(pairs, function(pair) {
    all(vapply(binding_fields, function(field) {
      identical(pair[[field]], first[[field]])
    }, logical(1L)))
  }, logical(1L))
  if (!all(same_release)) {
    stop("All signed categorical pairs must use the same signed vector release",
         call. = FALSE)
  }
  root_table <- first$joint_probabilities
  root_probabilities <- rowSums(root_table)
  valid_root <- is.numeric(root_probabilities) &&
    !is.null(names(root_probabilities)) && !anyNA(root_probabilities) &&
    all(is.finite(root_probabilities)) && all(root_probabilities >= 0) &&
    isTRUE(all.equal(sum(root_probabilities), 1, tolerance = 1e-12))
  if (!isTRUE(valid_root)) {
    stop("The signed root categorical probabilities are invalid", call. = FALSE)
  }
  root_levels <- names(root_probabilities)
  conditional_probabilities <- stats::setNames(lapply(children, function(child) {
    table <- pairs[[child]]$joint_probabilities
    valid_table <- is.matrix(table) && is.numeric(table) &&
      identical(rownames(table), root_levels) && !is.null(colnames(table)) &&
      !anyNA(table) && all(is.finite(table)) && all(table >= 0)
    totals <- if (isTRUE(valid_table)) rowSums(table) else numeric()
    if (!isTRUE(valid_table) || any(totals <= 0)) {
      stop("Each signed categorical pair must have positive root support",
           call. = FALSE)
    }
    conditional <- table / totals
    if (!isTRUE(all.equal(unname(rowSums(conditional)), rep(1, nrow(conditional)),
                         tolerance = 1e-12))) {
      stop("The signed categorical pair conditional probabilities are invalid",
           call. = FALSE)
    }
    conditional
  }), children)
  variables <- stats::setNames(vector("list", length(children)), children)
  if (isTRUE(include_root)) {
    variables <- c(stats::setNames(list(
      list(levels = root_levels, probabilities = root_probabilities)), root),
      variables)
  }
  for (child in children) {
    probabilities <- as.numeric(root_probabilities %*%
                                  conditional_probabilities[[child]])
    names(probabilities) <- colnames(conditional_probabilities[[child]])
    variables[[child]] <- list(
      levels = names(probabilities), probabilities = probabilities)
  }
  root_fields <- if (isTRUE(include_root)) {
    list(root_column = root, root_probabilities = root_probabilities)
  } else {
    list(conditioning_column = root,
         conditioning_probabilities = root_probabilities)
  }
  result <- c(list(
    status = if (all(vapply(pairs, `[[`, character(1L), "status") == "ok")) {
      "ok"
    } else "some_variables_dp_effective_count_zero",
    method = if (isTRUE(include_root)) {
      "signed_categorical_mcar_star_joint_v1"
    } else "signed_categorical_mcar_covariate_star_v1",
    joint_model = if (isTRUE(include_root)) {
      "strict_missing_signed_pairwise_star_completion_v1"
    } else "strict_missing_signed_conditional_star_completion_v1",
    formula = stats::as.formula(formula), impute_columns = outcomes,
    conditional_probabilities = conditional_probabilities,
    variables = variables, m = as.integer(m)), root_fields, first[binding_fields])
  result <- c(result, list(
    completed_draws_sha256 = digest::digest(
      list(version = if (isTRUE(include_root)) {
             "dsvert-mi-categorical-star-joint-v1"
           } else "dsvert-mi-categorical-covariate-star-v1",
           release_sha256 = first$release_sha256, root = root,
           root_probabilities = root_probabilities,
           conditional_probabilities = conditional_probabilities, m = as.integer(m)),
      algo = "sha256", serialize = TRUE, serializeVersion = 3L),
    sticky_replay = TRUE,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    assumption = if (isTRUE(include_root)) {
      paste("MCAR pair missingness under signed strict pair releases; non-root",
            "responses are conditionally independent given the first response")
    } else {
      paste("MCAR pair missingness under signed strict pair releases; responses",
            "are conditionally independent given the categorical conditioning column")
    },
    inference_scope = "No classical or Rubin sampling inference is provided",
    standard_errors = NULL, covariance = NULL, p_values = NULL,
    source_values_exposed = FALSE, intermediate_values_exposed = FALSE,
    production_ready = FALSE))
  class(result) <- c("ds.vertMI", "list")
  result
}

#' @keywords internal
.dsvert_mi_formula_spec_v1 <- function(formula) {
  valid_formula <- inherits(formula, "formula") && length(formula) == 3L
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
      "ds.vertMI supports only outcome ~ 1, outcome ~ conditioning,",
      "or cbind(outcome1, outcome2, ...) with one untransformed categorical",
      "conditioning column"),
         call. = FALSE)
  }
  rhs <- if (isTRUE(valid_formula)) formula[[3L]] else NULL
  intercept_only <- is.numeric(rhs) && length(rhs) == 1L && !is.na(rhs) &&
    identical(as.numeric(rhs), 1)
  conditioning <- if (isTRUE(intercept_only)) {
    NULL
  } else if (is.symbol(rhs)) {
    as.character(rhs)
  } else {
    NA_character_
  }
  if (!is.null(conditioning)) {
    if (!is.character(conditioning) || length(conditioning) != 1L ||
        is.na(conditioning) || !nzchar(conditioning) || conditioning %in% columns) {
      stop(paste(
        "ds.vertMI supports only outcome ~ 1, outcome ~ conditioning,",
        "or cbind(outcome1, outcome2, ...) with one untransformed categorical",
        "conditioning column"), call. = FALSE)
    }
  }
  list(outcomes = columns, conditioning = conditioning)
}

#' @keywords internal
.dsvert_mi_response_columns_v1 <- function(formula) {
  .dsvert_mi_formula_spec_v1(formula)$outcomes
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
    .frequency = .dsvert_dp_frequency_synopsis_result_v1,
    .contingency = .dsvert_dp_contingency_impl,
    .star_pair = .dsvert_mi_star_pair_release_v1,
    dependence = c("independent", "star")) {
  formula_spec <- .dsvert_mi_formula_spec_v1(formula)
  outcomes <- formula_spec$outcomes
  conditioning <- formula_spec$conditioning
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
  dependence <- match.arg(dependence)
  if (!is.null(conditioning)) {
    if (!identical(family, "auto") || !identical(dependence, "star") ||
        !is.function(.run) || !is.function(.star_pair)) {
      stop(paste(
        "categorical conditioning requires family = 'auto', dependence = 'star',",
        "and signed strict pair releases"), call. = FALSE)
    }
    run <- .run(datasources, .aggregate = .aggregate)
    return(.dsvert_mi_joint_star_result_v1(
      formula, data_name, outcomes, m, run, .star_pair,
      root = conditioning, include_root = FALSE))
  }
  if (length(outcomes) > 1L && !identical(family, "auto")) {
    stop(paste(
      "multivariable categorical MI requires family = 'auto' because",
      "each signed marginal determines its own binary or multinomial domain"),
         call. = FALSE)
  }
  if (length(outcomes) == 2L) {
    if (!is.function(.contingency)) {
      stop("Invalid categorical MI joint-pair dependency", call. = FALSE)
    }
    return(.dsvert_mi_joint_pair_result_v1(
      formula, data_name, outcomes, m, datasources, .aggregate, .contingency))
  }
  if (length(outcomes) >= 3L && identical(dependence, "star")) {
    if (!is.function(.run) || !is.function(.star_pair)) {
      stop("Invalid categorical MI signed-star dependency", call. = FALSE)
    }
    run <- .run(datasources, .aggregate = .aggregate)
    return(.dsvert_mi_joint_star_result_v1(
      formula, data_name, outcomes, m, run, .star_pair))
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
#' @description Returns categorical MCAR completions from signed sticky
#'   Synopsis releases. It accepts either \code{outcome ~ 1} or
#'   \code{cbind(outcome1, outcome2, ...) ~ 1}. The two-response form requires
#'   one strict-missing signed joint pair. Larger response sets default to
#'   independent signed marginals; \code{dependence = "star"} instead builds a
#'   conditional star model from signed pairs sharing one vector release. The
#'   same star route accepts one categorical conditioning column on the
#'   right-hand side. It never mutates source tables and its
#'   deterministic completion draws are post-processing of released vectors.
#' @param formula An intercept-only categorical response formula exactly of the
#'   form \code{outcome ~ 1} or \code{cbind(outcome1, outcome2, ...) ~ 1}; with
#'   \code{dependence = "star"}, it may instead have one untransformed
#'   categorical conditioning column on the right-hand side.
#' @param data Signed protected dataset name or federation.
#' @param impute_columns Must be omitted or exactly match the response columns
#'   in formula order.
#' @param m Number of deterministic categorical completion draws, from 2 to 100.
#' @param family One of \code{"auto"}, \code{"binomial"}, or
#'   \code{"multinomial"}. Multi-response routes require \code{"auto"}; the
#'   two-response route derives both domains from its signed pair.
#' @param dependence For three or more intercept-only responses,
#'   \code{"independent"} keeps the protected marginal completion.
#'   \code{"star"} requires one signed strict pair for each non-root response
#'   with the first response and assumes those responses are conditionally
#'   independent given that root. A formula with a conditioning column requires
#'   \code{"star"} and one signed strict pair from that column to every response.
#' @param max_iter,tol,lambda,intercept_only,verbose,seed Compatibility controls.
#'   Only \code{lambda = 0}, \code{intercept_only = "aggregate"}, the default
#'   iteration values, and \code{seed = NULL} are supported.
#' @param datasources DataSHIELD connections.
#' @return A \code{ds.vertMI} object with DP-projected completed-category
#'   probabilities and no classical or Rubin sampling inference. Two-response
#'   calls return a joint categorical probability table. Larger response sets
#'   return protected marginals by default or an explicitly scoped categorical
#'   star model. A conditioning-column call returns only signed conditional
#'   probability matrices and its conditioning marginal, never joint microdata.
#' @seealso \code{\link{ds.vertMethodStatus}}
#' @export
ds.vertMI <- function(formula, data = NULL, impute_columns = NULL,
                      m = 20L, family = c("auto", "binomial", "multinomial"),
                      dependence = c("independent", "star"),
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
  dependence <- match.arg(dependence)
  .dsvert_mi_synopsis_result_v1(
    formula, resolved$value, impute_columns, m, family,
    resolved$datasources, DSI::datashield.aggregate, dependence = dependence)
}

#' @export
print.ds.vertMI <- function(x, ...) {
  cat("dsVert signed categorical MCAR completion\n")
  if (identical(x$method, "signed_categorical_mcar_joint_pair_v3")) {
    cat("  Joint signed pair:", paste(x$impute_columns, collapse = " × "),
        "| missing (DP):", x$missing_count_dp, "\n")
    cat("  No joint microdata or Rubin inference is released.\n")
  } else if (identical(x$method, "signed_categorical_mcar_star_joint_v1")) {
    cat("  Signed categorical star model rooted at:", x$root_column, "\n")
    cat("  Conditional independence given the root; no joint microdata or Rubin inference.\n")
  } else if (identical(x$method, "signed_categorical_mcar_covariate_star_v1")) {
    cat("  Signed categorical conditional-star model given:",
        x$conditioning_column, "\n")
    cat("  No row-level imputations, joint microdata or Rubin inference.\n")
  } else if (is.list(x$variables)) {
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
