# Joint complete-case correlation from one signed sticky biomedical Synopsis.
# No protected column discovery or legacy exact/Ring correlation endpoint is
# reachable from this file.

.dsvert_dp_cor_identifier <- function(value, what) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", value)) {
    stop(what, " must be one non-empty capsule identifier", call. = FALSE)
  }
  enc2utf8(value)
}

.dsvert_dp_cor_reserved <- function(message) {
  stop("reserved_not_materialized: ", message, call. = FALSE)
}

.dsvert_signed_workload_condition <- function(
    class, message, method, argument, required_artifact_family,
    analysis_id = NULL, reason) {
  structure(
    list(
      message = message, call = NULL, method = method,
      argument = argument,
      required_artifact_family = required_artifact_family,
      analysis_id = analysis_id, reason = reason),
    class = c(class, "dsvert_signed_workload_error", "error", "condition"))
}

.dsvert_require_signed_analysis_id <- function(
    value, method, argument, required_artifact_family) {
  if (is.null(value)) {
    stop(.dsvert_signed_workload_condition(
      "dsvert_signed_analysis_required",
      paste0(
        method, " requires `", argument, "` for its signed ",
        required_artifact_family, " contract. Ask the server administrator ",
        "to declare the artifact in the capsule workload, then pass its ",
        "analysis_id explicitly; no legacy fallback was attempted."),
      method = method, argument = argument,
      required_artifact_family = required_artifact_family,
      reason = "signed_analysis_id_missing"))
  }
  .dsvert_dp_cor_identifier(value, argument)
}

.dsvert_stop_signed_workload_unavailable <- function(
    method, argument, required_artifact_family, analysis_id, reason,
    detail) {
  stop(.dsvert_signed_workload_condition(
    "dsvert_signed_workload_unavailable",
    paste0(
      method, " cannot use signed analysis_id '", analysis_id, "': ",
      detail, " Ask the server administrator to declare and materialize a ",
      "purpose-bound `", required_artifact_family,
      "` artifact; no Gaussian-correlation, pairwise, or legacy DSI fallback ",
      "was attempted."),
    method = method, argument = argument,
    required_artifact_family = required_artifact_family,
    analysis_id = analysis_id, reason = reason))
}

.dsvert_dp_cor_request <- function(variables, server, datasources) {
  requested_owner <- NULL
  owner_map <- NULL
  owner_count <- 0L
  if (is.list(variables)) {
    if (is.null(names(variables)) || anyNA(names(variables)) ||
        any(!nzchar(names(variables))) || anyDuplicated(names(variables)) ||
        any(!vapply(variables, is.character, logical(1L)))) {
      stop("variables must be a character vector or a named owner list",
           call. = FALSE)
    }
    nonempty <- names(variables)[lengths(variables) > 0L]
    if (!length(nonempty)) {
      stop("variables must request at least two signed variables",
           call. = FALSE)
    }
    owner_count <- length(nonempty)
    requested_owner <- if (owner_count == 1L) nonempty[[1L]] else NULL
    owner_map <- unlist(lapply(nonempty, function(owner) {
      stats::setNames(rep(owner, length(variables[[owner]])),
                      variables[[owner]])
    }), use.names = TRUE)
    variables <- unname(unlist(variables[nonempty], use.names = FALSE))
  } else if (!is.null(variables) && !is.character(variables)) {
    stop("variables must be NULL, a character vector, or a named owner list",
         call. = FALSE)
  }
  if (!is.null(server)) {
    server <- .dsvert_dp_server(server, datasources)
    if (owner_count == 1L && !is.null(requested_owner) &&
        !identical(server, requested_owner)) {
      stop("server and the named variables owner disagree", call. = FALSE)
    }
    requested_owner <- server
  }
  if (!is.null(variables)) {
    variables <- enc2utf8(unname(variables))
    if (length(variables) < 2L || anyNA(variables) ||
        any(!nzchar(variables)) || anyDuplicated(variables)) {
      stop("variables must contain at least two unique non-empty names",
           call. = FALSE)
    }
  }
  list(variables = variables, owner_peer = requested_owner,
       owner_map = owner_map, owner_count = owner_count,
       server = server)
}

.dsvert_dp_cor_artifact <- function(manifest, data_name, analysis_id,
                                    owner_peer = NULL) {
  artifacts <- tryCatch(
    manifest$workload$families$correlation_artifacts,
    error = function(error) NULL)
  artifact <- if (is.list(artifacts)) artifacts[[analysis_id]] else NULL
  required <- c(
    "version", "analysis_id", "dataset", "owner_peer", "variables",
    "pair_references", "pair_count", "coordinate_count",
    "implementation_state", "cross_owner_state")
  if (!is.list(artifact) || is.null(names(artifact)) || anyNA(names(artifact)) ||
      anyDuplicated(names(artifact)) || !setequal(names(artifact), required) ||
      !identical(artifact$version,
                 "same-owner-pairwise-correlation-artifact-v1") ||
      !identical(artifact$analysis_id, analysis_id) ||
      !identical(artifact$dataset, data_name) ||
      !identical(artifact$implementation_state, "same_owner_materialized") ||
      !identical(artifact$cross_owner_state, "reserved_not_materialized")) {
    .dsvert_dp_cor_reserved(paste0(
      "the signed Synopsis has no same-owner correlation artifact '",
      analysis_id, "' for dataset '", data_name, "'"))
  }
  variables <- .dsvert_dp_capsule_manifest_strings(
    artifact$variables, "correlation variables", sorted = TRUE)
  references <- .dsvert_dp_capsule_manifest_strings(
    artifact$pair_references, "correlation pair references", sorted = FALSE)
  pair_count <- length(variables) * (length(variables) - 1) / 2
  if (length(variables) < 2L || length(references) != pair_count ||
      anyDuplicated(references) ||
      !.dsvert_dp_is_integer(artifact$pair_count, pair_count, pair_count) ||
      !.dsvert_dp_is_integer(
        artifact$coordinate_count, 6 * pair_count, 6 * pair_count) ||
      !.dsvert_dp_is_string(artifact$owner_peer) ||
      (!is.null(owner_peer) && !identical(artifact$owner_peer, owner_peer))) {
    if (!is.null(owner_peer) && !identical(
          artifact$owner_peer, owner_peer)) {
      .dsvert_dp_cor_reserved(
        "the requested variables do not belong to the signed artifact owner")
    }
    stop("The signed correlation artifact is invalid", call. = FALSE)
  }
  artifact$variables <- variables
  artifact$pair_references <- references
  artifact
}

.dsvert_dp_cor_pair_key <- function(left, right) {
  pair <- sort(c(left, right), method = "radix")
  paste0(nchar(pair[[1L]], type = "bytes"), ":", pair[[1L]], "|",
         nchar(pair[[2L]], type = "bytes"), ":", pair[[2L]])
}

.dsvert_dp_cor_descriptor <- function(
    block, artifact, scale, capacity, adjacency) {
  descriptor <- block$descriptor
  required <- c(
    "version", "analysis_id", "dataset", "owner_peer", "left", "right",
    "numeric_grid_bits", "coordinate_count", "coordinate_order",
    "repeated_record_policy", "missingness_policy", "statistic_maximum",
    "source_raw_l1_sensitivity", "source_raw_l2_sensitivity",
    "natural_l1_sensitivity", "natural_l2_sensitivity", "adjacency",
    "adjacency_sensitivity_basis")
  sides_valid <- function(side) {
    is.list(side) && !is.null(names(side)) &&
      setequal(names(side), c("column", "lower", "upper")) &&
      .dsvert_dp_is_string(side$column) &&
      .dsvert_dp_is_number(side$lower) &&
      .dsvert_dp_is_number(side$upper) && side$lower < side$upper
  }
  maximum <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    descriptor$statistic_maximum, "numeric-pair maxima"),
    error = function(error) numeric())
  valid <- is.list(descriptor) && !is.null(names(descriptor)) &&
    !anyNA(names(descriptor)) && !anyDuplicated(names(descriptor)) &&
    setequal(names(descriptor), required) &&
    identical(descriptor$version,
              "pairwise-complete-normalized-moments-v1") &&
    identical(descriptor$analysis_id, artifact$analysis_id) &&
    identical(descriptor$dataset, artifact$dataset) &&
    identical(descriptor$owner_peer, artifact$owner_peer) &&
    sides_valid(descriptor$left) && sides_valid(descriptor$right) &&
    descriptor$left$column < descriptor$right$column &&
    all(c(descriptor$left$column, descriptor$right$column) %in%
          artifact$variables) &&
    .dsvert_dp_is_integer(descriptor$numeric_grid_bits, 8, 18) &&
    identical(as.numeric(2^as.numeric(descriptor$numeric_grid_bits)), scale) &&
    .dsvert_dp_is_integer(descriptor$coordinate_count, 6, 6) &&
    identical(
      descriptor$coordinate_order,
      paste0(
        "count,quantized_sum_left,quantized_sum_right,",
        "quantized_sumsq_left,quantized_sumsq_right,",
        "quantized_cross_product")) &&
    identical(
      descriptor$repeated_record_policy,
      paste0(
        "clip_finite_rows_then_mean_each_variable_once_per_",
        "admitted_unit_v1")) &&
    identical(
      descriptor$missingness_policy,
      "pairwise_complete_units_with_both_collapsed_values_v1") &&
    length(maximum) == 6L &&
    identical(maximum, c(capacity, rep(capacity * scale, 5L))) &&
    isTRUE(all.equal(
      as.numeric(descriptor$source_raw_l1_sensitivity),
      1 + 5 * scale, tolerance = 1e-12)) &&
    isTRUE(all.equal(
      as.numeric(descriptor$source_raw_l2_sensitivity),
      sqrt(1 + 5 * scale^2), tolerance = 1e-12)) &&
    identical(as.numeric(descriptor$natural_l1_sensitivity), 6) &&
    isTRUE(all.equal(as.numeric(descriptor$natural_l2_sensitivity),
                     sqrt(6), tolerance = 1e-12)) &&
    identical(descriptor$adjacency, adjacency) &&
    identical(
      descriptor$adjacency_sensitivity_basis,
      paste0(
        "zero_missing_pair_vs_complete_unit_is_worst_case_for_",
        "add_remove_and_replace_one"))
  if (!isTRUE(valid)) {
    stop("The signed numeric-pair descriptor is invalid", call. = FALSE)
  }
  descriptor
}

.dsvert_dp_cor_project_pair <- function(coordinates) {
  names(coordinates) <- c(
    "n", "sum_left", "sum_right", "sumsq_left", "sumsq_right", "cross")
  n <- coordinates[["n"]]
  if (!is.finite(n) || n <= 0) {
    return(list(status = "non_identifiable_nonpositive_dp_pair_count"))
  }
  sum_left <- min(n, max(0, coordinates[["sum_left"]]))
  sum_right <- min(n, max(0, coordinates[["sum_right"]]))
  sumsq_left <- min(
    sum_left, max(sum_left^2 / n, coordinates[["sumsq_left"]]))
  sumsq_right <- min(
    sum_right, max(sum_right^2 / n, coordinates[["sumsq_right"]]))
  cross <- min(
    sum_left, sum_right,
    max(0, sum_left + sum_right - n, coordinates[["cross"]]))
  centered_left <- max(0, sumsq_left - sum_left^2 / n)
  centered_right <- max(0, sumsq_right - sum_right^2 / n)
  denominator <- sqrt(centered_left * centered_right)
  tolerance <- 128 * .Machine$double.eps * max(1, n)
  if (!is.finite(denominator) || denominator <= tolerance) {
    return(list(status = "non_identifiable_zero_dp_pair_variance"))
  }
  centered_cross <- cross - sum_left * sum_right / n
  correlation <- min(1, max(-1, centered_cross / denominator))
  list(
    status = "ok", correlation = correlation, n = n,
    normalized_moments = c(
      sum_left = sum_left, sum_right = sum_right,
      sumsq_left = sumsq_left, sumsq_right = sumsq_right, cross = cross),
    centered_sums = c(
      left = centered_left, right = centered_right,
      cross = centered_cross), denominator = denominator,
    projection = paste(
      "coordinate clamp; bounded univariate moment projection; Frechet",
      "cross clamp; correlation clamp"))
}

.dsvert_dp_cor_interval <- function(
    coordinates, radius, capacity, scale, quantization = NULL) {
  if (!is.numeric(coordinates) || length(coordinates) != 6L ||
      anyNA(coordinates) || any(!is.finite(coordinates)) ||
      !is.numeric(radius) || length(radius) != 1L || is.na(radius) ||
      !is.finite(radius) || radius < 0 ||
      !is.numeric(capacity) || length(capacity) != 1L || is.na(capacity) ||
      !is.finite(capacity) || capacity <= 0 ||
      !is.numeric(scale) || length(scale) != 1L || is.na(scale) ||
      !is.finite(scale) || scale <= 0) {
    stop("Invalid signed correlation interval inputs", call. = FALSE)
  }
  down <- .dsvert_dp_vector_next_down
  up <- .dsvert_dp_vector_next_up
  n_interval <- c(
    lower = max(0, down(coordinates[[1L]] - radius)),
    upper = min(capacity, up(coordinates[[1L]] + radius)))
  if (is.null(quantization)) {
    quantization <- up(up(0.5 * n_interval[["upper"]]) / scale)
  }
  if (!is.numeric(quantization) || length(quantization) != 1L ||
      is.na(quantization) || !is.finite(quantization) || quantization < 0) {
    stop("Invalid signed correlation quantization enclosure",
         call. = FALSE)
  }
  if (quantization > 0) quantization <- up(quantization)
  coordinate_interval <- function(value) c(
    lower = max(0, down(down(value - radius) - quantization)),
    upper = min(
      n_interval[["upper"]], up(up(value + radius) + quantization)))
  x <- coordinate_interval(coordinates[[2L]])
  y <- coordinate_interval(coordinates[[3L]])
  xx <- coordinate_interval(coordinates[[4L]])
  yy <- coordinate_interval(coordinates[[5L]])
  xy <- coordinate_interval(coordinates[[6L]])
  product <- function(left, right) c(
    lower = down(left[["lower"]] * right[["lower"]]),
    upper = up(left[["upper"]] * right[["upper"]]))
  subtract <- function(left, right) c(
    lower = down(left[["lower"]] - right[["upper"]]),
    upper = up(left[["upper"]] - right[["lower"]]))
  numerator <- subtract(product(n_interval, xy), product(x, y))
  left_variance <- subtract(product(n_interval, xx), product(x, x))
  right_variance <- subtract(product(n_interval, yy), product(y, y))
  left_variance <- c(
    lower = max(0, left_variance[["lower"]]),
    upper = max(0, left_variance[["upper"]]))
  right_variance <- c(
    lower = max(0, right_variance[["lower"]]),
    upper = max(0, right_variance[["upper"]]))
  denominator_product_lower <- max(0, down(
    left_variance[["lower"]] * right_variance[["lower"]]))
  denominator_product_upper <- max(0, up(
    left_variance[["upper"]] * right_variance[["upper"]]))
  denominator <- c(
    lower = if (denominator_product_lower == 0) {
      0
    } else {
      down(sqrt(denominator_product_lower))
    },
    upper = up(sqrt(denominator_product_upper)))
  if (!is.finite(denominator[["lower"]]) ||
      denominator[["lower"]] <= 0) {
    correlation <- c(lower = -1, upper = 1)
    status <- "variance_box_reaches_zero_full_correlation_range"
  } else {
    candidates <- c(
      numerator[["lower"]] / denominator[["lower"]],
      numerator[["lower"]] / denominator[["upper"]],
      numerator[["upper"]] / denominator[["lower"]],
      numerator[["upper"]] / denominator[["upper"]])
    correlation <- c(
      lower = max(-1, down(min(candidates))),
      upper = min(1, up(max(candidates))))
    status <- "interval_arithmetic_with_positive_variance_lower_bound"
  }
  list(
    correlation = correlation, n = n_interval,
    normalized_coordinate_intervals = list(
      sum_left = x, sum_right = y, sumsq_left = xx,
      sumsq_right = yy, cross = xy),
    max_abs_quantization_per_sum = quantization,
    status = status)
}

.dsvert_dp_cor_psd <- function(matrix) {
  matrix <- (matrix + t(matrix)) / 2
  diag(matrix) <- 1
  decomposition <- eigen(matrix, symmetric = TRUE)
  clipped <- pmax(0, decomposition$values)
  weighted_vectors <- sweep(
    decomposition$vectors, 2L, clipped, `*`)
  projected <- weighted_vectors %*% t(decomposition$vectors)
  projected <- (projected + t(projected)) / 2
  diagonal <- diag(projected)
  if (anyNA(diagonal) || any(!is.finite(diagonal)) || any(diagonal <= 0)) {
    stop("The DP correlation matrix cannot be projected to a correlation matrix",
         call. = FALSE)
  }
  projected <- projected / sqrt(outer(diagonal, diagonal))
  projected <- (projected + t(projected)) / 2
  diag(projected) <- 1
  dimnames(projected) <- dimnames(matrix)
  final_eigen <- eigen(projected, symmetric = TRUE, only.values = TRUE)$values
  tolerance <- 256 * .Machine$double.eps * nrow(projected)
  if (min(final_eigen) < -tolerance) {
    stop("The DP correlation PSD projection is numerically uncertified",
         call. = FALSE)
  }
  list(
    matrix = projected,
    method = "eigenvalue_clipping_then_diagonal_rescaling_v1",
    exact_nearest_correlation_matrix = FALSE,
    input_min_eigenvalue = min(decomposition$values),
    output_min_eigenvalue = min(final_eigen),
    clipped_eigenvalues = as.integer(sum(decomposition$values < 0)),
    frobenius_distance = sqrt(sum((projected - matrix)^2)),
    numerical_tolerance = tolerance)
}

#' Differentially private same-owner Pearson correlation
#'
#' Reads pairwise-complete bounded moments from exactly one signed, sticky
#' biomedical Synopsis vector. The requested variables must all be co-located
#' in one dataset and owner. Cross-owner products remain deliberately
#' unavailable.
#' The returned raw matrix uses pairwise-complete DP moments and can therefore
#' be indefinite; `correlation` is an explicitly labelled eigenvalue-clipping
#' and diagonal-rescaling projection, not an exact correlation reconstruction
#' or an exact nearest-correlation solution.
#'
#' @param data_name Signed protected dataset name, or a reusable
#'   `ds.vertFederation` returned by `ds.vert.align()`.
#' @param analysis_id Signed correlation artifact id, normally
#'   `paste(data_name, owner, sep = "::")`.
#' @param variables Optional character subset, or a named list containing
#'   exactly one owner entry. At least two variables are required.
#' @param server Optional assertion of the single owner.
#' @param datasources DataSHIELD connections.
#' @return A `ds.vertDPCor`/`ds.cor` object containing raw and projected
#'   correlations, pair counts, simultaneous mechanism/quantization regions,
#'   and signed Synopsis provenance.
#' @export
ds.vertDPCor <- function(data_name, analysis_id, variables = NULL,
                         server = NULL, datasources = NULL) {
  resolved <- .dsvert_federation_argument(data_name, datasources)
  .dsvert_dp_cor_impl(
    resolved$value, analysis_id, variables, server, resolved$datasources,
    DSI::datashield.aggregate)
}

.dsvert_dp_cor_impl <- function(data_name, analysis_id, variables = NULL,
                                server = NULL, datasources = NULL,
                                .aggregate) {
  data_name <- .dsvert_dp_cor_identifier(data_name, "data_name")
  analysis_id <- .dsvert_dp_cor_identifier(analysis_id, "analysis_id")
  datasources <- .dsvert_dp_datasources(datasources)
  request <- .dsvert_dp_cor_request(variables, server, datasources)
  if (request$owner_count > 1L) {
    .dsvert_dp_cor_reserved(paste(
      "the explicit pairwise-complete artifact requires one owner;",
      "use ds.vertCor with a signed complete-case gaussian_models artifact",
      "for cross-owner correlation"))
  }
  run <- .dsvert_dp_synopsis_vector_run(
    datasources, .aggregate = .aggregate)
  context <- .dsvert_dp_vector_context(run, allow_synopsis = TRUE)
  metadata <- .dsvert_dp_vector_public_metadata(context)
  artifact <- .dsvert_dp_cor_artifact(
    context$manifest, data_name, analysis_id, request$owner_peer)
  variables <- if (is.null(request$variables)) {
    artifact$variables
  } else {
    request$variables
  }
  if (!all(variables %in% artifact$variables)) {
    .dsvert_dp_cor_reserved(paste(
      "one or more requested variables are not materialized in the signed",
      "same-owner correlation artifact"))
  }
  owner <- artifact$owner_peer
  scale <- as.numeric(context$lattice$output_lattice_scale)
  count_block <- .dsvert_dp_capsule_single_block(
    context$layout, "admitted_count",
    description = "signed admitted-count capacity block")
  capacity <- .dsvert_dp_vector_block_capacity(count_block)
  pair_blocks <- .dsvert_dp_capsule_vector_blocks(
    context$layout, "numeric_pair_moments", dataset = data_name,
    owner_peer = owner)
  pair_blocks <- pair_blocks[vapply(pair_blocks, function(block) {
    block$key %in% artifact$pair_references
  }, logical(1L))]
  if (length(pair_blocks) != artifact$pair_count ||
      !setequal(vapply(pair_blocks, `[[`, character(1L), "key"),
                artifact$pair_references)) {
    stop("The signed correlation artifact does not match its pair coordinates",
         call. = FALSE)
  }

  block_by_pair <- list()
  for (block in pair_blocks) {
    descriptor <- .dsvert_dp_cor_descriptor(
      block, artifact, scale, capacity, context$adjacency)
    key <- .dsvert_dp_cor_pair_key(
      descriptor$left$column, descriptor$right$column)
    if (!is.null(block_by_pair[[key]])) {
      stop("The signed correlation artifact contains a duplicate pair",
           call. = FALSE)
    }
    block$descriptor <- descriptor
    block_by_pair[[key]] <- block
  }

  p <- length(variables)
  raw <- pair_n <- matrix(
    NA_real_, p, p, dimnames = list(variables, variables))
  diag(raw) <- 1
  interval_lower <- interval_upper <- raw
  diag(interval_lower) <- diag(interval_upper) <- 1
  pair_details <- list()
  requested_pair_count <- p * (p - 1) / 2
  simultaneous <- .dsvert_dp_vector_accuracy_radius(
    context$release, context$manifest,
    coordinate_count = 6 * requested_pair_count,
    confidence = 0.95, maximum_error = capacity)
  for (left_index in seq_len(p - 1L)) {
    for (right_index in seq.int(left_index + 1L, p)) {
      left <- variables[[left_index]]
      right <- variables[[right_index]]
      key <- .dsvert_dp_cor_pair_key(left, right)
      block <- block_by_pair[[key]]
      if (is.null(block)) {
        stop("The signed correlation artifact is missing a requested pair",
             call. = FALSE)
      }
      coordinates <- .dsvert_dp_capsule_vector_values(
        context$release, block)
      if (length(coordinates) != 6L || any(coordinates < 0) ||
          any(coordinates > capacity)) {
        stop("A released numeric-pair block violates its signed bounds",
             call. = FALSE)
      }
      descriptor <- block$descriptor
      canonical <- c(descriptor$left$column, descriptor$right$column)
      projected <- .dsvert_dp_cor_project_pair(coordinates)
      if (!identical(projected$status, "ok")) {
        stop("non_identifiable: pair '", left, "'/'", right,
             "' has ", projected$status, call. = FALSE)
      }
      interval <- .dsvert_dp_cor_interval(
        coordinates, simultaneous$radius, capacity, scale)
      correlation <- projected$correlation
      if (!identical(canonical, c(left, right))) {
        if (!identical(canonical, c(right, left))) {
          stop("The signed numeric-pair orientation changed", call. = FALSE)
        }
      }
      raw[left_index, right_index] <- raw[right_index, left_index] <-
        correlation
      pair_n[left_index, right_index] <- pair_n[right_index, left_index] <-
        projected$n
      interval_lower[left_index, right_index] <-
        interval_lower[right_index, left_index] <-
        interval$correlation[["lower"]]
      interval_upper[left_index, right_index] <-
        interval_upper[right_index, left_index] <-
        interval$correlation[["upper"]]
      pair_details[[key]] <- list(
        variables = canonical, coordinates_dp = unname(coordinates),
        n_dp = projected$n,
        normalized_moments_projected = projected$normalized_moments,
        correlation_raw_pairwise = correlation,
        correlation_95_interval = interval$correlation,
        interval_status = interval$status,
        effective_count_95_interval = interval$n,
        max_abs_quantization_per_sum =
          interval$max_abs_quantization_per_sum,
        repeated_record_policy = descriptor$repeated_record_policy,
        missingness_policy = descriptor$missingness_policy)
    }
  }
  diag(pair_n) <- apply(pair_n, 1L, function(value) {
    values <- value[is.finite(value)]
    if (length(values)) min(values) else NA_real_
  })
  psd <- .dsvert_dp_cor_psd(raw)
  projected_lower <- projected_upper <- psd$matrix
  for (row in seq_len(p)) {
    for (column in seq_len(p)) {
      radius <- max(
        abs(psd$matrix[row, column] - interval_lower[row, column]),
        abs(psd$matrix[row, column] - interval_upper[row, column]))
      projected_lower[row, column] <- max(
        -1, psd$matrix[row, column] - radius)
      projected_upper[row, column] <- min(
        1, psd$matrix[row, column] + radius)
    }
  }
  n_obs <- min(pair_n[upper.tri(pair_n)])
  result <- c(metadata, list(
    status = "ok", analysis_id = analysis_id,
    estimand_missingness = "pairwise_complete",
    source_artifact_family = "correlation_artifacts",
    pca_eligible = FALSE,
    signed_artifact = artifact, server = owner, servers = owner,
    var_names = variables, n_obs = n_obs,
    n_obs_definition = "minimum_DP_noisy_pairwise_complete_count",
    pairwise_n = pair_n,
    correlation = psd$matrix,
    correlation_raw_pairwise = raw,
    local_correlations = stats::setNames(list(psd$matrix), owner),
    method = paste(
      "single-sticky-canonical-Synopsis pairwise-complete Pearson;",
      "explicit PSD post-processing"),
    pairwise_moment_postprocessing = paste(
      "bounded moment/Frechet projection and denominator validation;",
      "no legacy exact correlation route"),
    psd_projection = psd,
    psd_projection_applied = TRUE,
    psd_projection_changes_pairwise_estimand =
      psd$frobenius_distance > psd$numerical_tolerance,
    correlation_95_interval_raw_pairwise = list(
      lower = interval_lower, upper = interval_upper),
    correlation_95_enclosure_raw_estimand_around_projected_release = list(
      lower = projected_lower, upper = projected_upper),
    projected_enclosure_semantics = paste(
      "Simultaneous enclosure of the non-noised raw pairwise estimand,",
      "expressed around the released PSD point; it is not an interval for",
      "the projection of an unknown non-noised matrix"),
    pair_details = pair_details,
    artifact_l1_sensitivity_per_pair = 6,
    artifact_l2_sensitivity_per_pair = sqrt(6),
    accuracy_simultaneous_95_abs_raw_coordinates = simultaneous$radius,
    accuracy_simultaneous_coordinate_count =
      as.integer(6 * requested_pair_count),
    accuracy_simultaneous_confidence = simultaneous$confidence,
    accuracy_simultaneous_method = simultaneous$method,
    accuracy_implementation_tv_upper_bound =
      simultaneous$implementation_tv_upper_bound,
    accuracy_additional_privacy_cost =
      simultaneous$additional_privacy_cost,
    additional_server_calls_after_synopsis = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    quantization_rule =
      "each_normalized_first_second_cross_moment_rounded_to_signed_grid",
    uncertainty_scope = paste(
      "Simultaneous 95% regions cover Synopsis mechanism noise and a",
      "deterministic quantization enclosure; clipping and pairwise missingness",
      "define the finite-snapshot estimand; sampling uncertainty is excluded"),
    inferential_scope = paste(
      "DP-noised bounded pairwise correlations and an explicit PSD",
      "post-processing; no hypothesis test or population confidence interval"),
    cross_owner_state = "reserved_not_materialized",
    legacy_exact_route_called = FALSE))
  class(result) <- c("ds.vertDPCor", "ds.cor", "list")
  result
}

.dsvert_dp_cor_gaussian_coordinates <- function(coordinates, artifact) {
  if (!isTRUE(artifact$intercept)) {
    .dsvert_stop_non_identifiable(
      paste(
        "Pearson correlation is not identifiable from a signed Gaussian",
        "artifact without an intercept because marginal sums are absent."),
      reason = "pearson_requires_intercept_marginal_sums")
  }
  q <- length(artifact$design_terms)
  expected <- q * (q + 1) / 2 + q + 2
  if (!is.numeric(coordinates) || length(coordinates) != expected ||
      anyNA(coordinates) || any(!is.finite(coordinates))) {
    stop("The released complete-case Gaussian coordinates are invalid",
         call. = FALSE)
  }
  cursor <- 2L
  gram <- matrix(
    0, q, q, dimnames = list(artifact$design_terms,
                             artifact$design_terms))
  for (right in seq_len(q)) {
    for (left in seq_len(right)) {
      gram[left, right] <- gram[right, left] <- coordinates[[cursor]]
      cursor <- cursor + 1L
    }
  }
  cross <- coordinates[cursor:(cursor + q - 1L)]
  names(cross) <- artifact$design_terms
  cursor <- cursor + q
  variables <- c(artifact$predictor_order, artifact$outcome$column)
  sums <- c(
    as.numeric(gram[
      "(Intercept)", artifact$predictor_order, drop = TRUE]),
    as.numeric(cross[["(Intercept)"]]))
  names(sums) <- variables
  second <- matrix(
    0, length(variables), length(variables),
    dimnames = list(variables, variables))
  second[artifact$predictor_order, artifact$predictor_order] <-
    gram[artifact$predictor_order, artifact$predictor_order, drop = FALSE]
  second[artifact$predictor_order, artifact$outcome$column] <-
    cross[artifact$predictor_order]
  second[artifact$outcome$column, artifact$predictor_order] <-
    cross[artifact$predictor_order]
  second[artifact$outcome$column, artifact$outcome$column] <-
    coordinates[[cursor]]
  list(
    variables = variables, mass = gram["(Intercept)", "(Intercept)"],
    sums = sums, second = second)
}

.dsvert_dp_cor_gaussian_certificate_match <- function(result, verification) {
  artifact <- verification$artifact
  variables <- result$var_names
  all_variables <- if (is.list(artifact)) {
    c(artifact$predictor_order, artifact$outcome$column)
  } else character()
  roots <- c(
    "artifact_key", "execution_id", "contract_sha256", "attempt_sha256",
    "source_contract_sha256", "result_set_sha256", "final_vector_root")
  scalar_equal <- function(left, right) {
    isTRUE(all.equal(left, right, tolerance = 0, check.attributes = TRUE))
  }
  valid <- is.list(result) && is.list(verification) &&
    identical(verification$integrity_valid, TRUE) &&
    verification$authenticity %in% c(
      "caller_anchored", "session_transport_anchored") &&
    identical(result$analysis_id, verification$analysis_id) &&
    is.character(variables) && length(variables) >= 2L &&
    !anyNA(variables) && !anyDuplicated(variables) &&
    all(variables %in% all_variables) &&
    identical(
      .dsvert_joint_dp_client_json(result$signed_artifact),
      .dsvert_joint_dp_client_json(artifact)) &&
    all(vapply(roots, function(field) {
      identical(result[[field]], verification[[field]])
    }, logical(1L))) &&
    identical(result$coordinate_order_sha256,
              verification$coordinate_order_sha256) &&
    identical(result$provenance_certificate$certificate_sha256,
              verification$certificate$certificate_sha256) &&
    identical(result$cross_owner_state, "reserved_not_materialized") &&
    scalar_equal(as.numeric(result$epsilon),
                 as.numeric(verification$epsilon)) &&
    scalar_equal(as.numeric(result$delta),
                 as.numeric(verification$delta)) &&
    identical(result$mechanism, verification$mechanism) &&
    is.list(verification$validated_moment) &&
    is.numeric(verification$coordinates) &&
    .dsvert_dp_is_number(verification$coordinate_capacity, 0, 2^53 - 1,
                         lower_open = TRUE) &&
    .dsvert_dp_is_number(verification$output_lattice_scale, 0, 2^62,
                         lower_open = TRUE) &&
    is.list(verification$accuracy_simultaneous_95) &&
    .dsvert_dp_is_number(
      verification$accuracy_simultaneous_95$radius, 0, 2^53 - 1)
  if (!isTRUE(valid)) {
    stop("The DP correlation is detached from its Gaussian Synopsis certificate",
         call. = FALSE)
  }

  coordinate_summary <- .dsvert_dp_cor_gaussian_coordinates(
    verification$coordinates, artifact)
  moment <- verification$validated_moment
  augmented_names <- c(artifact$design_terms, artifact$outcome$column)
  augmented <- moment$augmented_projected
  dimnames(augmented) <- list(augmented_names, augmented_names)
  selected <- match(variables, augmented_names)
  intercept <- match("(Intercept)", augmented_names)
  mass <- augmented[intercept, intercept]
  sums <- augmented[intercept, selected]
  second <- augmented[selected, selected, drop = FALSE]
  centered <- second - outer(sums, sums) / mass
  centered <- (centered + t(centered)) / 2
  variances <- diag(centered)
  if (!is.finite(mass) || mass <= 0 || any(!is.finite(variances)) ||
      any(variances <= 0)) {
    stop("The certified Gaussian correlation moments are non-identifiable",
         call. = FALSE)
  }
  raw <- centered / sqrt(outer(variances, variances))
  raw <- (raw + t(raw)) / 2
  diag(raw) <- 1
  dimnames(raw) <- list(variables, variables)
  psd <- .dsvert_dp_cor_psd(raw)

  lower <- upper <- raw
  radius <- verification$accuracy_simultaneous_95$radius
  capacity <- as.numeric(verification$coordinate_capacity)
  scale <- as.numeric(verification$output_lattice_scale)
  for (left_index in seq_len(length(variables) - 1L)) {
    for (right_index in seq.int(left_index + 1L, length(variables))) {
      left <- variables[[left_index]]
      right <- variables[[right_index]]
      interval <- .dsvert_dp_cor_interval(c(
        coordinate_summary$mass,
        coordinate_summary$sums[[left]], coordinate_summary$sums[[right]],
        coordinate_summary$second[left, left],
        coordinate_summary$second[right, right],
        coordinate_summary$second[left, right]),
        radius, capacity, scale)
      lower[left_index, right_index] <- lower[right_index, left_index] <-
        interval$correlation[["lower"]]
      upper[left_index, right_index] <- upper[right_index, left_index] <-
        interval$correlation[["upper"]]
    }
  }
  projected_lower <- projected_upper <- psd$matrix
  for (row in seq_along(variables)) {
    for (column in seq_along(variables)) {
      enclosure_radius <- max(
        abs(psd$matrix[row, column] - lower[row, column]),
        abs(psd$matrix[row, column] - upper[row, column]))
      projected_lower[row, column] <- max(
        -1, psd$matrix[row, column] - enclosure_radius)
      projected_upper[row, column] <- min(
        1, psd$matrix[row, column] + enclosure_radius)
    }
  }
  complete_case_n <- matrix(
    mass, length(variables), length(variables),
    dimnames = list(variables, variables))
  bound <-
    scalar_equal(result$n_obs, mass) &&
    scalar_equal(result$complete_case_n, complete_case_n) &&
    scalar_equal(result$correlation_raw_complete_case, raw) &&
    scalar_equal(result$correlation, psd$matrix) &&
    scalar_equal(result$psd_projection, psd) &&
    scalar_equal(
      result$correlation_95_interval_complete_case,
      list(lower = lower, upper = upper)) &&
    scalar_equal(
      result$correlation_95_enclosure_raw_estimand_around_projected_release,
      list(lower = projected_lower, upper = projected_upper))
  if (!isTRUE(bound)) {
    stop("The DP correlation values do not match their Gaussian Synopsis certificate",
         call. = FALSE)
  }
  invisible(TRUE)
}

.dsvert_dp_cor_gaussian_impl <- function(
    data_name, analysis_id, variables = NULL, server = NULL,
    datasources = NULL, .aggregate) {
  data_name <- .dsvert_dp_cor_identifier(data_name, "data_name")
  analysis_id <- .dsvert_dp_cor_identifier(analysis_id, "analysis_id")
  datasources <- .dsvert_dp_datasources(datasources)
  request <- .dsvert_dp_cor_request(variables, server, datasources)
  source <- .dsvert_dp_gaussian_synopsis_release(
    data_name, analysis_id, request$server, datasources, .aggregate)
  context <- source$context
  scale <- source$scale
  capacity <- source$capacity
  artifact <- source$artifact
  all_variables <- c(artifact$predictor_order, artifact$outcome$column)
  variables <- if (is.null(request$variables)) {
    all_variables
  } else request$variables
  if (length(variables) < 2L || !all(variables %in% all_variables)) {
    stop(
      "Every requested correlation variable must belong to the signed ",
      "complete-case Gaussian artifact", call. = FALSE)
  }
  if (!is.null(request$owner_map)) {
    expected_owners <- stats::setNames(vapply(all_variables, function(name) {
      if (identical(name, artifact$outcome$column)) {
        artifact$outcome$owner_peer %||% artifact$owner_peer
      } else {
        artifact$predictors[[name]]$owner_peer %||% artifact$owner_peer
      }
    }, character(1L)), all_variables)
    if (!identical(unname(request$owner_map[variables]),
                   unname(expected_owners[variables]))) {
      stop("The named variable owners disagree with the signed artifact",
           call. = FALSE)
    }
  }
  coordinates <- source$coordinates
  coordinate_summary <- .dsvert_dp_cor_gaussian_coordinates(
    coordinates, artifact)
  moment <- source$moment
  augmented_names <- c(artifact$design_terms, artifact$outcome$column)
  augmented <- moment$augmented_projected
  dimnames(augmented) <- list(augmented_names, augmented_names)
  selected <- match(variables, augmented_names)
  intercept <- match("(Intercept)", augmented_names)
  mass <- augmented[intercept, intercept]
  tolerance <- 256 * .Machine$double.eps * max(1, abs(mass)) *
    length(variables)
  if (!is.finite(mass) || mass <= tolerance) {
    .dsvert_stop_non_identifiable(
      "The projected DP complete-case moment mass is non-positive.",
      reason = "non_positive_dp_complete_case_moment_mass")
  }
  sums <- augmented[intercept, selected]
  second <- augmented[selected, selected, drop = FALSE]
  centered <- second - outer(sums, sums) / mass
  centered <- (centered + t(centered)) / 2
  variances <- diag(centered)
  variance_tolerance <- 256 * .Machine$double.eps *
    max(1, max(abs(centered))) * length(variables)
  if (any(!is.finite(variances)) || any(variances <= variance_tolerance)) {
    .dsvert_stop_non_identifiable(
      "At least one projected DP complete-case variance is zero.",
      reason = "zero_dp_complete_case_variance")
  }
  raw <- centered / sqrt(outer(variances, variances))
  raw <- (raw + t(raw)) / 2
  diag(raw) <- 1
  dimnames(raw) <- list(variables, variables)
  psd <- .dsvert_dp_cor_psd(raw)

  simultaneous <- .dsvert_dp_vector_accuracy_radius(
    context$release, context$manifest,
    coordinate_count = artifact$coordinate_count,
    confidence = 0.95, maximum_error = capacity)
  quantization <- if (identical(
      artifact$version,
      .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_ARTIFACT_VERSION)) {
    as.numeric(
      artifact$quantization_contract$per_sum_max_abs_error_lattice_steps) /
      scale
  } else NULL
  lower <- upper <- raw
  pair_details <- list()
  for (left_index in seq_len(length(variables) - 1L)) {
    for (right_index in seq.int(left_index + 1L, length(variables))) {
      left <- variables[[left_index]]
      right <- variables[[right_index]]
      interval <- .dsvert_dp_cor_interval(c(
        coordinate_summary$mass,
        coordinate_summary$sums[[left]], coordinate_summary$sums[[right]],
        coordinate_summary$second[left, left],
        coordinate_summary$second[right, right],
        coordinate_summary$second[left, right]), simultaneous$radius,
        capacity, scale, quantization = quantization)
      lower[left_index, right_index] <- lower[right_index, left_index] <-
        interval$correlation[["lower"]]
      upper[left_index, right_index] <- upper[right_index, left_index] <-
        interval$correlation[["upper"]]
      pair_details[[.dsvert_dp_cor_pair_key(left, right)]] <- list(
        variables = c(left, right),
        complete_case_moment_mass_dp = coordinate_summary$mass,
        correlation_raw_complete_case = raw[left, right],
        correlation_95_interval = interval$correlation,
        interval_status = interval$status,
        missingness_policy = artifact$missingness_policy)
    }
  }
  projected_lower <- projected_upper <- psd$matrix
  for (row in seq_along(variables)) {
    for (column in seq_along(variables)) {
      radius <- max(
        abs(psd$matrix[row, column] - lower[row, column]),
        abs(psd$matrix[row, column] - upper[row, column]))
      projected_lower[row, column] <- max(
        -1, psd$matrix[row, column] - radius)
      projected_upper[row, column] <- min(
        1, psd$matrix[row, column] + radius)
    }
  }
  participants <- artifact$participating_peers %||% artifact$owner_peer
  computation <- artifact$computation_peers %||% character()
  complete_case_n <- matrix(
    mass, length(variables), length(variables),
    dimnames = list(variables, variables))
  result <- c(source$metadata, list(
    status = "ok", analysis_id = analysis_id,
    signed_artifact = artifact,
    source_artifact_family = "gaussian_models",
    estimand_missingness = "complete_case_joint",
    pca_eligible = TRUE,
    server = artifact$owner_peer, servers = participants,
    participating_servers = participants,
    computation_servers = computation,
    var_names = variables, n_obs = mass,
    n_obs_definition =
      "projected_DP_Gaussian_intercept_moment_mass_used_by_correlation",
    complete_case_n = complete_case_n,
    complete_case_moment_mass_dp = mass,
    complete_case_count_coordinate_dp = moment$n,
    complete_case_gram_intercept_released_dp = coordinate_summary$mass,
    pairwise_n = NULL,
    correlation = psd$matrix,
    correlation_raw_complete_case = raw,
    correlation_raw_pairwise = NULL,
    local_correlations = list(joint_complete_case = psd$matrix),
    method = paste(
      "single-sticky-canonical-Synopsis complete-case Pearson from signed",
      "gaussian_models sufficient statistics; explicit PSD post-processing"),
    psd_projection = psd, psd_projection_applied = TRUE,
    psd_projection_changes_complete_case_estimand =
      psd$frobenius_distance > psd$numerical_tolerance,
    correlation_95_interval_complete_case = list(
      lower = lower, upper = upper),
    correlation_95_enclosure_raw_estimand_around_projected_release = list(
      lower = projected_lower, upper = projected_upper),
    pair_details = pair_details,
    accuracy_simultaneous_95_abs_raw_coordinates = simultaneous$radius,
    accuracy_simultaneous_coordinate_count = artifact$coordinate_count,
    accuracy_simultaneous_confidence = simultaneous$confidence,
    accuracy_simultaneous_method = simultaneous$method,
    accuracy_implementation_tv_upper_bound =
      simultaneous$implementation_tv_upper_bound,
    accuracy_sampler_tv_upper_bound = simultaneous$sampler_tv_upper_bound,
    accuracy_additional_privacy_cost =
      simultaneous$additional_privacy_cost,
    uncertainty_scope = paste(
      "Simultaneous mechanism/quantization regions for the bounded",
      "finite-snapshot complete-case estimand; no population sampling",
      "confidence interval is claimed"),
    inferential_scope = paste(
      "DP-noised complete-case correlations only; no hypothesis tests,",
      "p-values, or sampling inference"),
    inference = list(
      p_values = NULL, confidence_intervals = NULL,
      sampling_inference_available = FALSE),
    additional_server_calls_after_synopsis = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    cross_owner_state = artifact$cross_owner_state,
    provenance_certificate = source$certificate,
    provenance_integrity = source$verification$integrity_valid,
    provenance_authenticity = source$verification$authenticity,
    legacy_exact_route_called = FALSE,
    disclosure_guard = list(
      satisfied = TRUE,
      basis = "formal_canonical_sticky_DP_synopsis_postprocessing")))
  class(result) <- c("ds.vertDPCor", "ds.cor", "list")
  result
}

#' Disclosure-safe compatibility adapter for correlation
#'
#' This existing name consumes only a signed `gaussian_models` sufficient-
#' statistic Synopsis artifact. A signed `analysis_id` is mandatory and every
#' requested variable must be same-owner. Cross-owner descriptors fail closed;
#' there is no silent fallback to either pairwise moments or a legacy capsule.
#'
#' @param data_name Signed protected dataset name.
#' @param variables Optional character subset, or a named list identifying the
#'   signed owner of each requested variable. At least two variables are
#'   required.
#' @param analysis_id Mandatory signed Gaussian artifact id.
#' @param verbose Logical progress flag retained for compatibility.
#' @param datasources DataSHIELD connections.
#' @return A `ds.vertDPCor`/`ds.cor` object.
#' @export
ds.vertCor <- function(data_name, variables = NULL, analysis_id = NULL,
                       verbose = TRUE, datasources = NULL) {
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be one non-missing logical", call. = FALSE)
  }
  resolved <- .dsvert_federation_argument(data_name, datasources)
  analysis_id <- .dsvert_dp_cor_identifier(analysis_id, "analysis_id")
  if (isTRUE(verbose)) {
    message("Reading signed complete-case DP Gaussian artifact '",
            analysis_id, "'...")
  }
  result <- .dsvert_dp_cor_gaussian_impl(
    data_name = resolved$value, analysis_id = analysis_id,
    variables = variables, datasources = resolved$datasources,
    .aggregate = DSI::datashield.aggregate)
  if (isTRUE(verbose)) message("DP correlation post-processing complete.")
  result
}

#' @title Print Method for ds.cor Objects
#' @description Prints a summary of the DP correlation result.
#' @param x A `ds.cor` object.
#' @param digits Number of displayed digits.
#' @param ... Additional arguments passed to `print`.
#' @export
print.ds.cor <- function(x, digits = 3, ...) {
  complete_case <- identical(x$estimand_missingness, "complete_case_joint")
  heading <- if (complete_case) {
    "Complete-case Pearson Correlation (Sticky DP Synopsis)"
  } else {
    "Pairwise Pearson Correlation (Sticky Joint DP Capsule)"
  }
  cat(heading, "\n", strrep("=", nchar(heading)), "\n\n", sep = "")
  cat("Method:", x$method, "\n")
  cat("Signed artifact:", x$analysis_id, "\n")
  if (complete_case) {
    cat("DP joint complete-case mass:", x$n_obs, "\n")
  } else {
    cat("Minimum DP pair count:", x$n_obs, "\n")
  }
  cat("Variables:", length(x$var_names), "\n")
  cat("Participating peers:", paste(x$servers, collapse = ", "), "\n")
  cat("PSD projection:", x$psd_projection$method,
      "(not exact nearest-correlation)\n\n")
  print(round(x$correlation, digits), ...)
  cat("\n", x$uncertainty_scope, "\n", sep = "")
  invisible(x)
}
