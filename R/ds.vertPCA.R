.dsvert_dp_pca_postprocess <- function(
    cor_result, n_components = NULL, verbose = FALSE,
    synopsis_read_performed = FALSE,
    verification) {
  if (!is.list(verification) ||
      !identical(verification$integrity_valid, TRUE) ||
      !.dsvert_dp_is_string(verification$authenticity)) {
    stop("PCA post-processing requires an explicit validated input",
         call. = FALSE)
  }
  matrix <- as.matrix(cor_result$correlation)
  raw_complete_case <- as.matrix(cor_result$correlation_raw_complete_case)
  complete_case_n <- as.matrix(cor_result$complete_case_n)
  variables <- cor_result$var_names
  p <- length(variables)
  tolerance <- 256 * .Machine$double.eps * max(1, p)
  if (!is.numeric(matrix) || !identical(dim(matrix), c(p, p)) || p < 2L ||
      anyNA(matrix) || any(!is.finite(matrix)) ||
      max(abs(matrix - t(matrix))) > tolerance ||
      max(abs(diag(matrix) - 1)) > tolerance ||
      !identical(rownames(matrix), variables) ||
      !identical(colnames(matrix), variables) ||
      !identical(dim(raw_complete_case), c(p, p)) ||
      !identical(dim(complete_case_n), c(p, p)) ||
      anyNA(raw_complete_case) || any(!is.finite(raw_complete_case)) ||
      anyNA(complete_case_n) || any(!is.finite(complete_case_n)) ||
      any(complete_case_n <= 0) ||
      max(abs(raw_complete_case - t(raw_complete_case))) > tolerance ||
      max(abs(diag(raw_complete_case) - 1)) > tolerance ||
      max(complete_case_n) - min(complete_case_n) >
        tolerance * max(1, max(complete_case_n)) ||
      !identical(dimnames(raw_complete_case), list(variables, variables)) ||
      !identical(dimnames(complete_case_n), list(variables, variables))) {
    stop("The signed DP correlation matrix is invalid", call. = FALSE)
  }
  if (!is.null(n_components) &&
      (!is.numeric(n_components) || length(n_components) != 1L ||
       is.na(n_components) || !is.finite(n_components) ||
       n_components != floor(n_components) ||
       n_components < 1L || n_components > p)) {
    stop("n_components must be an integer between 1 and the variable count",
         call. = FALSE)
  }
  if (is.null(n_components)) n_components <- p
  n_components <- as.integer(n_components)

  decomposition <- eigen((matrix + t(matrix)) / 2, symmetric = TRUE)
  raw_eigenvalues <- decomposition$values
  psd_tolerance <- max(
    tolerance, cor_result$psd_projection$numerical_tolerance)
  if (min(raw_eigenvalues) < -psd_tolerance) {
    stop("The signed DP correlation PSD certificate is inconsistent",
         call. = FALSE)
  }
  eigenvalues <- pmax(0, raw_eigenvalues)
  total <- sum(eigenvalues)
  if (!is.finite(total) || total <= 0) {
    stop("non_identifiable: PCA has no positive finite variance",
         call. = FALSE)
  }
  loadings <- decomposition$vectors
  rownames(loadings) <- variables
  colnames(loadings) <- paste0("PC", seq_len(p))
  names(eigenvalues) <- paste0("PC", seq_len(p))
  variance_pct <- 100 * eigenvalues / total
  cumulative_pct <- cumsum(variance_pct)

  intervals <-
    cor_result$correlation_95_enclosure_raw_estimand_around_projected_release
  lower <- tryCatch(as.matrix(intervals$lower), error = function(error) NULL)
  upper <- tryCatch(as.matrix(intervals$upper), error = function(error) NULL)
  if (!is.numeric(lower) || !is.numeric(upper) ||
      !identical(dim(lower), c(p, p)) || !identical(dim(upper), c(p, p)) ||
      anyNA(lower) || anyNA(upper) || any(!is.finite(lower)) ||
      any(!is.finite(upper)) || any(lower > upper) ||
      any(matrix < lower - tolerance) || any(matrix > upper + tolerance)) {
    stop("The DP correlation mechanism regions are invalid", call. = FALSE)
  }
  entry_radius <- matrix(0, p, p, dimnames = dimnames(matrix))
  for (row in seq_len(p)) {
    for (column in seq_len(p)) {
      below <- .dsvert_dp_vector_next_up(abs(
        matrix[row, column] - lower[row, column]))
      above <- .dsvert_dp_vector_next_up(abs(
        upper[row, column] - matrix[row, column]))
      entry_radius[row, column] <- max(below, above)
    }
  }
  # The Frobenius norm bounds the spectral norm, hence Weyl's eigenvalue
  # perturbation. This is deliberately conservative and does not turn a
  # mechanism region into a sampling confidence interval.
  spectral_squared_upper <- 0
  for (entry in entry_radius) {
    square_upper <- .dsvert_dp_vector_next_up(entry * entry)
    spectral_squared_upper <- .dsvert_dp_vector_next_up(
      spectral_squared_upper + square_upper)
  }
  spectral_radius <- if (spectral_squared_upper == 0) {
    0
  } else {
    .dsvert_dp_vector_next_up(sqrt(spectral_squared_upper))
  }
  eigenvalue_regions <- cbind(
    lower = pmax(2 - p, vapply(eigenvalues, function(value) {
      .dsvert_dp_vector_next_down(value - spectral_radius)
    }, numeric(1L))),
    upper = pmin(p, vapply(eigenvalues, function(value) {
      .dsvert_dp_vector_next_up(value + spectral_radius)
    }, numeric(1L))))
  rownames(eigenvalue_regions) <- names(eigenvalues)
  gaps <- if (p > 1L) eigenvalues[-p] - eigenvalues[-1L] else numeric()
  twice_radius_upper <- .dsvert_dp_vector_next_up(2 * spectral_radius)
  gap_lower <- pmax(0, vapply(gaps, function(value) {
    .dsvert_dp_vector_next_down(value - twice_radius_upper)
  }, numeric(1L)))
  gap_upper <- vapply(gaps, function(value) {
    .dsvert_dp_vector_next_up(value + twice_radius_upper)
  }, numeric(1L))
  names(gaps) <- names(gap_lower) <- names(gap_upper) <-
    paste0("PC", seq_len(length(gaps)), "-PC", seq.int(2L, p))
  stable <- gap_lower > 0
  angle_bound <- rep(NA_real_, length(gaps))
  angle_bound[stable] <- vapply(gap_lower[stable], function(gap) {
    ratio_upper <- .dsvert_dp_vector_next_up(spectral_radius / gap)
    .dsvert_dp_vector_next_up(asin(min(1, ratio_upper)))
  }, numeric(1L))
  names(angle_bound) <- names(gaps)

  inherited <- cor_result[c(
    "artifact_key", "execution_id", "manifest_sha256", "contract_sha256",
    "attempt_sha256", "source_contract_sha256", "result_set_sha256",
    "final_vector_root", "coordinate_order_sha256", "plan_sha256",
    "privacy", "unlimited_replay", "sticky_replay")]
  result <- c(inherited, list(
    loadings = loadings[, seq_len(n_components), drop = FALSE],
    eigenvalues = eigenvalues[seq_len(n_components)],
    variance_pct = variance_pct[seq_len(n_components)],
    cumulative_pct = cumulative_pct[seq_len(n_components)],
    var_names = variables, n_obs = cor_result$n_obs,
    complete_case_n = complete_case_n,
    correlation = matrix,
    correlation_raw_complete_case = raw_complete_case,
    correlation_raw_pairwise = NULL,
    correlation_release = cor_result,
    source_artifact_family = "gaussian_models",
    estimand_missingness = "complete_case_joint",
    pca_eligible = TRUE,
    analysis_id = cor_result$analysis_id,
    epsilon = cor_result$epsilon, delta = cor_result$delta,
    mechanism = cor_result$mechanism,
    method = paste(
      "client-only eigen decomposition of explicitly PSD-projected sticky",
      "DP joint complete-case correlation"),
    scores = NULL,
    scores_available = FALSE,
    source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    additional_server_calls = 0L,
    additional_server_calls_after_synopsis = 0L,
    synopsis_read_performed = synopsis_read_performed,
    synopsis_workflow_count = as.integer(synopsis_read_performed),
    psd_diagnostic = list(
      projection = cor_result$psd_projection,
      minimum_eigenvalue = min(raw_eigenvalues),
      tolerance = psd_tolerance,
      tiny_negative_eigenvalues_clamped =
        as.integer(sum(raw_eigenvalues < 0))),
    eigenvalue_95_mechanism_regions =
      eigenvalue_regions[seq_len(n_components), , drop = FALSE],
    spectral_mechanism_radius_95 = spectral_radius,
    spectral_radius_method = paste(
      "Frobenius upper bound on spectral norm from simultaneous projected",
      "point-versus-raw-estimand entry enclosures, followed by Weyl",
      "perturbation; the complete-case correlation target is joint"),
    loading_identifiability = list(
      sign_indeterminate = TRUE,
      adjacent_observed_gaps = gaps,
      adjacent_gap_95_lower = gap_lower,
      adjacent_gap_95_upper = gap_upper,
      separated_within_mechanism_region = stable,
      davis_kahan_angle_upper_radians = angle_bound,
      individual_directions_identifiable =
        length(stable) == 0L || all(stable[seq_len(min(
          length(stable), n_components))]),
      tied_or_unresolved_eigenspaces_identifiable = TRUE),
    uncertainty_scope = paste(
      "Spectral regions propagate the simultaneous DP-mechanism and",
      "quantization enclosure from the released PSD point to the non-noised",
      "raw complete-case estimand; they exclude sampling uncertainty and do",
      "not bound the projection of an unknown matrix"),
    inferential_scope = paste(
      "Finite-snapshot bounded DP PCA loadings; signs are arbitrary and",
      "directions without a certified eigengap are not individually stable"),
    inference = list(
      p_values = NULL, confidence_intervals = NULL,
      sampling_inference_available = FALSE),
    cross_owner_state = cor_result$cross_owner_state,
    provenance_certificate = cor_result$provenance_certificate,
    provenance_integrity = verification$integrity_valid,
    provenance_authenticity = verification$authenticity,
    legacy_exact_route_called = FALSE,
    disclosure_guard = list(
      satisfied = TRUE,
      basis = "formal_canonical_sticky_DP_synopsis_postprocessing")))
  class(result) <- c("ds.pca", "list")
  if (isTRUE(verbose)) {
    message("DP PCA complete: ", n_components, " components extracted.")
  }
  result
}

#' Principal components from a signed DP correlation artifact
#'
#' Performs client-only eigen decomposition of the explicitly PSD-projected
#' complete-case matrix returned by `ds.vertCor`. It never accepts the
#' pairwise-complete correlation artifact, invokes the former Ring63/exact
#' correlation protocol, or returns individual component scores.
#'
#' @param data_name Signed protected dataset name or reusable
#'   `ds.vertFederation`. Ignored when `cor_result` is supplied.
#' @param variables Optional signed variable subset. Named owner lists are
#'   supported when they agree with the Gaussian artifact.
#' @param n_components Number of components, or `NULL` for all.
#' @param analysis_id Mandatory signed Gaussian artifact id when
#'   `cor_result` is not supplied.
#' @param cor_result An existing complete-case `ds.vertCor` result from the
#'   same sticky Synopsis. Pairwise, arbitrary, and legacy `ds.cor` objects are
#'   rejected.
#' @param verbose Logical progress flag.
#' @param datasources DataSHIELD connections.
#' @return A `ds.pca` object with loadings, eigenvalues, explicit spectral
#'   mechanism diagnostics and inherited DP provenance. Scores are unavailable.
#' @export
ds.vertPCA <- function(data_name = NULL, variables = NULL,
                       n_components = NULL, analysis_id = NULL,
                       cor_result = NULL, verbose = TRUE,
                       datasources = NULL) {
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be one non-missing logical", call. = FALSE)
  }
  synopsis_read_performed <- is.null(cor_result)
  if (is.null(cor_result)) {
    if (inherits(data_name, "ds.vertFederation")) {
      resolved <- .dsvert_federation_argument(data_name, datasources)
      data_name <- resolved$value
      datasources <- resolved$datasources
    }
    data_name <- .dsvert_dp_cor_identifier(data_name, "data_name")
    analysis_id <- .dsvert_dp_cor_identifier(analysis_id, "analysis_id")
    if (isTRUE(verbose)) message("Computing signed sticky DP correlation...")
    cor_result <- ds.vertCor(
      data_name = data_name, variables = variables,
      analysis_id = analysis_id, verbose = verbose,
      datasources = datasources)
  } else {
    if (isTRUE(verbose)) {
      message("Using the provided signed DP correlation artifact...")
    }
  }

  basic_valid <- inherits(cor_result, "ds.vertDPCor") &&
    inherits(cor_result, "ds.cor") &&
    identical(cor_result$released, TRUE) &&
    identical(cor_result$source_values_exposed, FALSE) &&
    identical(cor_result$intermediate_values_exposed, FALSE) &&
    identical(cor_result$legacy_exact_route_called, FALSE) &&
    identical(cor_result$psd_projection_applied, TRUE) &&
    identical(cor_result$source_artifact_family, "gaussian_models") &&
    identical(cor_result$estimand_missingness, "complete_case_joint") &&
    identical(cor_result$pca_eligible, TRUE) &&
    is.null(cor_result$correlation_raw_pairwise) &&
    is.numeric(cor_result$correlation_raw_complete_case) &&
    is.numeric(cor_result$complete_case_n) &&
    .dsvert_dp_is_string(cor_result$analysis_id)
  if (!isTRUE(basic_valid)) {
    stop(paste(
      "cor_result must be an authenticated complete-case ds.vertCor",
      "Gaussian Synopsis artifact; pairwise, legacy, or arbitrary",
      "correlation matrices cannot enter the safe PCA adapter"),
      call. = FALSE)
  }
  verification <- tryCatch(
    ds.validateDPGaussianCertificate(cor_result$provenance_certificate),
    error = function(error) {
      stop("The PCA Gaussian Synopsis certificate is invalid: ",
           conditionMessage(error), call. = FALSE)
    })
  roots <- c(
    "artifact_key", "execution_id", "contract_sha256", "attempt_sha256",
    "source_contract_sha256", "result_set_sha256", "final_vector_root")
  valid <-
    identical(verification$integrity_valid, TRUE) &&
    verification$authenticity %in% c(
      "caller_anchored", "session_transport_anchored") &&
    all(vapply(roots, function(field) {
      identical(cor_result[[field]], verification[[field]])
    }, logical(1L))) &&
    identical(cor_result$coordinate_order_sha256,
              verification$coordinate_order_sha256)
  if (!isTRUE(valid)) {
    stop(paste(
      "cor_result must be an authenticated complete-case ds.vertCor",
      "Gaussian Synopsis artifact; pairwise, legacy, or arbitrary",
      "correlation matrices cannot enter the safe PCA adapter"),
      call. = FALSE)
  }
  .dsvert_dp_cor_gaussian_certificate_match(cor_result, verification)

  return(.dsvert_dp_pca_postprocess(
    cor_result = cor_result, n_components = n_components,
    verbose = verbose,
    synopsis_read_performed = synopsis_read_performed,
    verification = verification))
}

#' @title Print Method for ds.pca Objects
#' @description Prints a summary of signed DP PCA results.
#' @param x A `ds.pca` object.
#' @param ... Additional arguments passed to `print`.
#' @export
print.ds.pca <- function(x, ...) {
  dp_artifact <- is.character(x$analysis_id) &&
    length(x$analysis_id) == 1L && !is.na(x$analysis_id) &&
    nzchar(x$analysis_id)
  if (dp_artifact) {
    cat("Principal Component Analysis (Canonical DP Synopsis)\n")
    cat("====================================================\n\n")
    cat("Signed artifact:", x$analysis_id, "\n")
    cat("DP joint complete-case mass:", x$n_obs, "\n")
  } else {
    cat("Principal Component Analysis\n")
    cat("============================\n\n")
    cat("Observations:", x$n_obs, "\n")
  }
  cat("Variables:", length(x$var_names), "\n")
  cat("Components:", length(x$eigenvalues), "\n")
  if (dp_artifact) cat("Individual scores: unavailable\n")
  cat("\n")
  regions <- x$eigenvalue_95_mechanism_regions
  valid_regions <- is.matrix(regions) && is.numeric(regions) &&
    nrow(regions) == length(x$eigenvalues) &&
    all(c("lower", "upper") %in% colnames(regions))
  if (!valid_regions) {
    regions <- matrix(
      NA_real_, nrow = length(x$eigenvalues), ncol = 2L,
      dimnames = list(NULL, c("lower", "upper")))
  }
  table <- data.frame(
    Component = names(x$variance_pct),
    Eigenvalue = round(x$eigenvalues, 4),
    MechanismLower = round(regions[, "lower"], 4),
    MechanismUpper = round(regions[, "upper"], 4),
    Percent = round(x$variance_pct, 2),
    Cumulative = round(x$cumulative_pct, 2))
  print(table, row.names = FALSE, ...)
  cat("\nLoadings:\n")
  print(round(x$loadings, 3), ...)
  if (is.character(x$uncertainty_scope) &&
      length(x$uncertainty_scope) == 1L && !is.na(x$uncertainty_scope) &&
      nzchar(x$uncertainty_scope)) {
    cat("\n", x$uncertainty_scope, "\n", sep = "")
  }
  invisible(x)
}
