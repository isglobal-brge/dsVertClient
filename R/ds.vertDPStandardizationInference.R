#' Direct-standardisation inference with DP and sampling uncertainty
#'
#' Builds a conservative joint confidence region for a directly standardised
#' risk from one validated DP strata-by-binary-outcome release and fixed public
#' standard-population weights. It combines the signed simultaneous DP count
#' box with Bonferroni exact Clopper--Pearson intervals for every positive-weight
#' stratum. The final region is the public weighted sum of those stratum
#' intervals.
#'
#' The sampling statement assumes independent privacy units and conditionally
#' binomial outcomes within each public stratum. It conditions on the
#' confidential stratum sample sizes and treats the supplied standard weights
#' as fixed. No independence among stratum intervals is needed for the union
#' bound. A possible zero observed sample in a positive-weight stratum returns
#' an uninformative `[0, 1]` component rather than a fabricated estimate.
#'
#' @param x A released `ds.vertDPContingency` whose rows are strata and whose
#'   two columns are event/non-event outcomes.
#' @param standard_weights Fixed public non-negative weights, optionally named
#'   by the released stratum levels.
#' @param event Event column by name or index.
#' @param level Desired lower bound for joint sampling-plus-mechanism coverage.
#' @param mechanism_alpha_share Fraction of total non-coverage allocated to
#'   the DP mechanism event. The remainder is divided equally across all
#'   positive-weight stratum intervals.
#' @return A `ds.vertDPDirectStandardizationInference` object. It makes no
#'   server call and consumes no additional privacy.
#' @references Clopper, C. J. and Pearson, E. S. (1934). The use of
#'   confidence or fiducial limits illustrated in the case of the binomial.
#'   Biometrika 26(4), 404--413. \doi{10.1093/biomet/26.4.404}.
#' @export
ds.vertDPDirectStandardizationInference <- function(
    x, standard_weights, event = 2L, level = 0.95,
    mechanism_alpha_share = 0.5) {
  x <- .dsvert_dp_table_contract(x)
  if (ncol(x$table) != 2L || nrow(x$table) < 1L) {
    stop("x must be a strata-by-binary-outcome DP table", call. = FALSE)
  }
  values <- c(level = level, mechanism_alpha_share = mechanism_alpha_share)
  if (any(!is.finite(values)) || level <= 0 || level >= 1 ||
      mechanism_alpha_share <= 0 || mechanism_alpha_share >= 1) {
    stop("level and mechanism_alpha_share must be finite numbers in (0, 1)",
         call. = FALSE)
  }
  event_index <- .dsvert_dp_dimension_index(
    event, colnames(x$table), 2L, "event")
  nonevent_index <- setdiff(1:2, event_index)
  alpha <- 1 - level
  mechanism_alpha <- alpha * mechanism_alpha_share
  sampling_alpha <- alpha - mechanism_alpha
  mechanism_level <- 1 - mechanism_alpha
  mechanism <- ds.vertDPDirectStandardization(
    x, standard_weights = standard_weights,
    event = event_index, level = mechanism_level)
  integer_box <- .dsvert_dp_integer_count_box(
    mechanism$count_lower, mechanism$count_upper)
  positive_weight <- which(mechanism$weights > 0)
  if (!length(positive_weight)) {
    stop("At least one standardized stratum must have positive weight",
         call. = FALSE)
  }
  base_sampling_alpha <- sampling_alpha / length(positive_weight)

  stratum_regions <- matrix(
    c(0, 1), nrow = nrow(x$table), ncol = 2L, byrow = TRUE,
    dimnames = list(rownames(x$table), c("lower", "upper")))
  zero_sample_possible <- rep(TRUE, nrow(x$table))
  if (!integer_box$empty) {
    for (index in positive_weight) {
      stratum_regions[index, ] <- .dsvert_dp_cp_union_over_box(
        integer_box$lower[index, event_index],
        integer_box$upper[index, event_index],
        integer_box$lower[index, nonevent_index],
        integer_box$upper[index, nonevent_index], base_sampling_alpha)
    }
    zero_sample_possible <-
      rowSums(integer_box$lower) == 0
  }
  region <- c(
    lower = sum(mechanism$weights * stratum_regions[, "lower"]),
    upper = sum(mechanism$weights * stratum_regions[, "upper"]))
  positive_zero_possible <- any(zero_sample_possible[positive_weight])

  result <- list(
    status = mechanism$status,
    combined_region_status = if (integer_box$empty) {
      "vacuous_empty_integer_mechanism_box"
    } else if (positive_zero_possible) {
      "ok_includes_uninformative_zero_stratum_sample"
    } else {
      "ok"
    },
    estimate = mechanism$estimate,
    combined_region = region,
    mechanism_region = mechanism$mechanism_region,
    stratum_estimates = mechanism$stratum_estimates,
    stratum_combined_regions = stratum_regions,
    positive_weight_strata = positive_weight,
    positive_weight_zero_sample_possible =
      zero_sample_possible[positive_weight],
    confidential_count_integer_box = integer_box,
    weights = mechanism$weights,
    level = level,
    coverage_lower_bound = level,
    mechanism_level = mechanism_level,
    sampling_familywise_level = 1 - sampling_alpha,
    base_sampling_interval_level = 1 - base_sampling_alpha,
    alpha_allocation = c(
      total = alpha, mechanism = mechanism_alpha,
      sampling_familywise = sampling_alpha,
      each_positive_weight_stratum = base_sampling_alpha),
    positive_weight_stratum_count = length(positive_weight),
    coverage_method = paste(
      "union bound over the signed DP simultaneous count-box event and",
      length(positive_weight),
      "Bonferroni exact Clopper-Pearson positive-weight stratum intervals;",
      "the final region is fixed public weighted post-processing; an",
      "integer-empty mechanism box returns [0, 1]"),
    sampling_model = paste(
      "independent privacy units with conditionally binomial outcomes",
      "within each public stratum; confidential stratum sample sizes are",
      "conditioned on and standard-population weights are fixed public inputs"),
    uncertainty_scope =
      "joint DP-mechanism and superpopulation sampling uncertainty",
    inferential_scope = paste(
      "Conservative confidence region for one directly standardized risk;",
      "no hypothesis test, p-value, causal effect or transportability claim"),
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    additional_server_calls = 0L,
    epsilon = x$epsilon,
    delta = x$delta,
    server = x$server)
  class(result) <- c("ds.vertDPDirectStandardizationInference", "list")
  result
}

#' @export
print.ds.vertDPDirectStandardizationInference <- function(x, ...) {
  cat("dsVert DP direct-standardisation sampling inference:",
      x$status, "\n")
  cat("estimate: ", format(x$estimate), " | conservative joint region: [",
      format(x$combined_region[["lower"]]), ", ",
      format(x$combined_region[["upper"]]), "]\n", sep = "")
  cat(x$uncertainty_scope, "\n")
  invisible(x)
}

.dsvert_dp_garwood_count_interval <- function(observed, alpha) {
  values <- c(observed = observed, alpha = alpha)
  if (any(!is.finite(values)) || observed < 0 ||
      observed != floor(observed) || alpha <= 0 || alpha >= 1) {
    stop("Invalid exact-Poisson confidence-interval inputs", call. = FALSE)
  }
  interval <- c(
    lower = if (observed == 0) {
      0
    } else {
      0.5 * stats::qchisq(alpha / 2, 2 * observed)
    },
    upper = 0.5 * stats::qchisq(
      1 - alpha / 2, 2 * (observed + 1)))
  if (anyNA(interval) || interval[["lower"]] < 0 ||
      interval[["lower"]] > interval[["upper"]]) {
    stop("Exact-Poisson confidence limits are not representable",
         call. = FALSE)
  }
  interval
}

#' Indirect-standardisation inference with DP and sampling uncertainty
#'
#' Builds a conservative joint confidence region for an observed-to-expected
#' ratio from one validated DP strata-by-binary-outcome release and fixed public
#' expected rates. It combines the signed simultaneous DP count-box event with
#' one exact central Poisson Garwood interval and envelopes that interval over
#' every integer table in the box. This is client-only post-processing: it makes
#' no server call, draws no randomness, creates no release, and consumes no
#' additional privacy budget.
#'
#' The sampling statement assumes that the confidential total observed count is
#' Poisson with mean equal to the target ratio times a fixed, externally valid
#' expected total. Public stratum rates and confidential stratum populations
#' determine that expected total. If its integer-box range includes zero, or if
#' the mechanism box contains no integer table, the result is the vacuous
#' non-negative parameter region `[0, Inf]`. No causal or transportability claim
#' follows from the DP release or the supplied rates.
#'
#' @param x A released `ds.vertDPContingency` whose rows are strata and whose
#'   two columns are event/non-event outcomes.
#' @param expected_rates Fixed public expected event probabilities, optionally
#'   named by the released stratum levels.
#' @param event Event column by name or index.
#' @param level Desired lower bound for joint sampling-plus-mechanism coverage.
#' @param mechanism_alpha_share Fraction of total non-coverage allocated to the
#'   simultaneous DP mechanism event. The remainder is allocated to the one
#'   exact Poisson interval.
#' @return A `ds.vertDPIndirectStandardizationInference` object containing the
#'   DP point and mechanism region, conservative joint Garwood envelope, typed
#'   denominator states, alpha allocation, and source-release provenance.
#' @references Garwood, F. (1936). Fiducial limits for the Poisson
#'   distribution. \emph{Biometrika}, 28, 437--442.
#' @export
ds.vertDPIndirectStandardizationInference <- function(
    x, expected_rates, event = 2L, level = 0.95,
    mechanism_alpha_share = 0.5) {
  x <- .dsvert_dp_table_contract(x)
  if (ncol(x$table) != 2L || nrow(x$table) < 1L) {
    stop("x must be a strata-by-binary-outcome DP table", call. = FALSE)
  }
  values <- c(level = level, mechanism_alpha_share = mechanism_alpha_share)
  if (any(!is.finite(values)) || level <= 0 || level >= 1 ||
      mechanism_alpha_share <= 0 || mechanism_alpha_share >= 1) {
    stop("level and mechanism_alpha_share must be finite numbers in (0, 1)",
         call. = FALSE)
  }
  event_index <- .dsvert_dp_dimension_index(
    event, colnames(x$table), 2L, "event")
  alpha <- 1 - level
  mechanism_alpha <- alpha * mechanism_alpha_share
  sampling_alpha <- alpha - mechanism_alpha
  mechanism_level <- 1 - mechanism_alpha
  mechanism <- ds.vertDPIndirectStandardization(
    x, expected_rates = expected_rates,
    event = event_index, level = mechanism_level)
  integer_box <- .dsvert_dp_integer_count_box(
    mechanism$count_lower, mechanism$count_upper)

  observed_range <- expected_range <- c(
    lower = NA_real_, upper = NA_real_)
  observed_includes_zero <- expected_includes_zero <- NA
  if (integer_box$empty) {
    combined_region <- c(lower = 0, upper = Inf)
    combined_status <- "vacuous_empty_integer_mechanism_box"
  } else {
    observed_range <- c(
      lower = sum(integer_box$lower[, event_index]),
      upper = sum(integer_box$upper[, event_index]))
    expected_range <- c(
      lower = sum(mechanism$expected_rates * rowSums(integer_box$lower)),
      upper = sum(mechanism$expected_rates * rowSums(integer_box$upper)))
    observed_includes_zero <- observed_range[["lower"]] == 0
    expected_includes_zero <- expected_range[["lower"]] == 0
    if (expected_includes_zero) {
      combined_region <- c(lower = 0, upper = Inf)
      combined_status <- "vacuous_expected_denominator_includes_zero"
    } else {
      lower_count <- .dsvert_dp_garwood_count_interval(
        observed_range[["lower"]], sampling_alpha)
      upper_count <- .dsvert_dp_garwood_count_interval(
        observed_range[["upper"]], sampling_alpha)
      combined_region <- c(
        lower = lower_count[["lower"]] / expected_range[["upper"]],
        upper = upper_count[["upper"]] / expected_range[["lower"]])
      combined_status <- if (observed_includes_zero) {
        "ok_includes_zero_observed_events"
      } else {
        "ok"
      }
    }
  }

  provenance_fields <- c(
    "server", "servers", "datasets", "cross_owner", "capsule_id",
    "manifest_sha256", "final_vector_root", "coordinate_order_sha256",
    "privacy_epoch", "noise_key_id", "privacy_epochs", "noise_key_ids",
    "mechanism", "implementation", "sampler", "epsilon", "delta",
    "implementation_delta", "adjacency", "composition_rule",
    "security_claim")
  source_provenance <- c(
    list(source_class = "ds.vertDPContingency"),
    x[intersect(provenance_fields, names(x))])
  result <- list(
    status = mechanism$status,
    combined_region_status = combined_status,
    estimate = mechanism$estimate,
    observed_events_dp = mechanism$observed_events_dp,
    expected_events_dp = mechanism$expected_events_dp,
    expected_rates = mechanism$expected_rates,
    combined_region = combined_region,
    combined_region_type =
      .dsvert_dp_inference_region_type(combined_region),
    mechanism_region = mechanism$mechanism_region,
    mechanism_region_includes_non_estimable =
      mechanism$mechanism_region_includes_non_estimable,
    mechanism_region_includes_infinite =
      mechanism$mechanism_region_includes_infinite,
    confidential_count_integer_box = integer_box,
    confidential_observed_events_integer_range = observed_range,
    confidential_expected_events_range = expected_range,
    observed_count_box_includes_zero = observed_includes_zero,
    expected_denominator_box_includes_zero = expected_includes_zero,
    level = level,
    coverage_lower_bound = level,
    mechanism_level = mechanism_level,
    poisson_sampling_level = 1 - sampling_alpha,
    sampling_familywise_level = 1 - sampling_alpha,
    base_sampling_interval_level = 1 - sampling_alpha,
    alpha_allocation = c(
      total = alpha, mechanism = mechanism_alpha,
      sampling_familywise = sampling_alpha,
      one_poisson_garwood_interval = sampling_alpha),
    coverage_method = paste(
      "union bound over the signed DP simultaneous count-box event and one",
      "exact central Poisson Garwood interval for the confidential observed",
      "total; the returned region envelopes the Garwood family over every",
      "integer table in the box without assuming independence; an",
      "integer-empty box or possible zero expected denominator returns [0, Inf]"),
    sampling_model = paste(
      "the confidential total observed count is Poisson with mean equal to",
      "the observed-to-expected ratio times a fixed externally valid expected",
      "total formed from public stratum rates and confidential populations"),
    uncertainty_scope =
      "joint DP-mechanism and Poisson superpopulation sampling uncertainty",
    inferential_scope = paste(
      "Conservative confidence region for an observed-to-expected ratio;",
      "no p-value or hypothesis test and no causal effect or",
      "transportability claim"),
    p_values = NULL,
    hypothesis_tests = NULL,
    additional_server_calls = 0L,
    additional_random_draws = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    new_release = FALSE,
    source_release_provenance = source_provenance,
    epsilon = x$epsilon,
    delta = x$delta,
    server = x$server)
  class(result) <- c(
    "ds.vertDPIndirectStandardizationInference", "list")
  result
}

#' @export
print.ds.vertDPIndirectStandardizationInference <- function(x, ...) {
  cat("dsVert DP indirect-standardisation sampling inference:",
      x$status, "\n")
  cat("estimate: ", format(x$estimate), " | conservative joint region: [",
      format(x$combined_region[["lower"]]), ", ",
      format(x$combined_region[["upper"]]), "] (",
      x$combined_region_status, ")\n", sep = "")
  cat(x$uncertainty_scope, "\n")
  invisible(x)
}
