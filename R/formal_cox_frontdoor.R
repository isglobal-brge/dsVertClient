# Explicit client gate for read-only formal Cox releases.
#
# A formal_analysis_id names a custodian-configured, already completed sticky
# opening. The analyst supplies no privacy control, path, peer, seed, ring or
# optimisation setting, and this route cannot start a Cox computation.

.DSVERT_FORMAL_COX_FRONTDOOR_REQUEST_VERSION <-
  "dsvert-formal-cox-frontdoor-request-v1"

.dsvert_formal_cox_frontdoor_id <- function(value, what) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", value)) {
    stop(what, " must be one bounded ASCII identifier", call. = FALSE)
  }
  enc2utf8(value)
}

.dsvert_formal_cox_frontdoor_formula <- function(formula) {
  if (!inherits(formula, "formula")) {
    stop("formal_analysis_id requires an R formula", call. = FALSE)
  }
  canonical <- paste(deparse(formula, width.cutoff = 500L), collapse = " ")
  canonical <- trimws(gsub("[[:space:]]+", " ", canonical))
  if (!nzchar(canonical) || nchar(canonical, type = "bytes") > 4096L) {
    stop("formal_analysis_id requires a bounded formula", call. = FALSE)
  }
  list(
    canonical = enc2utf8(canonical),
    sha256 = digest::digest(
      paste0("dsVert/formal-cox/frontdoor-formula/v1|", canonical),
      algo = "sha256", serialize = FALSE))
}

.dsvert_formal_cox_frontdoor_request <- function(analysis_id, data, formula) {
  analysis_id <- .dsvert_formal_cox_frontdoor_id(
    analysis_id, "formal_analysis_id")
  data <- .dsvert_formal_cox_frontdoor_id(data, "data")
  formula <- .dsvert_formal_cox_frontdoor_formula(formula)
  request <- .dsvert_joint_dp_client_canonical(list(
    version = .DSVERT_FORMAL_COX_FRONTDOOR_REQUEST_VERSION,
    analysis_id = analysis_id, data_name = data,
    formula_sha256 = formula$sha256,
    public_selectors = as.list(c(
      "analysis_id", "data_name", "formula_sha256")),
    privacy_controls = "server_owned_not_transmitted_v1",
    role_selection = "pinned_server_policy_not_transmitted_v1"))
  list(
    value = request, json = .dsvert_joint_dp_client_json(request),
    sha256 = digest::digest(
      paste0("dsVert/formal-cox/frontdoor-request/v1|",
             .dsvert_joint_dp_client_json(request)),
      algo = "sha256", serialize = FALSE))
}

.dsvert_formal_cox_frontdoor_public_response <- function(value, request) {
  fields <- c("version", "analysis_id", "artifact_id", "certificate_sha256",
              "formula_sha256", "coefficients", "production_ready")
  invalid <- !is.list(value) || is.null(names(value)) || anyNA(names(value)) ||
    anyDuplicated(names(value)) || !setequal(names(value), fields) ||
    !identical(value$version, "dsvert-formal-cox-public-result-v1") ||
    !identical(value$analysis_id, request$analysis_id) ||
    !identical(value$formula_sha256, request$formula_sha256) ||
    !identical(value$production_ready, FALSE) ||
    !is.character(value$artifact_id) || length(value$artifact_id) != 1L ||
    !grepl("^[0-9a-f]{64}$", value$artifact_id) ||
    !is.character(value$certificate_sha256) ||
    length(value$certificate_sha256) != 1L ||
    !grepl("^[0-9a-f]{64}$", value$certificate_sha256) ||
    !is.list(value$coefficients) || !length(value$coefficients)
  if (isTRUE(invalid)) {
    stop("A server returned an invalid formal Cox public release.",
         call. = FALSE)
  }
  coefficients <- lapply(value$coefficients, function(coefficient) {
    expected <- c("coefficient", "beta_steps", "fraction_bits", "beta",
                  "hazard_ratio_lower", "hazard_ratio_upper",
                  "hazard_ratio_midpoint")
    if (!is.list(coefficient) || is.null(names(coefficient)) ||
        anyNA(names(coefficient)) || anyDuplicated(names(coefficient)) ||
        !setequal(names(coefficient), expected) ||
        !is.character(coefficient$coefficient) ||
        length(coefficient$coefficient) != 1L || is.na(coefficient$coefficient) ||
        !grepl("^[A-Za-z][A-Za-z0-9_.]{0,127}$", coefficient$coefficient) ||
        !is.character(coefficient$beta_steps) ||
        length(coefficient$beta_steps) != 1L || is.na(coefficient$beta_steps) ||
        !grepl("^(0|-[1-9][0-9]*|[1-9][0-9]*)$", coefficient$beta_steps) ||
        !is.numeric(coefficient$fraction_bits) ||
        length(coefficient$fraction_bits) != 1L ||
        !is.finite(coefficient$fraction_bits) ||
        coefficient$fraction_bits < 8L || coefficient$fraction_bits > 60L ||
        coefficient$fraction_bits != floor(coefficient$fraction_bits) ||
        !all(vapply(c("beta", "hazard_ratio_lower", "hazard_ratio_upper",
                      "hazard_ratio_midpoint"), function(field) {
          is.numeric(coefficient[[field]]) && length(coefficient[[field]]) == 1L &&
            !is.na(coefficient[[field]]) && is.finite(coefficient[[field]])
        }, logical(1L))) ||
        coefficient$hazard_ratio_lower <= 0 ||
        coefficient$hazard_ratio_upper < coefficient$hazard_ratio_lower ||
        coefficient$hazard_ratio_midpoint < coefficient$hazard_ratio_lower ||
        coefficient$hazard_ratio_midpoint > coefficient$hazard_ratio_upper) {
      stop("A server returned an invalid formal Cox coefficient.",
           call. = FALSE)
    }
    beta <- .dsvert_formal_lattice_value(
      coefficient$beta_steps, coefficient$fraction_bits, coefficient$beta,
      "formal Cox coefficient")
    list(
      coefficient = enc2utf8(coefficient$coefficient),
      beta_steps = coefficient$beta_steps,
      fraction_bits = as.integer(coefficient$fraction_bits),
      beta = beta,
      hazard_ratio_lower = as.numeric(coefficient$hazard_ratio_lower),
      hazard_ratio_upper = as.numeric(coefficient$hazard_ratio_upper),
      hazard_ratio_midpoint = as.numeric(coefficient$hazard_ratio_midpoint))
  })
  names <- vapply(coefficients, `[[`, character(1L), "coefficient")
  if (anyDuplicated(names)) {
    stop("A server returned duplicate formal Cox coefficients.", call. = FALSE)
  }
  list(
    version = value$version, analysis_id = value$analysis_id,
    artifact_id = value$artifact_id,
    certificate_sha256 = value$certificate_sha256,
    formula_sha256 = value$formula_sha256,
    coefficients = coefficients, production_ready = FALSE)
}

.dsvert_formal_cox_frontdoor_public_result <- function(
    request, datasources, .aggregate = DSI::datashield.aggregate) {
  responses <- .dsvert_aggregate_strict(
    conns = datasources,
    expr = call(
      name = "dsvertFormalCoxPublicResultDS",
      analysis_id = request$analysis_id, data_name = request$data_name,
      formula_sha256 = request$formula_sha256),
    operation = "formal Cox public-result retrieval", .aggregate = .aggregate)
  checked <- lapply(
    responses, .dsvert_formal_cox_frontdoor_public_response, request = request)
  canonical <- vapply(checked, .dsvert_joint_dp_client_json, character(1L))
  if (length(unique(canonical)) != 1L) {
    stop("Servers returned different formal Cox public releases.",
         call. = FALSE)
  }
  checked[[1L]]
}

.dsvert_formal_cox_frontdoor_fit <- function(release, request) {
  coefficient_names <- vapply(
    release$coefficients, `[[`, character(1L), "coefficient")
  coefficient_field <- function(field) stats::setNames(
    vapply(release$coefficients, `[[`, if (identical(field, "beta_steps")) {
      character(1L)
    } else if (identical(field, "fraction_bits")) {
      integer(1L)
    } else {
      numeric(1L)
    }, field), coefficient_names)
  result <- list(
    status = "public_certified_dp_point_estimate",
    coefficients = coefficient_field("beta"),
    coefficient_lattice_steps = coefficient_field("beta_steps"),
    fraction_bits = coefficient_field("fraction_bits"),
    hazard_ratio = coefficient_field("hazard_ratio_midpoint"),
    hazard_ratio_lower = coefficient_field("hazard_ratio_lower"),
    hazard_ratio_upper = coefficient_field("hazard_ratio_upper"),
    artifact_id = release$artifact_id,
    certificate_sha256 = release$certificate_sha256,
    formal_analysis_id = request$analysis_id,
    formula_sha256 = request$formula_sha256,
    converged = TRUE, n_obs = NA_real_, n_events = NA_real_,
    iterations = NA_integer_, covariance = NULL, std_errors = NULL,
    source_values_exposed = FALSE, intermediate_values_exposed = FALSE,
    production_ready = FALSE,
    inference = "unavailable_without_a_separate_attested_joint_artifact",
    called_via = "ds.vertCox_formal_analysis_id")
  class(result) <- c("dsvert_formal_dp_cox", "ds.vertCox", "list")
  result
}

.dsvert_formal_cox_frontdoor_adapter <- function(
    explicit_arguments, formula, data, verbose, datasources, analysis_id) {
  allowed <- c("formula", "data", "verbose", "datasources",
               "formal_analysis_id")
  unsupported <- setdiff(explicit_arguments, allowed)
  if (length(unsupported)) {
    stop(paste0(
      "formal_analysis_id does not accept analyst-controlled legacy Cox ",
      "arguments: ", paste(sort(unsupported, method = "radix"),
                            collapse = ", ")),
      call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose) ||
      is.null(formula) || is.null(data)) {
    stop(paste(
      "formal_analysis_id requires an explicit formula, aligned data name",
      "and one non-missing logical verbose value"), call. = FALSE)
  }
  request <- .dsvert_formal_cox_frontdoor_request(analysis_id, data, formula)
  datasources <- .dsvert_dp_datasources(datasources)
  resolved <- .dsvert_federation_argument(data, datasources)
  release <- .dsvert_formal_cox_frontdoor_public_result(
    request$value, resolved$datasources)
  .dsvert_formal_cox_frontdoor_fit(release, request$value)
}
