# Client adapter for read-only formal binomial/Poisson GLM releases.
#
# The public request contains selectors only.  It cannot carry epsilon, delta,
# bounds, contribution caps, random seeds, MPC roles or a numeric backend.  The
# route verifies an already completed public certificate or durable Phase21
# terminal remotely. It cannot request source materialisation, start the
# durable worker, or cause another joint-DP opening; those lifecycle stages
# remain outside this frontdoor.

.DSVERT_FORMAL_GLM_FRONTDOOR_REQUEST_VERSION <-
  "dsvert-formal-glm-frontdoor-request-v1"

.dsvert_formal_glm_frontdoor_id <- function(value, what) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", value)) {
    stop(what, " must be one bounded ASCII identifier", call. = FALSE)
  }
  enc2utf8(value)
}

.dsvert_formal_glm_frontdoor_formula <- function(formula) {
  parsed <- .dsvert_plain_formula(formula)
  identifiers <- c(parsed$response, parsed$predictors)
  if (any(!grepl("^[A-Za-z][A-Za-z0-9_.]{0,127}$", identifiers)) ||
      anyDuplicated(parsed$predictors)) {
    stop(paste(
      "The formal GLM formula must use unique bounded additive column",
      "identifiers"), call. = FALSE)
  }
  predictors <- sort(enc2utf8(parsed$predictors), method = "radix")
  canonical <- paste(
    enc2utf8(parsed$response), "~",
    paste(c(if (isTRUE(parsed$intercept)) "1" else "0", predictors),
          collapse = " + "))
  list(
    response = enc2utf8(parsed$response), predictors = predictors,
    intercept = isTRUE(parsed$intercept), canonical = canonical,
    sha256 = digest::digest(
      paste0("dsVert/formal-glm/frontdoor-formula/v1|", canonical),
      algo = "sha256", serialize = FALSE))
}

.dsvert_formal_glm_frontdoor_request <- function(
    analysis_id, data, family, formula) {
  analysis_id <- .dsvert_formal_glm_frontdoor_id(
    analysis_id, "formal_analysis_id")
  data <- .dsvert_formal_glm_frontdoor_id(data, "data")
  if (!is.character(family) || length(family) != 1L || is.na(family) ||
      !family %in% c("binomial", "poisson")) {
    stop(paste(
      "formal_analysis_id supports only family='binomial' or",
      "family='poisson'"), call. = FALSE)
  }
  formula <- .dsvert_formal_glm_frontdoor_formula(formula)
  request <- list(
    version = .DSVERT_FORMAL_GLM_FRONTDOOR_REQUEST_VERSION,
    analysis_id = analysis_id, data_name = data, family = family,
    formula_sha256 = formula$sha256,
    public_selectors = as.list(c(
      "analysis_id", "data_name", "family", "formula_sha256")),
    privacy_controls = "server_owned_not_transmitted_v1",
    role_selection = "pinned_server_policy_not_transmitted_v1")
  request <- .dsvert_joint_dp_client_canonical(request)
  list(
    value = request,
    json = .dsvert_joint_dp_client_json(request),
    sha256 = digest::digest(
      paste0("dsVert/formal-glm/frontdoor-request/v1|",
             .dsvert_joint_dp_client_json(request)),
      algo = "sha256", serialize = FALSE))
}

.dsvert_formal_glm_frontdoor_public_response <- function(value, request) {
  fields <- c(
    "version", "analysis_id", "artifact_id", "certificate_sha256",
    "family", "formula_sha256", "coefficients", "production_ready")
  invalid <- !is.list(value) || is.null(names(value)) || anyNA(names(value)) ||
    anyDuplicated(names(value)) || !setequal(names(value), fields) ||
    !identical(value$version, "dsvert-formal-glm-public-result-v1") ||
    !identical(value$analysis_id, request$analysis_id) ||
    !identical(value$family, request$family) ||
    !identical(value$formula_sha256, request$formula_sha256) ||
    !identical(value$production_ready, FALSE) ||
    !is.character(value$artifact_id) || length(value$artifact_id) != 1L ||
    !grepl("^[0-9a-f]{64}$", value$artifact_id) ||
    !is.character(value$certificate_sha256) ||
    length(value$certificate_sha256) != 1L ||
    !grepl("^[0-9a-f]{64}$", value$certificate_sha256) ||
    !is.list(value$coefficients) || !length(value$coefficients)
  if (isTRUE(invalid)) {
    stop("A server returned an invalid formal GLM public release.",
         call. = FALSE)
  }
  coefficients <- lapply(value$coefficients, function(coefficient) {
    expected <- c(
      "coefficient", "signed_steps", "output_lattice_bits", "value")
    if (!is.list(coefficient) || is.null(names(coefficient)) ||
        anyNA(names(coefficient)) || anyDuplicated(names(coefficient)) ||
        !setequal(names(coefficient), expected) ||
        !is.character(coefficient$coefficient) ||
        length(coefficient$coefficient) != 1L || is.na(coefficient$coefficient) ||
        !nzchar(coefficient$coefficient) ||
        nchar(coefficient$coefficient, type = "bytes") > 1024L ||
        !is.character(coefficient$signed_steps) ||
        length(coefficient$signed_steps) != 1L || is.na(coefficient$signed_steps) ||
        !grepl("^(0|-[1-9][0-9]*|[1-9][0-9]*)$", coefficient$signed_steps) ||
        !is.numeric(coefficient$output_lattice_bits) ||
        length(coefficient$output_lattice_bits) != 1L ||
        !is.finite(coefficient$output_lattice_bits) ||
        coefficient$output_lattice_bits < 0 ||
        coefficient$output_lattice_bits > 127 ||
        coefficient$output_lattice_bits != floor(coefficient$output_lattice_bits) ||
        !is.numeric(coefficient$value) || length(coefficient$value) != 1L ||
        !is.finite(coefficient$value)) {
      stop("A server returned an invalid formal GLM coefficient.",
           call. = FALSE)
    }
    value <- .dsvert_formal_lattice_value(
      coefficient$signed_steps, coefficient$output_lattice_bits,
      coefficient$value, "formal GLM coefficient")
    list(
      coefficient = enc2utf8(coefficient$coefficient),
      signed_steps = coefficient$signed_steps,
      output_lattice_bits = as.integer(coefficient$output_lattice_bits),
      value = value)
  })
  names(coefficients) <- NULL
  names <- vapply(coefficients, `[[`, character(1L), "coefficient")
  if (anyDuplicated(names)) {
    stop("A server returned duplicate formal GLM coefficients.", call. = FALSE)
  }
  list(
    version = value$version, analysis_id = value$analysis_id,
    artifact_id = value$artifact_id,
    certificate_sha256 = value$certificate_sha256,
    family = value$family, formula_sha256 = value$formula_sha256,
    coefficients = coefficients, production_ready = FALSE)
}

.dsvert_formal_glm_frontdoor_public_result <- function(
    request, datasources, .aggregate = DSI::datashield.aggregate) {
  responses <- .dsvert_aggregate_strict(
    conns = datasources,
    expr = call(
      name = "dsvertFormalGLMPublicResultDS",
      analysis_id = request$analysis_id, data_name = request$data_name,
      family = request$family, formula_sha256 = request$formula_sha256),
    operation = "formal GLM public-result retrieval", .aggregate = .aggregate)
  checked <- lapply(
    responses, .dsvert_formal_glm_frontdoor_public_response, request = request)
  canonical <- vapply(checked, .dsvert_joint_dp_client_json, character(1L))
  if (length(unique(canonical)) != 1L) {
    stop("Servers returned different formal GLM public releases.",
         call. = FALSE)
  }
  checked[[1L]]
}

.dsvert_formal_glm_frontdoor_fit <- function(
    release, request, called_via = "ds.vertGLM_formal_analysis_id") {
  coefficient_names <- vapply(
    release$coefficients, `[[`, character(1L), "coefficient")
  coefficients <- stats::setNames(
    vapply(release$coefficients, `[[`, numeric(1L), "value"),
    coefficient_names)
  steps <- stats::setNames(
    vapply(release$coefficients, `[[`, character(1L), "signed_steps"),
    coefficient_names)
  lattice_bits <- stats::setNames(
    vapply(release$coefficients, `[[`, integer(1L), "output_lattice_bits"),
    coefficient_names)
  result <- list(
    status = "public_certified_dp_point_estimate",
    family = release$family,
    coefficients = coefficients,
    coefficient_lattice_steps = steps,
    output_lattice_bits = lattice_bits,
    artifact_id = release$artifact_id,
    certificate_sha256 = release$certificate_sha256,
    formal_analysis_id = request$analysis_id,
    formula_sha256 = request$formula_sha256,
    converged = TRUE,
    deviance = NA_real_, log_likelihood = NA_real_, aic = NA_real_,
    n_obs = NA_real_, covariance = NULL, std_errors = NULL,
    fitted.values = NULL, residuals = NULL,
    source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    production_ready = FALSE,
    inference = "unavailable_without_a_separate_attested_joint_artifact",
    called_via = called_via)
  class(result) <- c("dsvert_formal_dp_glm", "ds.glm", "list")
  result
}

.dsvert_formal_glm_frontdoor_adapter <- function(
    explicit_arguments, formula, data, family, verbose, datasources,
    analysis_id) {
  allowed <- c(
    "formula", "data", "family", "verbose", "datasources",
    "formal_analysis_id")
  unsupported <- setdiff(explicit_arguments, allowed)
  if (length(unsupported)) {
    stop(paste0(
      "formal_analysis_id does not accept analyst-controlled legacy GLM ",
      "arguments: ", paste(sort(unsupported, method = "radix"),
                            collapse = ", ")),
      call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be one non-missing logical value", call. = FALSE)
  }
  if (is.null(formula) || is.null(data)) {
    stop(paste(
      "formal_analysis_id requires an explicit additive formula and",
      "aligned data name"), call. = FALSE)
  }
  request <- .dsvert_formal_glm_frontdoor_request(
    analysis_id = analysis_id, data = data, family = family,
    formula = formula)
  datasources <- .dsvert_dp_datasources(datasources)
  resolved <- .dsvert_federation_argument(data, datasources)
  release <- .dsvert_formal_glm_frontdoor_public_result(
    request$value, resolved$datasources)
  .dsvert_formal_glm_frontdoor_fit(release, request$value)
}

.dsvert_formal_glm_fresh_frontdoor_adapter <- function(
    explicit_arguments, formula, data, family, verbose, datasources,
    analysis_id) {
  .dsvert_block_retired_remote_route("formal_glm_fresh", .allow_test = FALSE)
}
