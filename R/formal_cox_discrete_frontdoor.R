# Read-only client boundary for a completed formal discrete-time hazard model.
#
# The public certificate is a binomial formal-GLM release over the server-owned
# person-period design.  The dedicated DS endpoint binds it to the caller's
# Surv() formula and fixed time grid, so this is not a Cox PH fallback and it
# never invokes the retained pooled-logistic protocol.

.DSVERT_FORMAL_COX_DISCRETE_FRONTDOOR_REQUEST_VERSION <-
  "dsvert-formal-cox-discrete-frontdoor-request-v1"

.dsvert_formal_cox_discrete_frontdoor_request <- function(
    analysis_id, data, formula) {
  analysis_id <- .dsvert_formal_cox_frontdoor_id(
    analysis_id, "formal_analysis_id")
  data <- .dsvert_formal_cox_frontdoor_id(data, "data")
  formula <- .dsvert_formal_cox_frontdoor_formula(formula)
  request <- .dsvert_joint_dp_client_canonical(list(
    version = .DSVERT_FORMAL_COX_DISCRETE_FRONTDOOR_REQUEST_VERSION,
    analysis_id = analysis_id, data_name = data,
    source_formula_sha256 = formula$sha256,
    public_selectors = as.list(c(
      "analysis_id", "data_name", "source_formula_sha256")),
    privacy_controls = "server_owned_not_transmitted_v1",
    role_selection = "pinned_server_policy_not_transmitted_v1"))
  list(
    value = request, json = .dsvert_joint_dp_client_json(request),
    sha256 = digest::digest(
      paste0("dsVert/formal-cox-discrete/frontdoor-request/v1|",
             .dsvert_joint_dp_client_json(request)),
      algo = "sha256", serialize = FALSE))
}

.dsvert_formal_cox_discrete_frontdoor_response <- function(value, request) {
  fields <- c(
    "version", "analysis_id", "artifact_id", "certificate_sha256",
    "target", "source_formula_sha256", "model_formula_sha256",
    "time_grid_sha256", "coefficients", "production_ready")
  invalid <- !is.list(value) || is.null(names(value)) || anyNA(names(value)) ||
    anyDuplicated(names(value)) || !setequal(names(value), fields) ||
    !identical(value$version,
               "dsvert-formal-cox-discrete-public-result-v1") ||
    !identical(value$analysis_id, request$analysis_id) ||
    !identical(value$target, "discrete_logit") ||
    !identical(value$source_formula_sha256, request$source_formula_sha256) ||
    !identical(value$production_ready, FALSE) ||
    !is.character(value$artifact_id) || length(value$artifact_id) != 1L ||
    !grepl("^[0-9a-f]{64}$", value$artifact_id) ||
    !is.character(value$certificate_sha256) ||
    length(value$certificate_sha256) != 1L ||
    !grepl("^[0-9a-f]{64}$", value$certificate_sha256) ||
    !is.character(value$model_formula_sha256) ||
    length(value$model_formula_sha256) != 1L ||
    !grepl("^[0-9a-f]{64}$", value$model_formula_sha256) ||
    !is.character(value$time_grid_sha256) ||
    length(value$time_grid_sha256) != 1L ||
    !grepl("^[0-9a-f]{64}$", value$time_grid_sha256) ||
    !is.list(value$coefficients) || !length(value$coefficients)
  if (isTRUE(invalid)) {
    stop("A server returned an invalid formal discrete-time public release.",
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
        !grepl("^(0|-[1-9][0-9]*|[1-9][0-9]*)$",
               coefficient$signed_steps) ||
        !is.numeric(coefficient$output_lattice_bits) ||
        length(coefficient$output_lattice_bits) != 1L ||
        !is.finite(coefficient$output_lattice_bits) ||
        coefficient$output_lattice_bits < 0 ||
        coefficient$output_lattice_bits > 127 ||
        coefficient$output_lattice_bits != floor(coefficient$output_lattice_bits) ||
        !is.numeric(coefficient$value) || length(coefficient$value) != 1L ||
        !is.finite(coefficient$value) ||
        abs(coefficient$value) > log(.Machine$double.xmax)) {
      stop("A server returned an invalid formal discrete-time coefficient.",
           call. = FALSE)
    }
    list(
      coefficient = enc2utf8(coefficient$coefficient),
      signed_steps = coefficient$signed_steps,
      output_lattice_bits = as.integer(coefficient$output_lattice_bits),
      value = as.numeric(coefficient$value))
  })
  names <- vapply(coefficients, `[[`, character(1L), "coefficient")
  if (anyDuplicated(names)) {
    stop("A server returned duplicate formal discrete-time coefficients.",
         call. = FALSE)
  }
  list(
    version = value$version, analysis_id = value$analysis_id,
    artifact_id = value$artifact_id,
    certificate_sha256 = value$certificate_sha256,
    target = value$target,
    source_formula_sha256 = value$source_formula_sha256,
    model_formula_sha256 = value$model_formula_sha256,
    time_grid_sha256 = value$time_grid_sha256,
    coefficients = coefficients, production_ready = FALSE)
}

.dsvert_formal_cox_discrete_frontdoor_public_result <- function(
    request, datasources, .aggregate = DSI::datashield.aggregate) {
  responses <- .dsvert_aggregate_strict(
    conns = datasources,
    expr = call(
      name = "dsvertFormalCoxDiscretePublicResultDS",
      analysis_id = request$analysis_id, data_name = request$data_name,
      source_formula_sha256 = request$source_formula_sha256),
    operation = "formal discrete-time public-result retrieval",
    .aggregate = .aggregate)
  checked <- lapply(
    responses, .dsvert_formal_cox_discrete_frontdoor_response,
    request = request)
  canonical <- vapply(checked, .dsvert_joint_dp_client_json, character(1L))
  if (length(unique(canonical)) != 1L) {
    stop("Servers returned different formal discrete-time public releases.",
         call. = FALSE)
  }
  checked[[1L]]
}

.dsvert_formal_cox_discrete_frontdoor_fit <- function(release, request) {
  coefficient_names <- vapply(
    release$coefficients, `[[`, character(1L), "coefficient")
  coefficients <- stats::setNames(
    vapply(release$coefficients, `[[`, numeric(1L), "value"),
    coefficient_names)
  result <- list(
    status = "public_certified_dp_point_estimate",
    target = "discrete_logit", coefficients = coefficients,
    hazard_odds_ratio = exp(coefficients),
    coefficient_lattice_steps = stats::setNames(
      vapply(release$coefficients, `[[`, character(1L), "signed_steps"),
      coefficient_names),
    output_lattice_bits = stats::setNames(
      vapply(release$coefficients, `[[`, integer(1L), "output_lattice_bits"),
      coefficient_names),
    artifact_id = release$artifact_id,
    certificate_sha256 = release$certificate_sha256,
    formal_analysis_id = request$analysis_id,
    source_formula_sha256 = request$source_formula_sha256,
    model_formula_sha256 = release$model_formula_sha256,
    time_grid_sha256 = release$time_grid_sha256,
    converged = TRUE, covariance = NULL, std_errors = NULL,
    source_values_exposed = FALSE, intermediate_values_exposed = FALSE,
    production_ready = FALSE,
    inference = "unavailable_without_a_separate_attested_joint_artifact",
    called_via = "ds.vertCoxDiscreteNonDisclosive_formal_analysis_id")
  class(result) <- c(
    "dsvert_formal_dp_cox_discrete",
    "ds.vertCoxDiscreteNonDisclosive", "list")
  result
}

.dsvert_formal_cox_discrete_frontdoor_adapter <- function(
    explicit_arguments, formula, data, verbose, datasources, analysis_id) {
  allowed <- c("formula", "data", "verbose", "datasources",
               "formal_analysis_id")
  unsupported <- setdiff(explicit_arguments, allowed)
  if (length(unsupported)) {
    stop(paste0(
      "formal_analysis_id does not accept analyst-controlled discrete-time ",
      "arguments: ", paste(sort(unsupported, method = "radix"),
                            collapse = ", ")),
      call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose) ||
      is.null(formula) || is.null(data)) {
    stop(paste(
      "formal_analysis_id requires an explicit Surv formula, aligned data ",
      "name and one non-missing logical verbose value"), call. = FALSE)
  }
  request <- .dsvert_formal_cox_discrete_frontdoor_request(
    analysis_id, data, formula)
  datasources <- .dsvert_dp_datasources(datasources)
  resolved <- .dsvert_federation_argument(data, datasources)
  release <- .dsvert_formal_cox_discrete_frontdoor_public_result(
    request$value, resolved$datasources)
  .dsvert_formal_cox_discrete_frontdoor_fit(release, request$value)
}
