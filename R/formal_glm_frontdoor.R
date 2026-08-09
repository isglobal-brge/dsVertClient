# Explicit client gate for the formal binomial/Poisson GLM route.
#
# The public request contains selectors only.  It cannot carry epsilon, delta,
# bounds, contribution caps, random seeds, MPC roles or a numeric backend.  The
# reviewed Go path is not yet registered remotely: the durable Phase-1.9
# schedule now has an internal R bridge, but Phase-1.8 source materialization
# and the common joint-DP finalizer still have no complete R/DSI lifecycle.
# Keep this gate before connection discovery so an incomplete route can never
# fall through to the legacy GLM surface.

.DSVERT_FORMAL_GLM_FRONTDOOR_REQUEST_VERSION <-
  "dsvert-formal-glm-frontdoor-request-v1"
.DSVERT_FORMAL_GLM_FRONTDOOR_BLOCKER <-
  "formal_glm_phase19_durable_r_dsi_release_bridge_not_promoted"
.DSVERT_FORMAL_GLM_FRONTDOOR_MISSING <- c(
  "registered_r_dsi_lifecycle_for_phase18_source_materialization_v1",
  "r_dsi_wrapper_for_phase19_full_durable_schedule_v1",
  "phase19_dp_shares_to_durable_common_one_draw_release_v1",
  "signed_public_release_adapter_bound_to_phase19_validity_v1",
  "local_multiprocess_dsi_e2e_restart_tamper_k2_k3_k4_k5_v1"
)

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

.dsvert_formal_glm_frontdoor_unavailable <- function(request) {
  stop(structure(list(
    message = paste(
      "The formal binomial/Poisson GLM remains sealed before DSI:",
      "its internal durable Phase-1.9 worker is not yet connected through",
      "the registered R/DSI lifecycle to its single common sticky joint-DP",
      "release."),
    call = NULL,
    code = .DSVERT_FORMAL_GLM_FRONTDOOR_BLOCKER,
    missing = .DSVERT_FORMAL_GLM_FRONTDOOR_MISSING,
    request_sha256 = request$sha256,
    dsi_calls = 0L, openings_performed = 0L,
    operation_limit = FALSE, request_limit = FALSE,
    history_can_deny_operation = FALSE,
    production_ready = FALSE),
    class = c("dsvert_formal_glm_frontdoor_unavailable",
              "error", "condition")))
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
  # `datasources` is deliberately not forced here.  In particular, neither
  # connection discovery nor a DataSHIELD aggregate call precedes this gate.
  request <- .dsvert_formal_glm_frontdoor_request(
    analysis_id = analysis_id, data = data, family = family,
    formula = formula)
  .dsvert_formal_glm_frontdoor_unavailable(request)
}
