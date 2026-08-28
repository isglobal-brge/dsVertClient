# Local fail-closed boundary for public names whose historical server surface
# is deliberately no longer registered.  The compatibility implementation is
# retained below the frontdoor for source-level migration tests only; an
# installed package has no option, environment variable or request flag that
# can re-enable it.

.dsvert_quarantine_test_mode <- function() FALSE

.DSVERT_RETIRED_REMOTE_ROUTES <- list(
  legacy_glm = list(
    method = "ds.vertGLM",
    state = "legacy_exact_glm_route_removed",
    replacement = paste0(
      "For Gaussian models use ds.vertGLM(..., dp_analysis_id=...) or ",
      "ds.vertDPGaussian(). For binomial/Poisson use ds.vertGLM(..., ",
      "analysis_id=...) for one signed finite grid or formal_analysis_id=...) ",
      "to read a completed certificate.")),
  formal_glm_fresh = list(
    method = "ds.vertGLM",
    state = "formal_glm_fresh_computation_sealed",
    replacement = paste0(
      "Use ds.vertGLM(..., analysis_id=...) for one signed finite ",
      "binomial/Poisson grid or formal_analysis_id=...) to read a completed ",
      "two-authority-signed certificate. Fresh formal GLM computation remains ",
      "sealed pending a production-attested protected runtime.")),
  cox = list(
    method = "ds.vertCox",
    state = "legacy_cox_route_quarantined",
    replacement = paste0(
      "Use ds.vertCox(..., formal_analysis_id=...) to read one completed ",
      "two-authority-signed certificate. Fresh Cox computation remains ",
      "sealed pending a production-attested protected runtime; classical ",
      "inference remains unavailable.")),
  negative_binomial = list(
    method = "ds.vertNBFullRegTheta",
    state = "legacy_nb2_route_quarantined",
    replacement = paste0(
      "For a bounded non-negative integer y ~ 1 domain, pass a validated ",
      "ds.vertDPFrequency object to ds.vertNBFullRegTheta(..., frequency = ",
      "...). For additive bare covariates, use one same-owner signed finite ",
      "beta/theta grid selected by analysis_id; arbitrary optimization and ",
      "sampling inference remain unavailable.")),
  multinomial = list(
    method = "ds.vertMultinomJointNewton",
    state = "multinomial_design_capsule_unavailable",
    replacement = paste0(
      "For y ~ 1, pass a validated ds.vertDPFrequency object to ",
      "ds.vertMultinom(..., frequency = ...). For additive bare covariates, ",
      "select one same-owner signed finite softmax grid by analysis_id; ",
      "interactions, transforms and sampling inference remain unavailable.")),
  ordinal = list(
    method = "ds.vertOrdinalJointNewton",
    state = "ordinal_score_capsule_unavailable",
    replacement = paste0(
      "For y ~ 1, pass a validated ds.vertDPFrequency object and a complete ",
      "clinical category order to ds.vertOrdinal(..., frequency = ...); ",
      "or select one same-owner signed finite cumulative-logit grid by ",
      "analysis_id for additive bare covariates. Interactions, transforms and ",
      "sampling inference remain unavailable.")),
  gee = list(
    method = "ds.vertGEE",
    state = "cluster_granular_gee_route_quarantined",
    replacement = paste0(
      "For an independence-working binomial/Poisson point estimate, use ",
      "ds.vertGEE(..., corstr='independence', analysis_id=...) for one ",
      "signed finite grid or formal_analysis_id=...) for one completed ",
      "certificate; clustered GEE remains unavailable.")),
  ipw = list(
    method = "ds.vertIPW",
    state = "legacy_ipw_route_quarantined",
    replacement = paste0(
      "Use ds.vertDPCausalStandardization() only when its saturated public-",
      "stratum identification contract matches the scientific question.")),
  mi = list(
    method = "ds.vertMI",
    state = "legacy_mutating_imputation_route_quarantined",
    replacement = paste0(
      "Use ds.vertMI(..., cbind(outcome1, outcome2, ...) ~ 1) for ",
      "signed categorical MCAR completions; dependence='star' is available ",
      "for three or more responses with signed root-to-response pairs, but ",
      "chained imputation remains unavailable.")),
  lasso_iter = list(
    method = "ds.vertLASSOIter",
    state = "binomial_poisson_lasso_path_unavailable",
    replacement = paste0(
      "Only the signed Gaussian Synopsis path is released; binomial and ",
      "Poisson require their purpose-bound score-design artifacts.")),
  legacy_joint_dp_capsule = list(
    method = "legacy joint-DP capsule lifecycle",
    state = "lifetime_admission_route_removed",
    replacement = paste0(
      "Use the per-artifact sticky Synopsis lifecycle; it has no lifetime ",
      "admission gate and replays one authenticated release.")),
  legacy_joint_dp_vector = list(
    method = "legacy joint-DP vector runner",
    state = "lifetime_admission_route_removed",
    replacement = paste0(
      "Use the per-artifact sticky Synopsis lifecycle; it has no lifetime ",
      "admission gate and replays one authenticated release.")),
  formal_finalizer_handoff = list(
    method = "formal finalizer handoff relay",
    state = "unregistered_source_route_removed",
    replacement = "No public handoff relay is released."),
  formal_glm_control = list(
    method = "formal GLM control relay",
    state = "unregistered_source_route_removed",
    replacement = "No public GLM control relay is released.")
)

.dsvert_block_retired_remote_route <- function(route, .allow_test = TRUE) {
  if (isTRUE(.allow_test) && isTRUE(.dsvert_quarantine_test_mode())) {
    return(invisible(FALSE))
  }
  contract <- .DSVERT_RETIRED_REMOTE_ROUTES[[route]]
  if (!is.list(contract) ||
      !identical(sort(names(contract)),
                 c("method", "replacement", "state"))) {
    stop("Invalid retired-route contract", call. = FALSE)
  }
  message <- paste0(
    "[dsvert_route_unavailable:v1] ", contract$method,
    " is unavailable before DSI (state=", contract$state, "). ",
    contract$replacement)
  stop(structure(
    list(
      message = message, call = NULL,
      code = "dsvert_route_unavailable",
      method = contract$method, state = contract$state,
      replacement = contract$replacement),
    class = c("dsvert_route_unavailable", "error", "condition")))
}
