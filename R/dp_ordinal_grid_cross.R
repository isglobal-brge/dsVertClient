# Stable cumulative-logit identity and probability-domain restrictions are
# signed in the categorical family contract. No production release is enabled
# by returning this descriptor to the future fused producer.
.dsvert_dp_ordinal_grid_cross_spec <- function(raw, policy, schema) {
  .dsvert_dp_categorical_grid_cross_spec(raw, policy, schema, "ordinal")
}
.dsvert_dp_ordinal_grid_cross_spec_validate <- function(value, policy, schema) {
  .dsvert_dp_categorical_grid_cross_spec_validate(value, policy, schema, "ordinal")
}
.dsvert_dp_ordinal_grid_cross_contract_validate <- function(
    value, policy, schema_manifest, .verifier = .dsvert_dp_glm_grid_cross_verify) {
  .dsvert_dp_categorical_grid_cross_contract_validate(value, policy,
    schema_manifest, "ordinal", .verifier)
}
.dsvert_dp_ordinal_grid_cross_register <- function() {
  list(family = "ordinal", spec_version = "ordinal_grid_cross_v1",
    operation = "dp.ordinal-grid-cross.v1", release_enabled = FALSE,
    spec = .dsvert_dp_ordinal_grid_cross_spec,
    validate = .dsvert_dp_ordinal_grid_cross_contract_validate,
    artifact = .dsvert_dp_categorical_grid_cross_artifact,
    source_contract = .dsvert_dp_glm_grid_cross_source_contract,
    postprocess = .dsvert_dp_ordinal_grid_cross_postprocess,
    release = .dsvert_dp_ordinal_grid_cross_release)
}

# Ordinal family adapter for the categorical cross-owner contract. The shared
# strict public validator is in dp_multinomial_grid_cross.R.
.dsvert_dp_ordinal_grid_cross_postprocess <- function(contract, noisy_losses) {
  .dsvert_dp_categorical_grid_cross_postprocess(contract, noisy_losses, "ordinal")
}

.dsvert_dp_ordinal_grid_cross_release <- function(contract, datasources) {
  .dsvert_dp_glm_grid_cross_fail()
}

#' Select a signed cross-owner ordinal cumulative-logit grid candidate
#'
#' Validates both custodians' signatures and the complete public numeric and
#' source contract. The release boundary fails closed before any DataSHIELD
#' call until Step 2 supplies its fused producer and authenticated joint-DP
#' release integration.
#' @inheritParams dp_multinomial_grid
#' @return After authenticated release integration, a finite-grid candidate
#'   with coefficients and thresholds and no standard errors. Currently raises
#'   the fixed transcript-safe unavailable error.
#' @details The ordinal domain has 2--8 ordered classes, 1--16 predictors and
#'   at most 256 signed candidates. Intercepts are fixed at zero; candidate
#'   slopes have L1 norm at most 8, thresholds have magnitude at most 8 and
#'   adjacent thresholds differ by at least 1/16. These restrictions bound
#'   each probability by the minimum of sigmoid(-T) and
#'   gap * exp(-T)/(1 + exp(-T))^2, where T is the slope L1 bound plus
#'   the largest threshold magnitude. This is more than 2^-28 for T <= 16;
#'   the two-class case uses the endpoint bound alone. The log loss is finite. Normalized predictor coefficients and thresholds are converted to
#'   original predictor units only after DP selection. No optimizer accesses
#'   protected data, and candidate ties use signed canonical order.
#'   The signed log-sigmoid profile uses 64 quadratic pieces and 32-bit words
#'   with 16 fractional bits. Its uniform loss error is at most 0.00012 per
#'   row, plus output rounding when fewer than 16 output fractional bits are
#'   used. Selection targets this certified profile; an exact-loss total
#'   differs by at most the admitted row count times the certified row error.
#' @note This integration entry point is namespace-internal until the shared
#'   public method registry and authenticated release path are wired together.
#' @keywords internal
dp_ordinal_grid <- function(formula, data, analysis_id, signed_contract,
                            policy, schema_manifest, datasources = NULL) {
  tryCatch({
    contract <- .dsvert_dp_ordinal_grid_cross_contract_validate(
      signed_contract, policy, schema_manifest)
    if (!identical(data, contract$spec$dataset) ||
        !identical(analysis_id, contract$spec$analysis_id)) {
      .dsvert_dp_glm_grid_cross_fail()
    }
    .dsvert_dp_categorical_grid_cross_formula(formula, contract$spec)
    .dsvert_dp_ordinal_grid_cross_release(contract, datasources)
  }, error = .dsvert_dp_glm_grid_cross_transcript_stop)
}
