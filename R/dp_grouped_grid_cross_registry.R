# Runtime availability and statistical promotion are tracked separately.
.dsvert_grouped_cross_inventory_register <- function(add) {
  for (family in c("lmm", "glmm", "gee")) {
    method <- paste0("dp_", family, "_grid")
    add(method, method, paste0(family, "_grid_cross_v1"),
        "authenticated_cross_owner_joint_dp_vector",
        "DP selection from a signed finite candidate grid; no standard errors.",
        c("two_owner_signatures", "certified_grouped_arithmetic",
          "authenticated_fusion", "measured_signed_capacity"),
        "requires_new_secure_protocol",
        current_route_status = "signed_workload_unavailable_quarantine")
  }
  add("dp_cox_grid", "dp_cox_grid", "cox_grid_cross_v1",
      "authenticated_cross_owner_joint_dp_vector",
      "DP selection from a signed Breslow finite coefficient grid; no standard errors.",
      c("all_source_owner_signatures", "certified_cox_arithmetic",
        "authenticated_staged_source", "measured_signed_capacity"),
      "synopsis_release_implemented",
      current_route_status = "formal_sticky_synopsis_artifact",
      artifact_implementation_state = "validated_synopsis_adapter_implemented",
      inference_implementation_state = "synopsis_postprocess_implemented")
  invisible(NULL)
}

.dsvert_grouped_cross_maturity_register <- function(out) {
  for (family in c("lmm", "glmm", "gee", "cox")) {
    row <- out[1L, , drop = FALSE]
    row$method <- row$canonical <- paste0("dp_", family, "_grid")
    row$status <- "quarantine"
    row$release_contract <- "disclosure_safe_protocol_no_statistic"
    row$numeric_contract <- "not_applicable_no_statistic"
    row$may_report_numerically_certified <- FALSE
    row$currently_numerically_certified <- FALSE
    row$safe_scope <- "Signed-contract validation only; protected invocation fails closed."
    row$principal_limitation <- "Authenticated grouped producer and measured capacity are not available."
    row$numeric_blocker <- "Scalar certificates do not certify a complete grouped release."
    if (family == "cox") {
      row$release_contract <- "formal_sticky_synopsis_artifact"
      row$numeric_contract <- "separate_integer_dp_contract"
      row$safe_scope <- "Signed Breslow finite-grid release at N<=400; promotion evidence pending."
      row$principal_limitation <- "Measured fleet capacity and promotion evidence are pending."
      row$numeric_blocker <- "Signed finite-grid arithmetic and DP provenance do not certify unrestricted Cox inference."
    }
    out <- rbind(out, row)
  }
  out
}
