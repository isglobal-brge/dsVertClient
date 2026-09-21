# Runtime availability and statistical promotion are tracked separately.
.dsvert_grouped_cross_inventory_register <- function(add) {
  for (family in c("lmm", "glmm", "gee")) {
    method <- paste0("dp_", family, "_grid")
    add(method, method, paste0(family, "_grid_cross_v1"),
        "authenticated_cross_owner_joint_dp_vector",
        "DP selection from a signed finite candidate grid; no standard errors.",
        c("two_owner_signatures", "certified_grouped_arithmetic",
          "authenticated_fusion", "measured_signed_capacity"),
        "synopsis_release_implemented",
        current_route_status = "formal_sticky_synopsis_artifact",
        artifact_implementation_state = "validated_synopsis_adapter_implemented",
        inference_implementation_state = "synopsis_postprocess_implemented")
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
    row$status <- "promoted"
    row$release_contract <- "formal_sticky_synopsis_artifact"
    row$numeric_contract <- "separate_integer_dp_contract"
    row$may_report_numerically_certified <- FALSE
    row$currently_numerically_certified <- FALSE
    row$safe_scope <- if (family == "cox")
      "Signed Breslow finite-grid sticky DP release at up to 400 aligned observations." else if (family == "gee")
      paste("Signed finite-grid sticky DP release for cross-owner binomial/Poisson GEE;",
            "independence (rho = 0) evaluated at 64 clusters of 4 participants.",
            "Server staged_fixed_rho_v1 admission permits signed cluster capacity up to 500",
            "and up to 8 slots per cluster; observation capacity must fit the signed domain.") else
      "Signed finite-grid sticky DP release for cross-owner LMM or binomial/Poisson GLMM."
    row$principal_limitation <- if (family == "gee")
      paste("Grid-resolution-limited selection; no standard errors or sandwich covariance.",
            "Working correlation is fixed and signed rather than estimated; independence evaluated.",
            "No unrestricted optimization. Capacity is descriptive only.") else
      "Grid-resolution-limited selection; no standard errors or unrestricted optimization. Capacity is descriptive only."
    row$numeric_blocker <- "Signed finite-grid arithmetic and DP provenance do not certify unrestricted inference."
    out <- rbind(out, row)
  }
  out
}
