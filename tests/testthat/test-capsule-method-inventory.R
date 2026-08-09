row_for <- function(inventory, method) {
  inventory[inventory$method == method, , drop = FALSE]
}

artifacts_for <- function(inventory, method) {
  row_for(inventory, method)$artifact_requirements[[1L]]
}

requirements_for <- function(inventory, method) {
  row_for(inventory, method)$inference_requirements[[1L]]
}

evidence_for <- function(inventory, method) {
  row_for(inventory, method)$legacy_remote_call_evidence[[1L]]
}

has_evidence <- function(inventory, method, call, release_class = NULL,
                         component = NULL) {
  item <- evidence_for(inventory, method)
  keep <- item$call == call
  if (!is.null(release_class)) {
    keep <- keep & item$release_class == release_class
  }
  if (!is.null(component)) keep <- keep & item$component == component
  any(keep)
}

test_that("capsule migration inventory covers the complete public surface", {
  inventory <- .dsvert_capsule_method_inventory()
  non_inference <- .dsvert_capsule_non_inference_exports()
  exports <- getNamespaceExports("dsVertClient")

  expect_equal(nrow(inventory), 95L)
  expect_equal(length(non_inference), 12L)
  expect_contains(non_inference, "ds.validateDPGaussianCertificate")
  expect_contains(non_inference, "ds.vertDPCapsulePlan")
  expect_identical(anyDuplicated(inventory$method), 0L)
  expect_length(intersect(inventory$method, non_inference), 0L)
  expect_setequal(c(inventory$method, non_inference), exports)
  expect_setequal(c(inventory$method, non_inference),
                  .dsvert_method_registry()$method)
})

test_that("inventory separates current, migration, artifact and inference state", {
  inventory <- .dsvert_capsule_method_inventory()
  expected_columns <- c(
    "method", "canonical_method", "canonical_family", "alias_of",
    "alias_kind", "artifact_requirements", "estimand",
    "inference_requirements", "current_route_status",
    "same_capsule_replay_history_can_deny",
    "new_capsule_reservation_history_can_deny", "migration_feasibility",
    "artifact_implementation_state", "inference_implementation_state",
    "legacy_remote_call_evidence")

  expect_identical(names(inventory), expected_columns)
  expect_identical(attr(inventory, "schema_version"),
                   "dsvert-capsule-method-inventory-v3")
  expect_true(all(nzchar(inventory$canonical_method)))
  expect_true(all(nzchar(inventory$canonical_family)))
  expect_true(all(nzchar(inventory$estimand)))
  expect_true(all(inventory$canonical_method %in%
                    getNamespaceExports("dsVertClient")))
  expect_true(all(inventory$alias_kind %in%
                    attr(inventory, "alias_kind_levels")))
  expect_true(all(inventory$current_route_status %in%
                    attr(inventory, "current_route_status_levels")))
  expect_true(all(inventory$migration_feasibility %in%
                    attr(inventory, "migration_feasibility_levels")))
  expect_true(all(inventory$artifact_implementation_state %in%
                    attr(inventory,
                         "artifact_implementation_state_levels")))
  expect_true(all(inventory$inference_implementation_state %in%
                    attr(inventory,
                         "inference_implementation_state_levels")))
  expect_true(all(lengths(inventory$artifact_requirements) > 0L))
  expect_true(all(lengths(inventory$inference_requirements) > 0L))
  expect_true(all(vapply(inventory$artifact_requirements, is.character,
                         logical(1L))))
  expect_true(all(vapply(inventory$inference_requirements, is.character,
                         logical(1L))))
  expect_false(any(inventory$same_capsule_replay_history_can_deny))
})

test_that("only new capsule reservation can be denied by lifetime history", {
  inventory <- .dsvert_capsule_method_inventory()
  capsule_producers <- c(
    "ds.vertDPCount", "ds.vertDPContingency", "ds.vertDPMeanVar",
    "ds.vertDPCor", "ds.vertDPDescribe", "ds.vertDPGaussian",
    "ds.vertDPSurvival")
  currently_deniable <- inventory$method[
    inventory$new_capsule_reservation_history_can_deny]

  expect_setequal(
    currently_deniable,
    inventory$method[
      inventory$current_route_status == "formal_joint_dp_capsule"])
  expect_true(all(inventory$current_route_status[
    inventory$method %in% capsule_producers] ==
      "formal_joint_dp_capsule"))
  expect_true(all(inventory$migration_feasibility[
    inventory$method %in% capsule_producers] ==
      "capsule_release_implemented"))
  expect_false(any(inventory$same_capsule_replay_history_can_deny))
  expect_true(all(inventory$new_capsule_reservation_history_can_deny[
    inventory$method %in% capsule_producers]))
})

test_that("aliases and wrappers retain honest routing semantics", {
  expected_aliases <- c(
    "ds.vert.chisq" = "ds.vertChisq",
    "ds.vert.chisq_cross" = "ds.vertChisqCross",
    "ds.vert.confint" = "ds.vertConfint",
    "ds.vert.contrast" = "ds.vertContrast",
    "ds.vert.cor" = "ds.vertCor",
    "ds.vert.coxph" = "ds.vertCoxProfileNonDisclosive",
    "ds.vert.desc" = "ds.vertDesc",
    "ds.vert.fisher" = "ds.vertFisher",
    "ds.vert.gee" = "ds.vertGEE",
    "ds.vert.glm" = "ds.vertGLM",
    "ds.vert.glmm" = "ds.vertGLMM",
    "ds.vert.ipw" = "ds.vertIPW",
    "ds.vert.lasso" = "ds.vertLASSO",
    "ds.vert.lasso_1step" = "ds.vertLASSO1Step",
    "ds.vert.lasso_cv" = "ds.vertLASSOCV",
    "ds.vert.lasso_iter" = "ds.vertLASSOIter",
    "ds.vert.lasso_proximal" = "ds.vertLASSOProximal",
    "ds.vert.lmm" = "ds.vertLMM",
    "ds.vert.lr" = "ds.vertLR",
    "ds.vert.mi" = "ds.vertMI",
    "ds.vert.multinom" = "ds.vertMultinomJointNewton",
    "ds.vert.nb" = "ds.vertNBFullRegTheta",
    "ds.vert.ordinal" = "ds.vertOrdinalJointNewton",
    "ds.vert.pca" = "ds.vertPCA",
    "ds.vert.wald" = "ds.vertWald",
    "ds.vertCox" = "ds.vertCoxProfileNonDisclosive",
    "ds.vertMultinom" = "ds.vertMultinomJointNewton",
    "ds.vertMultinomJoint" = "ds.vertMultinomJointNewton",
    "ds.vertOrdinal" = "ds.vertOrdinalJointNewton")
  wrappers <- c(
    "ds.vertCox", "ds.vert.coxph", "ds.vertMultinom", "ds.vertMultinomJoint",
    "ds.vertOrdinal")
  inventory <- .dsvert_capsule_method_inventory()
  actual <- inventory$alias_of
  names(actual) <- inventory$method
  actual <- actual[!is.na(actual)]

  expect_identical(actual[sort(names(actual), method = "radix")],
                   expected_aliases[sort(names(expected_aliases),
                                         method = "radix")])
  expect_true(all(inventory$alias_kind[inventory$method %in% wrappers] ==
                    "compatibility_wrapper"))
  expect_true(all(inventory$alias_kind[
    !is.na(inventory$alias_of) & !inventory$method %in% wrappers] ==
      "compatibility_alias"))

  k3 <- row_for(inventory, "ds.vertLMM.k3")
  expect_identical(k3$canonical_method, "ds.vertLMM")
  expect_true(is.na(k3$alias_of))
  expect_identical(k3$alias_kind, "deprecated_subroute")
  expect_contains(k3$inference_requirements[[1L]], "k_ge_3_only")

  contract_columns <- c(
    "canonical_method", "canonical_family", "artifact_requirements",
    "estimand", "inference_requirements", "current_route_status",
    "same_capsule_replay_history_can_deny",
    "new_capsule_reservation_history_can_deny", "migration_feasibility",
    "artifact_implementation_state", "inference_implementation_state",
    "legacy_remote_call_evidence")
  for (alias in names(actual)) {
    alias_row <- row_for(inventory, alias)
    target_row <- row_for(inventory, actual[[alias]])
    rownames(alias_row) <- NULL
    rownames(target_row) <- NULL
    expect_identical(alias_row[contract_columns], target_row[contract_columns],
                     info = paste("alias contract differs for", alias))
  }
})

test_that("remote-call evidence is classified and scoped", {
  inventory <- .dsvert_capsule_method_inventory()
  allowed_classes <- attr(inventory, "remote_release_class_levels")
  remote_counts <- integer(nrow(inventory))

  for (i in seq_len(nrow(inventory))) {
    item <- inventory$legacy_remote_call_evidence[[i]]
    expect_s3_class(item, "data.frame")
    expect_identical(names(item),
                     c("call", "component", "release_class", "evidence"))
    expect_identical(anyDuplicated(paste(item$call, item$component,
                                         sep = "::")), 0L)
    expect_true(all(nzchar(item$call)))
    expect_true(all(nzchar(item$component)))
    expect_true(all(nzchar(item$evidence)))
    expect_true(all(item$release_class %in% allowed_classes))
    remote_counts[[i]] <- nrow(item)
  }

  remote_migrations <- inventory$migration_feasibility %in% c(
    "requires_new_capsule_artifact", "requires_new_secure_protocol")
  local_migrations <- inventory$migration_feasibility ==
    "client_only_requires_attested_input"
  expect_true(all(remote_counts[remote_migrations] > 0L))
  expect_true(all(remote_counts[local_migrations] == 0L))
})

test_that("retired legacy DP endpoints are absent from public route evidence", {
  inventory <- .dsvert_capsule_method_inventory()
  retired <- c(
    "dsvertDPCountDS", "dsvertDPContingencyDS", "dsvertDPMeanVarDS",
    "dsvertDPDescribeDS", "dsvertDPSurvivalDS")
  evidence_calls <- unique(unlist(lapply(
    inventory$legacy_remote_call_evidence, `[[`, "call"),
    use.names = FALSE))

  expect_length(intersect(evidence_calls, retired), 0L)
  expect_false("legacy_dp_release" %in%
                 attr(inventory, "remote_release_class_levels"))
})

test_that("verified legacy disclosure evidence cannot regress to omission", {
  inventory <- .dsvert_capsule_method_inventory()
  standardized_families <- c(
    "ds.vertGLM", "ds.vertNBFullRegTheta",
    "ds.vertMultinomJointNewton", "ds.vertOrdinalJointNewton",
    "ds.vertLMM", "ds.vertGEE", "ds.vertGLMM", "ds.vertIPW",
    "ds.vertMI", "ds.vertLASSOIter")
  expect_true(all(vapply(standardized_families, function(method) {
    has_evidence(inventory, method, "glmStandardizeDS",
                 "plaintext_exact_aggregate")
  }, logical(1L))))

  expect_true(has_evidence(
    inventory, "ds.vertNBFullRegTheta", "dsvertLocalMomentsDS",
    "plaintext_exact_aggregate"))
  expect_true(has_evidence(
    inventory, "ds.vertNBFullRegTheta", "dsvertNBProfileSumsDS",
    "plaintext_exact_aggregate"))
  expect_true(has_evidence(
    inventory, "ds.vertNBFullRegTheta", "dsvertNBPsiAggregateDS",
    "plaintext_exact_aggregate"))
  expect_gte(length(unique(evidence_for(
    inventory, "ds.vertNBFullRegTheta")$call)), 20L)

  expect_true(has_evidence(
    inventory, "ds.vertLMM", "dsvertClusterZtZDS",
    "plaintext_exact_aggregate", "per_cluster_ztz_sizes_and_levels"))
  expect_true(has_evidence(
    inventory, "ds.vertMI", "dsvertImputeColumnDS",
    "plaintext_exact_aggregate", "imputation_counts"))
  expect_true(has_evidence(
    inventory, "ds.vertIPW", "k2ShareWeightsDS",
    "opaque_peer_ciphertext", "weight_share_blobs"))
  expect_true(has_evidence(
    inventory, "ds.vertGEE", "dsvertGEEAR1OrderBroadcastDS",
    "opaque_peer_ciphertext", "peer_blob"))
  expect_true(has_evidence(
    inventory, "ds.vertGLM", "k2GradientR1DS",
    "opaque_peer_ciphertext", "encrypted_peer_round"))
  expect_true(has_evidence(
    inventory, "ds.vertGLM", "k2GradientR2DS",
    "share_reconstructed_by_client", "gradient_and_residual_shares"))
})

test_that("planned artifacts describe the actual biomedical contracts", {
  inventory <- .dsvert_capsule_method_inventory()

  expect_contains(artifacts_for(inventory, "ds.vertDesc"),
                  "fixed_numeric_histograms")
  expect_false("categorical_marginals" %in%
                 artifacts_for(inventory, "ds.vertDesc"))
  expect_true(all(c(
    "gaussian_models", "signed_complete_case_gaussian_artifact") %in%
      artifacts_for(inventory, "ds.vertCor")))
  expect_false("numeric_pair_moments_same_owner" %in%
                 artifacts_for(inventory, "ds.vertCor"))
  expect_false("numeric_cross_products_same_and_cross_owner" %in%
                 artifacts_for(inventory, "ds.vertCor"))
  expect_contains(artifacts_for(inventory, "ds.vertChisqCross"),
                  "categorical_pairs_cross_owner")
  expect_true(all(c(
    "complete_case_patient_collapse",
    "gaussian_sufficient_statistics_same_and_cross_owner",
    "signed_gaussian_model_artifact") %in%
      artifacts_for(inventory, "ds.vertDPGaussian")))
  expect_contains(requirements_for(inventory, "ds.vertDPGaussian"),
                  "no_sampling_inference")

  expect_setequal(artifacts_for(inventory, "ds.vertLASSO1Step"), c(
    "authorized_fit_coefficients", "authorized_fit_covariance_or_fisher"))
  expect_true(all(c(
    "admitted_count", "gaussian_sufficient_statistics_same_and_cross_owner",
    "signed_gaussian_model_artifact",
    "validated_gaussian_provenance_certificate",
    "legacy_authorized_unpenalized_gaussian_fit") %in%
      artifacts_for(inventory, "ds.vertLASSOProximal")))
  expect_true(all(c(
    "admitted_count", "gaussian_sufficient_statistics_same_and_cross_owner",
    "signed_gaussian_model_artifact",
    "validated_gaussian_provenance_certificate",
    "legacy_authorized_fit_covariance_or_fisher") %in%
      artifacts_for(inventory, "ds.vertLASSOCV")))

  expect_true(all(c(
    "bounded_missingness_counts", "posterior_parameter_draws",
    "synthetic_imputation_draws") %in% artifacts_for(inventory, "ds.vertMI")))
  expect_true(all(c(
    "nb2_beta_score_information", "nb2_theta_score_information",
    "nb2_beta_theta_cross_information") %in%
      artifacts_for(inventory, "ds.vertNBFullRegTheta")))
  expect_true(all(c(
    "propensity_score_information", "outcome_score_information",
    "treatment_outcome_binding", "bounded_weight_distribution") %in%
      artifacts_for(inventory, "ds.vertIPW")))
  expect_true(all(c(
    "cox_partial_likelihood", "cox_score_hessian", "cox_baseline_hazard") %in%
      artifacts_for(inventory, "ds.vertCoxProfileNonDisclosive")))
  expect_contains(artifacts_for(inventory, "ds.vertDPDescribe"),
                  "fixed_numeric_histograms")
})

test_that("estimands and inference requirements match implemented semantics", {
  inventory <- .dsvert_capsule_method_inventory()

  expect_match(row_for(inventory, "ds.vertEpi2x2")$estimand,
               "among the exposed")
  expect_match(row_for(inventory, "ds.vertEpi2x2")$estimand,
               "number needed")
  expect_match(row_for(inventory, "ds.vertDPEpi2x2")$estimand,
               "attributable")
  expect_match(row_for(inventory, "ds.vertDPEpi2x2")$estimand,
               "number needed")
  expect_match(
    row_for(inventory, "ds.vertDPEpi2x2Inference")$estimand,
    "sampling uncertainty")
  expect_true(all(c(
    "binomial_sampling_model", "clopper_pearson_exact_intervals",
    "joint_mechanism_and_sampling_uncertainty") %in%
      requirements_for(inventory, "ds.vertDPEpi2x2Inference")))
  expect_match(
    row_for(inventory, "ds.vertDPPrevalenceRatio")$estimand,
    "exact aliases")
  expect_contains(
    requirements_for(inventory, "ds.vertDPPrevalenceRatio"),
    "cross_sectional_design_declared_not_inferred")
  expect_match(
    row_for(inventory, "ds.vertDPPrevalenceRatioInference")$estimand,
    "same conservative joint")
  expect_true(all(c(
    "cross_sectional_design_declared_not_inferred",
    "explicit_exposed_and_prevalent_orientation",
    "zero_call_numeric_identity") %in%
      requirements_for(
        inventory, "ds.vertDPPrevalenceRatioInference")))
  expect_match(row_for(inventory, "ds.vertMantelHaenszel")$estimand,
               "Common Mantel-Haenszel")
  expect_match(row_for(inventory, "ds.vertDPMantelHaenszel")$estimand,
               "no classical CMH p-value")
  expect_match(row_for(inventory, "ds.vertDPROC")$estimand,
               "AUC")
  expect_match(
    row_for(inventory, "ds.vertDPDiagnostic2x2Inference")$estimand,
    "sampling uncertainty")
  expect_true(all(c(
    "binomial_sampling_model", "clopper_pearson_exact_intervals",
    "joint_mechanism_and_sampling_uncertainty") %in%
      requirements_for(inventory, "ds.vertDPDiagnostic2x2Inference")))
  expect_match(row_for(inventory, "ds.vertDirectStandardization")$estimand,
               "cases and person-time")
  expect_match(
    row_for(inventory, "ds.vertDPDirectStandardizationInference")$estimand,
    "sampling uncertainty")
  expect_true(all(c(
    "clopper_pearson_exact_intervals", "fixed_public_standard_weights",
    "joint_mechanism_and_sampling_uncertainty") %in%
      requirements_for(
        inventory, "ds.vertDPDirectStandardizationInference")))

  expect_true(all(c(
    "binomial_or_poisson_only", "canonical_unweighted_deviance",
    "converged_unpenalized_fits", "same_cohort_missingness_and_offset") %in%
      requirements_for(inventory, "ds.vertLR")))
  expect_true(all(c(
    "ties_method", "strata_contract", "delayed_entry_contract") %in%
      requirements_for(inventory, "ds.vertCoxProfileNonDisclosive")))
  expect_true(all(c(
    "cryptographic_non_rerollable_draw_stream", "rubin_small_sample_df",
    "posterior_parameter_uncertainty") %in%
      requirements_for(inventory, "ds.vertMI")))
  expect_contains(requirements_for(inventory, "ds.vertIPW"),
                  "weight_provenance_binding")
})

test_that("implemented capsule producers and postprocessors are explicit", {
  inventory <- .dsvert_capsule_method_inventory()
  capsule_producers <- c(
    "ds.vertDPCount", "ds.vertDPContingency", "ds.vertDPMeanVar",
    "ds.vertDPCor", "ds.vertDPDescribe", "ds.vertDPGaussian",
    "ds.vertDPSurvival")
  postprocessors <- c(
    "ds.vertDPKaplanMeier", "ds.vertDPNelsonAalen",
    "ds.vertDPCumulativeIncidence", "ds.vertDPRMST", "ds.vertDPRMTL",
    "ds.vertDPSurvivalContrast", "ds.vertDPRMSTContrast",
    "ds.vertDPSurvivalQuantile", "ds.vertDPMedianSurvival",
    "ds.vertDPEpi2x2", "ds.vertDPEpi2x2Inference",
    "ds.vertDPPrevalenceRatio", "ds.vertDPPrevalenceRatioInference",
    "ds.vertDPMantelHaenszel", "ds.vertDPDiagnostic2x2",
    "ds.vertDPDiagnostic2x2Inference", "ds.vertDPROC",
    "ds.vertDPDirectStandardization",
    "ds.vertDPDirectStandardizationInference",
    "ds.vertDPIndirectStandardization",
    "ds.vertDPIndirectStandardizationInference")

  expect_true(all(inventory$migration_feasibility[
    inventory$method %in% capsule_producers] ==
      "capsule_release_implemented"))
  expect_true(all(inventory$current_route_status[
    inventory$method %in% capsule_producers] ==
      "formal_joint_dp_capsule"))
  expect_true(all(inventory$new_capsule_reservation_history_can_deny[
    inventory$method %in% capsule_producers]))
  expect_false(any(inventory$same_capsule_replay_history_can_deny[
    inventory$method %in% capsule_producers]))
  expect_true(all(inventory$artifact_implementation_state[
    inventory$method %in% setdiff(
      capsule_producers, c("ds.vertDPCor", "ds.vertDPGaussian"))] ==
      "joint_vector_release_implemented"))
  expect_identical(
    row_for(inventory, "ds.vertDPCor")$artifact_implementation_state,
    "validated_same_owner_capsule_adapter_implemented")
  expect_identical(
    row_for(inventory, "ds.vertDPGaussian")$
      artifact_implementation_state,
    "validated_same_and_cross_owner_capsule_adapter_implemented")
  chisq <- inventory$method %in% c("ds.vertChisq", "ds.vert.chisq")
  expect_true(all(inventory$current_route_status[chisq] ==
                    "formal_joint_dp_capsule"))
  expect_true(all(inventory$artifact_implementation_state[chisq] ==
                    "validated_capsule_adapter_implemented"))
  expect_true(all(inventory$inference_implementation_state[chisq] ==
                    "dp_aware_parametric_bootstrap_implemented"))
  fisher <- inventory$method %in% c("ds.vertFisher", "ds.vert.fisher")
  expect_true(all(inventory$artifact_implementation_state[fisher] ==
                    "validated_capsule_adapter_implemented"))
  expect_true(all(inventory$current_route_status[fisher] ==
                    "formal_joint_dp_capsule"))
  expect_true(all(inventory$migration_feasibility[fisher] ==
                    "capsule_release_implemented"))
  expect_true(all(inventory$inference_implementation_state[fisher] ==
                    "dp_aware_conditional_hypergeometric_bootstrap_implemented"))
  expect_true(all(vapply(inventory$legacy_remote_call_evidence[fisher],
                         nrow, integer(1L)) == 0L))
  lasso_proximal <- inventory$method %in% c(
    "ds.vertLASSOProximal", "ds.vert.lasso_proximal")
  expect_true(all(inventory$current_route_status[lasso_proximal] ==
                    "client_only_inherits_input"))
  expect_true(all(inventory$migration_feasibility[lasso_proximal] ==
                    "capsule_release_implemented"))
  expect_true(all(inventory$artifact_implementation_state[lasso_proximal] ==
                    "validated_same_and_cross_owner_capsule_adapter_implemented"))
  expect_true(all(inventory$inference_implementation_state[lasso_proximal] ==
                    "dp_gaussian_lasso_with_legacy_compatibility_implemented"))
  expect_true(all(vapply(
    inventory$legacy_remote_call_evidence[lasso_proximal],
    nrow, integer(1L)) == 0L))
  lasso_cv <- inventory$method %in% c(
    "ds.vertLASSOCV", "ds.vert.lasso_cv")
  expect_true(all(inventory$current_route_status[lasso_cv] ==
                    "client_only_inherits_input"))
  expect_true(all(inventory$migration_feasibility[lasso_cv] ==
                    "capsule_release_implemented"))
  expect_true(all(inventory$artifact_implementation_state[lasso_cv] ==
                    "validated_same_and_cross_owner_capsule_adapter_implemented"))
  expect_true(all(inventory$inference_implementation_state[lasso_cv] ==
                    "dp_gaussian_pseudo_ic_with_legacy_compatibility_implemented"))
  expect_true(all(vapply(inventory$legacy_remote_call_evidence[lasso_cv],
                         nrow, integer(1L)) == 0L))
  desc <- inventory$method %in% c("ds.vertDesc", "ds.vert.desc")
  expect_true(all(inventory$migration_feasibility[desc] ==
                    "capsule_release_implemented"))
  expect_true(all(inventory$current_route_status[desc] ==
                    "formal_joint_dp_capsule"))
  expect_true(all(inventory$artifact_implementation_state[desc] ==
                    "validated_capsule_adapter_implemented"))
  expect_true(all(inventory$inference_implementation_state[desc] ==
                    "capsule_postprocess_implemented"))
  expect_true(all(vapply(inventory$legacy_remote_call_evidence[desc],
                         nrow, integer(1L)) == 0L))
  expect_true(all(inventory$inference_implementation_state[
    inventory$method %in% setdiff(
      capsule_producers, c("ds.vertDPCor", "ds.vertDPGaussian"))] ==
      "formal_capsule_release_implemented"))
  expect_identical(
    row_for(inventory, "ds.vertDPCor")$inference_implementation_state,
    "pairwise_postprocess_and_explicit_psd_projection_implemented")
  expect_identical(
    row_for(inventory, "ds.vertDPGaussian")$
      inference_implementation_state,
    "capsule_postprocess_implemented")
  expect_true(all(vapply(inventory$legacy_remote_call_evidence[
    inventory$method %in% capsule_producers], nrow, integer(1L)) == 0L))

  expect_true(all(inventory$migration_feasibility[
    inventory$method %in% postprocessors] ==
      "capsule_release_implemented"))
  expect_true(all(inventory$current_route_status[
    inventory$method %in% postprocessors] ==
      "client_only_validated_capsule_postprocess"))
  expect_true(all(inventory$artifact_implementation_state[
    inventory$method %in% postprocessors] ==
      "validated_capsule_adapter_implemented"))
  expect_true(all(inventory$inference_implementation_state[
    inventory$method %in% postprocessors] ==
      "capsule_postprocess_implemented"))
  expect_true(all(vapply(inventory$legacy_remote_call_evidence[
    inventory$method %in% postprocessors], nrow, integer(1L)) == 0L))
  expect_false(any(inventory$same_capsule_replay_history_can_deny))
})

test_that("quarantine labels require secure redesign and concrete evidence", {
  inventory <- .dsvert_capsule_method_inventory()
  quarantined <- inventory[
    inventory$current_route_status == "legacy_granular_release_quarantine",
    , drop = FALSE]

  expect_gt(nrow(quarantined), 0L)
  expect_true(all(quarantined$migration_feasibility ==
                    "requires_new_secure_protocol"))
  expect_true(all(quarantined$artifact_implementation_state ==
                    "secure_artifact_not_implemented"))
  expect_true(all(lengths(quarantined$legacy_remote_call_evidence) > 0L))
  expect_true(all(quarantined$canonical_family %in% c(
    "linear_mixed_model", "gee", "glmm", "ipw", "negative_binomial")))
  expect_true(all(quarantined$method %in%
                    getNamespaceExports("dsVertClient")))
})

test_that("mixed variants and unavailable signed workloads cannot look promoted", {
  inventory <- .dsvert_capsule_method_inventory()
  mixed <- c("ds.vertGLM", "ds.vert.glm")
  broken <- c(
    "ds.vertMultinom", "ds.vertMultinomJoint",
    "ds.vertMultinomJointNewton", "ds.vert.multinom",
    "ds.vertLASSOIter", "ds.vert.lasso_iter")

  expect_true(all(inventory$current_route_status[
    inventory$method %in% mixed] ==
      "formal_capsule_variant_only_legacy_unavailable"))
  expect_true(all(inventory$current_route_status[
    inventory$method %in% c(
      "ds.vertChisqCross", "ds.vert.chisq_cross")] ==
      "formal_joint_dp_capsule"))
  expect_true(all(inventory$current_route_status[
    inventory$method %in% broken] ==
      "signed_workload_unavailable_quarantine"))
  expect_true(all(inventory$inference_implementation_state[
    inventory$method %in% c("ds.vertLASSOIter", "ds.vert.lasso_iter")] ==
      "signed_binomial_lasso_design_gram_not_materialized"))
  expect_true(all(inventory$inference_implementation_state[
    inventory$method %in% setdiff(broken, c(
      "ds.vertLASSOIter", "ds.vert.lasso_iter"))] ==
      "signed_multinomial_design_gram_not_materialized"))
  expect_length(intersect(
    ds.vertMethodStatus(status = "promoted")$method,
    c(mixed, broken)), 0L)
})

test_that("NB and mutating MI routes are explicitly quarantined", {
  inventory <- .dsvert_capsule_method_inventory()
  expect_true(all(inventory$current_route_status[
    inventory$method %in% c("ds.vertNBFullRegTheta", "ds.vert.nb")] ==
      "legacy_granular_release_quarantine"))
  expect_true(all(inventory$current_route_status[
    inventory$method %in% c("ds.vertMI", "ds.vert.mi")] ==
      "legacy_mutating_release_quarantine"))
  expect_true(all(ds.vertMethodStatus(c(
    "ds.vertNBFullRegTheta", "ds.vert.nb", "ds.vertMI", "ds.vert.mi"))$
      status == "quarantine"))
})
