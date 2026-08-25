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

  expect_equal(nrow(inventory), 96L)
  expect_equal(length(non_inference), 13L)
  expect_contains(non_inference, "ds.validateDPGaussianCertificate")
  expect_contains(non_inference, "ds.vertDPCapsulePlan")
  expect_identical(anyDuplicated(inventory$method), 0L)
  expect_length(intersect(inventory$method, non_inference), 0L)
  expect_setequal(c(inventory$method, non_inference), exports)
  expect_setequal(c(inventory$method, non_inference),
                  ds.vertMethodStatus()$method)
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
                   "dsvert-capsule-method-inventory-v4")
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

test_that("sticky Synopsis artifacts never deny a canonical replay", {
  inventory <- .dsvert_capsule_method_inventory()
  synopsis_producers <- c(
    "ds.vertDPContingency", "ds.vertDPMeanVar",
    "ds.vertDPCor", "ds.vertDPDescribe", "ds.vertDPGaussian",
    "ds.vertDPSurvival")
  currently_deniable <- inventory$method[
    inventory$new_capsule_reservation_history_can_deny]

  expect_length(currently_deniable, 0L)
  expect_true(all(inventory$current_route_status[
    inventory$method %in% synopsis_producers] ==
      "formal_sticky_synopsis_artifact"))
  expect_true(all(inventory$migration_feasibility[
    inventory$method %in% synopsis_producers] ==
      "synopsis_release_implemented"))
  expect_false(any(inventory$same_capsule_replay_history_can_deny))
  expect_false(any(inventory$new_capsule_reservation_history_can_deny[
    inventory$method %in% synopsis_producers]))
  count <- row_for(inventory, "ds.vertDPCount")
  expect_identical(count$current_route_status,
                   "formal_sticky_count_artifact")
  expect_identical(count$migration_feasibility,
                   "count_operation_implemented")
  expect_identical(count$artifact_implementation_state,
                   "signed_count_artifact_implemented")
  expect_false(count$new_capsule_reservation_history_can_deny)

  frequency <- row_for(inventory, "ds.vertDPFrequency")
  expect_identical(frequency$current_route_status,
                   "formal_sticky_frequency_artifact")
  expect_identical(frequency$migration_feasibility,
                   "frequency_operation_implemented")
  expect_identical(frequency$artifact_implementation_state,
                   "signed_frequency_artifact_implemented")
  expect_identical(frequency$inference_implementation_state,
                   "formal_frequency_release_implemented")
  expect_false(frequency$new_capsule_reservation_history_can_deny)
  expect_false(any(grepl("capsule", c(
    artifacts_for(inventory, "ds.vertDPFrequency"),
    requirements_for(inventory, "ds.vertDPFrequency")),
    ignore.case = TRUE)))

  frequency_inference <-
    row_for(inventory, "ds.vertDPFrequencyInference")
  expect_identical(frequency_inference$current_route_status,
                   "client_only_inherits_input")
  expect_identical(frequency_inference$migration_feasibility,
                   "frequency_operation_implemented")
  expect_identical(frequency_inference$artifact_implementation_state,
                   "validated_frequency_artifact_adapter_implemented")
  expect_identical(frequency_inference$inference_implementation_state,
                   "frequency_postprocess_implemented")
  expect_false(frequency_inference$new_capsule_reservation_history_can_deny)
})

test_that("aliases and wrappers retain honest routing semantics", {
  expected_aliases <- c(
    "ds.vert.chisq" = "ds.vertChisq",
    "ds.vert.chisq_cross" = "ds.vertChisqCross",
    "ds.vert.confint" = "ds.vertConfint",
    "ds.vert.contrast" = "ds.vertContrast",
    "ds.vert.cor" = "ds.vertCor",
    "ds.vert.cox" = "ds.vertCox",
    "ds.vert.coxph" = "ds.vertCox",
    "ds.vertCoxProfileNonDisclosive" = "ds.vertCox",
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
    "ds.vertLMM" = "ds.vertDPLMM",
    "ds.vert.lmm" = "ds.vertDPLMM",
    "ds.vertLMM.k3" = "ds.vertDPLMM",
    "ds.vert.lr" = "ds.vertLR",
    "ds.vert.mi" = "ds.vertMI",
    "ds.vert.multinom" = "ds.vertMultinom",
    "ds.vertMultinomJoint" = "ds.vertMultinom",
    "ds.vertMultinomJointNewton" = "ds.vertMultinom",
    "ds.vert.nb" = "ds.vertNBFullRegTheta",
    "ds.vert.ordinal" = "ds.vertOrdinal",
    "ds.vertOrdinalJointNewton" = "ds.vertOrdinal",
    "ds.vert.pca" = "ds.vertPCA",
    "ds.vert.wald" = "ds.vertWald")
  wrappers <- c(
    "ds.vert.cox", "ds.vert.coxph", "ds.vertCoxProfileNonDisclosive",
    "ds.vertMultinomJoint", "ds.vertMultinomJointNewton",
    "ds.vertOrdinalJointNewton", "ds.vertLMM.k3")
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
  expect_identical(k3$canonical_method, "ds.vertDPLMM")
  expect_identical(k3$alias_of, "ds.vertDPLMM")
  expect_identical(k3$alias_kind, "compatibility_wrapper")

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
  standardized_families <- "ds.vertGLM"
  expect_true(all(vapply(standardized_families, function(method) {
    has_evidence(inventory, method, "glmStandardizeDS",
                 "plaintext_exact_aggregate")
  }, logical(1L))))

  lmm <- row_for(inventory, "ds.vertDPLMM")
  expect_identical(lmm$current_route_status,
                   "formal_sticky_synopsis_artifact")
  expect_identical(lmm$migration_feasibility,
                   "synopsis_release_implemented")
  expect_identical(nrow(lmm$legacy_remote_call_evidence[[1L]]), 0L)
  mi <- row_for(inventory, "ds.vertMI")
  expect_identical(mi$current_route_status,
                   "client_only_validated_synopsis_postprocess")
  expect_identical(mi$migration_feasibility,
                   "synopsis_release_implemented")
  expect_identical(nrow(mi$legacy_remote_call_evidence[[1L]]), 0L)
  glmm <- row_for(inventory, "ds.vertGLMM")
  expect_identical(glmm$current_route_status,
                   "client_only_validated_synopsis_postprocess")
  expect_identical(glmm$migration_feasibility,
                   "synopsis_release_implemented")
  expect_identical(nrow(glmm$legacy_remote_call_evidence[[1L]]), 0L)
  ipw <- row_for(inventory, "ds.vertIPW")
  expect_identical(ipw$current_route_status,
                   "client_only_validated_synopsis_postprocess")
  expect_identical(ipw$migration_feasibility,
                   "synopsis_release_implemented")
  expect_identical(nrow(ipw$legacy_remote_call_evidence[[1L]]), 0L)
  gee <- row_for(inventory, "ds.vertGEE")
  expect_identical(gee$current_route_status,
                   "formal_completed_public_certificate_only")
  expect_identical(gee$migration_feasibility,
                   "formal_public_certificate_implemented")
  expect_true(all(c(
    "formal_glm_public_certificate", "formal_glm_sticky_opening",
    "formal_glm_two_authority_signatures") %in% gee$artifact_requirements[[1L]]))
  # The retired implementation remains inventoried as historical disclosure
  # evidence, but it is not a current route for this narrow certificate reader.
  expect_gte(length(gee$legacy_remote_call_evidence[[1L]]), 1L)
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
    "gaussian_sufficient_statistics_same_owner",
    "signed_gaussian_model_artifact") %in%
      artifacts_for(inventory, "ds.vertDPGaussian")))
  expect_contains(requirements_for(inventory, "ds.vertDPGaussian"),
                  "no_sampling_inference")

  lasso_paths <- c("ds.vertLASSO", "ds.vertLASSO1Step",
                   "ds.vertLASSOProximal", "ds.vertLASSOCV")
  expect_true(all(vapply(lasso_paths, function(method) {
    all(c("admitted_count", "gaussian_sufficient_statistics_same_owner",
          "signed_gaussian_model_artifact",
          "validated_gaussian_provenance_certificate") %in%
          artifacts_for(inventory, method))
  }, logical(1L))))
  expect_true(all(c(
    "admitted_count", "strict_missing_categorical_marginal",
    "signed_synopsis_release_provenance",
    "deterministic_release_bound_completion_draws") %in%
      artifacts_for(inventory, "ds.vertMI")))
  expect_true(all(c(
    "validated_sticky_frequency_artifact",
    "bounded_nonnegative_integer_domain", "zero_call_postprocessing") %in%
      artifacts_for(inventory, "ds.vertNBFullRegTheta")))
  expect_true(all(c(
    "validated_sticky_categorical_pair_artifact",
    "binary_treatment_outcome_domain",
    "intercept_only_propensity_identity", "zero_call_postprocessing") %in%
      artifacts_for(inventory, "ds.vertIPW")))
  expect_true(all(c(
    "formal_cox_public_certificate", "formal_cox_sticky_opening",
    "formal_cox_two_authority_signatures") %in%
      artifacts_for(inventory, "ds.vertCox")))
  expect_true(all(c(
    "formal_cox_public_certificate", "formal_cox_sticky_opening",
    "formal_cox_two_authority_signatures") %in%
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
    "certificate_revalidation", "stable_full_and_reduced_gram_inverse",
    "explicit_non_sampling_interpretation") %in%
      requirements_for(inventory, "ds.vertLR")))
  expect_true(all(c(
    "canonical_certificate_validation", "completed_release_only",
    "two_authority_signatures") %in%
      requirements_for(inventory, "ds.vertCoxProfileNonDisclosive")))
  expect_true(all(c(
    "intercept_only_response_columns", "strict_missingness_contract",
    "mcar_assumption", "independent_marginals_not_joint_imputation",
    "no_rubin_sampling_inference") %in%
      requirements_for(inventory, "ds.vertMI")))
  expect_true(all(c(
    "intercept_only_propensity_model", "treated_level_binding",
    "mechanism_and_sampling_regions", "no_individual_weights") %in%
      requirements_for(inventory, "ds.vertIPW")))
})

test_that("implemented Synopsis producers and postprocessors are explicit", {
  inventory <- .dsvert_capsule_method_inventory()
  synopsis_producers <- c(
    "ds.vertDPContingency", "ds.vertDPMeanVar",
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
    "ds.vertDPIndirectStandardizationInference",
    "ds.vertIPW", "ds.vert.ipw")

  expect_true(all(inventory$migration_feasibility[
    inventory$method %in% synopsis_producers] ==
      "synopsis_release_implemented"))
  expect_true(all(inventory$current_route_status[
    inventory$method %in% synopsis_producers] ==
      "formal_sticky_synopsis_artifact"))
  expect_false(any(inventory$new_capsule_reservation_history_can_deny[
    inventory$method %in% synopsis_producers]))
  expect_false(any(inventory$same_capsule_replay_history_can_deny[
    inventory$method %in% synopsis_producers]))
  expect_true(all(inventory$artifact_implementation_state[
    inventory$method %in% setdiff(
      synopsis_producers, c("ds.vertDPCor", "ds.vertDPDescribe",
                             "ds.vertDPGaussian"))] ==
      "validated_synopsis_adapter_implemented"))
  expect_identical(
    row_for(inventory, "ds.vertDPCor")$artifact_implementation_state,
    "validated_synopsis_adapter_implemented")
  expect_identical(
    row_for(inventory, "ds.vertDPGaussian")$
      artifact_implementation_state,
    "validated_same_owner_synopsis_adapter_implemented")
  count <- row_for(inventory, "ds.vertDPCount")
  expect_identical(count$current_route_status,
                   "formal_sticky_count_artifact")
  expect_identical(count$inference_implementation_state,
                   "formal_count_release_implemented")
  expect_false(grepl(
    "capsule|lifetime|budget|remaining|limit",
    paste(unlist(count), collapse = " "), ignore.case = TRUE))
  chisq <- inventory$method %in% c("ds.vertChisq", "ds.vert.chisq")
  expect_true(all(inventory$current_route_status[chisq] ==
                    "formal_sticky_synopsis_artifact"))
  expect_true(all(inventory$artifact_implementation_state[chisq] ==
                    "validated_same_owner_synopsis_adapter_implemented"))
  expect_true(all(inventory$inference_implementation_state[chisq] ==
                    "dp_aware_parametric_bootstrap_implemented"))
  fisher <- inventory$method %in% c("ds.vertFisher", "ds.vert.fisher")
  expect_true(all(inventory$artifact_implementation_state[fisher] ==
                    "validated_same_owner_synopsis_adapter_implemented"))
  expect_true(all(inventory$current_route_status[fisher] ==
                    "formal_sticky_synopsis_artifact"))
  expect_true(all(inventory$migration_feasibility[fisher] ==
                    "synopsis_release_implemented"))
  expect_true(all(inventory$inference_implementation_state[fisher] ==
                    "dp_aware_conditional_hypergeometric_bootstrap_implemented"))
  expect_true(all(vapply(inventory$legacy_remote_call_evidence[fisher],
                         nrow, integer(1L)) == 0L))
  lasso_proximal <- inventory$method %in% c(
    "ds.vertLASSOProximal", "ds.vert.lasso_proximal")
  expect_true(all(inventory$current_route_status[lasso_proximal] ==
                    "client_only_validated_synopsis_postprocess"))
  expect_true(all(inventory$migration_feasibility[lasso_proximal] ==
                    "synopsis_release_implemented"))
  expect_true(all(inventory$artifact_implementation_state[lasso_proximal] ==
                    "validated_same_owner_synopsis_adapter_implemented"))
  expect_true(all(inventory$inference_implementation_state[lasso_proximal] ==
                    "dp_gaussian_lasso_path_implemented"))
  expect_true(all(vapply(
    inventory$legacy_remote_call_evidence[lasso_proximal],
    nrow, integer(1L)) == 0L))
  lasso_cv <- inventory$method %in% c(
    "ds.vertLASSOCV", "ds.vert.lasso_cv")
  expect_true(all(inventory$current_route_status[lasso_cv] ==
                    "client_only_validated_synopsis_postprocess"))
  expect_true(all(inventory$migration_feasibility[lasso_cv] ==
                    "synopsis_release_implemented"))
  expect_true(all(inventory$artifact_implementation_state[lasso_cv] ==
                    "validated_same_owner_synopsis_adapter_implemented"))
  expect_true(all(inventory$inference_implementation_state[lasso_cv] ==
                    "dp_gaussian_pseudo_ic_implemented"))
  expect_true(all(vapply(inventory$legacy_remote_call_evidence[lasso_cv],
                         nrow, integer(1L)) == 0L))
  desc <- inventory$method %in% c("ds.vertDesc", "ds.vert.desc")
  expect_true(all(inventory$migration_feasibility[desc] ==
                    "synopsis_release_implemented"))
  expect_true(all(inventory$current_route_status[desc] ==
                    "formal_sticky_synopsis_artifact"))
  expect_true(all(inventory$artifact_implementation_state[desc] ==
                    "validated_synopsis_adapter_implemented"))
  expect_true(all(inventory$inference_implementation_state[desc] ==
                    "synopsis_postprocess_implemented"))
  expect_true(all(vapply(inventory$legacy_remote_call_evidence[desc],
                         nrow, integer(1L)) == 0L))
  expect_true(all(inventory$inference_implementation_state[
    inventory$method %in% setdiff(
      synopsis_producers, c("ds.vertDPCor", "ds.vertDPDescribe",
                             "ds.vertDPGaussian"))] ==
      "synopsis_postprocess_implemented"))
  expect_identical(
    row_for(inventory, "ds.vertDPCor")$inference_implementation_state,
    "pairwise_postprocess_and_explicit_psd_projection_implemented")
  expect_identical(
    row_for(inventory, "ds.vertDPGaussian")$
      inference_implementation_state,
    "synopsis_postprocess_implemented")
  expect_true(all(vapply(inventory$legacy_remote_call_evidence[
    inventory$method %in% synopsis_producers], nrow, integer(1L)) == 0L))

  expect_true(all(inventory$migration_feasibility[
    inventory$method %in% postprocessors] ==
      "synopsis_release_implemented"))
  expect_true(all(inventory$current_route_status[
    inventory$method %in% postprocessors] ==
      "client_only_validated_synopsis_postprocess"))
  expect_true(all(inventory$artifact_implementation_state[
    inventory$method %in% postprocessors] ==
      "validated_synopsis_adapter_implemented"))
  expect_true(all(inventory$inference_implementation_state[
    inventory$method %in% postprocessors] ==
      "synopsis_postprocess_implemented"))
  expect_true(all(vapply(inventory$legacy_remote_call_evidence[
    inventory$method %in% postprocessors], nrow, integer(1L)) == 0L))
  expect_false(any(inventory$same_capsule_replay_history_can_deny))
})

test_that("MI is a scoped Synopsis postprocessor, not the retired mutating route", {
  inventory <- .dsvert_capsule_method_inventory()
  mi <- inventory[inventory$method %in% c("ds.vertMI", "ds.vert.mi"), ]
  expect_true(all(mi$current_route_status ==
                    "client_only_validated_synopsis_postprocess"))
  expect_true(all(mi$migration_feasibility ==
                    "synopsis_release_implemented"))
  expect_true(all(mi$artifact_implementation_state ==
                    "validated_synopsis_adapter_implemented"))
  expect_true(all(vapply(mi$legacy_remote_call_evidence,
                         nrow, integer(1L)) == 0L))
})

test_that("mixed variants and Frequency compatibility names have explicit scopes", {
  inventory <- .dsvert_capsule_method_inventory()
  mixed <- c("ds.vertGLM", "ds.vert.glm")
  frequency_multinom <- c(
    "ds.vertMultinom", "ds.vert.multinom", "ds.vertMultinomJoint",
    "ds.vertMultinomJointNewton")
  frequency_ordinal <- c(
    "ds.vertOrdinal", "ds.vert.ordinal", "ds.vertOrdinalJointNewton")
  frequency_nb2 <- c("ds.vertNBFullRegTheta", "ds.vert.nb")

  expect_true(all(inventory$current_route_status[
    inventory$method %in% mixed] ==
      "formal_completed_public_certificate_or_registered_fresh"))
  expect_true(all(inventory$migration_feasibility[
    inventory$method %in% mixed] ==
      "formal_public_certificate_implemented"))
  expect_true(all(vapply(inventory$artifact_requirements[
    inventory$method %in% mixed], function(requirements) {
      all(c("formal_glm_public_certificate",
            "formal_glm_registered_fresh_terminal",
            "formal_glm_two_authority_signatures") %in% requirements)
    }, logical(1L))))
  expect_true(all(inventory$current_route_status[
    inventory$method %in% c(
      "ds.vertChisqCross", "ds.vert.chisq_cross")] ==
      "formal_sticky_synopsis_artifact"))
  expect_true(all(inventory$current_route_status[
    inventory$method %in% frequency_multinom] ==
      "client_only_validated_capsule_postprocess"))
  expect_true(all(inventory$inference_implementation_state[
    inventory$method %in% frequency_multinom] ==
      "frequency_postprocess_implemented"))
  expect_true(all(inventory$current_route_status[
    inventory$method %in% frequency_ordinal] ==
      "client_only_validated_capsule_postprocess"))
  expect_true(all(inventory$inference_implementation_state[
    inventory$method %in% frequency_ordinal] ==
      "frequency_postprocess_implemented"))
  expect_true(all(inventory$current_route_status[
    inventory$method %in% frequency_nb2] ==
      "client_only_validated_capsule_postprocess"))
  expect_true(all(inventory$inference_implementation_state[
    inventory$method %in% frequency_nb2] ==
      "frequency_postprocess_implemented"))
  lasso_iter <- inventory[inventory$method %in%
    c("ds.vertLASSOIter", "ds.vert.lasso_iter"), , drop = FALSE]
  expect_true(all(lasso_iter$current_route_status ==
                    "formal_same_owner_synopsis_variant_only_legacy_unavailable"))
  expect_true(all(lasso_iter$artifact_implementation_state ==
                    "validated_same_owner_synopsis_adapter_implemented"))
  expect_true(all(lasso_iter$inference_implementation_state ==
                    "dp_gaussian_lasso_path_implemented"))
  expect_true(all(ds.vertMethodStatus(mixed)$status == "promoted"))
  expect_true(all(ds.vertMethodStatus(c(
    "ds.vertLASSOIter", "ds.vert.lasso_iter"))$status == "promoted"))
  expect_true(all(ds.vertMethodStatus(frequency_multinom)$status == "promoted"))
  expect_true(all(ds.vertMethodStatus(frequency_ordinal)$status == "promoted"))
  expect_true(all(ds.vertMethodStatus(frequency_nb2)$status == "promoted"))
})

test_that("NB2 slope and categorical MI compatibility routes retain explicit scopes", {
  inventory <- .dsvert_capsule_method_inventory()
  expect_true(all(inventory$current_route_status[
    inventory$method %in% c("ds.vertNBFullRegTheta", "ds.vert.nb")] ==
      "client_only_validated_capsule_postprocess"))
  expect_true(all(inventory$current_route_status[
    inventory$method %in% c("ds.vertMI", "ds.vert.mi")] ==
      "client_only_validated_synopsis_postprocess"))
  expect_true(all(ds.vertMethodStatus(c(
    "ds.vertNBFullRegTheta", "ds.vert.nb"))$status == "promoted"))
  expect_true(all(ds.vertMethodStatus(c("ds.vertMI", "ds.vert.mi"))$
      status == "promoted"))
})
