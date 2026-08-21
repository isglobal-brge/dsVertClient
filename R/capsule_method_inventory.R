# Internal migration inventory for the disclosure-safe capsule surface.
#
# This file is deliberately independent of ds.vertMethodStatus().  The latter
# is the public maturity registry; this inventory is a migration contract that
# must be able to say that a public name still needs work without changing the
# advertised product status prematurely.

.dsvert_capsule_non_inference_exports <- function() {
  sort(c(
    "ds.getIdentityPks",
    "ds.isPsiAligned",
    "ds.psiAlign",
    "ds.vert.align",
    "ds.vert.is_aligned",
    "ds.vertDPCalibrate",
    "ds.vertDPCapsulePlan",
    "ds.vertDPStatus",
    "ds.validateDPGaussianCertificate",
    "ds.vertMethodStatus",
    "ds.vertNumericPreflight",
    "ds.vertSecurityStatus"
  ), method = "radix")
}

.dsvert_capsule_method_inventory <- function() {
  rows <- list()

  evidence <- function(call, component, release_class, source) {
    call <- as.character(call)
    component <- as.character(component)
    release_class <- as.character(release_class)
    source <- as.character(source)
    n <- max(length(call), length(component), length(release_class),
             length(source))
    data.frame(
      call = rep_len(call, n),
      component = rep_len(component, n),
      release_class = rep_len(release_class, n),
      evidence = rep_len(source, n),
      stringsAsFactors = FALSE)
  }

  remote_catalog <- do.call(rbind, list(
    evidence("dsvertColNamesDS", "column_names", "public_metadata",
             paste(
               "dsVert/R/mpcUtils.R returns only custodian-published policy",
               "schema without resolving the protected object.")),
    evidence("getObsCountDS", "observation_counts",
             "plaintext_exact_aggregate",
             "dsVert/R/getObsCountDS.R returns exact complete-case counts."),
    evidence("glmStandardizeDS", "means_sds_and_response_moments",
             "plaintext_exact_aggregate",
             paste(
               "dsVert/R/glmStandardizeDS.R returns exact predictor means/SDs",
               "and, when requested, response mean/SD.")),
    evidence("dsvertLocalMomentsDS", "local_numeric_moments",
             "plaintext_exact_aggregate",
             "dsVert/R/localMomentsDS.R returns exact local moments."),
    evidence("dsvertOneHotDS", "levels_and_marginal_counts",
             "plaintext_exact_aggregate",
             paste(
               "dsVert/R/oneHotDS.R exposes the fixed domain and exact local",
               "marginal counts used by the cross-site route.")),
    evidence("k2ChisqCrossAccumulateCountDS", "joint_count_share",
             "share_reconstructed_by_client",
             paste(
               "dsVert/R/exactGCChisqGuardDS.R returns additive count shares;",
               "the current client route combines the peer outputs.")),
    evidence("exactGCChisqJointReleaseDS", "guarded_joint_table",
             "client_reconstructed_exact_statistic",
             paste(
               "dsVert/R/exactGCChisqGuardDS.R releases the exact guarded",
               "joint table after the secure protocol.")),
    evidence("k2GradientR1DS", "encrypted_peer_round",
             "opaque_peer_ciphertext",
             paste(
               "dsVert/R/k2InputSharingDS.R returns encrypted_r1 sealed to",
               "the pinned peer.")),
    evidence("k2GradientR1DS", "residual_sum_share",
             "share_reconstructed_by_client",
             paste(
               "dsVert/R/k2InputSharingDS.R also returns a residual-sum share",
               "that the current client combines across parties.")),
    evidence("k2GradientR2DS", "gradient_and_residual_shares",
             "share_reconstructed_by_client",
             paste(
               "dsVert/R/k2InputSharingDS.R returns gradient/residual shares",
               "that the current client combines across parties.")),
    evidence("glmRing63DevianceSumsDS", "deviance_sum_share",
             "share_reconstructed_by_client",
             paste(
               "dsVert/R/glmRing63DS.R returns a fixed-point deviance-sum",
               "share consumed by the client route.")),
    evidence("glmRing63DevianceSumsDS", "poisson_null_term",
             "plaintext_exact_aggregate",
             "dsVert/R/glmRing63DS.R returns the exact Poisson null term."),
    evidence("dsvertCoxEventTimeShareMaskDS", "risk_and_event_shares",
             "opaque_peer_ciphertext",
             paste(
               "dsVert/R/coxDiscreteShareDS.R seals patient-by-time mask and",
               "event shares to the pinned peer.")),
    evidence("dsvertCoxEventTimeShareMaskDS", "guarded_dimensions",
             "public_metadata",
             paste(
               "dsVert/R/coxDiscreteShareDS.R returns only guarded scalar",
               "dimensions alongside the sealed blobs.")),
    evidence("dsvertCoxDiscreteExpandXDS", "expanded_frame_dimensions",
             "public_metadata",
             paste(
               "dsVert/R/coxDiscreteShareDS.R mutates server state and returns",
               "stored/n/p metadata, not the expanded rows.")),
    evidence("k2BeaverSumShareDS", "sum_share",
             "share_reconstructed_by_client",
             paste(
               "dsVert/R/coxDiscreteShareDS.R returns an additive sum share",
               "that the current client combines across parties.")),
    evidence("dsvertNBProfileSumsDS", "profile_sums_n_and_mean",
             "plaintext_exact_aggregate",
             paste(
               "dsVert/R/nbProfileDS.R returns exact digamma/trigamma sums,",
               "n and outcome mean for each queried theta.")),
    evidence("dsvertNBPsiAggregateDS", "psi_sums_and_n",
             "plaintext_exact_aggregate",
             paste(
               "dsVert/R/nbFullRegShareDS.R returns exact digamma/trigamma",
               "sums and n for each queried theta.")),
    evidence("dsvertGEEInterceptShareDS", "intercept_share",
             "share_reconstructed_by_client",
             "dsVert/R/nbFullRegShareDS.R returns a ring share of the intercept term."),
    evidence("dsvertNBEtaShareDS", "eta_share_blob",
             "opaque_peer_ciphertext",
             "dsVert/R/nbFullRegShareDS.R seals an eta share to a pinned peer."),
    evidence("dsvertNBEtaShareConfirmDS", "receipt",
             "public_metadata",
             "dsVert/R/nbFullRegShareDS.R returns only a share receipt."),
    evidence("dsvertNBEtaTotalReceiveDS", "receipt",
             "public_metadata",
             "dsVert/R/nbFullRegShareDS.R stores peer material and returns a receipt."),
    evidence("dsvertNBYThetaShareDS", "y_theta_share_blob",
             "opaque_peer_ciphertext",
             "dsVert/R/nbFullRegShareDS.R seals y/theta share material to a pinned peer."),
    evidence("dsvertNBYThetaShareReceiveDS", "receipt",
             "public_metadata",
             "dsVert/R/nbFullRegShareDS.R stores peer material and returns a receipt."),
    evidence("glmRing63ExportOwnShareDS", "x_share_blob",
             "opaque_peer_ciphertext",
             "dsVert/R/glmRing63DS.R exports an X share sealed to a pinned peer."),
    evidence("glmRing63ReceiveExtraShareDS", "receipt",
             "public_metadata",
             "dsVert/R/glmRing63DS.R stores a sealed peer share and returns a receipt."),
    evidence("glmRing63ReorderXFullDS", "receipt",
             "public_metadata",
             "dsVert/R/glmRing63DS.R mutates session layout and returns metadata."),
    evidence("glmRing63ShareExtraInputDS", "extra_input_blob",
             "opaque_peer_ciphertext",
             "dsVert/R/glmRing63DS.R seals extra input shares to pinned peers."),
    evidence("glmRing63TransportInitDS", "transport_public_key",
             "public_metadata",
             "dsVert/R/glmRing63DS.R returns transport identity metadata."),
    evidence("k2BeaverExtractColumnDS", "receipt",
             "public_metadata",
             "dsVert/R/k2InputSharingDS.R stores a share-domain column and returns metadata."),
    evidence("k2ComputeEtaShareDS", "receipt",
             "public_metadata",
             "dsVert/R/k2InputSharingDS.R computes session state and returns metadata."),
    evidence("k2ReceiveShareDS", "receipt",
             "public_metadata",
             "dsVert/R/k2InputSharingDS.R consumes a sealed peer share and returns metadata."),
    evidence("k2Ring127AffineCombineDS", "receipt",
             "public_metadata",
             "dsVert/R/k2InputSharingDS.R stores an affine share result and returns metadata."),
    evidence("k2ShareInputDS", "input_share_blob",
             "opaque_peer_ciphertext",
             "dsVert/R/k2InputSharingDS.R seals input shares to pinned peers."),
    evidence("mpcStoreTransportKeysDS", "receipt",
             "public_metadata",
             "dsVert/R/mpcSessionDS.R stores pinned transport keys and returns metadata."),
    evidence("dsvertOutcomeLevelsDS", "outcome_domain_and_counts",
             "plaintext_exact_aggregate",
             paste(
               "dsVert/R/multinomJointDS.R returns the outcome domain and",
               "guarded exact class counts.")),
    evidence("dsvertOrdinalShareClassMasksDS", "class_mask_blobs",
             "opaque_peer_ciphertext",
             paste(
               "dsVert/R/ordinalJointDS.R returns class-mask shares sealed",
               "to pinned peers.")),
    evidence("dsvertClusterSizesDS", "per_cluster_sizes",
             "plaintext_exact_aggregate",
             "dsVert/R/dsvertClusterDS.R returns exact guarded cluster sizes."),
    evidence("dsvertClusterResidualsDS", "per_cluster_residual_sums",
             "plaintext_exact_aggregate",
             paste(
               "dsVert/R/dsvertClusterDS.R returns exact residual and squared",
               "residual sums for every guarded cluster.")),
    evidence("dsvertClusterZtZDS", "per_cluster_ztz_sizes_and_levels",
             "plaintext_exact_aggregate",
             paste(
               "dsVert/R/dsvertClusterDS.R returns each guarded cluster's",
               "Z'Z matrix, size and level.")),
    evidence("dsvertLMMVarianceComponentsDS", "variance_component_summaries",
             "plaintext_exact_aggregate",
             paste(
               "dsVert/R/dsvertClusterDS.R returns exact SSW/SSB, per-cluster",
               "SSW and sizes.")),
    evidence("dsvertPerClusterSumShareDS", "per_cluster_sum_shares",
             "share_reconstructed_by_client",
             paste(
               "dsVert/R/clusterShareDS.R returns one additive share per",
               "cluster; its own documentation says the client combines parties.")),
    evidence("dsvertPerClusterSumShareDS", "cluster_sizes",
             "plaintext_exact_aggregate",
             "dsVert/R/clusterShareDS.R also returns exact guarded cluster sizes."),
    evidence("dsvertLMMPerClusterSumDS", "per_cluster_sum_shares",
             "share_reconstructed_by_client",
             paste(
               "dsVert/R/dsvertLMMClusterBroadcastDS.R returns additive",
               "per-cluster shares that the current K=2 client combines.")),
    evidence("dsvertGEEAR1OrderBroadcastDS", "peer_blob",
             "opaque_peer_ciphertext",
             paste(
               "dsVert/R/clusterShareDS.R encrypts predecessor/successor",
               "indices to the pinned peer; the relay sees only ciphertext.")),
    evidence("dsvertGEEAR1OrderBroadcastDS", "cluster_pair_dimensions",
             "public_metadata",
             paste(
               "dsVert/R/clusterShareDS.R returns guarded cluster/pair/lag",
               "dimensions alongside the ciphertext.")),
    evidence("dsvertPearsonR2ColDS", "missingness_counts",
             "plaintext_exact_aggregate",
             paste(
               "dsVert/R/dsvertPearsonR2ColDS.R mutates server state and",
               "returns exact observed/missing counts.")),
    evidence("k2ShareWeightsDS", "weight_share_blobs",
             "opaque_peer_ciphertext",
             paste(
               "dsVert/R/weightsDS.R pins both recipients and returns only",
               "sealed complementary weight shares.")),
    evidence("k2ShareWeightsDS", "dimensions_and_attestation",
             "public_metadata",
             "dsVert/R/weightsDS.R also returns n and numeric attestation metadata."),
    evidence("dsvertImputeColumnDS", "imputation_counts",
             "plaintext_exact_aggregate",
             paste(
               "dsVert/R/dsvertImputeColumnDS.R mutates the server data and",
               "returns exact n_imputed and n_observed on every round.")),
    evidence("dsvertImputeColumnDS", "model_diagnostics",
             "public_metadata",
             paste(
               "dsVert/R/dsvertImputeColumnDS.R returns imputation method,",
               "regularization and predictor-count metadata."))
  ))
  rownames(remote_catalog) <- NULL

  remote_evidence <- function(calls = character()) {
    calls <- unique(as.character(calls))
    if (!length(calls)) return(remote_catalog[FALSE, , drop = FALSE])
    unknown <- setdiff(calls, remote_catalog$call)
    if (length(unknown)) {
      stop("remote evidence is missing for: ", paste(unknown, collapse = ", "),
           call. = FALSE)
    }
    out <- remote_catalog[remote_catalog$call %in% calls, , drop = FALSE]
    out <- out[order(out$call, out$component, method = "radix"), ,
               drop = FALSE]
    rownames(out) <- NULL
    out
  }

  add <- function(methods, canonical_method, canonical_family,
                  artifact_requirements, estimand, inference_requirements,
                  migration_feasibility,
                  legacy_remote_calls = character(), aliases = character(),
                  alias_kinds = character(), current_route_status = NULL,
                  artifact_implementation_state = NULL,
                  inference_implementation_state = NULL) {
    methods <- as.character(methods)
    alias_of <- rep(NA_character_, length(methods))
    names(alias_of) <- methods
    alias_kind <- rep("canonical", length(methods))
    names(alias_kind) <- methods
    if (length(aliases)) {
      alias_names <- names(aliases)
      if (is.null(alias_names) || any(!alias_names %in% methods)) {
        stop("aliases must be a named subset of methods", call. = FALSE)
      }
      alias_of[alias_names] <- unname(aliases)
      alias_kind[alias_names] <- "compatibility_alias"
    }
    if (length(alias_kinds)) {
      alias_kind_names <- names(alias_kinds)
      if (is.null(alias_kind_names) ||
          any(!alias_kind_names %in% methods)) {
        stop("alias_kinds must be a named subset of methods", call. = FALSE)
      }
      alias_kind[alias_kind_names] <- unname(alias_kinds)
    }
    defaults <- switch(migration_feasibility,
      count_operation_implemented = list(
        route = "formal_sticky_count_artifact",
        artifact = "signed_count_artifact_implemented",
        inference = "formal_count_release_implemented"),
      frequency_operation_implemented = list(
        route = "formal_sticky_frequency_artifact",
        artifact = "signed_frequency_artifact_implemented",
        inference = "formal_frequency_release_implemented"),
      synopsis_release_implemented = list(
        route = "formal_sticky_synopsis_artifact",
        artifact = "signed_synopsis_describe_artifact_implemented",
        inference = "formal_synopsis_describe_release_implemented"),
      capsule_release_implemented = list(
        route = "legacy_joint_dp_capsule_incompatible",
        artifact = "joint_vector_release_implemented",
        inference = "legacy_capsule_release_incompatible"),
      requires_new_capsule_artifact = list(
        route = "legacy_exact_release_not_capsule_safe",
        artifact = "planned_no_materializer",
        inference = "existing_inference_requires_capsule_backend"),
      requires_new_secure_protocol = list(
        route = "legacy_granular_release_requires_secure_redesign",
        artifact = "secure_artifact_not_implemented",
        inference = "legacy_inference_requires_secure_redesign"),
      client_only_requires_attested_input = list(
        route = "client_only_inherits_input",
        artifact = "client_input_not_capsule_attested",
        inference = "implemented_client_algebra_inherits_input"),
      stop("unknown migration_feasibility", call. = FALSE))
    if (is.null(current_route_status)) {
      current_route_status <- defaults$route
    }
    if (is.null(artifact_implementation_state)) {
      artifact_implementation_state <- defaults$artifact
    }
    if (is.null(inference_implementation_state)) {
      inference_implementation_state <- defaults$inference
    }
    remote <- remote_evidence(legacy_remote_calls)
    n <- length(methods)
    rows[[length(rows) + 1L]] <<- data.frame(
      method = methods,
      canonical_method = rep(canonical_method, n),
      canonical_family = rep(canonical_family, n),
      alias_of = unname(alias_of),
      alias_kind = unname(alias_kind),
      artifact_requirements = I(rep(list(sort(unique(
        as.character(artifact_requirements)), method = "radix")), n)),
      estimand = rep(estimand, n),
      inference_requirements = I(rep(list(sort(unique(
        as.character(inference_requirements)), method = "radix")), n)),
      current_route_status = rep(current_route_status, n),
      same_capsule_replay_history_can_deny = rep(FALSE, n),
      # Historical capsule validators remain available for old certificates,
      # but no current public route may be denied by lifetime or reservation
      # history.
      new_capsule_reservation_history_can_deny = rep(FALSE, n),
      migration_feasibility = rep(migration_feasibility, n),
      artifact_implementation_state = rep(
        artifact_implementation_state, n),
      inference_implementation_state = rep(
        inference_implementation_state, n),
      legacy_remote_call_evidence = I(rep(list(remote), n)),
      stringsAsFactors = FALSE
    )
  }

  add(
    c("ds.vertDesc", "ds.vert.desc"), "ds.vertDesc", "descriptive",
    c("admitted_count", "fixed_numeric_histograms", "numeric_moments",
      "signed_synopsis_describe_artifact"),
    "Bounded univariate moments, fixed-grid distributions and quantiles.",
    c("finite_public_bounds", "fixed_histogram_grid", "fixed_workload",
      "mechanism_uncertainty"),
    "synopsis_release_implemented",
    character(),
    c("ds.vert.desc" = "ds.vertDesc"),
    artifact_implementation_state =
      "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    c("ds.vertCor", "ds.vert.cor"), "ds.vertCor", "correlation",
    c("joint_dp_noise", "gaussian_models",
      "signed_complete_case_gaussian_artifact",
      "validated_synopsis_provenance"),
    paste(
      "Bounded joint complete-case Pearson correlations for signed",
      "same-owner Gaussian Synopsis artifacts; cross-owner requests are",
      "quarantined without a capsule fallback."),
    c("finite_public_bounds", "identifiability", "mechanism_uncertainty",
      "positive_semidefinite_projection"),
    "synopsis_release_implemented", character(),
    c("ds.vert.cor" = "ds.vertCor"),
    artifact_implementation_state =
      "validated_complete_case_gaussian_synopsis_adapter_implemented",
    inference_implementation_state =
      "complete_case_postprocess_and_explicit_psd_projection_implemented")

  add(
    c("ds.vertPCA", "ds.vert.pca"), "ds.vertPCA", "pca",
    c("joint_dp_noise", "gaussian_models",
      "signed_complete_case_gaussian_artifact",
      "validated_synopsis_provenance"),
    paste(
      "Eigenstructure of the explicitly PSD-projected bounded DP joint",
      "complete-case correlation matrix."),
    c("finite_public_bounds", "mechanism_uncertainty",
      "positive_semidefinite_projection", "spectral_stability"),
    "synopsis_release_implemented", character(),
    c("ds.vert.pca" = "ds.vertPCA"),
    artifact_implementation_state =
      "validated_complete_case_gaussian_synopsis_adapter_implemented",
    inference_implementation_state =
      "client_only_spectral_postprocess_with_eigengap_regions")

  add(
    c("ds.vertChisq", "ds.vert.chisq"), "ds.vertChisq",
    "categorical_independence",
    c("admitted_count", "categorical_marginals",
      "categorical_pairs_same_owner", "joint_dp_noise"),
    "Pearson/Yates independence statistic for a fixed-domain contingency table.",
    c("dp_aware_null_distribution", "fixed_category_domain",
      "minimum_expected_cell_handling"),
    "synopsis_release_implemented", character(),
    c("ds.vert.chisq" = "ds.vertChisq"),
    artifact_implementation_state =
      "validated_same_owner_synopsis_adapter_implemented",
    inference_implementation_state =
      "dp_aware_parametric_bootstrap_implemented")

  add(
    c("ds.vertFisher", "ds.vert.fisher"), "ds.vertFisher",
    "categorical_fisher",
    c("admitted_count", "categorical_marginals",
      "categorical_pairs_same_owner", "joint_dp_noise"),
    paste(
      "DP-aware conditional hypergeometric plug-in test for a signed",
      "fixed-domain 2x2 release; the displayed cross-product odds ratio",
      "is descriptive and is not Fisher's conditional maximum-likelihood",
      "estimate."),
    c("dp_aware_null_distribution", "fixed_category_domain",
      "signed_sampler_total_variation", "two_by_two_domain"),
    "synopsis_release_implemented", character(),
    c("ds.vert.fisher" = "ds.vertFisher"),
    artifact_implementation_state =
      "validated_same_owner_synopsis_adapter_implemented",
    inference_implementation_state =
      "dp_aware_conditional_hypergeometric_bootstrap_implemented")

  add(
    c("ds.vertChisqCross", "ds.vert.chisq_cross"), "ds.vertChisqCross",
    "cross_vertical_categorical",
    c("categorical_pairs_cross_owner", "joint_dp_noise",
      "private_psi_alignment", "signed_no_lifetime_projection_protocol",
      "validated_synopsis_provenance"),
    paste(
      "DP-aware Pearson/Yates or conditional hypergeometric plug-in",
      "inference for fixed-domain variables held by different peers."),
    c("dp_aware_null_distribution", "fixed_category_domain",
      "numeric_certificate", "signed_sampler_total_variation",
      "cross_signed_allocation_before_source_access"),
    "synopsis_release_implemented", character(),
    c("ds.vert.chisq_cross" = "ds.vertChisqCross"),
    current_route_status = "formal_sticky_synopsis_artifact",
    artifact_implementation_state =
      "validated_synopsis_adapter_implemented",
    inference_implementation_state =
      "dp_aware_parametric_bootstrap_implemented")

  add(
    "ds.vertDPGaussian", "ds.vertDPGaussian", "gaussian_regression",
    c("admitted_count", "complete_case_patient_collapse",
      "gaussian_sufficient_statistics_same_owner", "joint_dp_noise",
      "public_numeric_bounds", "signed_gaussian_model_artifact",
      "validated_synopsis_provenance"),
    paste(
      "Bounded, clipped same-owner complete-case Gaussian",
      "least-squares coefficients for one signed public model; a positive",
      "ridge requests an explicitly different normalized-design estimand;",
      "cross-owner descriptors fail closed without capsule fallback."),
    c("explicit_regularization_estimand", "finite_public_bounds",
      "fixed_complete_case_rule", "identifiability", "mechanism_uncertainty",
      "no_sampling_inference", "positive_semidefinite_projection"),
    "synopsis_release_implemented", character(),
    artifact_implementation_state =
      "validated_same_owner_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    c("ds.vertGLM", "ds.vert.glm"), "ds.vertGLM", "glm",
    c("admitted_count", "model_deviance", "model_score_hessian",
      "numeric_cross_products_same_and_cross_owner", "numeric_moments"),
    paste(
      "Bounded Gaussian coefficients only when explicit dp_analysis_id selects",
      "the same-owner Synopsis adapter; cross-owner, identifier-free,",
      "binomial, and Poisson calls fail closed without capsule fallback;",
      "the unsupported variants fail before protected analysis."),
    c("convergence", "dp_aware_covariance", "identifiability",
      "numeric_certificate"),
    "requires_new_capsule_artifact",
    c("getObsCountDS", "glmRing63DevianceSumsDS", "glmStandardizeDS",
      "k2GradientR1DS", "k2GradientR2DS"),
    c("ds.vert.glm" = "ds.vertGLM"),
    current_route_status =
      "formal_same_owner_synopsis_variant_only_legacy_unavailable")

  add(
    c("ds.vertCox", "ds.vertCoxProfileNonDisclosive", "ds.vert.coxph"),
    "ds.vertCoxProfileNonDisclosive", "cox_profile",
    c("admitted_count", "cox_baseline_hazard", "cox_partial_likelihood",
      "cox_risk_sets", "cox_score_hessian", "numeric_cross_products"),
    paste(
      "Cox proportional-hazards log hazard ratios under an explicit",
      "ties, strata and delayed-entry contract."),
    c("baseline_hazard_contract", "convergence", "delayed_entry_contract",
      "dp_aware_covariance", "identifiability", "strata_contract",
      "ties_method"),
    "requires_new_capsule_artifact",
    c("dsvertCoxEventTimeShareMaskDS", "k2BeaverSumShareDS"),
    c("ds.vertCox" = "ds.vertCoxProfileNonDisclosive",
      "ds.vert.coxph" = "ds.vertCoxProfileNonDisclosive"),
    alias_kinds = c(
      "ds.vertCox" = "compatibility_wrapper",
      "ds.vert.coxph" = "compatibility_wrapper"))

  add(
    "ds.vert.cox", "ds.vert.cox", "cox_dispatch",
    c("admitted_count", "cox_baseline_hazard", "cox_partial_likelihood",
      "cox_risk_sets", "cox_score_hessian", "discrete_hazard_grid",
      "numeric_cross_products"),
    "Route-selected Cox PH or explicitly labelled discrete-time hazard coefficients.",
    c("convergence", "delayed_entry_contract", "dp_aware_covariance",
      "estimand_label", "identifiability", "strata_contract",
      "ties_method"),
    "requires_new_capsule_artifact",
    c("dsvertCoxDiscreteExpandXDS", "dsvertCoxEventTimeShareMaskDS",
      "k2BeaverSumShareDS"))

  add(
    "ds.vertCoxDiscreteNonDisclosive", "ds.vertCoxDiscreteNonDisclosive",
    "discrete_time_hazard",
    c("admitted_count", "discrete_hazard_grid", "discrete_risk_masks",
      "model_score_hessian", "numeric_cross_products"),
    "Discrete-time pooled-logistic hazard coefficients (not a Cox estimand).",
    c("convergence", "dp_aware_covariance", "estimand_label",
      "fixed_time_grid", "identifiability"),
    "requires_new_capsule_artifact",
    c("dsvertCoxDiscreteExpandXDS", "dsvertCoxEventTimeShareMaskDS",
      "k2BeaverSumShareDS"))

  add(
    c("ds.vertNBFullRegTheta", "ds.vert.nb"), "ds.vertNBFullRegTheta",
    "negative_binomial",
    c("admitted_count", "nb2_beta_score_information",
      "nb2_beta_theta_cross_information", "nb2_deviance",
      "nb2_theta_score_information", "numeric_cross_products",
      "numeric_moments"),
    "NB2 log-link regression coefficients and size/dispersion theta.",
    c("convergence", "dp_aware_joint_beta_theta_covariance",
      "identifiability", "nb2_parameterization", "theta_link",
      "theta_uncertainty"),
    "requires_new_secure_protocol",
    c("dsvertGEEInterceptShareDS", "dsvertLocalMomentsDS",
      "dsvertNBEtaShareConfirmDS", "dsvertNBEtaShareDS",
      "dsvertNBEtaTotalReceiveDS", "dsvertNBProfileSumsDS",
      "dsvertNBPsiAggregateDS", "dsvertNBYThetaShareDS",
      "dsvertNBYThetaShareReceiveDS", "getObsCountDS",
      "glmRing63DevianceSumsDS", "glmRing63ExportOwnShareDS",
      "glmRing63ReceiveExtraShareDS", "glmRing63ReorderXFullDS",
      "glmRing63ShareExtraInputDS", "glmRing63TransportInitDS",
      "glmStandardizeDS", "k2BeaverExtractColumnDS",
      "k2BeaverSumShareDS", "k2ComputeEtaShareDS", "k2GradientR1DS",
      "k2GradientR2DS", "k2ReceiveShareDS",
      "k2Ring127AffineCombineDS", "k2ShareInputDS",
      "mpcStoreTransportKeysDS"),
    c("ds.vert.nb" = "ds.vertNBFullRegTheta"),
    current_route_status = "legacy_granular_release_quarantine",
    inference_implementation_state =
      "nb2_joint_beta_theta_inference_not_implemented")

  add(
    c("ds.vertMultinom", "ds.vertMultinomJoint",
      "ds.vertMultinomJointNewton", "ds.vert.multinom"),
    "ds.vertMultinomJointNewton", "multinomial_regression",
    c("admitted_count", "categorical_outcome_domain",
      "model_score_hessian", "numeric_cross_products_same_and_cross_owner",
      "numeric_moments", "signed_multinomial_design_gram"),
    "Joint softmax log-odds coefficients relative to the declared reference class.",
    c("convergence", "dp_aware_covariance", "identifiability",
      "outcome_domain"),
    "requires_new_capsule_artifact",
    c("dsvertLocalMomentsDS", "dsvertOutcomeLevelsDS", "getObsCountDS",
      "glmStandardizeDS", "k2GradientR1DS", "k2GradientR2DS"),
    c("ds.vertMultinom" = "ds.vertMultinomJointNewton",
      "ds.vertMultinomJoint" = "ds.vertMultinomJointNewton",
      "ds.vert.multinom" = "ds.vertMultinomJointNewton"),
    alias_kinds = c(
      "ds.vertMultinom" = "compatibility_wrapper",
      "ds.vertMultinomJoint" = "compatibility_wrapper"),
    current_route_status = "signed_workload_unavailable_quarantine",
    inference_implementation_state =
      "signed_multinomial_design_gram_not_materialized")

  add(
    c("ds.vertOrdinal", "ds.vertOrdinalJointNewton", "ds.vert.ordinal"),
    "ds.vertOrdinalJointNewton", "ordinal_regression",
    c("admitted_count", "categorical_outcome_domain",
      "model_score_hessian", "numeric_cross_products_same_and_cross_owner"),
    "Joint proportional-odds thresholds and slopes.",
    c("convergence", "dp_aware_covariance", "identifiability",
      "ordered_outcome_domain", "proportional_odds_diagnostics"),
    "requires_new_capsule_artifact",
    c("dsvertOrdinalShareClassMasksDS", "getObsCountDS",
      "glmStandardizeDS", "k2GradientR1DS", "k2GradientR2DS"),
    c("ds.vertOrdinal" = "ds.vertOrdinalJointNewton",
      "ds.vert.ordinal" = "ds.vertOrdinalJointNewton"),
    alias_kinds = c("ds.vertOrdinal" = "compatibility_wrapper"))

  add(
    c("ds.vertLMM", "ds.vert.lmm"), "ds.vertLMM",
    "linear_mixed_model",
    c("bounded_cluster_cross_products", "bounded_cluster_moments",
      "cluster_random_effect_design", "lmm_score_information"),
    "Random-intercept or declared random-slope LMM fixed effects and variance components.",
    c("cluster_contribution_bounds", "dp_aware_covariance",
      "identifiability", "ml_reml_contract", "random_effects_contract"),
    "requires_new_secure_protocol",
    c("dsvertClusterResidualsDS", "dsvertClusterSizesDS",
      "dsvertClusterZtZDS", "dsvertLMMPerClusterSumDS",
      "dsvertLMMVarianceComponentsDS", "dsvertPerClusterSumShareDS",
      "getObsCountDS", "glmStandardizeDS", "k2GradientR1DS",
      "k2GradientR2DS"),
    c("ds.vert.lmm" = "ds.vertLMM"),
    current_route_status = "legacy_granular_release_quarantine")

  add(
    "ds.vertLMM.k3", "ds.vertLMM", "linear_mixed_model",
    c("bounded_cluster_cross_products", "bounded_cluster_moments",
      "cluster_random_effect_design", "lmm_score_information"),
    paste(
      "Deprecated K>=3 subroute for random-intercept or declared",
      "random-slope LMM fixed effects and variance components."),
    c("cluster_contribution_bounds", "dp_aware_covariance",
      "identifiability", "k_ge_3_only", "ml_reml_contract",
      "random_effects_contract"),
    "requires_new_secure_protocol",
    c("dsvertClusterResidualsDS", "dsvertClusterSizesDS",
      "dsvertClusterZtZDS", "dsvertLMMPerClusterSumDS",
      "dsvertLMMVarianceComponentsDS", "dsvertPerClusterSumShareDS",
      "getObsCountDS", "glmStandardizeDS", "k2GradientR1DS",
      "k2GradientR2DS"),
    alias_kinds = c("ds.vertLMM.k3" = "deprecated_subroute"),
    current_route_status = "legacy_granular_release_quarantine")

  add(
    c("ds.vertGEE", "ds.vert.gee"), "ds.vertGEE", "gee",
    c("bounded_cluster_score_meat", "cluster_working_correlation",
      "model_score_hessian", "numeric_cross_products_same_and_cross_owner"),
    "Population-average GEE coefficients and robust covariance.",
    c("cluster_contribution_bounds", "convergence", "dp_aware_covariance",
      "identifiability", "working_correlation_contract"),
    "requires_new_secure_protocol",
    c("dsvertGEEAR1OrderBroadcastDS", "dsvertPearsonR2ColDS",
      "dsvertPerClusterSumShareDS", "getObsCountDS", "glmStandardizeDS",
      "k2GradientR1DS", "k2GradientR2DS"),
    c("ds.vert.gee" = "ds.vertGEE"),
    current_route_status = "legacy_granular_release_quarantine")

  add(
    c("ds.vertGLMM", "ds.vert.glmm"), "ds.vertGLMM", "glmm",
    c("bounded_cluster_moments", "glmm_working_score_information",
      "numeric_cross_products_same_and_cross_owner"),
    "Conditional binomial random-intercept PQL coefficients and variance.",
    c("cluster_contribution_bounds", "convergence", "dp_aware_covariance",
      "identifiability", "target_estimator_validation"),
    "requires_new_secure_protocol",
    c("dsvertClusterSizesDS", "dsvertPerClusterSumShareDS",
      "getObsCountDS", "glmStandardizeDS", "k2GradientR1DS",
      "k2GradientR2DS"),
    c("ds.vert.glmm" = "ds.vertGLMM"),
    current_route_status = "legacy_granular_release_quarantine")

  add(
    c("ds.vertIPW", "ds.vert.ipw"), "ds.vertIPW", "ipw",
    c("bounded_weight_distribution", "outcome_score_information",
      "propensity_score_information", "treatment_outcome_binding"),
    paste(
      "Declared ATE, ATT or ATC marginal effect under a bound propensity",
      "and outcome-model contract."),
    c("consistency_assumption", "dp_aware_covariance",
      "exchangeability_assumption", "positivity", "stabilization_contract",
      "treatment_outcome_contract", "weight_clipping_and_bounds",
      "weight_provenance_binding"),
    "requires_new_secure_protocol",
    c("getObsCountDS", "glmRing63DevianceSumsDS", "glmStandardizeDS",
      "k2GradientR1DS", "k2GradientR2DS", "k2ShareWeightsDS"),
    c("ds.vert.ipw" = "ds.vertIPW"),
    current_route_status = "legacy_granular_release_quarantine",
    inference_implementation_state =
      "propensity_to_weight_estimand_pipeline_not_implemented")

  add(
    c("ds.vertMI", "ds.vert.mi"), "ds.vertMI", "multiple_imputation",
    c("bounded_missingness_counts", "imputation_model_sufficient_statistics",
      "posterior_parameter_draws", "synthetic_imputation_draws"),
    "Rubin-pooled regression coefficients under the declared imputation model.",
    c("between_imputation_variance", "chained_equation_predictor_contract",
      "cryptographic_non_rerollable_draw_stream",
      "imputation_model_validity", "missingness_mechanism_assumption",
      "posterior_parameter_uncertainty", "rubin_small_sample_df",
      "synthetic_release_provenance"),
    "requires_new_secure_protocol",
    c("dsvertColNamesDS", "dsvertImputeColumnDS", "getObsCountDS",
      "glmRing63DevianceSumsDS", "glmStandardizeDS", "k2GradientR1DS",
      "k2GradientR2DS"),
    c("ds.vert.mi" = "ds.vertMI"),
    current_route_status = "legacy_mutating_release_quarantine",
    inference_implementation_state =
      "multiple_imputation_protocol_and_inference_require_redesign")

  add(
    c("ds.vertLASSO", "ds.vert.lasso"), "ds.vertLASSO",
    "lasso_gaussian",
    c("admitted_count", "gaussian_sufficient_statistics_same_owner",
      "signed_gaussian_model_artifact",
      "validated_gaussian_provenance_certificate"),
    "Bound-normalized Gaussian LASSO path over one signed projected DP moment artifact.",
    c("certificate_integrity_validation", "gaussian_only",
      "identifiability", "kkt_validation", "no_sampling_inference",
      "objective_scale_contract", "authentic_federation_e2e_validation"),
    "synopsis_release_implemented", character(),
    aliases = c("ds.vert.lasso" = "ds.vertLASSO"),
    current_route_status = "client_only_validated_synopsis_postprocess",
    artifact_implementation_state =
      "validated_same_owner_synopsis_adapter_implemented",
    inference_implementation_state = "dp_gaussian_lasso_path_implemented")

  add(
    c("ds.vertLASSO1Step", "ds.vert.lasso_1step"), "ds.vertLASSO1Step",
    "lasso_gaussian",
    c("admitted_count", "gaussian_sufficient_statistics_same_owner",
      "signed_gaussian_model_artifact",
      "validated_gaussian_provenance_certificate"),
    "Historical one-step frontdoor for the bound-normalized Gaussian LASSO path over one signed projected DP moment artifact.",
    c("certificate_integrity_validation", "gaussian_only",
      "identifiability", "kkt_validation", "no_sampling_inference",
      "objective_scale_contract", "authentic_federation_e2e_validation"),
    "synopsis_release_implemented", character(),
    aliases = c("ds.vert.lasso_1step" = "ds.vertLASSO1Step"),
    current_route_status = "client_only_validated_synopsis_postprocess",
    artifact_implementation_state =
      "validated_same_owner_synopsis_adapter_implemented",
    inference_implementation_state = "dp_gaussian_lasso_path_implemented")

  add(
    c("ds.vertLASSOProximal", "ds.vert.lasso_proximal"),
    "ds.vertLASSOProximal", "lasso_gaussian",
    c("admitted_count", "gaussian_sufficient_statistics_same_owner",
      "signed_gaussian_model_artifact",
      "validated_gaussian_provenance_certificate"),
    "Bound-normalized Gaussian LASSO on one signed projected DP moment artifact.",
    c("certificate_integrity_validation", "gaussian_only",
      "identifiability", "kkt_validation", "no_sampling_inference",
      "objective_scale_contract", "authentic_federation_e2e_validation"),
    "synopsis_release_implemented", character(),
    aliases = c("ds.vert.lasso_proximal" = "ds.vertLASSOProximal"),
    current_route_status = "client_only_validated_synopsis_postprocess",
    artifact_implementation_state =
      "validated_same_owner_synopsis_adapter_implemented",
    inference_implementation_state = "dp_gaussian_lasso_path_implemented")

  add(
    c("ds.vertLASSOIter", "ds.vert.lasso_iter"), "ds.vertLASSOIter",
    "lasso_iterative",
    c("admitted_count", "gaussian_sufficient_statistics_same_owner",
      "signed_gaussian_model_artifact",
      "validated_gaussian_provenance_certificate"),
    "Gaussian L1 coefficient paths from one signed Synopsis artifact.",
    c("gaussian_only", "identifiability", "kkt_validation",
      "no_sampling_inference", "objective_scale_contract"),
    "synopsis_release_implemented", character(),
    c("ds.vert.lasso_iter" = "ds.vertLASSOIter"),
    current_route_status =
      "formal_same_owner_synopsis_variant_only_legacy_unavailable",
    artifact_implementation_state =
      "validated_same_owner_synopsis_adapter_implemented",
    inference_implementation_state =
      "dp_gaussian_lasso_path_implemented")

  add(
    c("ds.vertLASSOCV", "ds.vert.lasso_cv"), "ds.vertLASSOCV",
    "lasso_information_criterion",
    c("admitted_count", "gaussian_sufficient_statistics_same_owner",
      "signed_gaussian_model_artifact",
      "validated_gaussian_provenance_certificate"),
    paste(
      "DP-projected pseudo-AIC/BIC/EBIC selection on one signed Gaussian",
      "moment artifact; it is not cross-validation."),
    c("certificate_integrity_validation", "estimand_label",
      "no_sampling_inference", "pseudo_information_criterion_label",
      "selection_uncertainty", "authentic_federation_e2e_validation"),
    "synopsis_release_implemented", character(),
    aliases = c("ds.vert.lasso_cv" = "ds.vertLASSOCV"),
    current_route_status = "client_only_validated_synopsis_postprocess",
    artifact_implementation_state =
      "validated_same_owner_synopsis_adapter_implemented",
    inference_implementation_state = "dp_gaussian_pseudo_ic_implemented")

  add(
    c("ds.vertLR", "ds.vert.lr"), "ds.vertLR", "likelihood_ratio",
    c("authorized_full_fit", "authorized_reduced_fit",
      "canonical_unweighted_deviance"),
    paste(
      "Likelihood-ratio comparison of compatible nested, converged,",
      "unpenalized binomial or Poisson models."),
    c("authorized_fit_provenance", "binomial_or_poisson_only",
      "canonical_unweighted_deviance", "converged_unpenalized_fits",
      "nested_model_contract", "same_cohort_missingness_and_offset"),
    "client_only_requires_attested_input",
    aliases = c("ds.vert.lr" = "ds.vertLR"))

  add(
    c("ds.vertConfint", "ds.vert.confint"), "ds.vertConfint",
    "confidence_interval",
    c("authorized_fit_coefficients", "authorized_fit_covariance"),
    "Wald confidence intervals for an authorized fitted model.",
    c("authorized_fit_provenance", "valid_covariance"),
    "client_only_requires_attested_input",
    aliases = c("ds.vert.confint" = "ds.vertConfint"))

  add(
    c("ds.vertWald", "ds.vert.wald"), "ds.vertWald", "wald_test",
    c("authorized_fit_coefficients", "authorized_fit_covariance"),
    "Wald test for declared fitted-model parameters.",
    c("authorized_fit_provenance", "valid_covariance"),
    "client_only_requires_attested_input",
    aliases = c("ds.vert.wald" = "ds.vertWald"))

  add(
    c("ds.vertContrast", "ds.vert.contrast"), "ds.vertContrast",
    "linear_contrast",
    c("authorized_fit_coefficients", "authorized_fit_covariance"),
    "Linear contrasts of authorized fitted-model coefficients.",
    c("authorized_fit_provenance", "estimable_contrast",
      "valid_covariance"),
    "client_only_requires_attested_input",
    aliases = c("ds.vert.contrast" = "ds.vertContrast"))

  add(
    "ds.vertEpi2x2", "ds.vertEpi2x2", "epi_2x2",
    "authorized_table",
    paste(
      "Risk ratio, odds ratio, risk difference, attributable fraction",
      "among the exposed and in the population, and number needed to",
      "benefit/harm from a 2x2 table."),
    c("authorized_input_provenance", "sampling_uncertainty"),
    "client_only_requires_attested_input")

  add(
    "ds.vertMantelHaenszel", "ds.vertMantelHaenszel",
    "mantel_haenszel",
    "authorized_stratified_2x2_tables",
    paste(
      "Common Mantel-Haenszel odds ratio and conditional classical",
      "inference from consistently oriented strata."),
    c("authorized_input_provenance", "common_odds_ratio_model",
      "public_cell_mapping", "public_strata", "sampling_uncertainty"),
    "client_only_requires_attested_input")

  add(
    "ds.vertDirectStandardization", "ds.vertDirectStandardization",
    "direct_standardization",
    "authorized_strata_table",
    paste(
      "Directly standardized rate from stratum-specific cases and",
      "person-time under public reference weights."),
    c("authorized_input_provenance", "compatible_strata",
      "sampling_uncertainty"),
    "client_only_requires_attested_input")

  add(
    "ds.vertIndirectStandardization", "ds.vertIndirectStandardization",
    "indirect_standardization",
    "authorized_strata_table",
    "Observed-to-expected standardized mortality/incidence ratio.",
    c("authorized_input_provenance", "compatible_populations",
      "sampling_uncertainty"),
    "client_only_requires_attested_input")

  add(
    "ds.vertDPCount", "ds.vertDPCount", "dp_count",
    c("current_aligned_snapshot", "signed_count_analysis_contract",
      "joint_dp_noise_or_fixed_public_consensus"),
    "Privacy-unit count for the canonical signed current snapshot.",
    c("adjacency_contract", "mechanism_uncertainty",
      "signed_K_peer_provenance"),
    "count_operation_implemented")

  add(
    "ds.vertDPContingency", "ds.vertDPContingency", "dp_contingency",
    c("admitted_count", "categorical_marginals",
      "categorical_pairs_same_owner", "categorical_pairs_cross_owner",
      "joint_dp_noise", "signed_categorical_pair_projection",
      "validated_synopsis_provenance"),
    paste(
      "Fixed-domain noisy contingency table for one signed same- or",
      "cross-owner categorical pair."),
    c("adjacency_contract", "fixed_category_domain",
      "mechanism_uncertainty", "validated_synopsis_provenance"),
    "synopsis_release_implemented",
    artifact_implementation_state =
      "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPFrequency", "ds.vertDPFrequency", "dp_frequency",
    c("current_aligned_snapshot", "fixed_category_domain",
      "signed_frequency_analysis_contract",
      "sticky_two_authority_ring128_vector"),
    paste(
      "Fixed-domain finite-snapshot frequency distribution for one declared",
      "categorical variable after the signed repeated-record collapse."),
    c("adjacency_contract", "explicit_source_owner", "fixed_category_domain",
      "mechanism_uncertainty", "signed_K_peer_provenance", "sticky_replay"),
    "frequency_operation_implemented")

  add(
    "ds.vertDPFrequencyInference", "ds.vertDPFrequencyInference",
    "dp_frequency_sampling_inference",
    "validated_sticky_frequency_artifact",
    paste(
      "Conservative simultaneous population-proportion regions combining",
      "the signed DP count-box event with exact multinomial sampling",
      "uncertainty."),
    c("iid_multinomial_sampling_model", "clopper_pearson_exact_intervals",
      "joint_mechanism_and_sampling_uncertainty",
      "validated_frequency_provenance", "zero_call_postprocessing"),
    "frequency_operation_implemented",
    current_route_status = "client_only_inherits_input",
    artifact_implementation_state =
      "validated_frequency_artifact_adapter_implemented",
    inference_implementation_state = "frequency_postprocess_implemented")

  add(
    "ds.vertDPMeanVar", "ds.vertDPMeanVar", "dp_mean_variance",
    c("admitted_count", "joint_dp_noise", "numeric_moments",
      "signed_synopsis_artifact"),
    "Clipped bounded mean and variance for one declared numeric variable.",
    c("adjacency_contract", "finite_public_bounds",
      "mechanism_uncertainty", "validated_synopsis_provenance"),
    "synopsis_release_implemented",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPCor", "ds.vertDPCor", "dp_correlation",
    c("joint_dp_noise", "numeric_pair_moments_same_owner",
      "signed_correlation_artifact"),
    paste(
      "Pairwise-complete bounded correlations and explicit PSD",
      "post-processing for one signed same-owner variable set."),
    c("adjacency_contract", "finite_public_bounds", "identifiability",
      "mechanism_uncertainty", "validated_synopsis_provenance"),
    "synopsis_release_implemented",
    artifact_implementation_state =
      "validated_synopsis_adapter_implemented",
    inference_implementation_state =
      "pairwise_postprocess_and_explicit_psd_projection_implemented")

  add(
    "ds.vertDPDescribe", "ds.vertDPDescribe", "dp_descriptive",
    c("admitted_count", "categorical_marginals",
      "fixed_numeric_histograms", "joint_dp_noise", "numeric_moments",
      "signed_synopsis_describe_artifact"),
    "Bounded moments and fixed-grid quantiles for a custodian-declared variable set.",
    c("adjacency_contract", "finite_public_bounds", "fixed_workload",
      "mechanism_uncertainty", "validated_synopsis_provenance"),
    "synopsis_release_implemented")

  add(
    "ds.vertDPQuantile", "ds.vertDPQuantile",
    "dp_descriptive_postprocess",
    "validated_synopsis_describe_artifact",
    paste(
      "Fixed-grid binned quantiles from the coordinatewise-",
      "nonnegative histogram in one validated DP describe release."),
    c("fixed_public_grid", "mechanism_uncertainty",
      "sampling_uncertainty_if_claimed", "validated_synopsis_provenance"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPMedian", "ds.vertDPMedian", "dp_descriptive_postprocess",
    "validated_synopsis_describe_artifact",
    paste(
      "Fixed-grid binned medians from the coordinatewise-nonnegative",
      "histogram in one validated DP describe release."),
    c("fixed_public_grid", "mechanism_uncertainty",
      "sampling_uncertainty_if_claimed", "validated_synopsis_provenance"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPSurvival", "ds.vertDPSurvival", "dp_survival_release",
    c("admitted_count", "joint_dp_noise", "survival_fixed_grid"),
    "Fixed-grid entry, censoring and cause-specific event release.",
    c("adjacency_contract", "fixed_time_grid", "mechanism_uncertainty",
      "validated_synopsis_provenance"),
    "synopsis_release_implemented",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPKaplanMeier", "ds.vertDPKaplanMeier",
    "dp_survival_postprocess", "validated_synopsis_survival_artifact",
    "Kaplan-Meier survival curve from one validated DP survival artifact.",
    c("mechanism_uncertainty", "sampling_uncertainty_if_claimed",
      "validated_synopsis_provenance"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPNelsonAalen", "ds.vertDPNelsonAalen",
    "dp_survival_postprocess", "validated_synopsis_survival_artifact",
    "Nelson-Aalen cumulative hazard from one validated DP survival artifact.",
    c("mechanism_uncertainty", "sampling_uncertainty_if_claimed",
      "validated_synopsis_provenance"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPCumulativeIncidence", "ds.vertDPCumulativeIncidence",
    "dp_survival_postprocess", "validated_synopsis_survival_artifact",
    "Cause-specific cumulative incidence from one validated DP survival artifact.",
    c("competing_risk_contract", "mechanism_uncertainty",
      "sampling_uncertainty_if_claimed", "validated_synopsis_provenance"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPRMST", "ds.vertDPRMST", "dp_survival_postprocess",
    "validated_synopsis_survival_artifact",
    "Fixed-grid restricted mean survival time through public tau.",
    c("fixed_time_grid", "mechanism_uncertainty",
      "sampling_uncertainty_if_claimed", "validated_synopsis_provenance"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPRMTL", "ds.vertDPRMTL", "dp_survival_postprocess",
    "validated_synopsis_survival_artifact",
    paste(
      "Fixed-grid restricted mean time lost over the public interval as the",
      "exact restriction-width complement of released RMST."),
    c("fixed_time_grid", "mechanism_uncertainty",
      "sampling_uncertainty_if_claimed", "validated_synopsis_provenance",
      "zero_call_numeric_identity"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPSurvivalContrast", "ds.vertDPSurvivalContrast",
    "dp_survival_postprocess",
    "two_compatible_validated_synopsis_survival_artifacts",
    paste(
      "Fixed-grid comparison-minus-reference survival difference and",
      "comparison/reference survival ratio."),
    c("identical_signed_public_time_grid", "mechanism_uncertainty",
      "bonferroni_joint_event_for_distinct_releases",
      "typed_zero_denominator", "sampling_uncertainty_if_claimed",
      "validated_synopsis_provenance", "zero_call_postprocessing"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPRMSTContrast", "ds.vertDPRMSTContrast",
    "dp_survival_postprocess",
    "two_compatible_validated_synopsis_survival_artifacts",
    paste(
      "Fixed-grid comparison-minus-reference RMST difference and",
      "comparison/reference RMST ratio through common tau."),
    c("identical_signed_public_time_grid", "fixed_time_grid",
      "mechanism_uncertainty",
      "bonferroni_joint_event_for_distinct_releases",
      "typed_zero_denominator", "sampling_uncertainty_if_claimed",
      "validated_synopsis_provenance", "zero_call_postprocessing"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPSurvivalQuantile", "ds.vertDPSurvivalQuantile",
    "dp_survival_postprocess", "validated_synopsis_survival_artifact",
    paste(
      "Fixed-grid survival quantiles as first public endpoints crossing",
      "requested event-distribution probabilities, with inverted",
      "simultaneous mechanism limits and explicit beyond-grid states."),
    c("fixed_time_grid", "mechanism_uncertainty",
      "sampling_uncertainty_if_claimed", "validated_synopsis_provenance",
      "beyond_grid_estimability_state"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPMedianSurvival", "ds.vertDPMedianSurvival",
    "dp_survival_postprocess", "validated_synopsis_survival_artifact",
    paste(
      "Median survival as the exact probability-one-half view of the",
      "validated fixed-grid survival-quantile result."),
    c("fixed_time_grid", "mechanism_uncertainty",
      "sampling_uncertainty_if_claimed", "validated_synopsis_provenance",
      "beyond_grid_estimability_state", "zero_call_numeric_identity"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPEpi2x2", "ds.vertDPEpi2x2", "dp_epi_2x2",
    "validated_synopsis_2x2_artifact",
    paste(
      "Group and population risks, risk difference, risk/odds ratios,",
      "attributable fractions and number needed from one DP 2x2 table."),
    c("mechanism_uncertainty", "sampling_uncertainty_if_claimed",
      "validated_synopsis_provenance"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPEpi2x2Inference", "ds.vertDPEpi2x2Inference",
    "dp_epi_2x2_sampling_inference", "validated_synopsis_2x2_artifact",
    paste(
      "Conservative joint confidence regions for group/population risks,",
      "risk difference, risk/odds ratios, attributable fractions and number",
      "needed, combining signed DP-mechanism and binomial sampling uncertainty."),
    c("binomial_sampling_model", "clopper_pearson_exact_intervals",
      "joint_mechanism_and_sampling_uncertainty",
      "validated_synopsis_provenance"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPPrevalenceRatio", "ds.vertDPPrevalenceRatio",
    "dp_prevalence_ratio", "validated_synopsis_2x2_artifact",
    paste(
      "Caller-declared cross-sectional exposed/unexposed prevalences,",
      "prevalence difference, prevalence/odds ratios, attributable",
      "prevalence fractions and number needed, as exact aliases of one",
      "validated DP epidemiology result."),
    c("cross_sectional_design_declared_not_inferred",
      "explicit_exposed_and_prevalent_orientation",
      "mechanism_uncertainty", "sampling_uncertainty_if_claimed",
      "validated_synopsis_provenance", "zero_call_numeric_identity"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPPrevalenceRatioInference",
    "ds.vertDPPrevalenceRatioInference",
    "dp_prevalence_ratio_sampling_inference",
    "validated_synopsis_2x2_artifact",
    paste(
      "Caller-declared cross-sectional prevalence effects with the exact",
      "same conservative joint DP-mechanism and binomial sampling regions",
      "as the validated DP 2x2 inference result."),
    c("binomial_sampling_model", "clopper_pearson_exact_intervals",
      "cross_sectional_design_declared_not_inferred",
      "explicit_exposed_and_prevalent_orientation",
      "joint_mechanism_and_sampling_uncertainty",
      "validated_synopsis_provenance", "zero_call_numeric_identity"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPMantelHaenszel", "ds.vertDPMantelHaenszel",
    "dp_mantel_haenszel", "validated_synopsis_strata_by_four_cells_artifact",
    paste(
      "Finite-snapshot common Mantel-Haenszel odds ratio with a conservative",
      "simultaneous mechanism region and no classical CMH p-value."),
    c("mechanism_uncertainty", "public_cell_mapping", "public_strata",
      "sampling_uncertainty_if_claimed", "validated_synopsis_provenance"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPDiagnostic2x2", "ds.vertDPDiagnostic2x2",
    "dp_diagnostic_2x2", "validated_synopsis_2x2_artifact",
    paste(
      "Sensitivity, specificity, predictive values, balanced accuracy, F1,",
      "likelihood ratios and diagnostic odds ratio from one DP 2x2 table."),
    c("diagnostic_axis_contract", "mechanism_uncertainty",
      "sampling_uncertainty_if_claimed", "validated_synopsis_provenance"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPDiagnostic2x2Inference",
    "ds.vertDPDiagnostic2x2Inference", "dp_diagnostic_2x2_sampling_inference",
    "validated_synopsis_2x2_artifact",
    paste(
      "Conservative joint confidence regions for diagnostic accuracy,",
      "combining signed DP-mechanism and exact-binomial sampling uncertainty."),
    c("binomial_sampling_model", "clopper_pearson_exact_intervals",
      "diagnostic_axis_contract",
      "joint_mechanism_and_sampling_uncertainty",
      "validated_synopsis_provenance"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPROC", "ds.vertDPROC", "dp_diagnostic_curve",
    "validated_synopsis_ordered_diagnostic_artifact",
    paste(
      "Threshold ROC curve and tie-adjusted finite-snapshot AUC from one",
      "ordered disease-status by score-bin DP table."),
    c("diagnostic_axis_contract", "fixed_public_score_order",
      "mechanism_uncertainty", "sampling_uncertainty_if_claimed",
      "validated_synopsis_provenance"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPDirectStandardization", "ds.vertDPDirectStandardization",
    "dp_direct_standardization", "validated_synopsis_strata_artifact",
    "Directly standardized risk from one DP strata-by-outcome table.",
    c("compatible_strata", "mechanism_uncertainty",
      "sampling_uncertainty_if_claimed", "validated_synopsis_provenance"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPDirectStandardizationInference",
    "ds.vertDPDirectStandardizationInference",
    "dp_direct_standardization_sampling_inference",
    "validated_synopsis_strata_artifact",
    paste(
      "Conservative joint confidence region for one directly standardized",
      "risk, combining signed DP-mechanism and exact-binomial stratum",
      "sampling uncertainty with fixed public weights."),
    c("binomial_sampling_model", "clopper_pearson_exact_intervals",
      "compatible_strata", "fixed_public_standard_weights",
      "joint_mechanism_and_sampling_uncertainty",
      "validated_synopsis_provenance"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPCausalStandardization", "ds.vertDPCausalStandardization",
    "dp_causal_stratified_standardization",
    "validated_synopsis_stratum_treatment_by_outcome_artifact",
    paste(
      "Saturated stratum-standardised treated/control risks and derived",
      "contrasts from one DP table and fixed public target weights."),
    c("binary_treatment_contract", "compatible_strata",
      "fixed_public_standard_weights", "causal_identification_assumptions",
      "mechanism_uncertainty", "validated_synopsis_provenance"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPCausalStandardizationInference",
    "ds.vertDPCausalStandardizationInference",
    "dp_causal_stratified_standardization_sampling_inference",
    "validated_synopsis_stratum_treatment_by_outcome_artifact",
    paste(
      "Conservative joint regions for saturated stratum-standardised",
      "treated/control risks and contrasts, combining signed DP-mechanism",
      "and exact-binomial sampling uncertainty."),
    c("binary_treatment_contract", "binomial_sampling_model",
      "clopper_pearson_exact_intervals", "compatible_strata",
      "fixed_public_standard_weights", "causal_identification_assumptions",
      "joint_mechanism_and_sampling_uncertainty",
      "validated_synopsis_provenance"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPIndirectStandardization", "ds.vertDPIndirectStandardization",
    "dp_indirect_standardization", "validated_synopsis_strata_artifact",
    "Observed-to-expected standardized ratio from one DP strata table.",
    c("compatible_populations", "mechanism_uncertainty",
      "sampling_uncertainty_if_claimed", "validated_synopsis_provenance"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  add(
    "ds.vertDPIndirectStandardizationInference",
    "ds.vertDPIndirectStandardizationInference",
    "dp_indirect_standardization_sampling_inference",
    "validated_synopsis_strata_artifact",
    paste(
      "Conservative joint confidence region for an observed-to-expected",
      "ratio, combining signed DP-mechanism and exact Poisson Garwood",
      "sampling uncertainty with fixed public expected rates."),
    c("compatible_populations", "fixed_public_expected_rates",
      "poisson_total_count_model", "garwood_exact_interval",
      "joint_mechanism_and_sampling_uncertainty",
      "validated_synopsis_provenance", "zero_call_postprocessing"),
    "synopsis_release_implemented",
    current_route_status =
      "client_only_validated_synopsis_postprocess",
    artifact_implementation_state = "validated_synopsis_adapter_implemented",
    inference_implementation_state = "synopsis_postprocess_implemented")

  out <- do.call(rbind, rows)
  out <- out[order(out$method, method = "radix"), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "schema_version") <- "dsvert-capsule-method-inventory-v4"
  attr(out, "alias_kind_levels") <- c(
    "canonical", "compatibility_alias", "compatibility_wrapper",
    "deprecated_subroute")
  attr(out, "current_route_status_levels") <- c(
    "client_only_inherits_input",
    "client_only_validated_capsule_postprocess",
    "client_only_validated_synopsis_postprocess",
    "formal_sticky_count_artifact",
    "formal_sticky_frequency_artifact",
    "formal_sticky_synopsis_artifact",
    "legacy_joint_dp_capsule_incompatible", "known_broken_route_quarantine",
    "legacy_exact_release_not_capsule_safe",
    "legacy_granular_release_not_capsule_safe",
    "legacy_granular_release_quarantine",
    "legacy_mutating_release_not_capsule_safe",
    "legacy_mutating_release_quarantine",
    "signed_workload_unavailable_quarantine",
    "formal_capsule_variant_only_legacy_unavailable",
    "formal_same_owner_synopsis_variant_only_legacy_unavailable")
  attr(out, "migration_feasibility_levels") <- c(
    "count_operation_implemented",
    "frequency_operation_implemented",
    "synopsis_release_implemented",
    "capsule_release_implemented",
    "client_only_requires_attested_input",
    "requires_new_capsule_artifact", "requires_new_secure_protocol")
  attr(out, "artifact_implementation_state_levels") <- c(
    "client_input_not_capsule_attested",
    "joint_vector_release_implemented",
    "signed_count_artifact_implemented",
    "signed_frequency_artifact_implemented",
    "signed_synopsis_describe_artifact_implemented",
    "planned_no_materializer", "reserved_not_materialized",
    "secure_artifact_not_implemented",
    "validated_capsule_adapter_implemented",
    "validated_frequency_artifact_adapter_implemented",
    "validated_synopsis_adapter_implemented",
    "validated_complete_case_gaussian_synopsis_adapter_implemented",
    "validated_same_owner_synopsis_adapter_implemented",
    "validated_complete_case_gaussian_capsule_adapter_implemented",
    "validated_same_and_cross_owner_capsule_adapter_implemented",
    "validated_same_owner_capsule_adapter_implemented")
  attr(out, "inference_implementation_state_levels") <- c(
    "dp_gaussian_lasso_path_implemented",
    "dp_gaussian_pseudo_ic_implemented",
    "dp_aware_conditional_hypergeometric_bootstrap_implemented",
    "dp_aware_conditional_null_not_implemented",
    "dp_aware_null_distribution_not_implemented",
    "dp_aware_parametric_bootstrap_implemented",
    "existing_inference_requires_capsule_backend",
    "legacy_capsule_release_incompatible",
    "formal_count_release_implemented",
    "formal_frequency_release_implemented",
    "formal_synopsis_describe_release_implemented",
    "frequency_postprocess_implemented",
    "implemented_client_algebra_inherits_input",
    "capsule_postprocess_implemented",
    "synopsis_postprocess_implemented",
    "client_only_spectral_postprocess_with_eigengap_regions",
    "complete_case_postprocess_and_explicit_psd_projection_implemented",
    "pairwise_postprocess_and_explicit_psd_projection_implemented",
    "legacy_inference_requires_secure_redesign",
    "signed_binomial_lasso_design_gram_not_materialized",
    "signed_multinomial_design_gram_not_materialized",
    "multiple_imputation_protocol_and_inference_require_redesign",
    "nb2_joint_beta_theta_inference_not_implemented",
    "propensity_to_weight_estimand_pipeline_not_implemented",
    "target_capsule_release_not_implemented")
  attr(out, "remote_release_class_levels") <- c(
    "client_reconstructed_exact_statistic",
    "opaque_peer_ciphertext", "plaintext_exact_aggregate",
    "public_metadata", "share_reconstructed_by_client")
  attr(out, "scope") <- paste(
    "All public inferential methods. Administrative, identity, alignment,",
    "calibration and status exports are tracked separately by",
    ".dsvert_capsule_non_inference_exports().")
  out
}
