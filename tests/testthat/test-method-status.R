test_that("method maturity registry covers every public analysis entry point", {
  registry <- ds.vertMethodStatus()
  expect_s3_class(registry, "ds.vertMethodStatus")
  expect_identical(anyDuplicated(registry$method), 0L)
  expect_true(all(registry$status %in%
                    c("promoted", "provisional", "compatibility", "quarantine")))
  expect_true(all(registry$release_contract %in% c(
    "formal_joint_dp_capsule",
    "formal_joint_dp_capsule_only_legacy_unavailable",
    "formal_sticky_count_artifact",
    "formal_sticky_frequency_artifact",
    "formal_sticky_synopsis_artifact",
    "formal_completed_public_certificate",
    "disclosure_safe_protocol_no_statistic",
    "postprocessing_inherits_input",
    "legacy_exact_release_not_capsule_safe")))
  expect_false(any(grepl("strict|budget|exhaust|remaining",
                         registry$release_contract, ignore.case = TRUE)))
  expect_match(attr(registry, "threat_model"),
               "distinct Count artifacts compose", fixed = TRUE)
  expect_match(attr(registry, "threat_model"),
               "distinct Synopsis artifacts compose", fixed = TRUE)
  expect_true(all(registry$numeric_contract %in% c(
    "not_applicable_no_statistic", "separate_integer_dp_contract",
    "inherits_input_contract", "data_free_preflight_only",
    "glm_preflight_runtime_incomplete", "input_guard_or_inherited_only",
    "unattested_result")))
  expect_false(any(registry$currently_numerically_certified))
  expect_length(
    registry$method[registry$may_report_numerically_certified], 0L)

  namespace <- readLines(system.file("NAMESPACE", package = "dsVertClient"),
                         warn = FALSE)
  exports <- sub("^export\\((.*)\\)$", "\\1",
                 grep("^export\\(", namespace, value = TRUE))
  public <- exports[grepl(paste0(
    "^(ds[.]vert|ds[.]psiAlign$|ds[.]isPsiAligned$|",
    "ds[.]getIdentityPks$|ds[.]validateDP(Gaussian|LMM)Certificate$)"),
                          exports)]
  expect_setequal(registry$method, public)
})

test_that("no public route may report a result numeric certificate", {
  registry <- ds.vertMethodStatus()
  non_reporters <- registry$method[
    !registry$may_report_numerically_certified &
      registry$release_contract !=
        "disclosure_safe_protocol_no_statistic"]
  namespace <- asNamespace("dsVertClient")
  offenders <- vapply(non_reporters, function(method) {
    if (!exists(method, envir = namespace, mode = "function",
                inherits = FALSE)) {
      return(FALSE)
    }
    body_text <- paste(deparse(body(get(method, envir = namespace))),
                       collapse = "\n")
    grepl("numerically_certified", body_text, fixed = TRUE)
  }, logical(1L))
  expect_false(any(offenders),
               info = paste(names(offenders)[offenders], collapse = ", "))
})

test_that("readmitted MI and existing Synopsis routes are presented honestly", {
  quarantined <- ds.vertMethodStatus(status = "quarantine")$method
  readmitted_mi <- c("ds.vertMI", "ds.vert.mi")
  expect_false(any(readmitted_mi %in% quarantined))
  expect_true(all(ds.vertMethodStatus(readmitted_mi)$status == "promoted"))
  expect_true(all(ds.vertMethodStatus(readmitted_mi)$release_contract ==
                  "formal_sticky_synopsis_artifact"))
  expect_true(all(ds.vertMethodStatus(c("ds.vertGLMM", "ds.vert.glmm"))$
                  status == "promoted"))
  expect_true(all(ds.vertMethodStatus(c("ds.vertGLMM", "ds.vert.glmm"))$
                  release_contract == "formal_sticky_synopsis_artifact"))
  expect_true(all(ds.vertMethodStatus(c("ds.vertLMM", "ds.vert.lmm"))$
                  status == "promoted"))
  expect_match(
    ds.vertMethodStatus("ds.vertLMM")$safe_scope,
    "fixed-effect formula", fixed = TRUE)
  expect_true(all(ds.vertMethodStatus(c("ds.vertIPW", "ds.vert.ipw"))$
                  status == "promoted"))
  expect_true(all(ds.vertMethodStatus(c("ds.vertIPW", "ds.vert.ipw"))$
                  release_contract == "postprocessing_inherits_input"))
  expect_match(
    ds.vertMethodStatus("ds.vertIPW")$safe_scope,
    "one categorical saturated-stratum", fixed = TRUE)
  expect_true(all(ds.vertMethodStatus(c("ds.vertGEE", "ds.vert.gee"))$
                  status == "promoted"))
  expect_true(all(ds.vertMethodStatus(c("ds.vertGEE", "ds.vert.gee"))$
                  release_contract == "formal_completed_public_certificate"))
  frequency_multinom <- c(
    "ds.vertMultinom", "ds.vert.multinom", "ds.vertMultinomJoint",
    "ds.vertMultinomJointNewton")
  expect_true(all(ds.vertMethodStatus(frequency_multinom)$status == "promoted"))
  expect_true(all(ds.vertMethodStatus(frequency_multinom)$release_contract ==
      "postprocessing_inherits_input"))
  frequency_ordinal <- c(
    "ds.vertOrdinal", "ds.vert.ordinal", "ds.vertOrdinalJointNewton")
  expect_true(all(ds.vertMethodStatus(frequency_ordinal)$status == "promoted"))
  expect_true(all(ds.vertMethodStatus(frequency_ordinal)$release_contract ==
      "postprocessing_inherits_input"))
  expect_true(all(ds.vertMethodStatus(c(
    "ds.vertNBFullRegTheta", "ds.vert.nb"))$status == "promoted"))
  expect_true(all(ds.vertMethodStatus(c(
    "ds.vertNBFullRegTheta", "ds.vert.nb"))$release_contract ==
      "postprocessing_inherits_input"))
  expect_true(all(readmitted_mi %in%
                    ds.vertMethodStatus(status = "promoted")$method))
  psi <- ds.vertMethodStatus(c(
    "ds.psiAlign", "ds.vert.align",
    "ds.isPsiAligned", "ds.vert.is_aligned"))
  expect_true(all(psi$status == "promoted"))
  expect_true(all(psi$release_contract ==
                    "disclosure_safe_protocol_no_statistic"))
  expect_match(
    ds.vertMethodStatus("ds.psiAlign")$principal_limitation,
    "compute peers")
  expect_match(
    ds.vertMethodStatus("ds.isPsiAligned")$safe_scope,
    "HMAC-bound")
  expect_identical(ds.vertMethodStatus("ds.vertEpi2x2")$status, "promoted")
  expect_identical(
    ds.vertMethodStatus("ds.vertMantelHaenszel")$status, "promoted")
  expect_true(all(ds.vertMethodStatus(c(
    "ds.vertChisq", "ds.vert.chisq"))$release_contract ==
      "formal_sticky_synopsis_artifact"))
  expect_true(all(ds.vertMethodStatus(c(
    "ds.vertChisq", "ds.vert.chisq",
    "ds.vertFisher", "ds.vert.fisher"))$status == "promoted"))
  expect_true(all(ds.vertMethodStatus(c(
    "ds.vertFisher", "ds.vert.fisher"))$release_contract ==
      "formal_sticky_synopsis_artifact"))
  expect_match(ds.vertMethodStatus("ds.vertFisher")$principal_limitation,
               "not certified")
  expect_identical(
    ds.vertMethodStatus("ds.vertDPCount")$release_contract,
    "formal_sticky_count_artifact")
  expect_true(all(ds.vertMethodStatus(c(
    "ds.vertDesc", "ds.vert.desc"))$release_contract ==
      "formal_sticky_synopsis_artifact"))
  expect_identical(ds.vertMethodStatus("ds.vertDesc")$status, "promoted")
  expect_match(ds.vertMethodStatus("ds.vertDesc")$principal_limitation,
               "explicit analysis_id")
  count <- ds.vertMethodStatus("ds.vertDPCount")
  expect_identical(count$status, "promoted")
  expect_match(count$safe_scope, "Ring128 convolution", fixed = TRUE)
  expect_false(grepl("exact MPC", count$safe_scope, fixed = TRUE))
  frequency <- ds.vertMethodStatus("ds.vertDPFrequency")
  expect_identical(frequency$status, "promoted")
  expect_identical(frequency$release_contract,
                   "formal_sticky_frequency_artifact")
  expect_identical(frequency$numeric_contract,
                   "separate_integer_dp_contract")
  expect_match(frequency$numeric_blocker, "not an accountant", fixed = TRUE)
  expect_false(any(grepl("capsule", unlist(frequency),
                         ignore.case = TRUE)))
  expect_match(attr(frequency, "threat_model"),
               "distinct Frequency artifacts compose", fixed = TRUE)
  frequency_inference <-
    ds.vertMethodStatus("ds.vertDPFrequencyInference")
  expect_identical(frequency_inference$status, "promoted")
  expect_identical(frequency_inference$release_contract,
                   "postprocessing_inherits_input")
  expect_match(frequency_inference$safe_scope, "Zero-call")
  combined_epi <- ds.vertMethodStatus("ds.vertDPEpi2x2Inference")
  expect_identical(combined_epi$status, "promoted")
  expect_identical(
    combined_epi$release_contract, "postprocessing_inherits_input")
  expect_match(combined_epi$safe_scope, "sampling uncertainty")
  prevalence <- ds.vertMethodStatus(c(
    "ds.vertDPPrevalenceRatio",
    "ds.vertDPPrevalenceRatioInference"))
  expect_true(all(prevalence$status == "promoted"))
  expect_true(all(prevalence$release_contract ==
                    "postprocessing_inherits_input"))
  expect_true(all(grepl(
    "not inferred", prevalence$principal_limitation, fixed = TRUE)))
  combined_diagnostic <-
    ds.vertMethodStatus("ds.vertDPDiagnostic2x2Inference")
  expect_identical(combined_diagnostic$status, "promoted")
  expect_identical(
    combined_diagnostic$release_contract, "postprocessing_inherits_input")
  expect_match(combined_diagnostic$safe_scope, "sampling uncertainty")
  combined_standardization <-
    ds.vertMethodStatus("ds.vertDPDirectStandardizationInference")
  expect_identical(combined_standardization$status, "promoted")
  expect_identical(
    combined_standardization$release_contract,
    "postprocessing_inherits_input")
  expect_match(combined_standardization$safe_scope, "sampling uncertainty")
  expect_identical(
    ds.vertMethodStatus("ds.vertDPContingency")$status,
    "promoted")
  expect_identical(
    ds.vertMethodStatus("ds.vertDPMeanVar")$status,
    "promoted")
  expect_identical(
    ds.vertMethodStatus("ds.vertDPDescribe")$status,
    "promoted")
  expect_true(all(ds.vertMethodStatus(c(
    "ds.vertDPQuantile", "ds.vertDPMedian"))$status == "promoted"))
  expect_identical(
    ds.vertMethodStatus("ds.vertDPDescribe")$release_contract,
    "formal_sticky_synopsis_artifact")
  expect_identical(
    ds.vertMethodStatus("ds.vertDPSurvival")$status,
    "promoted")
  expect_identical(
    ds.vertMethodStatus("ds.vertDPCor")$status,
    "promoted")
  gaussian_synopsis <- ds.vertMethodStatus(c(
    "ds.vertDPGaussian", "ds.vertCor", "ds.vertPCA"))
  expect_true(all(gaussian_synopsis$status == "promoted"))
  expect_true(all(gaussian_synopsis$release_contract ==
                    "formal_sticky_synopsis_artifact"))
  synopsis <- ds.vertMethodStatus(c(
    "ds.vertDesc", "ds.vert.desc", "ds.vertDPDescribe",
    "ds.vertDPMeanVar", "ds.vertDPCor", "ds.vertDPSurvival",
    "ds.vertDPContingency", "ds.vertChisq", "ds.vert.chisq",
    "ds.vertFisher", "ds.vert.fisher",
    "ds.vertChisqCross", "ds.vert.chisq_cross"))
  expect_true(all(synopsis$release_contract ==
                    "formal_sticky_synopsis_artifact"))
  expect_false(any(grepl(
    "capsule|lifetime|accountant",
    apply(synopsis, 1L, paste, collapse = " "), ignore.case = TRUE)))
  survival_views <- ds.vertMethodStatus(c(
    "ds.vertDPKaplanMeier", "ds.vertDPNelsonAalen",
    "ds.vertDPCumulativeIncidence", "ds.vertDPRMST", "ds.vertDPRMTL",
    "ds.vertDPSurvivalContrast", "ds.vertDPRMSTContrast",
    "ds.vertDPLogRank",
    "ds.vertDPSurvivalQuantile", "ds.vertDPMedianSurvival"))
  expect_true(all(survival_views$release_contract ==
                    "postprocessing_inherits_input"))
  expect_true(all(grepl("Synopsis", survival_views$safe_scope, fixed = TRUE)))
  expect_false(any(grepl(
    "capsule|lifetime|accountant",
    apply(survival_views, 1L, paste, collapse = " "), ignore.case = TRUE)))
  one_release_survival_views <- ds.vertMethodStatus(c(
    "ds.vertDPKaplanMeier", "ds.vertDPNelsonAalen",
    "ds.vertDPCumulativeIncidence", "ds.vertDPRMST", "ds.vertDPRMTL",
    "ds.vertDPSurvivalQuantile", "ds.vertDPMedianSurvival"))
  expect_true(all(one_release_survival_views$status == "promoted"))
  expect_match(
    ds.vertMethodStatus("ds.vertDPGaussian")$principal_limitation,
    "sampling inference")
  validator <- ds.vertMethodStatus("ds.validateDPGaussianCertificate")
  expect_identical(validator$status, "provisional")
  expect_identical(
    validator$release_contract,
    "disclosure_safe_protocol_no_statistic")
  expect_match(validator$safe_scope, "without DSI")
  expect_match(validator$principal_limitation, "trusted pinset")
  mechanism_intervals <- ds.vertMethodStatus(c(
    "ds.vertConfint", "ds.vert.confint"))
  expect_true(all(mechanism_intervals$status == "promoted"))
  expect_true(all(grepl("mechanism", mechanism_intervals$safe_scope,
                        fixed = TRUE)))
  inference_helpers <- ds.vertMethodStatus(c(
    "ds.vertWald", "ds.vert.wald", "ds.vertContrast", "ds.vert.contrast",
    "ds.vertLR", "ds.vert.lr"))
  expect_true(all(inference_helpers$status == "promoted"))
  expect_true(all(grepl("not ", inference_helpers$principal_limitation,
                        fixed = TRUE)))
  expect_true(all(grepl("statistic", inference_helpers$principal_limitation,
                        fixed = TRUE)))
  lasso <- ds.vertMethodStatus(c(
    "ds.vertLASSOProximal", "ds.vert.lasso_proximal",
    "ds.vertLASSOCV", "ds.vert.lasso_cv"))
  proximal <- lasso$canonical == "ds.vertLASSOProximal"
  expect_true(all(lasso$status[proximal] == "promoted"))
  expect_true(all(lasso$status[!proximal] == "promoted"))
  expect_true(all(lasso$release_contract ==
                    "postprocessing_inherits_input"))
  expect_match(lasso$principal_limitation[proximal][[1L]],
               "cross-owner", fixed = TRUE)
  expect_true(all(grepl(
    "cross-validation", lasso$principal_limitation[!proximal], fixed = TRUE)))
  expect_match(
    ds.vertMethodStatus("ds.vertLASSOProximal")$safe_scope,
    "KKT", fixed = TRUE)
  expect_match(ds.vertMethodStatus("ds.vertLASSOCV")$safe_scope,
               "pseudo-AIC/BIC/EBIC")
  lasso_path <- ds.vertMethodStatus(c("ds.vertLASSO", "ds.vert.lasso"))
  expect_true(all(lasso_path$status == "promoted"))
  expect_match(lasso_path$safe_scope[[1L]], "KKT", fixed = TRUE)
  lasso_1step <- ds.vertMethodStatus(c(
    "ds.vertLASSO1Step", "ds.vert.lasso_1step"))
  expect_true(all(lasso_1step$status == "promoted"))
  expect_match(lasso_1step$safe_scope[[1L]], "KKT", fixed = TRUE)
  expect_true(all(c(lasso_path$release_contract,
                    lasso_1step$release_contract) ==
                    "postprocessing_inherits_input"))
  expect_true(all(c(lasso_path$numeric_contract,
                    lasso_1step$numeric_contract) ==
                    "inherits_input_contract"))
  expect_match(
    ds.vertMethodStatus("ds.vertDPCor")$principal_limitation,
    "Cross-owner")
  expect_match(
    ds.vertMethodStatus("ds.vertDPContingency")$principal_limitation,
    "threat boundary")
  expect_match(
    ds.vertMethodStatus("ds.vertDPCount")$principal_limitation,
    "per canonical signed")
  count_status <- ds.vertMethodStatus("ds.vertDPCount")
  count_public_contract <- paste(c(
    unlist(count_status), attr(count_status, "threat_model")), collapse = " ")
  expect_false(grepl(
    paste(
      "formal composition applies when each capsule",
      "operations over that capsule",
      "DP mechanism/accountant contract", sep = "|"),
    count_public_contract, ignore.case = TRUE))
  expect_match(attr(count_status, "threat_model"),
               "no finite global Count composition claim", fixed = TRUE)
  expect_match(
    ds.vertMethodStatus("ds.vertDPCount")$numeric_blocker,
    "not an accountant", fixed = TRUE)
  expect_identical(
    ds.vertMethodStatus("ds.vertGLM")$release_contract,
    "formal_sticky_synopsis_artifact")
  cox_public <- ds.vertMethodStatus(c(
    "ds.vertCox", "ds.vert.cox", "ds.vert.coxph",
    "ds.vertCoxProfileNonDisclosive"))
  expect_true(all(cox_public$status == "promoted"))
  expect_true(all(cox_public$release_contract ==
                    "formal_sticky_synopsis_artifact"))
  expect_true(all(grepl("analysis_id", cox_public$safe_scope,
                        fixed = TRUE)))
  expect_false(any(grepl("fresh_formal_analysis_id", cox_public$safe_scope,
                         fixed = TRUE)))
  expect_true(all(grepl("no covariance", cox_public$principal_limitation,
                        fixed = TRUE)))
  cox_discrete <- ds.vertMethodStatus("ds.vertCoxDiscreteNonDisclosive")
  expect_identical(cox_discrete$status, "promoted")
  expect_identical(cox_discrete$release_contract,
                   "formal_completed_public_certificate")
  expect_match(cox_discrete$principal_limitation,
               "not a Cox proportional-hazards estimand", fixed = TRUE)
  glm <- ds.vertMethodStatus(c("ds.vertGLM", "ds.vert.glm"))
  expect_true(all(glm$status == "promoted"))
  expect_match(glm$safe_scope[[1L]], "explicit dp_analysis_id", fixed = TRUE)
  expect_match(glm$principal_limitation[[1L]], "Binomial and Poisson",
               fixed = TRUE)
  expect_false(any(grepl("fresh_formal_analysis_id", glm$safe_scope,
                         fixed = TRUE)))
  expect_true(any(grepl("fresh computation is sealed",
                        glm$principal_limitation, fixed = TRUE)))
  expect_true(all(ds.vertMethodStatus(c(
    "ds.vertChisqCross", "ds.vert.chisq_cross"))$release_contract ==
      "formal_sticky_synopsis_artifact"))
  expect_true(all(ds.vertMethodStatus(c(
    "ds.vertChisqCross", "ds.vert.chisq_cross"))$status == "promoted"))
  expect_match(ds.vertMethodStatus("ds.vertChisqCross")$principal_limitation,
               "finite-sample")
  joint_frequency <- ds.vertMethodStatus(c(
    "ds.vertMultinomJoint", "ds.vertMultinomJointNewton",
    "ds.vertOrdinalJointNewton"))
  expect_true(all(joint_frequency$status == "promoted"))
  expect_true(all(grepl("intercept", joint_frequency$safe_scope,
                        ignore.case = TRUE)))
  lasso_iter <- ds.vertMethodStatus(c(
    "ds.vertLASSOIter", "ds.vert.lasso_iter"))
  expect_true(all(lasso_iter$status == "promoted"))
  expect_true(all(lasso_iter$release_contract ==
                    "formal_sticky_synopsis_artifact"))
  expect_match(lasso_iter$principal_limitation[[1L]],
               "binomial and Poisson", ignore.case = TRUE)
  expect_match(ds.vertMethodStatus("ds.vertMI")$principal_limitation,
               "intercept-only categorical", fixed = TRUE)
  expect_match(ds.vertMethodStatus("ds.vertNBFullRegTheta")$
                 principal_limitation,
               "finite beta/theta candidates", fixed = TRUE)
  expect_identical(
    ds.vertMethodStatus("ds.psiAlign")$release_contract,
    "disclosure_safe_protocol_no_statistic")
  expect_identical(
    ds.vertMethodStatus("ds.vertDPEpi2x2")$release_contract,
    "postprocessing_inherits_input")
  expect_identical(
    ds.vertMethodStatus("ds.vertDPMantelHaenszel")$release_contract,
    "postprocessing_inherits_input")
  expect_match(
    ds.vertMethodStatus("ds.vertDPMantelHaenszel")$principal_limitation,
    "no classical CMH p-value", fixed = TRUE)
  expect_identical(
    ds.vertMethodStatus("ds.vertDPDiagnostic2x2")$release_contract,
    "postprocessing_inherits_input")
  expect_identical(
    ds.vertMethodStatus("ds.vertDPROC")$release_contract,
    "postprocessing_inherits_input")
  expect_identical(
    ds.vertMethodStatus("ds.vertDPRMST")$release_contract,
    "postprocessing_inherits_input")
  expect_identical(
    ds.vertMethodStatus("ds.vertDPRMTL")$release_contract,
    "postprocessing_inherits_input")
  expect_match(
    ds.vertMethodStatus("ds.vertDPRMTL")$principal_limitation,
    "sampling uncertainty", fixed = TRUE)
  survival_contrasts <- ds.vertMethodStatus(c(
    "ds.vertDPSurvivalContrast", "ds.vertDPRMSTContrast",
    "ds.vertDPLogRank"))
  expect_true(all(survival_contrasts$status == "promoted"))
  expect_true(all(survival_contrasts$release_contract ==
                    "postprocessing_inherits_input"))
  expect_true(all(grepl(
    "Bonferroni", survival_contrasts$principal_limitation, fixed = TRUE)))
  survival_quantiles <- ds.vertMethodStatus(c(
    "ds.vertDPSurvivalQuantile", "ds.vertDPMedianSurvival"))
  expect_true(all(survival_quantiles$status == "promoted"))
  expect_true(all(survival_quantiles$release_contract ==
                    "postprocessing_inherits_input"))
  expect_true(all(grepl(
    "beyond-grid", survival_quantiles$principal_limitation, fixed = TRUE)))
  expect_true(all(ds.vertMethodStatus(c(
    "ds.vertLR", "ds.vert.lr"))$release_contract ==
      "postprocessing_inherits_input"))
  expect_match(ds.vertMethodStatus("ds.vertLR")$principal_limitation,
               "not a likelihood-ratio statistic", fixed = TRUE)
  indirect <- ds.vertMethodStatus(c(
    "ds.vertDPIndirectStandardization",
    "ds.vertDPIndirectStandardizationInference"))
  expect_true(all(indirect$status == "promoted"))
  expect_true(all(indirect$release_contract ==
                    "postprocessing_inherits_input"))
  expect_match(
    indirect$principal_limitation[[2L]], "Poisson", fixed = TRUE)
  causal <- ds.vertMethodStatus(c(
    "ds.vertDPCausalStandardization",
    "ds.vertDPCausalStandardizationInference"))
  expect_true(all(causal$status == "promoted"))
  expect_true(all(causal$release_contract ==
                    "postprocessing_inherits_input"))
  expect_true(all(grepl(
    "causal", causal$principal_limitation, ignore.case = TRUE)))
  dp_status <- ds.vertMethodStatus("ds.vertDPStatus")
  expect_match(dp_status$safe_scope, "signed Synopsis")
  expect_match(dp_status$principal_limitation,
               "No request, rate, or catalog limit")
  expect_match(dp_status$principal_limitation, "distinct artifacts compose")
  capsule_plan <- ds.vertMethodStatus("ds.vertDPCapsulePlan")
  expect_identical(
    capsule_plan$release_contract,
    "disclosure_safe_protocol_no_statistic")
  expect_match(capsule_plan$safe_scope, "Data-free")
  expect_error(ds.vertMethodStatus("not.a.method"), "Unknown method")
})
