#' @title Audited maturity status of dsVert client methods
#' @description Returns the conservative method-maturity registry used to keep
#'   experimental or compatibility routes out of the promoted path. The status
#'   reflects both statistical validity and disclosure behavior under a pinned,
#'   honest-but-curious peer model with an untrusted analyst/relay.
#'
#' @param method Optional character vector of exported method names.
#' @param status Optional subset of `"promoted"`, `"provisional"`,
#'   `"compatibility"`, or `"quarantine"`.
#'
#' @return A data frame with one row per public method, its canonical route,
#'   maturity status, release contract, defensible scope, and principal
#'   limitation. This registry is an auditable release contract, not a
#'   cryptographic certification.
#' @export
ds.vertMethodStatus <- function(method = NULL, status = NULL) {
  registry <- .dsvert_method_registry()
  if (!is.null(method)) {
    if (!is.character(method) || anyNA(method) || any(!nzchar(method))) {
      stop("method must be NULL or non-empty method names", call. = FALSE)
    }
    unknown <- setdiff(method, registry$method)
    if (length(unknown)) {
      stop("Unknown method(s): ", paste(unknown, collapse = ", "),
           call. = FALSE)
    }
    registry <- registry[registry$method %in% method, , drop = FALSE]
  }
  if (!is.null(status)) {
    allowed <- c("promoted", "provisional", "compatibility", "quarantine")
    if (!is.character(status) || anyNA(status) ||
        any(!status %in% allowed)) {
      stop("status must use promoted, provisional, compatibility, or quarantine",
           call. = FALSE)
    }
    registry <- registry[registry$status %in% status, , drop = FALSE]
  }
  rownames(registry) <- NULL
  class(registry) <- c("ds.vertMethodStatus", class(registry))
  attr(registry, "audit_date") <- "2026-08-24"
  attr(registry, "threat_model") <- paste(
    "pinned honest-but-curious peers; untrusted analyst/relay;",
    "privacy and composition follow each row's release contract; Count is",
    "sticky per canonical signed artifact, distinct Count artifacts compose,",
    "and no finite global Count composition claim is made; capsule-backed",
    "routes retain their separately documented capsule contract; Frequency",
    "is independently sticky per canonical signed artifact, distinct",
    "Frequency artifacts compose, and no finite global Frequency composition",
    "claim is made; Synopsis releases are independently sticky per canonical",
    "signed artifact, distinct Synopsis artifacts compose, and no finite",
    "global Synopsis composition claim is made;",
    "no malicious-peer security or unlimited exact-output non-reconstruction claim"
  )
  registry
}

.dsvert_method_registry <- function() {
  rows <- list()
  add <- function(methods, canonical, status, scope, limitation) {
    rows[[length(rows) + 1L]] <<- data.frame(
      method = methods,
      canonical = canonical,
      status = status,
      safe_scope = scope,
      principal_limitation = limitation,
      stringsAsFactors = FALSE
    )
  }

  add("ds.getIdentityPks", "ds.getIdentityPks", "promoted",
      "Public Ed25519 identity discovery.",
      "Deployment must still test rotation, restart persistence and pin mismatch.")
  add(c("ds.psiAlign", "ds.vert.align"), "ds.psiAlign", "promoted",
      paste(
        "Pinned ECDH-PSI with a server-owned fixed-capacity bucket, signed",
        "contract/receipts, encrypted fixed-shape relay and exact GC/OT",
        "all-server membership; it releases no statistic or cardinality."),
      paste(
        "The public bucket is an input upper bound; availability attacks",
        "remain possible, and collusion of both designated compute peers is",
        "outside the membership-confidentiality claim."))
  add(c("ds.isPsiAligned", "ds.vert.is_aligned"), "ds.isPsiAligned",
      "promoted",
      paste(
        "Validates the persistent HMAC-bound padded-PSI manifest locally and",
        "accepts only one identical fixed-schema attestation from all peers."),
      paste(
        "It intentionally returns no count; externally pre-aligned data",
        "without the authenticated manifest fail closed."))
  add(c("ds.vertDesc", "ds.vert.desc"), "ds.vertDesc", "promoted",
      paste(
        "Compatibility data frame over one custodian-owned sticky DP describe",
        "artifact with bounded moments and fixed-grid quantiles."),
      paste(
        "It requires an explicit analysis_id; observed extrema and adaptive",
        "ranges/bins are unavailable, and mechanism regions exclude sampling",
        "uncertainty."))
  add("ds.vertDPStatus", "ds.vertDPStatus", "promoted",
      paste(
        "Read-only signed Synopsis bootstrap for the immutable workload,",
        "exact pinned-peer set and deterministic sticky-release continuity."),
      paste(
        "No request, rate, or catalog limit is imposed; distinct artifacts",
        "compose and exact replay is deterministic. The status cannot remotely attest a",
        "malicious host, simultaneous rollback of",
        "every root/state copy, or all-peer collusion."))
  add("ds.vertDPCapsulePlan", "ds.vertDPCapsulePlan", "promoted",
      paste(
        "Data-free, signed dry-run of the server-authoritative Synopsis",
        "workload, coordinate layout, sensitivity and selected mechanism."),
      paste(
        "It creates no release and accesses no protected snapshot; utility",
        "still depends on the immutable custodian workload that it reports."))
  add("ds.vertDPCalibrate", "ds.vertDPCalibrate", "promoted",
      paste(
        "Data-free utility calibration for granular Laplace and the signed",
        "fixed-work, TV-accounted discrete-Gaussian mechanism."),
      paste(
        "It quantifies mechanism noise, not sampling error or a universally",
        "optimal epsilon; feasibility still depends on the signed mechanism",
        "shape and the declared implementation-delta allowance."))
  add("ds.vertDPCount", "ds.vertDPCount", "promoted",
      paste(
        "A canonical current-snapshot Count artifact signed by all K peers:",
        "two pinned authorities contribute independent sticky",
        "discrete-Laplace noise through the authenticated Ring128",
        "convolution path for add/remove adjacency, while fixed-cohort",
        "adjacency uses a signed public K-consensus value with zero",
        "sensitivity."),
      paste(
        "The formal (epsilon, delta)-DP claim is per canonical signed",
        "artifact; distinct artifacts compose and no finite global",
        "composition claim is made. It assumes published bounds, secret",
        "persistent seeds and at least one non-colluding pinned authority."))
  add("ds.vertDPContingency", "ds.vertDPContingency", "promoted",
      paste(
        "Fixed-domain, one-contribution-per-unit table from one canonical",
        "sticky Synopsis projection selected from signed same- or cross-owner",
        "categorical metadata."),
      paste(
        "ordinary chi-square/Fisher laws are invalid for the noisy table, and",
        "the concordant-unit estimand and cross-owner threat boundary must be",
        "reported."))
  add("ds.vertDPFrequency", "ds.vertDPFrequency", "promoted",
      paste(
        "One canonical signed fixed-domain categorical vector from the",
        "explicit PSI-aligned source owner; K peers compile and sign, while",
        "the source and one pinned secondary authority execute one sticky",
        "Ring128 release with simultaneous mechanism-only regions."),
      paste(
        "The formal DP claim is per canonical signed artifact; distinct",
        "artifacts compose and no finite global composition claim is made.",
        "It assumes persistent secret seeds and one non-colluding authority;",
        "base regions exclude population-sampling uncertainty."))
  add("ds.vertDPFrequencyInference", "ds.vertDPFrequencyInference", "promoted",
      paste(
        "Zero-call conservative joint regions combining the validated",
        "frequency release's simultaneous DP event with exact multinomial",
        "level-wise binomial sampling intervals."),
      paste(
        "Requires iid privacy units and scientifically ignorable exclusions;",
        "Bonferroni/Clopper-Pearson regions can be wide and provide no p-value."))
  add("ds.vertDPMeanVar", "ds.vertDPMeanVar", "promoted",
      paste(
        "Bounded per-unit count, normalized sum and sum-of-squares from the",
        "one canonical sticky Synopsis artifact with consistent",
        "post-processing and unlimited replay."),
      paste(
        "Clipping, quantisation and privacy noise change finite-sample",
        "utility; the custodian bounds must be scientifically defensible."))
  add("ds.vertDPDescribe", "ds.vertDPDescribe", "promoted",
      paste(
        "One custodian-specified sticky release of bounded quantized moments,",
        "fixed-grid histograms and quantile post-processing."),
      paste(
        "Moment/quantile regions cover mechanism noise and public",
        "quantisation/grid resolution only; sampling uncertainty is excluded."))
  add(c("ds.vertDPQuantile", "ds.vertDPMedian"),
      "ds.vertDPDescribe post-processing", "promoted",
      paste(
        "Release-only binned quantiles or medians from one validated sticky",
        "DP Synopsis Describe artifact, with zero additional DSI calls and",
        "privacy cost."),
      paste(
        "The estimand is identified only to a fixed public histogram bin;",
        "mechanism regions exclude sampling uncertainty and no exact sample",
        "quantile or within-bin interpolation is claimed."))
  add("ds.vertDPCor", "ds.vertDPCor", "promoted",
      paste(
        "Same-owner pairwise-complete bounded correlations from one sticky",
        "canonical Synopsis artifact, with unlimited exact replay and explicit",
        "PSD post-processing."),
      paste(
        "Cross-owner products remain reserved_not_materialized; mechanism",
        "regions exclude sampling uncertainty and PSD projection changes the",
        "raw pairwise estimand when the matrix is indefinite."))
  add("ds.vertDPSurvival", "ds.vertDPSurvival", "promoted",
      paste(
        "One canonical fixed-grid sticky Synopsis artifact for Kaplan-Meier,",
        "Nelson-Aalen, competing-risk post-processing and unlimited replay."),
      paste(
        "Current accuracy metadata covers DP histogram noise and public-grid",
        "discretisation separately; sampling confidence bands are not yet provided."))
  add(c("ds.vertDPKaplanMeier", "ds.vertDPNelsonAalen",
        "ds.vertDPCumulativeIncidence"),
      "ds.vertDPSurvival post-processing", "promoted",
      paste(
        "Pure curve post-processing with conservative simultaneous",
        "mechanism bands from one validated fixed-grid Synopsis artifact."),
      "No sampling confidence bands or hypothesis tests are provided.")
  add("ds.vertDPRMST", "ds.vertDPSurvival post-processing", "promoted",
      "Zero-cost fixed-grid RMST and simultaneous mechanism limits from one validated Synopsis survival artifact.",
      paste(
        "The limits exclude sampling uncertainty and continuous-time grid",
        "error; tau must lie within the public release bounds."))
  add("ds.vertDPRMTL", "ds.vertDPSurvival post-processing", "promoted",
      paste(
        "Exact zero-cost fixed-grid RMTL complement of RMST over the public",
        "Synopsis interval, with unchanged simultaneous mechanism coverage."),
      paste(
        "The limits exclude sampling uncertainty and continuous-time grid",
        "error; RMTL is tau-RMST only when the public lower bound is zero."))
  add("ds.vertDPSurvivalContrast", "ds.vertDPSurvival post-processing",
      "promoted",
      paste(
        "Zero-call fixed-grid survival differences and ratios from two",
        "validated Synopsis survival artifacts, with conservative joint",
        "mechanism-only regions and typed zero denominators."),
      paste(
        "Distinct artifacts retain only the Bonferroni joint-confidence",
        "lower bound; sampling uncertainty, hypothesis tests and",
        "continuous-time grid error are excluded."))
  add("ds.vertDPRMSTContrast", "ds.vertDPSurvival post-processing",
      "promoted",
      paste(
        "Zero-call fixed-grid RMST differences and ratios from two",
        "validated Synopsis survival artifacts through a common public tau."),
      paste(
        "Distinct artifacts retain only the Bonferroni joint-confidence",
        "lower bound; sampling uncertainty and continuous-time grid error",
        "are excluded, and denominator-zero ratios can be unbounded."))
  add(c("ds.vertDPSurvivalQuantile", "ds.vertDPMedianSurvival"),
      "ds.vertDPSurvival post-processing", "promoted",
      paste(
        "Zero-cost fixed-grid survival quantiles obtained by inverting one",
        "validated Synopsis Kaplan-Meier release and its mechanism band."),
      paste(
        "Targets not reached by the signed public horizon are explicitly",
        "typed as beyond-grid; limits exclude sampling uncertainty and",
        "continuous-time grid error."))
  add("ds.vertDPEpi2x2", "ds.vertDPEpi2x2", "promoted",
      "Client-only epidemiologic effects and simultaneous DP-mechanism regions from one authorized DP 2x2 release.",
      "The regions cover release noise, not sampling uncertainty; zero and near-zero denominators can make effects unbounded.")
  add("ds.vertDPEpi2x2Inference", "ds.vertDPEpi2x2Inference", "promoted",
      paste(
        "Client-only conservative joint confidence regions combining the",
        "signed DP count-box event with exact binomial sampling uncertainty",
        "from one authorized DP 2x2 release."),
      paste(
        "Requires iid privacy units under the declared joint",
        "exposure/outcome sampling model; Bonferroni/Clopper-Pearson regions",
        "can be wide or unbounded and provide no p-value."))
  add("ds.vertDPPrevalenceRatio", "ds.vertDPPrevalenceRatio", "promoted",
      paste(
        "Zero-call cross-sectional prevalence naming view of one validated",
        "DP 2-by-2 epidemiology release, with numerically identical points,",
        "mechanism regions, attributable fractions and number-needed output."),
      paste(
        "The caller must declare exposed and prevalent orientations; study",
        "design and temporality are not inferred from the table, and the",
        "regions exclude population-sampling uncertainty."))
  add("ds.vertDPPrevalenceRatioInference",
      "ds.vertDPPrevalenceRatioInference", "promoted",
      paste(
        "Zero-call cross-sectional prevalence naming view of the conservative",
        "joint DP-mechanism and exact-binomial sampling regions from one",
        "validated DP 2-by-2 release."),
      paste(
        "The caller-declared cross-sectional iid sampling model must be",
        "scientifically valid; design is not inferred from the table and",
        "regions can be wide or unbounded."))
  add("ds.vertDPMantelHaenszel", "ds.vertDPMantelHaenszel", "promoted",
      paste(
        "Zero-cost common Mantel-Haenszel odds ratio and conservative",
        "simultaneous mechanism region from one validated DP strata-by-cells",
        "table under the one-global-cell add/remove contribution contract."),
      paste(
        "The public cell/stratum contract must be correct; the region excludes",
        "sampling uncertainty, can be unbounded, and no classical CMH p-value",
        "is provided."))
  add("ds.vertDPDiagnostic2x2", "ds.vertDPDiagnostic2x2", "promoted",
      paste(
        "Client-only diagnostic-accuracy measures, including balanced accuracy",
        "and F1, with exhaustive-box-validated",
        "simultaneous DP-mechanism regions from one authorized DP 2x2 release."),
      paste(
        "Rows must be disease status and columns test result; regions exclude",
        "sampling uncertainty and can include infinite or non-estimable boundaries."))
  add("ds.vertDPDiagnostic2x2Inference",
      "ds.vertDPDiagnostic2x2Inference", "promoted",
      paste(
        "Client-only conservative joint confidence regions combining the",
        "signed DP count-box event with exact binomial sampling uncertainty",
        "for diagnostic accuracy from one authorized DP 2x2 release."),
      paste(
        "Requires iid privacy units under one declared joint disease/test",
        "population; Bonferroni/Clopper-Pearson regions can be wide or",
        "unbounded and provide no p-value."))
  add("ds.vertDPROC", "ds.vertDPROC", "promoted",
      paste(
        "Client-only threshold ROC curve and tie-adjusted finite-snapshot AUC",
        "with simultaneous DP-mechanism regions from one ordered DP table."),
      paste(
        "Score-bin order and direction must be public and explicit; the AUC",
        "region is conservative and excludes sampling uncertainty."))
  add("ds.vertDPIndirectStandardization",
      "ds.vertDPIndirectStandardization", "promoted",
      "Zero-cost observed-to-expected ratio and simultaneous mechanism region from one authorised DP strata table and public expected rates.",
      paste(
        "The region excludes sampling uncertainty; denominator-zero boxes",
        "can be non-estimable or unbounded and no Garwood interval is used."))
  add("ds.vertDPIndirectStandardizationInference",
      "ds.vertDPIndirectStandardizationInference", "promoted",
      paste(
        "Zero-call conservative joint region combining one signed DP count",
        "box with an exact Poisson Garwood family for the observed-to-expected",
        "ratio under fixed public expected rates."),
      paste(
        "Requires a valid Poisson total-count model and externally valid fixed",
        "expected rates; zero expected denominators are vacuous, and no",
        "p-value, causal effect or transportability claim is provided."))
  add("ds.vertDPDirectStandardization", "ds.vertDPDirectStandardization",
      "promoted",
      "Client-only directly standardized risk from one fixed-domain DP strata-by-outcome release.",
      "The interval covers release noise, not sampling uncertainty, and public reference weights must match the declared strata.")
  add("ds.vertDPDirectStandardizationInference",
      "ds.vertDPDirectStandardizationInference", "promoted",
      paste(
        "Client-only conservative joint confidence region combining the",
        "signed DP count-box event with exact binomial stratum sampling",
        "uncertainty and fixed public standard-population weights."),
      paste(
        "Requires conditionally binomial outcomes within public strata;",
        "Clopper-Pearson/Bonferroni limits can be wide and add no causal or",
        "population-transportability guarantee."))
  add("ds.vertDPCausalStandardization",
      "ds.vertDPCausalStandardization", "promoted",
      paste(
        "Zero-cost saturated stratum-standardised g-formula and simultaneous",
        "DP-mechanism regions from one validated stratum-treatment-by-outcome",
        "release with fixed public target weights."),
      paste(
        "A causal interpretation additionally requires consistency,",
        "conditional exchangeability, positivity, no interference, a correct",
        "public row mapping and scientifically valid target weights; no",
        "propensity score or doubly robust nuisance model is estimated."))
  add("ds.vertDPCausalStandardizationInference",
      "ds.vertDPCausalStandardizationInference", "promoted",
      paste(
        "Conservative joint DP-mechanism and exact-binomial sampling regions",
        "for the saturated public-stratum g-formula, with zero additional",
        "server calls or privacy cost."),
      paste(
        "All causal identification assumptions remain external scientific",
        "requirements; Bonferroni/Clopper-Pearson regions can be wide or",
        "unbounded and no p-value is provided."))
  add("ds.vertSecurityStatus", "ds.vertSecurityStatus", "promoted",
      paste(
        "Read-only consortium handshake for the unanimous release profile",
        "and coherent joint-capsule policy, pinset and privacy epoch."),
      "It validates reported contracts; it cannot remotely attest a malicious custodian implementation.")
  add("ds.vertNumericPreflight", "ds.vertNumericPreflight", "promoted",
      "Data-free fail-closed selection and certification against custodian-owned numeric bounds.",
      "A backend is usable only after its advertised exact primitives pass end-to-end attestation.")
  add(c("ds.vertCor", "ds.vert.cor"), "ds.vertCor", "promoted",
      paste(
        "Joint complete-case Pearson correlation from one signed same-owner",
        "Gaussian Synopsis artifact."),
      paste(
        "A signed analysis id and intercept are mandatory; mechanism regions",
        "exclude population-sampling uncertainty. Cross-owner descriptors",
        "are quarantined without a capsule fallback."))
  add(c("ds.vertPCA", "ds.vert.pca"), "ds.vertPCA", "promoted",
      paste(
        "Client-only eigenstructure of the explicitly PSD-projected sticky DP",
        "joint complete-case Gaussian correlation artifact."),
      paste(
        "Pairwise artifacts and individual scores are rejected; signs are",
        "arbitrary and directions without a mechanism-certified eigengap",
        "are not individually stable."))
  add(c("ds.vertChisq", "ds.vert.chisq"),
      "ds.vertChisq", "promoted",
      paste(
        "DP-aware parametric-bootstrap inference over one signed sticky",
        "same-owner contingency Synopsis artifact."),
      paste(
        "The plug-in null is asymptotic rather than finite-sample conditional;",
        "mechanism and Monte Carlo uncertainty are reported separately."))
  add(c("ds.vertFisher", "ds.vert.fisher"),
      "ds.vertFisher", "promoted",
      paste(
        "DP-aware conditional hypergeometric plug-in bootstrap over one",
        "signed sticky same-owner 2-by-2 contingency Synopsis artifact."),
      paste(
        "It is asymptotic rather than Fisher-exact for confidential data;",
        "Gaussian-mechanism calibration and a conditional odds-ratio",
        "confidence interval are not certified."))
  add(c("ds.vertChisqCross", "ds.vert.chisq_cross"), "ds.vertChisqCross",
      "promoted",
      paste(
        "One canonical signed cross-owner categorical Synopsis release, then",
        "client-only DP-aware Pearson/Yates and optional conditional",
        "hypergeometric plug-in post-processing."),
      paste(
        "The plug-in reference laws are asymptotic rather than finite-sample",
        "exact; mechanism and Monte Carlo uncertainty are reported",
        "separately, and no malicious-peer security is claimed."))
  add("ds.vertDPGaussian", "ds.vertDPGaussian", "promoted",
      paste(
        "Bounded, clipped complete-case Gaussian coefficients from one signed",
        "no-lifetime Synopsis sufficient-statistic artifact for same-owner",
        "predictors."),
      paste(
        "Binomial/Poisson links and population-sampling inference are not",
        "implemented on this route; singular released designs require an",
        "explicit ridge that changes the estimand. Cross-owner descriptors",
        "fail closed without capsule fallback."))
  add(c("ds.vertDPLMM", "ds.vertLMM", "ds.vert.lmm"),
      "ds.vertDPLMM", "promoted",
      paste(
        "Bounded random-intercept method-of-moments components from one",
        "signed same-owner sticky Synopsis artifact; the compatibility names",
        "require a signed analysis_id and matching cluster column. The",
        "artifact may define outcome ~ 1 or one canonical additive",
        "fixed-effect formula with its signed finite variance-ratio grid."),
      paste(
        "This is not ML/REML: random slopes, covariance estimates, standard",
        "errors, p-values and population-sampling inference are unavailable."))
  add("ds.validateDPGaussianCertificate",
      "ds.validateDPGaussianCertificate", "provisional",
      paste(
        "Client-only version-dispatched revalidation of the parallel Synopsis",
        "provenance certificate v1 and byte-compatible legacy capsule v3",
        "without DSI."),
      paste(
        "A self-contained certificate establishes internal integrity only;",
        "peer authenticity requires a caller/session-anchored trusted pinset,",
        "and it certifies the release mechanism rather than population-",
        "sampling validity."))
  add("ds.validateDPLMMCertificate", "ds.validateDPLMMCertificate",
      "provisional",
      "Client-only revalidation of a signed random-intercept LMM Synopsis certificate without DSI.",
      "It verifies artifact and release integrity, not ML/REML or population-sampling validity.")
  add(c("ds.vertGLM", "ds.vert.glm"), "ds.vertGLM", "promoted",
      paste(
        "Gaussian GLM uses an explicit dp_analysis_id selecting a bounded",
        "signed Synopsis artifact.",
        "Binomial/Poisson formal_analysis_id reads an already completed",
        "two-authority-signed public certificate or durable Phase21 terminal;",
        "fresh_formal_analysis_id may only resume the matching",
        "custodian-configured durable formal analysis.",
        "Same- and cross-owner",
        "artifacts must match the requested formula; the legacy estimator",
        "is unreachable."),
      paste(
        "Gaussian needs a matching complete-case artifact. Binomial and Poisson",
        "are completed releases, or the one configured fresh durable analysis;",
        "they cannot choose source work or a new DP opening. Neither route",
        "returns sampling inference or",
        "individual fitted values."))
  add("ds.vertCox", "ds.vertCox", "promoted",
      paste(
        "Read-only coefficients and hazard-ratio ranges from one completed",
        "two-authority-signed sticky formal Cox opening selected by",
        "formal_analysis_id, or the matching custodian-configured durable",
        "fresh analysis selected by fresh_formal_analysis_id."),
      paste(
        "The fresh route cannot choose source work or a new opening and",
        "returns no covariance,",
        "standard errors, p-values, baseline hazard, prediction, residual or",
        "sampling inference artifact."))
  add("ds.vert.coxph", "ds.vertCox", "promoted",
      paste(
        "Read-only coefficients and hazard-ratio ranges from one completed",
        "two-authority-signed sticky formal Cox opening selected by",
        "formal_analysis_id; no legacy Cox computation is reachable."),
      paste(
        "This compatibility route cannot start a Cox analysis and returns no",
        "covariance, standard errors, p-values, baseline hazard, prediction,",
        "residual or sampling inference artifact."))
  add(c("ds.vertCoxProfileNonDisclosive"),
      "ds.vertCox", "promoted",
      paste(
        "Read-only coefficients and hazard-ratio ranges from one completed",
        "two-authority-signed sticky formal Cox opening selected by",
        "formal_analysis_id; no legacy Cox computation is reachable."),
      paste(
        "This route cannot start a Cox analysis and returns no covariance,",
        "standard errors, p-values, baseline hazard, prediction, residual or",
        "sampling inference artifact."))
  add(c("ds.vertCoxDiscreteNonDisclosive"),
      "ds.vertCoxDiscreteNonDisclosive", "promoted",
      paste(
        "Read-only discrete-time pooled-logistic coefficients from one",
        "completed two-authority-signed binomial formal release selected by",
        "formal_analysis_id and bound to a fixed custodian time grid."),
      paste(
        "This is not a Cox proportional-hazards estimand and cannot start",
        "person-period expansion or a new DP release. Covariance, p-values,",
        "baseline hazards and predictions are unavailable."))
  add("ds.vert.cox", "ds.vertCox", "promoted",
      paste(
        "Alias for the completed formal Cox public certificate when",
        "formal_analysis_id is supplied, or the matching fresh durable Cox",
        "analysis when fresh_formal_analysis_id is supplied."),
      paste(
        "Fresh analysis is available only for method='profile'; no covariance,",
        "p-values or analyst-selected source work or opening is available."))
  add(c("ds.vertNBFullRegTheta", "ds.vert.nb"), "ds.vertNBFullRegTheta",
      "promoted",
      paste(
        "Read-only intercept-only NB2 method-of-moments fit from one",
        "validated sticky DP Frequency release; it performs no new DSI call",
        "or opening."),
      paste(
        "Only a bounded non-negative integer y ~ 1 domain is available.",
        "Covariates, likelihood optimization, covariance, standard errors",
        "and inference remain unavailable until a protected NB2",
        "score-and-information artifact exists."))
  add(c("ds.vertMultinom", "ds.vert.multinom", "ds.vertMultinomJoint",
        "ds.vertMultinomJointNewton"),
      "ds.vertMultinom", "promoted",
      paste(
        "Read-only intercept-only multinomial log-odds from one validated",
        "sticky DP Frequency release, including the historical joint names;",
        "it performs no new DSI call or opening."),
      paste(
        "Only y ~ 1 is available, with deterministic Jeffreys smoothing.",
        "Covariates, joint softmax, covariance, standard errors and inference",
        "remain unavailable until a purpose-bound protected score/design",
        "artifact exists."))
  add(c("ds.vertOrdinal", "ds.vert.ordinal", "ds.vertOrdinalJointNewton"),
      "ds.vertOrdinal", "promoted",
      paste(
        "Read-only intercept-only cumulative-logit thresholds from one",
        "validated sticky DP Frequency release, including the historical",
        "Newton name; it performs no new DSI call or opening."),
      paste(
        "The caller must give the complete signed category domain in clinical",
        "order. Covariates, protected score/information, covariance, standard",
        "errors and inference remain unavailable."))
  add("ds.vertLMM.k3", "ds.vertDPLMM", "promoted",
      "K>=3 compatibility wrapper for the signed random-intercept method-of-moments Synopsis.",
      "It accepts no K>=3 profile/REML, random slopes, covariance, standard errors or classical inference.")
  add(c("ds.vertGEE", "ds.vert.gee"), "ds.vertGEE", "promoted",
      paste(
        "A completed two-authority formal binomial/Poisson GLM certificate",
        "or the matching configured fresh durable GLM analysis can be read",
        "as an independence-working GEE point estimate; the",
        "score equation is the same and GEE adds no private computation",
        "beyond that configured GLM analysis."),
      paste(
        "Only corstr='independence' is available. Cluster identifiers,",
        "exchangeable/AR(1) correlations, sandwich covariance, standard",
        "errors and inference remain unavailable until a protected cluster",
        "score-and-meat artifact is implemented."))
  add(c("ds.vertGLMM", "ds.vert.glmm"), "ds.vertGLMM", "promoted",
      paste(
        "Binary outcome ~ 1 population-average log-odds and observed-scale",
        "ICC from one validated sticky random-intercept Synopsis; no new",
        "private computation occurs after the signed release."),
      paste(
        "Only a same-owner binary random-intercept moment approximation is",
        "available. It is not PQL, Laplace, ML or REML and has no covariates,",
        "random slopes, standard errors, p-values or sampling inference."))
  add(c("ds.vertIPW", "ds.vert.ipw"), "ds.vertIPW", "promoted",
      paste(
        "Exact binary intercept-only ATE IPW identity from one validated",
        "sticky treatment-by-outcome DP contingency artifact; it releases",
        "no individual weights or propensity fit."),
      paste(
        "Only outcome ~ treatment with treatment ~ 1 is available.",
        "Covariates, ATT/ATC, weighting, outcome regression, standard errors",
        "and p-values remain unavailable pending a purpose-bound causal",
        "artifact."))
  add(c("ds.vertMI", "ds.vert.mi"), "ds.vertMI", "promoted",
      paste(
        "Categorical MCAR completion for one or more intercept-only signed",
        "marginals from one strict-missing sticky Synopsis release. A",
        "three-or-more response call can instead use signed root-to-response",
        "pairs to build an explicit conditional star model; all draws are",
        "deterministic post-processing and add no privacy cost."),
      paste(
        "Only intercept-only categorical marginals, a two-response signed",
        "joint pair, or an opt-in star model with conditional independence",
        "given its first response are available. Covariates, chained",
        "equations, analyst seeds, standard errors, p-values, and Rubin",
        "sampling inference remain unavailable."))
  add(c("ds.vertLASSO", "ds.vert.lasso"), "ds.vertLASSO", "promoted",
      paste(
        "Gaussian L1 penalty path over one validated same-owner signed",
        "Synopsis moment artifact; each member has a deterministic KKT",
        "certificate and incurs no new private computation."),
      paste(
        "Gaussian same-owner designs only; no sampling inference or",
        "coefficient regions are available, and the retired ds.glm",
        "thresholding sketch is rejected."))
  add(c("ds.vertLASSO1Step", "ds.vert.lasso_1step"), "ds.vertLASSO1Step",
      "promoted",
      paste(
        "Historical one-step name for the Gaussian L1 penalty path over one",
        "validated same-owner signed Synopsis moment artifact; each member",
        "has a deterministic KKT certificate and incurs no new private",
        "computation."),
      paste(
        "Gaussian same-owner designs only; no sampling inference or",
        "coefficient regions are available, and the retired ds.glm",
        "quadratic-surrogate route is rejected."))
  add(c("ds.vertLASSOProximal", "ds.vert.lasso_proximal"),
      "ds.vertLASSOProximal", "promoted",
      paste(
        "Preferred zero-call Gaussian LASSO post-processing of one validated",
        "same-owner signed Synopsis moment artifact; each result carries a",
        "deterministic KKT certificate."),
      paste(
        "Gaussian only; no sampling inference or coefficient regions;",
        "cross-owner designs and unauthenticated legacy ds.glm inputs are",
        "unavailable."))
  add(c("ds.vertLASSOIter", "ds.vert.lasso_iter"), "ds.vertLASSOIter",
      "promoted",
      paste(
        "Gaussian L1 paths from one explicit signed same-owner Synopsis",
        "artifact; each path member carries a deterministic KKT certificate",
        "and creates no second private release."),
      paste(
        "Only family='gaussian' with dp_analysis_id is available. Binomial",
        "and Poisson fail closed: their score designs lack a purpose-bound",
        "signed Gram and their whole-path DP contracts are not implemented."))
  add(c("ds.vertLASSOCV", "ds.vert.lasso_cv"), "ds.vertLASSOCV",
      "promoted",
      paste(
        "Preferred zero-call DP-projected pseudo-AIC/BIC/EBIC selection over",
        "one validated same-owner signed Gaussian Synopsis LASSO path."),
      paste(
        "It is not cross-validation or a one-standard-error rule;",
        "selection uncertainty and sampling inference are unavailable;",
        "Gaussian same-owner designs only and cross-owner designs are",
        "unavailable; unauthenticated legacy ds.glm inputs are rejected."))
  add(c("ds.vertConfint", "ds.vert.confint"), "ds.vertConfint", "promoted",
      paste(
        "For a current signed same-owner Gaussian Synopsis,",
        "type='mechanism' returns a certificate-revalidated simultaneous",
        "mechanism region for the bounded complete-case ridge estimator",
        "without another server call or release."),
      paste(
        "This is not a population-sampling confidence interval and it can",
        "fail closed when the certified perturbation is too large relative",
        "to the released design. type='sampling' still requires an attested",
        "sampling covariance."))
  add(c("ds.vertWald", "ds.vert.wald", "ds.vertContrast", "ds.vert.contrast"),
      "ds.vert inference helpers", "promoted",
      paste(
        "For a current signed same-owner Gaussian Synopsis, type='mechanism'",
        "returns deterministic simultaneous DP-mechanism null regions for",
        "a coefficient or declared linear contrasts without another release."),
      paste(
        "These are not Wald statistics, p-values, sampling covariance or",
        "population inference. type='sampling' still requires an attested",
        "joint covariance artifact."))
  add(c("ds.vertLR", "ds.vert.lr"), "ds.vertLR", "promoted",
      paste(
        "For one current signed unpenalized Gaussian Synopsis,",
        "type='mechanism' returns a simultaneous outer region for the",
        "nested bounded normalized RSS reduction without another release."),
      paste(
        "This is not a likelihood-ratio statistic, chi-square reference,",
        "p-value or population inference. type='sampling' still requires",
        "two attested binomial/Poisson fits with canonical deviance."))
  add("ds.vertEpi2x2", "ds.vertEpi2x2", "promoted",
      "Client-only measures from an already-authorized 2x2 table.",
      "Arbitrary matrices carry no DataSHIELD provenance; small-sample Wald intervals are approximate.")
  add("ds.vertMantelHaenszel", "ds.vertMantelHaenszel", "promoted",
      paste(
        "Client-only common odds ratio and classical MH inference from",
        "already-authorized stratified 2x2 aggregate tables."),
      paste(
        "Bare tables carry caller-attested rather than cryptographically",
        "verifiable provenance; inference is asymptotic and conditional on",
        "valid strata, independent units and a common-odds-ratio model."))
  add("ds.vertDirectStandardization", "ds.vertDirectStandardization", "promoted",
      "Client-only direct rates from authorized stratum aggregates.",
      "Correctness depends on compatible strata and reference weights.")
  add("ds.vertIndirectStandardization", "ds.vertIndirectStandardization", "promoted",
      "Client-only SMR/SIR and Garwood interval from authorized aggregates.",
      "Correctness depends on valid expected counts and population definitions.")
  add("ds.vertMethodStatus", "ds.vertMethodStatus", "promoted",
      "Read-only audited release registry.",
      "Status is conservative release metadata, not formal certification.")

  out <- do.call(rbind, rows)
  out$release_contract <- "legacy_exact_release_not_capsule_safe"
  synopsis_releases <- c(
    "ds.vertDesc", "ds.vert.desc", "ds.vertDPDescribe",
    "ds.vertDPMeanVar", "ds.vertDPCor", "ds.vertDPSurvival",
    "ds.vertDPGaussian", "ds.vertDPLMM", "ds.vertLMM", "ds.vert.lmm",
    "ds.vertLMM.k3", "ds.vertGLMM", "ds.vert.glmm",
    "ds.vertMI", "ds.vert.mi",
    "ds.vertCor", "ds.vert.cor",
    "ds.vertPCA", "ds.vert.pca",
    "ds.vertDPContingency", "ds.vertChisq", "ds.vert.chisq",
    "ds.vertFisher", "ds.vert.fisher",
    "ds.vertChisqCross", "ds.vert.chisq_cross")
  safe_nonreleases <- c(
    "ds.getIdentityPks",
    "ds.psiAlign", "ds.vert.align",
    "ds.isPsiAligned", "ds.vert.is_aligned", "ds.vertDPStatus",
    "ds.vertDPCapsulePlan",
    "ds.vertDPCalibrate", "ds.validateDPGaussianCertificate",
    "ds.validateDPLMMCertificate",
    "ds.vertSecurityStatus", "ds.vertNumericPreflight",
    "ds.vertMethodStatus")
  inherited_postprocessing <- c(
    "ds.vertDPQuantile", "ds.vertDPMedian",
    "ds.vertDPFrequencyInference",
    "ds.vertDPEpi2x2", "ds.vertDPEpi2x2Inference",
    "ds.vertDPPrevalenceRatio", "ds.vertDPPrevalenceRatioInference",
    "ds.vertDPMantelHaenszel", "ds.vertDPDiagnostic2x2",
    "ds.vertDPDiagnostic2x2Inference", "ds.vertDPROC",
    "ds.vertDPDirectStandardization",
    "ds.vertDPDirectStandardizationInference",
    "ds.vertDPCausalStandardization",
    "ds.vertDPCausalStandardizationInference",
    "ds.vertIPW", "ds.vert.ipw",
    "ds.vertDPIndirectStandardization",
    "ds.vertDPIndirectStandardizationInference",
    "ds.vertDPKaplanMeier", "ds.vertDPNelsonAalen",
    "ds.vertDPCumulativeIncidence", "ds.vertDPRMST", "ds.vertDPRMTL",
    "ds.vertDPSurvivalContrast", "ds.vertDPRMSTContrast",
    "ds.vertDPSurvivalQuantile", "ds.vertDPMedianSurvival",
    "ds.vertConfint", "ds.vert.confint", "ds.vertWald", "ds.vert.wald",
    "ds.vertContrast", "ds.vert.contrast", "ds.vertEpi2x2",
    "ds.vertMantelHaenszel",
    "ds.vertDirectStandardization", "ds.vertIndirectStandardization",
    "ds.vertLR", "ds.vert.lr",
    "ds.vertLASSO", "ds.vert.lasso",
    "ds.vertLASSO1Step", "ds.vert.lasso_1step",
    "ds.vertLASSOProximal", "ds.vert.lasso_proximal",
    "ds.vertLASSOCV", "ds.vert.lasso_cv")
  inherited_postprocessing <- c(
    inherited_postprocessing, "ds.vertMultinom", "ds.vert.multinom",
    "ds.vertMultinomJoint", "ds.vertMultinomJointNewton",
    "ds.vertOrdinal", "ds.vert.ordinal", "ds.vertOrdinalJointNewton",
    "ds.vertNBFullRegTheta", "ds.vert.nb")
  glm_synopsis_wrappers <- c(
    "ds.vertGLM", "ds.vert.glm",
    "ds.vertLASSOIter", "ds.vert.lasso_iter")
  out$release_contract[out$method %in% synopsis_releases] <-
    "formal_sticky_synopsis_artifact"
  out$release_contract[out$method == "ds.vertDPCount"] <-
    "formal_sticky_count_artifact"
  out$release_contract[out$method == "ds.vertDPFrequency"] <-
    "formal_sticky_frequency_artifact"
  out$release_contract[out$method %in% safe_nonreleases] <-
    "disclosure_safe_protocol_no_statistic"
  out$release_contract[out$method %in% inherited_postprocessing] <-
    "postprocessing_inherits_input"
  out$release_contract[out$method %in% glm_synopsis_wrappers] <-
    "formal_sticky_synopsis_artifact"
  formal_cox_public <- c(
    "ds.vertCox", "ds.vert.cox", "ds.vert.coxph",
    "ds.vertCoxProfileNonDisclosive",
    "ds.vertCoxDiscreteNonDisclosive",
    "ds.vertGEE", "ds.vert.gee")
  out$release_contract[out$method %in% formal_cox_public] <-
    "formal_completed_public_certificate"

  # Numeric certification is a separate axis from statistical maturity and
  # disclosure release status.  Keep it explicit so a newly promoted method
  # cannot inherit the GLM certificate merely because it reuses part of a GLM
  # session or object.
  out$numeric_contract <- "unattested_result"
  out$numeric_blocker <- paste(
    "No operation-complete workload adapter, intermediate-bound attestation,",
    "and end-to-end estimator error certificate.")
  out$may_report_numerically_certified <- FALSE
  out$currently_numerically_certified <- FALSE

  protocol_only <-
    out$release_contract == "disclosure_safe_protocol_no_statistic" |
    out$canonical == "ds.psiAlign"
  out$numeric_contract[protocol_only] <- "not_applicable_no_statistic"
  out$numeric_blocker[protocol_only] <-
    "This endpoint does not return a statistical numeric result."

  dp_release <- out$release_contract %in% c(
    "formal_joint_dp_capsule",
    "formal_joint_dp_capsule_only_legacy_unavailable",
    "formal_sticky_count_artifact",
    "formal_sticky_frequency_artifact",
    "formal_sticky_synopsis_artifact",
    "formal_completed_public_certificate")
  out$numeric_contract[dp_release] <- "separate_integer_dp_contract"
  out$numeric_blocker[dp_release] <- paste(
    "Covered by the DP mechanism/accountant contract, not the Ring numeric",
    "certificate.")
  count_release <- out$method == "ds.vertDPCount"
  out$numeric_blocker[count_release] <- paste(
    "Covered by the signed per-artifact integer DP or exact public",
    "K-consensus contract, not an accountant or Ring numeric certificate.")
  frequency_release <- out$method == "ds.vertDPFrequency"
  out$numeric_blocker[frequency_release] <- paste(
    "Covered by the signed per-artifact integer DP contract, not an",
    "accountant or Ring numeric certificate.")
  synopsis_release <- out$method %in% synopsis_releases
  out$numeric_blocker[synopsis_release] <- paste(
    "Covered by the signed per-artifact integer DP contract, not the Ring",
    "numeric certificate.")

  inherited <- out$release_contract == "postprocessing_inherits_input"
  out$numeric_contract[inherited] <- "inherits_input_contract"
  out$numeric_blocker[inherited] <- paste(
    "Client post-processing cannot strengthen the numeric or disclosure",
    "contract of its supplied aggregate/fit.")

  preflight <- out$method == "ds.vertNumericPreflight"
  out$numeric_contract[preflight] <- "data_free_preflight_only"
  out$numeric_blocker[preflight] <- paste(
    "Preflight eligibility is not execution certification and returns no",
    "estimator.")

  out <- out[c(
    "method", "canonical", "status", "release_contract",
    "numeric_contract", "may_report_numerically_certified",
    "currently_numerically_certified", "safe_scope",
    "principal_limitation", "numeric_blocker")]
  out <- out[order(out$method), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @export
print.ds.vertMethodStatus <- function(x, ...) {
  NextMethod("print", x, row.names = FALSE, ...)
  invisible(x)
}
