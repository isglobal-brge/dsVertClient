# Package index

## Public frontdoor API

- [`ds.vert.align()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.is_aligned()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.desc()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.cor()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.pca()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.chisq()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.fisher()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.chisq_cross()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.glm()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.cox()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.coxph()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.nb()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.multinom()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.ordinal()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.lmm()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.gee()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.glmm()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.ipw()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.mi()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.lasso()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.lasso_iter()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.lasso_proximal()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.lasso_1step()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.lasso_cv()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.lr()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.confint()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.wald()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  [`ds.vert.contrast()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md)
  : Public ds.vert.\* API aliases

## PSI alignment and identity

- [`ds.psiAlign()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md)
  : Align vertically partitioned records with pinned, fixed-capacity PSI
- [`ds.isPsiAligned()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.isPsiAligned.md)
  : Check whether a data symbol is pinned-padded-PSI aligned
- [`ds.getIdentityPks()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.getIdentityPks.md)
  : Query Server Identity Public Keys

## DP policy, security and numeric contracts

- [`ds.vertMethodStatus()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMethodStatus.md)
  : Audited maturity status of dsVert client methods
- [`ds.vertSecurityStatus()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertSecurityStatus.md)
  : Verify the consortium disclosure profile
- [`ds.vertNumericPreflight()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertNumericPreflight.md)
  : Preflight a vertically partitioned GLM numeric backend
- [`ds.validateDPGaussianCertificate()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.validateDPGaussianCertificate.md)
  : Verify a bounded Gaussian capsule certificate without DSI

## Canonical sticky Count

- [`ds.vertDPCount()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPCount.md)
  [`print(`*`<ds.vertDPCount>`*`)`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPCount.md)
  : Differentially private privacy-unit count

## Sticky differentially private biomedical capsule

- [`ds.vertDPCalibrate()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPCalibrate.md)
  : Calibrate deployed DP mechanism candidates
- [`ds.vertDPCapsulePlan()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPCapsulePlan.md)
  : Dry-run one server-authoritative DP capsule
- [`ds.vertDPCausalStandardization()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPCausalStandardization.md)
  : DP stratified causal standardisation
- [`ds.vertDPCausalStandardizationInference()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPCausalStandardizationInference.md)
  : Joint DP and sampling inference for causal standardisation
- [`ds.vertDPContingency()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPContingency.md)
  [`print(`*`<ds.vertDPContingency>`*`)`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPContingency.md)
  : Differentially private fixed-domain contingency table
- [`ds.vertDPCor()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPCor.md)
  : Differentially private same-owner Pearson correlation
- [`ds.vertDPCumulativeIncidence()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPCumulativeIncidence.md)
  : Competing-Risks Cumulative Incidence from One DP Release
- [`ds.vertDPDescribe()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPDescribe.md)
  [`print(`*`<ds.vertDPDescribe>`*`)`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPDescribe.md)
  : Differentially Private Fixed-Grid Descriptive Statistics
- [`ds.vertDPDiagnostic2x2()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPDiagnostic2x2.md)
  : Diagnostic-accuracy measures with simultaneous DP-noise regions
- [`ds.vertDPDiagnostic2x2Inference()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPDiagnostic2x2Inference.md)
  : Diagnostic-accuracy inference with DP and sampling uncertainty
- [`ds.vertDPDirectStandardization()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPDirectStandardization.md)
  : Direct standardisation with simultaneous DP-noise bounds
- [`ds.vertDPDirectStandardizationInference()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPDirectStandardizationInference.md)
  : Direct-standardisation inference with DP and sampling uncertainty
- [`ds.vertDPEpi2x2()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPEpi2x2.md)
  : Epidemiological measures with simultaneous DP-noise regions
- [`ds.vertDPEpi2x2Inference()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPEpi2x2Inference.md)
  : Epidemiological 2-by-2 inference with DP and sampling uncertainty
- [`ds.vertDPFrequency()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPFrequency.md)
  : Differentially private one-way frequency distribution
- [`ds.vertDPFrequencyInference()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPFrequencyInference.md)
  : Conservative sampling regions for a DP frequency distribution
- [`ds.vertDPGaussian()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPGaussian.md)
  : Bounded Gaussian regression from the sticky DP capsule
- [`ds.vertDPIndirectStandardization()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPIndirectStandardization.md)
  : Differentially private indirect standardisation
- [`ds.vertDPIndirectStandardizationInference()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPIndirectStandardizationInference.md)
  : Indirect-standardisation inference with DP and sampling uncertainty
- [`ds.vertDPKaplanMeier()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPKaplanMeier.md)
  : Kaplan–Meier Curve from One DP Survival Release
- [`ds.vertDPMantelHaenszel()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPMantelHaenszel.md)
  : DP common Mantel-Haenszel odds ratio from one fixed capsule table
- [`ds.vertDPMeanVar()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPMeanVar.md)
  [`print(`*`<ds.vertDPMeanVar>`*`)`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPMeanVar.md)
  : Differentially private bounded mean and variance
- [`ds.vertDPMedian()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPMedian.md)
  : Binned DP median from one validated describe release
- [`ds.vertDPMedianSurvival()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPMedianSurvival.md)
  : Median survival from one DP survival release
- [`ds.vertDPNelsonAalen()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPNelsonAalen.md)
  : Nelson–Aalen Curve from One DP Survival Release
- [`ds.vertDPPrevalenceRatio()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPPrevalenceRatio.md)
  [`ds.vertDPPrevalenceRatioInference()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPPrevalenceRatio.md)
  : Cross-sectional prevalence effects from one DP 2-by-2 release
- [`ds.vertDPQuantile()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPQuantile.md)
  : Binned DP quantiles from one validated describe release
- [`ds.vertDPRMST()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPRMST.md)
  : Restricted mean survival time from one DP survival release
- [`ds.vertDPRMSTContrast()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPRMSTContrast.md)
  : Contrast restricted mean survival time from two DP releases
- [`ds.vertDPRMTL()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPRMTL.md)
  : Restricted mean time lost from one DP survival release
- [`ds.vertDPROC()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPROC.md)
  : Diagnostic ROC curve and AUC from one ordered DP table
- [`ds.vertDPStatus()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPStatus.md)
  : Differential-privacy policy status
- [`ds.vertDPSurvival()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPSurvival.md)
  [`print(`*`<ds.vertDPSurvival>`*`)`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPSurvival.md)
  : Differentially Private Non-Parametric Survival Curves
- [`ds.vertDPSurvivalContrast()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPSurvivalContrast.md)
  : Contrast two differentially private survival releases
- [`ds.vertDPSurvivalQuantile()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPSurvivalQuantile.md)
  : Fixed-grid survival quantiles from one DP survival release

## Epidemiology from authorised inputs

- [`ds.vertEpi2x2()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertEpi2x2.md)
  : Epidemiological effect measures from a protected 2-by-2 table
- [`ds.vertMantelHaenszel()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMantelHaenszel.md)
  : Common Mantel-Haenszel odds ratio from authorised stratified tables
- [`ds.vertDirectStandardization()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDirectStandardization.md)
  : Direct standardization of aggregate stratum-specific rates
- [`ds.vertIndirectStandardization()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertIndirectStandardization.md)
  : Indirect standardization from aggregate observed and expected events

## Conditional and compatibility analyses

- [`ds.vertDesc()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDesc.md)
  : Disclosure-safe descriptive statistics compatibility adapter
- [`ds.vertCor()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCor.md)
  : Disclosure-safe compatibility adapter for correlation
- [`ds.vertPCA()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertPCA.md)
  : Principal components from a signed DP correlation artifact
- [`ds.vertChisq()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertChisq.md)
  : DP-aware independence test on a 2-way contingency table
- [`ds.vertChisqCross()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertChisqCross.md)
  : DP-aware inference for a signed cross-owner categorical table
- [`ds.vertFisher()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertFisher.md)
  : DP-aware conditional test for a 2-by-2 contingency release
- [`ds.vertGLM()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGLM.md)
  : DP-capsule GLM compatibility frontdoor
- [`ds.vertLASSO()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSO.md)
  : Post-hoc soft-thresholded GLM coefficients (naive LASSO)
- [`ds.vertLASSOProximal()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSOProximal.md)
  : Gaussian LASSO via client-side coordinate descent
- [`ds.vertLASSO1Step()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSO1Step.md)
  : One-step LASSO via quadratic-surrogate proximal gradient
- [`ds.vertLASSOCV()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSOCV.md)
  : Client-side information-criterion selection for Gaussian LASSO
- [`ds.vertLR()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLR.md)
  : Likelihood-ratio test on two nested ds.vertGLM fits
- [`ds.vertConfint()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertConfint.md)
  : Wald confidence intervals for ds.vertGLM coefficients
- [`ds.vertWald()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertWald.md)
  : Univariate Wald test for a single ds.vertGLM coefficient
- [`ds.vertContrast()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertContrast.md)
  : Multi-coefficient Wald test via linear contrast K\*beta

## Quarantined compatibility frontdoors

Retained public names that fail locally before DSI; consult
ds.vertMethodStatus().

- [`ds.vertCox()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCox.md)
  : Quarantined Cox compatibility frontdoor
- [`ds.vertCoxProfileNonDisclosive()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCoxProfileNonDisclosive.md)
  : Quarantined Cox-profile compatibility frontdoor
- [`ds.vertCoxDiscreteNonDisclosive()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCoxDiscreteNonDisclosive.md)
  : Quarantined discrete-hazard compatibility frontdoor
- [`ds.vertNBFullRegTheta()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertNBFullRegTheta.md)
  : Quarantined negative-binomial compatibility frontdoor
- [`ds.vertMultinom()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMultinom.md)
  : Quarantined multinomial compatibility frontdoor
- [`ds.vertMultinomJoint()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMultinomJoint.md)
  : Quarantined joint-multinomial compatibility frontdoor
- [`ds.vertMultinomJointNewton()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMultinomJointNewton.md)
  : Quarantined joint-softmax Newton compatibility frontdoor
- [`ds.vertOrdinal()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertOrdinal.md)
  : Quarantined ordinal-regression compatibility frontdoor
- [`ds.vertOrdinalJointNewton()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertOrdinalJointNewton.md)
  : Quarantined proportional-odds Newton compatibility frontdoor
- [`ds.vertLMM()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLMM.md)
  : Quarantined linear-mixed-model compatibility frontdoor
- [`ds.vertLMM.k3()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLMM.k3.md)
  : Quarantined deprecated K\>=3 LMM frontdoor
- [`ds.vertGEE()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGEE.md)
  : Quarantined GEE compatibility frontdoor
- [`ds.vertGLMM()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGLMM.md)
  : Quarantined generalized-linear mixed-model frontdoor
- [`ds.vertIPW()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertIPW.md)
  : Quarantined IPW compatibility frontdoor
- [`ds.vertMI()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMI.md)
  : Quarantined multiple-imputation compatibility frontdoor
- [`ds.vertLASSOIter()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSOIter.md)
  : Quarantined iterative-LASSO compatibility frontdoor

## Helpers and S3 methods

- [`k2-beaver-lbfgs-client`](https://isglobal-brge.github.io/dsVertClient/reference/k2-beaver-lbfgs-client.md)
  : K=2 Beaver L-BFGS Pipeline
- [`k3-ring63-dcf-loop`](https://isglobal-brge.github.io/dsVertClient/reference/k3-ring63-dcf-loop.md)
  : K\>=3 Ring63 DCF + Beaver Gradient Loop
- [`chunk-utils`](https://isglobal-brge.github.io/dsVertClient/reference/chunk-utils.md)
  : Acknowledged Chunking for DataSHIELD Transport
- [`coef(`*`<ds.glm>`*`)`](https://isglobal-brge.github.io/dsVertClient/reference/coef.ds.glm.md)
  : Coefficients Method for ds.glm Objects
- [`ds.vertDPContingency()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPContingency.md)
  [`print(`*`<ds.vertDPContingency>`*`)`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPContingency.md)
  : Differentially private fixed-domain contingency table
- [`ds.vertDPCount()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPCount.md)
  [`print(`*`<ds.vertDPCount>`*`)`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPCount.md)
  : Differentially private privacy-unit count
- [`ds.vertDPDescribe()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPDescribe.md)
  [`print(`*`<ds.vertDPDescribe>`*`)`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPDescribe.md)
  : Differentially Private Fixed-Grid Descriptive Statistics
- [`ds.vertDPMeanVar()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPMeanVar.md)
  [`print(`*`<ds.vertDPMeanVar>`*`)`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPMeanVar.md)
  : Differentially private bounded mean and variance
- [`ds.vertDPSurvival()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPSurvival.md)
  [`print(`*`<ds.vertDPSurvival>`*`)`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPSurvival.md)
  : Differentially Private Non-Parametric Survival Curves
- [`print(`*`<ds.cor>`*`)`](https://isglobal-brge.github.io/dsVertClient/reference/print.ds.cor.md)
  : Print Method for ds.cor Objects
- [`print(`*`<ds.glm>`*`)`](https://isglobal-brge.github.io/dsVertClient/reference/print.ds.glm.md)
  : Print Method for ds.glm Objects
- [`print(`*`<ds.pca>`*`)`](https://isglobal-brge.github.io/dsVertClient/reference/print.ds.pca.md)
  : Print Method for ds.pca Objects
- [`summary(`*`<ds.glm>`*`)`](https://isglobal-brge.github.io/dsVertClient/reference/summary.ds.glm.md)
  : Summary Method for ds.glm Objects
