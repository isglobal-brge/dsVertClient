# Package index

## Record alignment

Privacy-preserving record alignment using ECDH-PSI

- [`ds.psiAlign()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md)
  : ECDH-PSI Record Alignment (Blind Relay)
- [`ds.isPsiAligned()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.isPsiAligned.md)
  : Check whether a data symbol is already PSI-aligned on every server
- [`ds.getIdentityPks()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.getIdentityPks.md)
  : Query Server Identity Public Keys

## Correlation, PCA, descriptive

Distributed second-order summaries

- [`ds.vertCor()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCor.md)
  : Ring63 Privacy-Preserving Correlation for Vertically Partitioned
  Data
- [`ds.vertPCA()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertPCA.md)
  : Principal Component Analysis for Vertically Partitioned Data
- [`ds.vertDesc()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDesc.md)
  : Federated descriptive statistics with approximate quantiles
- [`ds.vertChisq()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertChisq.md)
  : Federated chi-square test on a 2-way contingency table
- [`ds.vertChisqCross()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertChisqCross.md)
  : Cross-server chi-square (one-hot + Beaver cross-products)
- [`ds.vertFisher()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertFisher.md)
  : Federated Fisher exact test (same-server case)

## Generalised linear models (GLM)

K=2 / K=3 federated GLM fits via Ring63 Beaver MPC

- [`ds.vertGLM()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGLM.md)
  : Generalized Linear Model for Vertically Partitioned Data
- [`ds.vertConfint()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertConfint.md)
  : Wald confidence intervals for ds.vertGLM coefficients
- [`ds.vertContrast()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertContrast.md)
  : Multi-coefficient Wald test via linear contrast K\*beta
- [`ds.vertWald()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertWald.md)
  : Univariate Wald test for a single ds.vertGLM coefficient
- [`ds.vertLR()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLR.md)
  : Likelihood-ratio test on two nested ds.vertGLM fits

## Survival

Cox proportional-hazards (K=2 OT-Beaver) and K=3 Poisson trick

- [`ds.vertCox()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCox.md)
  : Federated Cox proportional-hazards regression
- [`ds.vertCox.k3()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCox.k3.md)
  : Federated Cox proportional-hazards via Allison-1982 Poisson trick
  (K=3)
- [`ds.vertCoxDiscreteNonDisclosive()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCoxDiscreteNonDisclosive.md)
  : Cox K=2 discrete-time pooled-logistic – non-disclosive option B
  (#D')

## Negative-binomial

NB regression with iid-mu / MoM / full-reg theta pivots

- [`ds.vertNB()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertNB.md)
  : Federated negative binomial regression with dispersion estimate
- [`ds.vertNBMoMTheta()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertNBMoMTheta.md)
  : NB regression with Method-of-Moments theta-estimator (K=2-safe)
- [`ds.vertNBFullRegTheta()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertNBFullRegTheta.md)
  : Federated NB regression with full-regression theta refinement

## Multinomial / ordinal

OVR + softmax-anchor (warm) and joint Newton (Tutz / Bohning)

- [`ds.vertMultinom()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMultinom.md)
  : Federated multinomial logistic regression via K-1 one-vs-rest fits
- [`ds.vertMultinomJoint()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMultinomJoint.md)
  : Federated joint-softmax multinomial logistic regression
- [`ds.vertMultinomJointNewton()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMultinomJointNewton.md)
  : Federated joint-softmax multinomial logistic regression via Ring127
  MPC-orchestrated Newton iteration
- [`ds.vertOrdinal()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertOrdinal.md)
  : Federated ordinal logistic regression (naive cumulative binomials)
- [`ds.vertOrdinalJointNewton()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertOrdinalJointNewton.md)
  : Federated joint proportional-odds ordinal regression via Ring127
  MPC-orchestrated Newton iteration

## Mixed models / longitudinal

REML LMM (K=2 / K=3), GEE sandwich, binomial Laplace GLMM

- [`ds.vertLMM()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLMM.md)
  : Federated linear mixed model with a single random intercept
- [`ds.vertLMM.k3()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLMM.k3.md)
  : Federated linear mixed model via REML 1-D profile (K=3)
- [`ds.vertGEE()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGEE.md)
  : Federated generalised estimating equations
- [`ds.vertGLMM()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGLMM.md)
  : Federated binomial GLMM via Laplace approximation

## Causal inference / robustness

Inverse-probability-weighting and multiple-imputation pooling

- [`ds.vertIPW()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertIPW.md)
  : Federated inverse-probability-weighted GLM (two-stage)
- [`ds.vertMI()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMI.md)
  : Federated multiple imputation with Rubin pooling

## LASSO / penalised regression

Proximal-gradient / FISTA on the federated GLM normal-equations

- [`ds.vertLASSO()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSO.md)
  : Post-hoc soft-thresholded GLM coefficients (naive LASSO)
- [`ds.vertLASSO1Step()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSO1Step.md)
  : One-step LASSO via quadratic-surrogate proximal gradient
- [`ds.vertLASSOIter()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSOIter.md)
  : Iterative proximal-gradient LASSO over the MPC GLM gradient
- [`ds.vertLASSOCV()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSOCV.md)
  : Information-criterion lambda selection for one-step LASSO
- [`ds.vertLASSOProximal()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSOProximal.md)
  : Proper LASSO via client-side proximal gradient on normal equations

## Utility / printers

Helper and administrative methods

- [`coef(`*`<ds.glm>`*`)`](https://isglobal-brge.github.io/dsVertClient/reference/coef.ds.glm.md)
  : Coefficients Method for ds.glm Objects
- [`print(`*`<ds.cor>`*`)`](https://isglobal-brge.github.io/dsVertClient/reference/print.ds.cor.md)
  : Print Method for ds.cor Objects
- [`print(`*`<ds.pca>`*`)`](https://isglobal-brge.github.io/dsVertClient/reference/print.ds.pca.md)
  : Print Method for ds.pca Objects
- [`print(`*`<ds.glm>`*`)`](https://isglobal-brge.github.io/dsVertClient/reference/print.ds.glm.md)
  : Print Method for ds.glm Objects
- [`summary(`*`<ds.glm>`*`)`](https://isglobal-brge.github.io/dsVertClient/reference/summary.ds.glm.md)
  : Summary Method for ds.glm Objects

## Internal protocol implementations

- [`chunk-utils`](https://isglobal-brge.github.io/dsVertClient/reference/chunk-utils.md)
  : Adaptive Chunking for DataSHIELD Transport
- [`k2-beaver-lbfgs-client`](https://isglobal-brge.github.io/dsVertClient/reference/k2-beaver-lbfgs-client.md)
  : K=2 Beaver L-BFGS Pipeline
- [`k3-ring63-dcf-loop`](https://isglobal-brge.github.io/dsVertClient/reference/k3-ring63-dcf-loop.md)
  : K\>=3 Ring63 DCF + Beaver Gradient Loop
