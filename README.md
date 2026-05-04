# dsVertClient - DataSHIELD Client for Vertically Partitioned Data

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version](https://img.shields.io/badge/version-1.1.0-blue.svg)](NEWS.md)

## Overview

**dsVertClient** provides user-friendly R functions for privacy-preserving analysis on vertically partitioned data across DataSHIELD servers. The analyst calls simple functions; all cryptographic protocols (ECDH-PSI, Ring63 / Ring127 Beaver MPC, DCF wide-spline, OT-Beaver dishonest-majority triples, X25519 + AES-256-GCM transport, Ed25519 identity verification) run transparently.

Pair with the server-side companion package [**dsVert**](https://github.com/isglobal-brge/dsVert).

## Quick Start

```r
library(DSI); library(DSOpal); library(dsVertClient)

# Connect to servers
conns <- DSI::datashield.login(logins = builder$build())
DSI::datashield.assign.table(conns, "D", "project.data")

# 1. Align records (NAs removed automatically, like glm na.omit)
ds.psiAlign("D", "patient_id", "DA", datasources = conns)

# 2. Fit a federated GLM (auto-detects which server has each variable)
fit <- ds.vertGLM(diabetes ~ age + bmi + glu + bp,
                  data = "DA", family = "binomial",
                  datasources = conns)
print(fit)

# 3. Cox PH via the non-disclosive Breslow profile route
cox <- ds.vertCoxProfileNonDisclosive(
  survival::Surv(time, event) ~ age + bmi + bp,
  data = "DA",
  datasources = conns
)

# 4. Correlation + PCA
cor <- ds.vertCor("DA", datasources = conns)
pca <- ds.vertPCA(cor_result = cor)

DSI::datashield.logout(conns)
```

## Functions (v1.1.0)

| Family | Functions |
|---|---|
| **Record alignment** | `ds.psiAlign()`, `ds.isPsiAligned()`, `ds.getIdentityPks()` |
| **Descriptive / 2nd-order** | `ds.vertDesc()`, `ds.vertCor()`, `ds.vertPCA()`, `ds.vertChisq()`, `ds.vertChisqCross()`, `ds.vertFisher()` |
| **GLM** (gaussian / binomial / poisson) | `ds.vertGLM()` with `offset`, `weights`, `ring`, `binomial_sigmoid_intervals`, `keep_session`, `no_intercept`, `std_mode` |
| **Inference helpers** | `ds.vertConfint()`, `ds.vertContrast()`, `ds.vertWald()`, `ds.vertLR()` |
| **Survival** | `ds.vertCox()` / `ds.vertCoxProfileNonDisclosive()` (non-disclosive Breslow Cox PH), `ds.vertCoxDiscreteNonDisclosive()` (pooled-logistic discrete survival); legacy rank/person-time routes removed |
| **Negative binomial** | `ds.vertNBFullRegTheta(variant = "full_reg_nd")` (default share-domain full-reg θ), `ds.vertNB()` / `ds.vertNBMoMTheta()` for lighter scalar-theta variants |
| **Multinomial** | `ds.vertMultinom()` / `ds.vertMultinomJointNewton()` (joint softmax Newton); warm OVR is internal initialisation only |
| **Ordinal (proportional odds)** | `ds.vertOrdinal()` / `ds.vertOrdinalJointNewton()` (joint proportional-odds Newton); warm cumulative-binomial is internal initialisation only |
| **Mixed models** | `ds.vertLMM()` (REML closed-form, K=2; random intercept + slopes), `ds.vertLMM.k3()` (REML 1-D profile, K>=3), `ds.vertGEE()` (sandwich SE, `binomial_sigmoid_intervals`), `ds.vertGLMM()` (binomial GLMM-PQL) |
| **Causal / robustness** | `ds.vertIPW()` (two-stage propensity + weighted GLM), `ds.vertMI()` (multiple imputation + Rubin pooling) |
| **Penalised regression** | `ds.vertLASSO()`, `ds.vertLASSO1Step()`, `ds.vertLASSOIter()` (Gaussian/binomial/Poisson standardized L1), `ds.vertLASSOCV()` (AIC / BIC / EBIC selector), `ds.vertLASSOProximal()` |

## Output sample (`ds.vertGLM`)

```
Coefficients:
                Estimate Std.Error z value  Pr(>|z|)
(Intercept)     -9.2784    2.0981  -4.425  0.000010 ***
age              0.0733    0.0223   3.301  0.000993 ***
bmi              0.0798    0.0522   1.534  0.126024
glu              0.0290    0.0101   2.878  0.004046 **

Converged: TRUE (14 iterations)
Deviance: 22.09
```

## K=2 vs K≥3 support

| Family | K=2 product route | K≥3 product route | Legacy / not offered |
|---|---|---|---|
| GLM (gauss / binom / poisson) | `ds.vertGLM()` K=2 Beaver MPC | `ds.vertGLM()` secure-agg / DCF-pair route | manual mapping aliases only |
| Cox PH | `ds.vertCox()` / `ds.vertCoxProfileNonDisclosive()`; discrete-time ND route also available | same profile ND route | `legacy_rank`, `ds.vertCox.k3()` removed |
| Negative binomial | `ds.vertNBFullRegTheta(variant = "full_reg_nd")`; iid/MoM scalar-theta variants | same full-reg ND route; iid/MoM scalar-theta variants | disclosive `variant = "full_reg"` removed |
| Multinomial | `ds.vertMultinom()` / `ds.vertMultinomJointNewton()` | same joint softmax route | warm / OVR final-estimator route removed from the exported API |
| Ordinal | `ds.vertOrdinal()` / `ds.vertOrdinalJointNewton()` | same proportional-odds route | warm final-estimator and patient-level joint reconstruction routes removed from the exported API |
| LMM | `ds.vertLMM()` K=2 closed form | `ds.vertLMM.k3()` REML 1-D profile | direct client-supplied cluster-vector helper is not product |
| GEE | `ds.vertGEE()` exchangeable / guarded AR1 | same route | unguarded order metadata not accepted |
| GLMM | `ds.vertGLMM()` | same PQL aggregate route | `method = "em"` removed |
| IPW / MI / LASSO / Cor / PCA / Chisq / Desc | product wrappers | same product wrappers | small-cell / high-dimensional diagnostics gated |

See `inst/docs/product_surface.md` for the disclosure/accuracy status table
used by the cleanup audit.

## Validation evidence

| Method | max\|Δβ\| | Reference | Theoretical bound |
|---|---|---|---|
| GLM binomial (Pima, K=2) | 3.18 × 10⁻⁴ | `lm()` / `glm()` | Catrina-Saxena 2010 fp20 |
| Cox PH (NCCTG, Ring127) | 1.33 × 10⁻⁵ | `survival::coxph` | Catrina-Saxena 2010 fp50 |
| Cox PH (Pima, Ring127) | 8.09 × 10⁻⁵ | `survival::coxph` | Catrina-Saxena 2010 fp50 |
| Cox PH in LMM+Cox harness | 8.32 × 10⁻⁶ | `survival::coxph` | Catrina-Saxena 2010 fp50 |
| NB iid θ | 6.83 × 10⁻¹² | `MASS::theta.ml` | Newton precision |
| NB full-reg θ | 4.44 × 10⁻³ | `MASS::glm.nb` | Catrina-Saxena fp50 + 5-iter NR (Goldschmidt 1964 + Pugh 2004) — σ-ratio 105×, **SUB-NOISE** |
| Multinomial joint (Bohning) | 4.70 × 10⁻² | `nnet::multinom` | Catrina-Saxena Ring63 + Bohning |
| Ordinal joint PO | 6.12 × 10⁻² | `MASS::polr` | 3× McCullagh-Agresti L1 floor |
| LASSO proximal | 1-2 × 10⁻³ | `glmnet::cv.glmnet` | FHT 2010 + Beck-Teboulle 2009 |

All methods inside their theoretical floors (paper §V.A). **Sub-noise margin** (paper §V.B) at ≥ 100× the per-fit Wald SE for all four federated K=2 estimators (ord_joint β/θ, Cox β, NB θ).

## Security

- **No product observation-level disclosure**: client sees only model-scale
  aggregates returned by the registered server methods
- **Server-generated Beaver triples**: client never sees cryptographic material
- **Dealer rotation**: different server generates triples each iteration (K ≥ 3); fixed dealer with OT-Beaver for K = 2
- **Transport encryption**: X25519 + AES-256-GCM between servers
- **Identity verification**: Ed25519 signed peer transport keys (`dsvert.require_trusted_peers`)
- **Ring**: Ring63 (frac_bits = 20) and Ring127 (frac_bits = 50), selected by
  method precision needs
- **Discarded estimators removed**: routes that revealed rank metadata,
  patient-level working vectors, or old approximations are not part of the
  exported product API

## Installation

```r
# From GitHub
devtools::install_github("isglobal-brge/dsVertClient")

# Or from source
R CMD build --no-build-vignettes .
R CMD INSTALL "$(ls -t dsVertClient_*.tar.gz | head -1)"
```

Requires the server-side package [dsVert](https://github.com/isglobal-brge/dsVert) installed on all DataSHIELD servers.

## Documentation

- [Reference index](https://isglobal-brge.github.io/dsVertClient/) (pkgdown site)
- [Method validation evidence](https://isglobal-brge.github.io/dsVertClient/articles/vert_validation_evidence.html)
- [Product surface audit](inst/docs/product_surface.md)
- [NEWS](NEWS.md)

## License

MIT - see [LICENSE](LICENSE.md).

## Citation

See `paper/jbhi_dsvert.tex` (IEEE J-BHI submission r2.5) for the full validation table, theoretical bounds, and disclosure ledger.
