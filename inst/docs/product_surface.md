# dsVertClient Product Surface Audit

Status: development cleanup, 2026-05-04.

This table records the routes intended for users and paper validation. It is
not a list of every internal helper in the package. Routes that were
disclosive or materially less accurate than an available product route are
removed from the exported API.

| Method family | K=2 product route | K>=3 product route | Disclosure status | Accuracy status | Not offered by default |
|---|---|---|---|---|---|
| PSI / alignment | `ds.psiAlign()` | `ds.psiAlign()` | Non-disclosive PSI intersection; row indices are not product-returned | Functional prerequisite, not an estimator | `psiGetMatchedIndicesDS` style diagnostics |
| GLM | `ds.vertGLM()` K=2 Beaver MPC | `ds.vertGLM()` secure aggregation with selected DCF pair | Model-scale gradients, Fisher/Hessian summaries, deviance; no patient-level vectors | DSLite-validated against central GLM within accepted numerical envelope | Manual server mapping aliases are compatibility only |
| Cox PH | `ds.vertCox()` / `ds.vertCoxProfileNonDisclosive()`; `ds.vertCoxDiscreteNonDisclosive()` for discrete-time target | Same profile ND route; discrete-time ND route uses a selected DCF pair when requested | Profile route avoids event-rank metadata exposure; discrete route uses share-mask gating | DSLite evidence supports profile ND vs `survival::coxph`; discrete route is a safe target alternative | `legacy_rank` and `ds.vertCox.k3()` removed |
| Negative binomial | `ds.vertNBFullRegTheta(variant = "full_reg_nd")`; `ds.vertNB()` and `ds.vertNBMoMTheta()` for scalar-theta alternatives | Same full-reg ND route, with selected DCF pair | Share-domain theta refinement; no per-patient non-label eta reveal | Full-reg ND is the preferred accuracy route; iid/MoM are lighter approximations | Disclosive `variant = "full_reg"` eta-transport removed |
| Multinomial | `ds.vertMultinom()` / `ds.vertMultinomJointNewton()` | Same joint softmax route | Aggregate score/Hessian route over secret shares | Joint route is the paper route; OVR gap is intrinsic and not used as product accuracy target | Warm/OVR final estimator removed from exported API; internal warm start only |
| Ordinal | `ds.vertOrdinal()` / `ds.vertOrdinalJointNewton()` | Same proportional-odds route | Strict path avoids patient-level eta/F reveal | Joint route stays within accepted practical envelope in DSLite evidence | Patient-level ordinal joint reconstruction and warm-only final estimator removed from exported API |
| LMM | `ds.vertLMM()` K=2 closed form | `ds.vertLMM.k3()` REML 1-D profile | Guarded cluster-size/cluster-summary disclosure only | DSLite-validated against centralized mixed-model targets in the current evidence set | Client-supplied cluster-vector helper as product route |
| GEE | `ds.vertGEE()` exchangeable / guarded AR1 | Same route | Guarded order metadata for AR1; sandwich summaries only | Current accepted route uses precision controls for binomial sigmoid | Unguarded AR1/order routes |
| GLMM | `ds.vertGLMM()` | Same PQL aggregate route | Share-domain cluster sufficient statistics; no BLUP vector returned | PQL aggregate is the supported approximate mixed-model route | Legacy EM/Laplace-style route removed |
| IPW | `ds.vertIPW()` | `ds.vertIPW()` | Propensity and outcome GLM aggregates only | Covered by DSLite method evidence | None known beyond base GLM gates |
| MI | `ds.vertMI()` | `ds.vertMI()` | Server-local imputations plus Rubin pooling summaries; returned object audited for no row vectors | Covered by returned-object audit and DSLite evidence | Returning imputed patient-level data |
| LASSO | `ds.vertLASSOIter()`, `ds.vertLASSOProximal()`, related selectors | Same wrappers | Built from non-disclosive GLM/cross-product aggregates | Product route is accepted for coefficient/support agreement | Post-hoc-only shortcuts are not preferred when exact/proximal route applies |
| Cor / PCA | `ds.vertCor()`, `ds.vertPCA()` | Same wrappers | Standard second-order aggregate disclosure tier | Product baseline for acceptable aggregate disclosure | High-dimensional release blocked by guards |
| Chisq / Fisher / Desc | `ds.vertChisq*()`, `ds.vertFisher()`, `ds.vertDesc()` | Same wrappers | Small cells and high-risk releases are guarded | Descriptive/statistical target, not iterative estimator | Small-cell diagnostics unless explicitly enabled |

Operational rule for cleanup: if a route is disclosive or materially less
accurate and a product route exists for the same statistical target, the route
must not be user-invokable. Warm starts needed by product Newton methods may
remain as non-exported helpers, but not as final estimators.
