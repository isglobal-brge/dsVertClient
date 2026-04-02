# K=2 Secure GLM Implementation Plan

## Overview

Implement family-specific two-party (K=2) secure GLM pathways that prevent
the label server from seeing the non-label server's η contribution.

Current state: K=2 falls to "transport" mode where η_nonlabel is visible
to the coordinator. This is the primary blind spot.

## Architecture

Three new modes, sharing common infrastructure:

```
K=2 auto-selection:
  gaussian  → k2_gaussian_exact   (encrypted residual, exact block update)
  binomial  → k2_binomial_helink  (polynomial sigmoid under CKKS)
  poisson   → k2_poisson_helink   (polynomial exp under CKKS)
```

All three share:
- Owner-targeted decryption (each server only sees its own gradient)
- Encrypted η aggregation (ct_η = ct_ηL + ct_ηN, never plaintext)
- Encrypted residual (ct_r = ct_y - ct_μ, never plaintext)

## Implementation order

### Phase 0: Shared infrastructure
1. Owner-targeted decryption in mheFullProtocol.R
2. mheEncryptEtaDS (encrypt η_k = X_k β_k under CPK)
3. mheAddCiphertextsDS (ct_η = ct_ηL + ct_ηN)
4. mheGradientFromResidualDS (ct_gk = X_k^T ct_r, owner-tagged)
5. K2 mode selector in ds.vertGLM.R
6. Policy enforcement in manifest/FSM

### Phase 1: K=2 Gaussian (exact, no polynomial needed)
1. mheResidualGaussianDS (ct_r = ct_y - ct_ηL - ct_ηN)
2. K2 Gaussian BCD loop in ds.vertGLM.R
3. Tests: correctness vs centralized, security invariants
4. Vignette: k2-gaussian-pima.Rmd

### Phase 2: K=2 Binomial (HE-Link with polynomial sigmoid)
1. Generate sigmoid polynomial coefficients (Remez/Chebyshev, degree 7, [-8,8])
2. mheEvalPolynomialDS (Horner evaluation under CKKS)
3. K2 Binomial gradient descent loop in ds.vertGLM.R
4. Tests: correctness, polynomial accuracy, security invariants
5. Vignette: k2-binomial-pima-helink.Rmd

### Phase 3: K=2 Poisson (HE-Link with polynomial exp)
1. Generate exp polynomial coefficients (degree 7, [-4,4])
2. Reuse mheEvalPolynomialDS with different coefficients
3. K2 Poisson gradient descent loop in ds.vertGLM.R
4. Tests: correctness, stability, security invariants
5. Vignette: k2-poisson-pima-helink.Rmd

## Key design decisions

### Why not IRLS for K=2 binomial/poisson?
IRLS needs μ, w in plaintext for the coordinator step. In K=2, exposing μ
to the label server leaks η (invertible link). Gradient descent avoids this
by keeping μ encrypted throughout.

### Why exact block update for Gaussian but gradient for binomial/poisson?
Gaussian link is identity (μ = η, w = 1), so there's no non-linear function
to approximate. The residual r = y - η can be computed exactly under HE,
and the block update (X^T X + λI)^-1 (X^T X β + g) is exact.

For binomial/poisson, the non-linear link (sigmoid/exp) must be approximated
by a polynomial, which introduces approximation error. Using exact IRLS
would require μ in plaintext, defeating the purpose.

### Owner-targeted decryption
Instead of a fixed fusion server that sees all decrypted values, each
gradient ciphertext is tagged with an owner_server. Only the owner can
fuse the partial shares and see the plaintext gradient. The other party
contributes its partial share but never sees the result.

### Clipping for polynomial stability
Before encrypting η_k, clip to [-B_loc, B_loc]:
- Gaussian: no clipping needed (identity link, no polynomial)
- Binomial: B_loc = 4 (total η in [-8,8] for sigmoid)
- Poisson: B_loc = 2 (total η in [-4,4] for exp)

Features are standardized + clipped to [-5, 5] before analysis.

## Files to create/modify

### New files
- dsVert/R/glmK2.R — K=2 server-side functions
- dsVert/inst/extdata/poly_sigmoid_deg7.json — sigmoid coefficients
- dsVert/inst/extdata/poly_exp_deg7.json — exp coefficients
- dsVertClient/R/ds.vertGLM.k2.R — K=2 client orchestration (or integrate into ds.vertGLM.R)
- dsVertClient/vignettes/k2-gaussian-pima.Rmd
- dsVertClient/vignettes/k2-binomial-pima-helink.Rmd
- dsVertClient/vignettes/k2-poisson-pima-helink.Rmd

### Modified files
- dsVert/R/mheFullProtocol.R — owner-targeted decryption
- dsVert/DESCRIPTION — new AggregateMethods
- dsVert/NAMESPACE — new exports
- dsVertClient/R/ds.vertGLM.R — K2 mode selector + loops
- Paper: npjdigitalmed_main_final.tex — K=2 Methods subsection

## New DataSHIELD methods needed

### Server-side (dsVert)
```
# Shared K=2 infrastructure
mheEncryptEtaDS(data_name, x_vars, beta, clip_radius, session_id)
mheAddCiphertextsDS(ct_key_a, ct_key_b, result_key, session_id)
mheGradientFromResidualDS(data_name, x_vars, residual_key, owner_server, session_id)
mheFuseOwnedDS(ct_key, session_id)  # owner-only fusion

# Gaussian K=2
mheResidualGaussianDS(ct_y_key, ct_eta_label_key, ct_eta_nonlabel_key, session_id)

# Binomial/Poisson K=2
mheEvalPolynomialDS(input_key, coefficients, result_key, session_id)
mheResidualFromMuDS(ct_y_key, ct_mu_key, result_key, session_id)
```

### Go binary (mhe-tool) new commands
```
ckks-add          # ct_a + ct_b → ct_result
ckks-poly-eval    # evaluate polynomial on ciphertext (Horner scheme)
partial-decrypt-owned  # partial decrypt tagged with owner
fuse-owned        # fuse shares, only if caller is owner
```

## Test plan

### Pima K=2 split
- Server 1 (non-label): patient_id, age, bmi, ped (180 patients)
- Server 2 (label): patient_id, glu, bp, skin, npreg, diabetes (175 patients)

### Test models
1. Gaussian: y = glu, predictors = age, bmi, ped (s1) + bp, skin, npreg, diabetes (s2)
2. Binomial: y = diabetes, predictors = age, bmi, ped (s1) + glu, bp, skin, npreg (s2)
3. Poisson: y = npreg, predictors = age, bmi, ped (s1) + glu, bp, skin, diabetes (s2)

### Correctness criteria
- Gaussian: max_abs_coef_diff ≤ 1e-4 (exact block update, no polynomial error)
- Binomial: max_abs_coef_diff ≤ 5e-3 (polynomial + gradient descent approximation)
- Poisson: max_abs_coef_diff ≤ 1e-2 (polynomial + gradient descent + narrower range)

### Security invariants (automated tests)
- Label server never stores/receives η_nonlabel in plaintext
- Non-label server never stores/receives y in plaintext
- Owner-targeted decrypt: label cannot fuse non-label's gradient CT
- No n-length plaintext vectors cross the session boundary

### Stress test
- Single-predictor non-label: Server 1 has only `age`
- Verify that even with p=1, the label cannot recover X_1 from gradient history

## References for paper
- Kim et al. 2018 (logistic regression over approximate HE)
- Kelkar et al. 2022 (secure Poisson regression)
- Bonawitz et al. 2016 (secure aggregation — why it fails for K=2)
- Jiang et al. 2022 (VFL prediction leakage)
- Ghavamipour et al. 2022 (secret-sharing logistic regression)

## Timeline estimate
- Phase 0 (infrastructure): 1-2 sessions
- Phase 1 (Gaussian): 1 session
- Phase 2 (Binomial): 1-2 sessions
- Phase 3 (Poisson): 1 session
- Paper updates + vignettes: 1 session
- Total: ~5-7 sessions of focused work
