# Federated linear mixed model via REML 1-D profile (K=3)

Random-intercept LMM \$\$y\_{ij} = X\_{ij} \beta + b_i +
\varepsilon\_{ij}, \quad b_i \sim \mathcal{N}(0, \sigma_b^2),\\
\varepsilon\_{ij} \sim \mathcal{N}(0, \sigma^2)\$\$ estimated by REML
1-D profile over the intra-cluster correlation \\\rho = \sigma_b^2 /
(\sigma_b^2 + \sigma^2/n_i)\\. For balanced designs the GLS-transformed
score equations are solved exactly by a weighted Gaussian `ds.vertGLM`
call; the profile likelihood at each candidate \\\rho\\ is evaluated
from the inner fit's `$deviance` plus the per-cluster log-determinant
closed form (Christensen 2019 Sec.A.3, Pinheiro & Bates 2000 Sec.2.4).

K=3 implementation reuses the secure-aggregation `ds.vertGLM` pipeline
at each profile evaluation – no new MPC primitives. The per-cluster GLS
weights are a public function of \\\rho\\ and the cluster sizes (already
disclosed under `datashield.privacyLevel`), so D-INV-1..5 are preserved
identically to the K=2 `ds.vertLMM` path.

Outer optimisation is golden-section search over \\\rho \in
(\rho\_{\min}, \rho\_{\max})\\; default \\(0.001, 0.999)\\ so the search
never enters the boundary singularities at \\\rho = 0\\ (no random
effect – degenerate GLS) or \\\rho = 1\\ (infinite \\\sigma_b^2\\).

## Usage

``` r
ds.vertLMM.k3(
  formula,
  data,
  cluster_col,
  rho_lo = 0.001,
  rho_hi = 0.999,
  tol = 0.001,
  max_outer = 30L,
  verbose = TRUE,
  datasources = NULL
)
```

## Arguments

- formula:

  Fixed-effects formula `y ~ X`.

- data:

  Aligned data-frame name on each server.

- cluster_col:

  Cluster id column on the outcome server.

- rho_lo, rho_hi:

  Profile search interval (default 0.001, 0.999).

- tol:

  Profile convergence tolerance on \\\rho\\ (default 1e-3).

- max_outer:

  Maximum golden-section iterations.

- verbose:

  Print progress.

- datasources:

  DataSHIELD K=3 connections.

## Value

list of class `ds.vertLMM.k3` with `coefficients`, `sigma_b2`, `sigma2`,
`rho_hat`, `n_clusters`, `cluster_sizes`, and the inner ds.glm fit.

## References

Christensen, R. H. B. (2019). *Linear Models*.
[`ordinal::clm.fit`](https://rdrr.io/pkg/ordinal/man/clm.fit.html)
Sec.A.3 – REML profile likelihood for variance components. Pinheiro, J.
C. & Bates, D. M. (2000). *Mixed-Effects Models in S/S-PLUS*, Sec.2.4.
Laird, N. M. & Ware, J. H. (1982). Random-effects models for
longitudinal data. *Biometrics*, 38(4), 963-974. Lindstrom, M. J. &
Bates, D. M. (1990). Nonlinear mixed effects models for repeated
measures data. *Biometrics*, 46(3), 673-687.

## See also

[`ds.vertGLM`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGLM.md),
[`ds.vertLMM`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLMM.md)
(K=2 exact closed-form path).
