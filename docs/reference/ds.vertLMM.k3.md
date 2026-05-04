# Federated linear mixed model via share-domain residual GLS (K\>=3)

Random-intercept LMM \$\$y\_{ij} = X\_{ij} \beta + b_i +
\varepsilon\_{ij}, \quad b_i \sim \mathcal{N}(0, \sigma_b^2),\\
\varepsilon\_{ij} \sim \mathcal{N}(0, \sigma^2)\$\$ estimated by an
initial REML 1-D profile followed by share-domain residual sufficient
statistics and the closed-form cluster-mean GLS transform. The final
fixed-effect step applies \\\tilde v\_{ij} = v\_{ij} - \lambda_i \bar
v_i\\ locally on each server and fits the transformed Gaussian system
with an explicit transformed intercept \\1-\lambda_i\\.

K\>=3 implementation reuses the secure-aggregation `ds.vertGLM`
pipeline. Residual sums and squared-residual sums are evaluated at a
supplied beta inside the DCF share domain; the analyst client never
receives row-level residuals, fitted values, eta, or transformed
columns. Cluster membership is sent as encrypted server-to-server
integer metadata, with original labels withheld and small clusters
blocked by `datashield.privacyLevel`.

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
  ring = c("ring127", "ring63"),
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

- ring:

  Character (`"ring127"` or `"ring63"`). MPC ring used for the K\>=3
  Gaussian GLM, residual squared sums, and GLS refits. Ring127 is the
  default because variance components are sensitive to fixed-point noise
  in the secure \\r^2\\ pass.

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
