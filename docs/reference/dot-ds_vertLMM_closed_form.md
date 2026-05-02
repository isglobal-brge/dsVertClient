# LMM closed-form GLS driver

Computes the exact Laird-Ware generalised-least-squares estimate for a
random-intercept LMM without relying on `ds.vertGLM`'s L-BFGS iteration.
Directly assembles the Gram matrix X'X and right-hand side X'y over the
cluster-mean- centred design via: - local per-server block
computations - Beaver vecmul + FP-sum for cross-server entries and
solves `beta = solve(XtX, Xty)` client-side. Matches
[`lme4::lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) to FP precision
(~1e-5 on small n).

Assumptions: - Cluster IDs have been broadcast via
`dsvertLMMBroadcastClusterIDsDS` / `dsvertLMMReceiveClusterIDsDS`. -
Transport keys are initialised in the session. - Exactly K=2 servers
(outcome + one peer).

## Usage

``` r
.ds_vertLMM_closed_form(
  conns,
  server_names,
  y_srv,
  peer_srv,
  data,
  y_var,
  x_ysrv,
  x_peer,
  lambda_i,
  transport_pks,
  session_id,
  verbose = FALSE,
  share_scale = 1,
  column_scales = NULL,
  standardize = FALSE,
  ring = "ring63"
)
```

## Arguments

- conns:

  Named list of DS connections.

- server_names:

  Character vector of server names (names of `conns`).

- y_srv, peer_srv:

  Character. Outcome server name and peer server name.

- data:

  Aligned data-frame name on both servers.

- y_var:

  Response column (on `y_srv`).

- x_ysrv, x_peer:

  Character vectors of predictor names on each server (after formula
  parsing).

- lambda_i:

  Numeric vector of length n_clusters.

- transport_pks:

  Named list of public keys per server.

- session_id:

  Active MPC session (cluster IDs stored within).

## Value

list(coefficients, XtX, Xty, yty, n)
