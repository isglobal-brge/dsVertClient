# One-step Newton Cox at beta=0 (client orchestrator)

Called by ds.vertCox with one_step_newton = TRUE. Pre-conditions (caller
has already:) \* standardised X via .glm_mpc_setup() \* shared X
(k2ShareInputDS / k2ReceiveShareDS) \* registered Cox times
(k2SetCoxTimesDS / k2ReceiveCoxMetaDS) \* applied permutation
(k2ApplyCoxPermutationDS) so ss\$k2_x_share_fp, ss\$k2_peer_x_share_fp,
ss\$k2_cox_delta_fp are all sorted by time on both parties.

Does: 1. dsvertCoxNewtonPrepDS on both parties (local extract + cumsum).
2. dsvertCoxNewtonGradDS on both parties -\> aggregate grad vector. 3.
For each Fisher pair (j \<= k), run two Beaver vecmul protocols (one for
the X*X term, one for the S*S term), then scalar aggregate. 4. Solve
beta = solve(Fisher, grad). Returns std-scale beta.

## Usage

``` r
.ds_vertCox_newton_one_step(
  datasources,
  server_names,
  server_list,
  y_server,
  nl,
  session_id,
  n_obs,
  transport_pks,
  p_coord,
  p_nl,
  .dsAgg = NULL,
  .sendBlob = NULL,
  verbose = TRUE,
  ring = 63L
)
```
