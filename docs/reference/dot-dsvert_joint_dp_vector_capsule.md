# Materialize or replay the one server-authoritative DP capsule vector

Internal orchestration. There is no analyst-selected epsilon, bound,
workload, peer, source, seed, allocation index, request counter or
fallback mechanism. An authenticated stale-root signal emitted before
the first valid START refreshes signed roots and retries the immutable
handshake at most once. After START has claimed an instance, another
instance fails closed.

## Usage

``` r
.dsvert_joint_dp_vector_capsule(
  datasources,
  status = NULL,
  manifest_bundle = NULL,
  .aggregate = DSI::datashield.aggregate,
  .source_transport = .dsvert_dp_capsule_source_transport_context,
  .cross_orchestrate = .dsvert_dp_cross_orchestrate,
  .setup_transport = .dsvert_setup_peer_transport,
  .setup_exact = .dsvert_setup_exact_gc_transport,
  .run_exact = .dsvert_joint_dp_vector_exact_gc_run,
  .store_typed = .dsvert_store_typed_blob,
  .allocation = .dsvert_joint_dp_vector_allocation,
  .retry_refresh = NULL
)
```
