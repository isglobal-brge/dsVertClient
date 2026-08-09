# Run one exact server-to-server operation over the DSI connections

Internal orchestration primitive. It performs one fan-out exchange per
cycle and retries the exact same request after transient DSI failures.
There is deliberately no request counter or rate limit; the optional
monotonic timeout is an inactivity lease renewed only by verified bytes,
offsets or a terminal state transition.

## Usage

``` r
.dsvert_exact_gc_run(
  datasources,
  server_names = names(datasources),
  servers = seq_along(datasources),
  session_id,
  operation_id,
  source_key,
  output_key,
  operation,
  ring,
  frac_bits = 0L,
  vector_len,
  purpose,
  transport_ready = FALSE,
  timeout_seconds = getOption("dsvert.exact_gc.timeout_seconds", 900),
  initialized = NULL,
  mul_contract = NULL,
  alignment_source_count = NULL,
  .aggregate = DSI::datashield.aggregate
)
```
