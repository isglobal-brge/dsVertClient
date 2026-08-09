# Acknowledged Chunking for DataSHIELD Transport

Utilities for sending large payloads through DataSHIELD's R expression
parser with fixed transfer geometry and typed acknowledgements.

## Details

DataSHIELD passes function arguments as inline R expressions via HTTP.
Large string arguments (base64-encoded ciphertexts, EC points, encrypted
vectors) can exceed the parser or HTTP body size limit. This module
implements acknowledged chunking:

1.  Start with a portable initial chunk size (default 640 KiB of Base64
    text / at most 480 KiB raw, configurable via
    `options(dsvert.chunk_size = N)`)

2.  Freeze the chunk size, chunk count, indices, and bytes before the
    first request

3.  Require an exact logical `TRUE` acknowledgement from the intended
    DataSHIELD target before advancing

4.  If an idempotent store loses its response, replay only the exact
    same request until a monotonic availability deadline; never silently
    change transfer geometry

Automatic fallback to a smaller geometry is deliberately not attempted
after transmission starts: the server may already have committed the
first chunk even when its response is lost. Changing `n_chunks` at that
point would turn a recoverable response loss into a conflicting
transfer.
