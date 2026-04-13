# Adaptive Chunking for DataSHIELD Transport

Utilities for sending large payloads through DataSHIELD's R expression
parser with automatic chunk size detection and fallback.

## Details

DataSHIELD passes function arguments as inline R expressions via HTTP.
Large string arguments (base64-encoded ciphertexts, EC points, encrypted
vectors) can exceed the parser or HTTP body size limit. This module
implements adaptive chunking:

1.  Start with a large initial chunk size (default 100KB, configurable
    via `options(dsvert.chunk_size = N)`)

2.  On first large send, probe whether the chunk size is accepted

3.  If the probe fails, reduce the chunk size by 25\\

4.  Cache the successful chunk size for all subsequent sends

The probe is only conclusive when a chunk of at least half the target
chunk size is successfully sent. Small payloads that succeed trivially
do not confirm that the full chunk size is accepted.
