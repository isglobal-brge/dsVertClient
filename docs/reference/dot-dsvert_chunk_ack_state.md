# Classify a one-target DataSHIELD acknowledgement (internal)

Classify a one-target DataSHIELD acknowledgement (internal)

## Usage

``` r
.dsvert_chunk_ack_state(response, target)
```

## Arguments

- response:

  Result returned by
  [`DSI::datashield.aggregate`](https://datashield.github.io/DSI/reference/datashield.aggregate.html).

- target:

  Expected logical DataSHIELD connection name.

## Value

One of `"ack"`, `"missing"`, or `"invalid"`.
