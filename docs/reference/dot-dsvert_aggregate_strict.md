# Strict DataSHIELD aggregate result contract (internal)

Strict DataSHIELD aggregate result contract (internal)

## Usage

``` r
.dsvert_aggregate_strict(
  conns,
  expr,
  operation = "DataSHIELD aggregate",
  allow_null = character(),
  result_contract = c("non_null", "logical_true"),
  .aggregate = DSI::datashield.aggregate
)
```

## Arguments

- conns:

  Named DataSHIELD connection list.

- expr:

  One aggregate call or a named list of calls.

- operation:

  Public, non-sensitive phase label used in errors.

- allow_null:

  Named sites whose documented protocol contract permits a NULL result.
  Empty by default.

- result_contract:

  Either `"non_null"` or `"logical_true"`.

- .aggregate:

  Aggregate implementation, injectable for tests.

## Value

A result list ordered exactly like `conns`.
