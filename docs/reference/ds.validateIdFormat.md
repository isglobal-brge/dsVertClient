# Validate Identifier Format Across Servers

Client-side function that validates identifier format consistency across
all DataSHIELD servers before record alignment.

## Usage

``` r
ds.validateIdFormat(data_name, id_col, pattern = NULL, datasources = NULL)
```

## Arguments

- data_name:

  Character string. Name of the data frame on each server.

- id_col:

  Character string. Name of the identifier column.

- pattern:

  Character string (optional). Regular expression pattern that IDs
  should match on all servers.

- datasources:

  DataSHIELD connection object(s). If NULL, uses all available
  connections.

## Value

A list with class "ds.id.validation" containing:

- `valid`: Logical, TRUE if formats are consistent across servers

- `servers`: Data frame with validation results per server

- `format_match`: Logical, TRUE if all servers have same format
  signature

- `pattern_match`: Logical, TRUE if all servers match the pattern (only
  if pattern provided)

- `warnings`: Character vector of any warnings detected

## Details

This function should be called before `ds.alignRecords` to ensure that
identifier formats are consistent across all data partitions. It helps
catch common issues like:

- Different ID formats (e.g., "001" vs "1" vs "ID001")

- Missing or duplicate identifiers

- Type mismatches (numeric vs character)

The function uses format signatures (hashed format characteristics) to
compare formats without revealing actual identifier values.

## See also

[`ds.hashId`](https://isglobal-brge.github.io/dsVertClient/reference/ds.hashId.md),
[`ds.alignRecords`](https://isglobal-brge.github.io/dsVertClient/reference/ds.alignRecords.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Validate ID format before alignment
validation <- ds.validateIdFormat("D", "patient_id", datasources = conns)
print(validation)

if (validation$valid) {
  # Proceed with alignment
  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes, "D_aligned", conns)
}

# With pattern validation
validation <- ds.validateIdFormat("D", "patient_id",
                                  pattern = "^[A-Z]{2}[0-9]{6}$",
                                  datasources = conns)
} # }
```
