# PSI alignment validation

## What is validated

Functions:
[`ds.psiAlign()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md),
[`ds.isPsiAligned()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.isPsiAligned.md),
[`ds.getIdentityPks()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.getIdentityPks.md)

ECDH-PSI uses the commutativity of elliptic-curve scalar multiplication.
For an identifier hash point $`H(id)`$, the reference server computes
$`\alpha H(id)`$, the target server computes $`\beta H(id)`$, and both
sides can compare $`\alpha\beta H(id)`$ without revealing $`id`$ or
either scalar. The multi-server intersection keeps only rows present at
every site.

## Centralized reference

The centralized reference is built from the same deterministic rows and
formula as the DSLite run. The long method harness stores the exact
seed, split and reference object in the cache listed below.

``` r


ids <- list(
  s1 = sprintf("P%03d", 1:80),
  s2 = sprintf("P%03d", 11:90),
  s3 = sprintf("P%03d", c(11:80, 101:105))
)
Reduce(intersect, ids)
```

## Vertical DSLite split

The validation split creates independent server tables with the same
`patient_id` key and then aligns them with
[`ds.psiAlign()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md).

``` r


tables <- list(
  s1 = data.frame(patient_id = ids$s1, x1 = rnorm(length(ids$s1))),
  s2 = data.frame(patient_id = ids$s2, x2 = rnorm(length(ids$s2))),
  s3 = data.frame(patient_id = ids$s3, y = rnorm(length(ids$s3)))
)
```

``` r


server <- DSLite::newDSLiteServer(
  tables = tables,
  config = DSLite::defaultDSConfiguration(include = c("dsBase", "dsVert")))
builder <- DSI::newDSLoginBuilder()
for (nm in names(tables)) {
  builder$append(server = nm, url = "server", table = nm,
                 driver = "DSLiteDriver")
}
conns <- DSI::datashield.login(builder$build(), assign = TRUE,
                               symbol = "D", opts = list(server = server))
psi <- dsVertClient::ds.psiAlign("D", "patient_id", "DA",
                                 datasources = conns, verbose = FALSE)
```

To reproduce the cache from the repository root:

``` r
Rscript scripts/validate_method_psi.R
```

## Disclosure review

The client sees aggregate counts and opaque encrypted blobs. Raw
identifiers, matched row-index vectors and per-row alignment positions
are not returned. Pinned peers use signed transport keys so a server can
reject untrusted peers before encrypted PSI payloads are accepted.

## Executed evidence check

This chunk is evaluated when the vignette renders. It fails if either K
mode is not marked non-disclosive, is not `PASS`, or exceeds its
accepted tolerance.

``` r

rows <- validation_rows("psi")
assert_validation(rows)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | cache |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|:---|
| K=2 | ds.psiAlign | synthetic ID intersection | deterministic set intersection | intersection_count_delta | 0 | 0 | strict-precise | PASS | psi_dslite_20260503-003316.rds |
| K\>=3 | ds.psiAlign | synthetic ID intersection | deterministic set intersection | intersection_count_delta | 0 | 0 | strict-precise | PASS | psi_dslite_20260503-003316.rds |

## Verdict

Both K=2 and K\>=3 validation rows are inside their accepted numerical
envelope and use the current non-disclosive product route. Any legacy
route mentioned in the package is excluded from this evidence path.
