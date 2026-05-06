# PSI alignment validation (K\>=3)

## Scope

Functions:
[`ds.vert.align()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md).

This vignette validates the `K>=3` modality. It creates the local pooled
fixture, partitions it vertically, starts an in-memory DSLite server,
aligns records with PSI, runs the distributed dsVert route, runs the
centralized R reference on the same pooled fixture, and computes the
numerical delta.

No external RDS, cache, or pre-loaded server table is used. Rendering
the vignette executes the workflow from dataset generation to assertion.

## Method

ECDH-PSI aligns vertically split tables on a shared identifier. Each
server masks identifiers, the protocol computes the common set, and the
aligned object keeps only matched rows in a common order.

## Mathematical target

The target is the set intersection I = intersection_k I_k. Validation
checks matched counts and that deterministic id-derived columns align
perfectly after PSI.

## Disclosure review

The analyst receives counts and status fields only. Patient IDs and
matched row-index vectors are not returned; the legacy matched-index
aggregate is checked as blocked.

The DSLite validation uses `datashield.privacyLevel = 5`. Trusted-peer
pinning is disabled only because DSLite is an in-memory test backend,
not an Opal/Rock deployment. That does not weaken the method-level check
of what the product route returns to the analyst.

## Executed validation

``` r

K <- 3L

contains_ids <- function(x) {
  if (is.character(x)) return(any(grepl("^P[0-9]{4}$", x)))
  if (is.list(x)) return(any(vapply(x, contains_ids, logical(1))))
  FALSE
}

make_table <- function(ids, server_index) {
  id_num <- as.integer(sub("^P", "", ids))
  ord <- sample(seq_along(ids))
  data.frame(patient_id = ids[ord], id_num = id_num[ord],
             value = id_num[ord] + server_index / 10)
}

set.seed(9100 + K)
id_sets <- if (K == 2L) {
  list(s1 = sprintf("P%04d", 1:80),
       s2 = sprintf("P%04d", 11:90))
} else {
  list(s1 = sprintf("P%04d", 1:85),
       s2 = sprintf("P%04d", 11:95),
       s3 = sprintf("P%04d", 6:80))
}
common_ids <- Reduce(intersect, id_sets)
tables <- Map(make_table, id_sets, seq_along(id_sets))

knitr::kable(partition_summary(tables))
```

|     | server |   n | columns                   |
|:----|:-------|----:|:--------------------------|
| s1  | s1     |  85 | patient_id, id_num, value |
| s2  | s2     |  85 | patient_id, id_num, value |
| s3  | s3     |  75 | patient_id, id_num, value |

``` r


result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi <- psi_align(conns)
    counts <- DSI::datashield.aggregate(
      conns, call(name = "getObsCountDS", data_name = "DA"))
    n_by_server <- vapply(counts, function(x) as.integer(x$n_obs), integer(1))
    variables <- stats::setNames(lapply(names(conns), function(.x) "id_num"),
                                 names(conns))
    cor_fit <- dsVertClient::ds.vert.cor(
      "DA", variables = variables,
      verbose = validation_demo_verbose(), datasources = conns)
    legacy_blocked <- tryCatch({
      DSI::datashield.aggregate(conns[1],
        call(name = "psiGetMatchedIndicesDS"))
      FALSE
    }, error = function(e) TRUE)
    list(psi = psi, n_by_server = n_by_server, cor_fit = cor_fit,
         legacy_blocked = legacy_blocked)
  }, finally = try(DSI::datashield.logout(conns), silent = TRUE))
})
#> 
#> Logging into the collaborating servers
#> [NA handling] Per-server na.omit (before PSI)...
#> === ECDH-PSI Record Alignment (Blind Relay) ===
#> Reference: s1, Targets: s2, s3
#> [Phase 0] Transport key exchange...
#>   Transport keys exchanged.
#> [Phase 1] Reference server masking IDs...
#>   s1: 85 IDs masked (stored server-side)
#> [Phase 2] s1: exporting encrypted points for s2...
#> [Phase 2] s1: exporting encrypted points for s3...
#> [Phase 3] s2, s3: processing (parallel)...
#>   s2: 85 IDs masked
#>   s3: 75 IDs masked
#> [Phase 4-5] s2: double-masking via s1 (blind relay)...
#> [Phase 6-7] s2: matching and aligning (blind relay)...
#> [Phase 4-5] s3: double-masking via s1 (blind relay)...
#> [Phase 6-7] s3: matching and aligning (blind relay)...
#> [Phase 7] s1: self-aligning...
#> [Phase 8] Computing encrypted multi-server intersection...
#>   Common records: 70
#> Server 's1': 70 of 85 records matched (82.4%)
#> Server 's2': 70 of 85 records matched (82.4%)
#> Server 's3': 70 of 75 records matched (93.3%)
#> PSI alignment complete.
#> === Ring63 Correlation (p=3, 3 servers) ===
#> [Phase 0] Transport keys...
#> [Phase 1] Standardizing...
#> [Phase 2] Local correlations...
#> [Phase 3] Input sharing...
#> [Phase 4] Beaver cross-products (3 columns)...
#>   Column 1/3 done
#>   Column 2/3 done
#>   Column 3/3 done
#> Correlation complete (1s)

offdiag <- result$value$cor_fit$correlation[
  upper.tri(result$value$cor_fit$correlation)]
observed <- max(abs(result$value$n_by_server - length(common_ids)),
                abs((result$value$psi$n_common %||%
                       result$value$psi[[1]]$n_matched) - length(common_ids)),
                max(abs(offdiag - 1)))
audit_ok <- !contains_ids(result$value$psi) &&
  isTRUE(result$value$legacy_blocked)

rows <- row_result(
  "psi", "PSI alignment", K, "ds.vert.align",
  "synthetic ID intersection", "deterministic set intersection",
  "max(count_delta, correlation_delta)", observed, 1e-8,
  "strict-precise",
  "Returns match counts/status only; matched IDs and row indices are not returned.",
  result$runtime_s,
  if (audit_ok) "index reveal blocked" else "audit warning")

display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K\>=3 | ds.vert.align | synthetic ID intersection | deterministic set intersection | max(count_delta, correlation_delta) | 0 | 0 | strict-precise | PASS | 2.9 |

``` r

assert_validation(rows)
```

## Verdict

Rendering fails if the distributed result leaves the accepted numerical
envelope or if the method row is marked disclosive.
