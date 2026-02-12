# Getting Started with dsVertClient

## Vertical Partitioning

In a vertically partitioned setting, different institutions hold
different variables for the same set of patients. A shared patient
identifier links records across sites, but no single institution has
access to the complete variable set. dsVertClient extends DataSHIELD to
run privacy-preserving statistical methods across these disjoint columns
without ever centralising the raw data.

| Server  | Institution | Variables                                     |
|---------|-------------|-----------------------------------------------|
| server1 | Hospital A  | patient_id, age, bmi                          |
| server2 | Hospital B  | patient_id, glucose, bp, hypertension, visits |
| server3 | Lab C       | patient_id, cholesterol, heart_rate           |

The goal is to analyse the combined variable space across all three
institutions without sharing any raw data.

## Connect to Opal Servers

Build a login data frame containing the URL, credentials, and table
reference for each Opal server, then call
[`datashield.login()`](https://datashield.github.io/DSI/reference/datashield.login.html)
to open connections and assign each table to a server-side symbol called
`"D"`.

``` r

library(dsVertClient)
library(DSI)
library(DSOpal)

builder <- DSI::newDSLoginBuilder()
builder$append(server = "server1", url = "https://opal1.example.org",
               table = "project.server1_data", user = "analyst",
               password = "password", driver = "OpalDriver")
builder$append(server = "server2", url = "https://opal2.example.org",
               table = "project.server2_data", user = "analyst",
               password = "password", driver = "OpalDriver")
builder$append(server = "server3", url = "https://opal3.example.org",
               table = "project.server3_data", user = "analyst",
               password = "password", driver = "OpalDriver")

connections <- datashield.login(builder$build(), assign = TRUE, symbol = "D")
```

    Logging into the collaborating servers

## Validate Identifiers

Before aligning records, confirm that the patient identifiers are
present, unique, and consistently formatted at every server.
[`ds.validateIdFormat()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.validateIdFormat.md)
performs these checks without revealing the actual identifier values.

``` r

validation <- ds.validateIdFormat("D", "patient_id", datasources = connections)
validation
```

``` r

## # Identifier Validation Summary
##
## | Server  | Valid IDs | Total | Format        | Example       |
## |---------|----------|-------|---------------|---------------|
## | server1 | 50       | 50    | PATIENT_NNNNN | PATIENT_00001 |
## | server2 | 50       | 50    | PATIENT_NNNNN | PATIENT_00001 |
## | server3 | 50       | 50    | PATIENT_NNNNN | PATIENT_00001 |
```

## Align Records

Records across servers are stored in different row orders.
[`ds.psiAlign()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md)
uses an ECDH-based Private Set Intersection (PSI) protocol to determine
the common set of patient identifiers and reorder every server so that
the same row position corresponds to the same patient. Each server masks
its identifiers with a secret elliptic-curve scalar; the client sees
only opaque curve points and cannot recover the original identifiers.

``` r

ds.psiAlign("D", "patient_id", "D_aligned", datasources = connections)
```

    Server 'server1': 50 of 50 records matched (100.0%)
    Server 'server2': 50 of 50 records matched (100.0%)
    Server 'server3': 50 of 50 records matched (100.0%)

## Quick Preview

With records aligned, compute a cross-server correlation matrix to
verify that the data is accessible and the alignment is correct. The
analyst specifies which variables reside on which server.

``` r

x_vars <- list(
  server1 = c("age", "bmi"),
  server2 = c("glucose", "bp"),
  server3 = c("cholesterol", "heart_rate")
)

cor_result <- ds.vertCor("D_aligned", x_vars, datasources = connections)
round(cor_result$correlation, 2)
```

``` r

##             age   bmi glucose    bp cholesterol heart_rate
## age        1.00 -0.32   -0.10  0.01       0.22      -0.09
## bmi       -0.32  1.00   -0.09  0.23      -0.02       0.31
## glucose   -0.10 -0.09    1.00 -0.04      -0.00       0.13
## bp         0.01  0.23   -0.04  1.00      -0.05       0.24
## cholesterol 0.22 -0.02  -0.00 -0.05       1.00      -0.24
## heart_rate -0.09  0.31    0.13  0.24      -0.24       1.00
```

## Disconnect

Always close the server connections when the analysis session is
finished.

``` r

datashield.logout(connections)
```

## Next Steps

The companion vignettes cover the full workflow in detail:

- **Statistical Analysis** – correlation, PCA, and GLM examples with
  interpretation guidance:
  [`vignette("b-statistical-analysis")`](https://isglobal-brge.github.io/dsVertClient/articles/b-statistical-analysis.md)
- **Methodology** – mathematical foundations of Block Coordinate
  Descent, MHE threshold decryption, and the security model:
  [`vignette("c-methodology")`](https://isglobal-brge.github.io/dsVertClient/articles/c-methodology.md)
- **Validation** – empirical validation of accuracy guarantees:
  [`vignette("d-validation")`](https://isglobal-brge.github.io/dsVertClient/articles/d-validation.md)
- **Security** – threat model and privacy analysis:
  [`vignette("e-security")`](https://isglobal-brge.github.io/dsVertClient/articles/e-security.md)
