# dsVert Validation: Pima Diabetes Dataset

This vignette validates the full **dsVert** pipeline against centralized
R computations on the **Pima Indian women diabetes dataset** (`Pima.tr`
from the MASS package, $`n = 200`$). The dataset is vertically
partitioned across three DataSHIELD (Opal/Rock) servers, and all six
protocols are exercised: PSI alignment, correlation, PCA, and three GLM
families (Gaussian, binomial, Poisson).

The local ground truth is computed here at build time. The dsVert
results shown below were obtained on a three-server Opal deployment and
are embedded as pre-computed values; to reproduce them, run
`demo/run_pima_demo.R` with three running Opal servers.

## 1. Dataset and vertical partitioning

``` r

library(MASS)
data("Pima.tr", package = "MASS")
df <- as.data.frame(Pima.tr)
df$patient_id <- sprintf("P%03d", seq_len(nrow(df)))
df$diabetes   <- ifelse(df$type == "Yes", 1L, 0L)

cat("Dataset:", nrow(df), "rows x", ncol(df) - 1, "numeric columns\n")
#> Dataset: 200 rows x 9 numeric columns
```

The data are vertically split across three servers. Row order is
independently shuffled per server so that PSI alignment has real work to
do.

| Server  | Role                  | Variables                         |
|---------|-----------------------|-----------------------------------|
| server1 | Anthropometric / risk | `patient_id`, `age`, `bmi`, `ped` |
| server2 | Clinical              | `patient_id`, `glu`, `bp`, `skin` |
| server3 | Outcomes              | `patient_id`, `npreg`, `diabetes` |

## 2. Local ground truth

We compute the centralized results that dsVert should approximate.

### Correlation and PCA

``` r

num_vars  <- c("age", "bmi", "ped", "glu", "bp", "skin", "npreg")
local_cor <- cor(df[, num_vars])
local_pca <- eigen(local_cor)

cat("Local correlation matrix (rounded to 4 decimals):\n")
#> Local correlation matrix (rounded to 4 decimals):
print(round(local_cor, 4))
#>           age    bmi     ped    glu      bp   skin   npreg
#> age    1.0000 0.1319 -0.0714 0.3434  0.3911 0.2519  0.5989
#> bmi    0.1319 1.0000  0.1906 0.2168  0.2388 0.6590  0.0583
#> ped   -0.0714 0.1906  1.0000 0.0607 -0.0474 0.0954 -0.1195
#> glu    0.3434 0.2168  0.0607 1.0000  0.2694 0.2176  0.1705
#> bp     0.3911 0.2388 -0.0474 0.2694  1.0000 0.2650  0.2521
#> skin   0.2519 0.6590  0.0954 0.2176  0.2650 1.0000  0.1090
#> npreg  0.5989 0.0583 -0.1195 0.1705  0.2521 0.1090  1.0000

cat("\nLocal PCA eigenvalues:\n")
#> 
#> Local PCA eigenvalues:
cat(round(local_pca$values, 6), sep = ", ")
#> 2.409261, 1.496447, 0.911868, 0.799971, 0.690135, 0.389819, 0.302499
cat("\n")
```

### GLM Gaussian (y = `glu`)

``` r

local_gauss <- glm(glu ~ age + bmi + ped + bp + skin + npreg + diabetes,
                   data = df, family = gaussian)
cat("Local Gaussian GLM coefficients:\n")
#> Local Gaussian GLM coefficients:
print(round(coef(local_gauss), 6))
#> (Intercept)         age         bmi         ped          bp        skin 
#>   68.961590    0.558322    0.233011   -2.369034    0.317864    0.070458 
#>       npreg    diabetes 
#>   -0.842445   26.299284
```

### GLM Binomial (y = `diabetes`)

``` r

local_binom <- glm(diabetes ~ age + bmi + ped + glu + bp + skin + npreg,
                   data = df, family = binomial)
cat("Local Binomial GLM coefficients:\n")
#> Local Binomial GLM coefficients:
print(round(coef(local_binom), 6))
#> (Intercept)         age         bmi         ped         glu          bp 
#>   -9.773062    0.041184    0.083624    1.820410    0.032117   -0.004768 
#>        skin       npreg 
#>   -0.001917    0.103183
```

### GLM Poisson (y = `npreg`)

``` r

local_pois <- glm(npreg ~ age + bmi + ped + glu + bp + skin + diabetes,
                  data = df, family = poisson)
cat("Local Poisson GLM coefficients:\n")
#> Local Poisson GLM coefficients:
print(round(coef(local_pois), 6))
#> (Intercept)         age         bmi         ped         glu          bp 
#>   -0.026629    0.039114    0.002999   -0.378371   -0.002727    0.004165 
#>        skin    diabetes 
#>   -0.002472    0.286849
```

## 3. dsVert pipeline (DataSHIELD code)

The following code blocks show the dsVert calls used on the three-server
Opal deployment. They are not executed during vignette build (they
require live servers); the results are shown in Section 4.

### Connect to DataSHIELD

``` r

library(DSI)
library(DSOpal)
library(dsVertClient)

builder <- DSI::newDSLoginBuilder()
builder$append(server = "server1", url = "https://opal1:8443",
               user = "administrator", password = "password",
               driver = "OpalDriver", table = "project.data")
builder$append(server = "server2", url = "https://opal2:8443",
               user = "administrator", password = "password",
               driver = "OpalDriver", table = "project.data")
builder$append(server = "server3", url = "https://opal3:8443",
               user = "administrator", password = "password",
               driver = "OpalDriver", table = "project.data")

conns <- DSI::datashield.login(builder$build(), assign = TRUE, symbol = "D")
```

### A) PSI alignment

``` r

psi_result <- ds.psiAlign(
  data_name   = "D",
  id_col      = "patient_id",
  newobj      = "D_aligned",
  datasources = conns
)
```

### B) Correlation

``` r

variables <- list(
  server1 = c("age", "bmi", "ped"),
  server2 = c("glu", "bp", "skin"),
  server3 = c("npreg")
)

cor_result <- ds.vertCor(
  data_name   = "D_aligned",
  variables   = variables,
  log_n       = 12,
  log_scale   = 40,
  datasources = conns
)
```

### C) PCA

``` r

pca_result <- ds.vertPCA(cor_result = cor_result, n_components = 7)
```

### D) GLM Gaussian

``` r

glm_gauss <- ds.vertGLM(
  data_name = "D_aligned", y_var = "glu",
  x_vars    = list(server1 = c("age", "bmi", "ped"),
                   server2 = c("bp", "skin"),
                   server3 = c("npreg", "diabetes")),
  y_server  = "server2", family = "gaussian",
  lambda = 1e-4, datasources = conns
)
```

### E) GLM Binomial

``` r

glm_binom <- ds.vertGLM(
  data_name = "D_aligned", y_var = "diabetes",
  x_vars    = list(server1 = c("age", "bmi", "ped"),
                   server2 = c("glu", "bp", "skin"),
                   server3 = c("npreg")),
  y_server  = "server3", family = "binomial",
  lambda = 1e-4, datasources = conns
)
```

### F) GLM Poisson

``` r

glm_pois <- ds.vertGLM(
  data_name = "D_aligned", y_var = "npreg",
  x_vars    = list(server1 = c("age", "bmi", "ped"),
                   server2 = c("glu", "bp", "skin"),
                   server3 = c("diabetes")),
  y_server  = "server3", family = "poisson",
  lambda = 1e-4, datasources = conns
)
```

## 4. Validation: dsVert vs. local R

The results below were obtained on a three-server Opal deployment (Opal
4.8 + Rock, Docker, macOS host) running dsVert 2.0.0 and dsVertClient
2.0.0.

### A) PSI alignment

All 200 patients were matched across all three servers (100% match
rate). This confirms that the ECDH-PSI protocol correctly identifies the
intersection despite independent row shuffling.

``` r

cat("PSI result: 200 / 200 matched (100%)\n")
#> PSI result: 200 / 200 matched (100%)
```

### B) Correlation

``` r

# Pre-computed dsVert correlation max absolute error
dsvert_cor_max_error <- 7.37e-06

cat(sprintf("Max |r_dsVert - r_local|: %.2e\n", dsvert_cor_max_error))
#> Max |r_dsVert - r_local|: 7.37e-06
cat("All 21 unique off-diagonal entries match to at least 5 decimal places.\n")
#> All 21 unique off-diagonal entries match to at least 5 decimal places.
```

### C) PCA

``` r

dsvert_eig_max_error <- 3.18e-06

cat("Local eigenvalues: ", paste(round(local_pca$values, 4), collapse = ", "), "\n")
#> Local eigenvalues:  2.4093, 1.4964, 0.9119, 0.8, 0.6901, 0.3898, 0.3025
cat(sprintf("Max |eigenvalue_dsVert - eigenvalue_local|: %.2e\n", dsvert_eig_max_error))
#> Max |eigenvalue_dsVert - eigenvalue_local|: 3.18e-06
```

### D) GLM coefficient comparison

``` r

# Pre-computed dsVert max coefficient errors (L-infinity norm)
dsvert_gauss_err <- 6.23e-03
dsvert_binom_err <- 4.02e-04
dsvert_pois_err  <- 1.84e-04

comparison <- data.frame(
  Family   = c("Gaussian (y=glu)", "Binomial (y=diabetes)", "Poisson (y=npreg)"),
  Local_n_coefs = c(length(coef(local_gauss)),
                    length(coef(local_binom)),
                    length(coef(local_pois))),
  Max_coef_error = c(dsvert_gauss_err, dsvert_binom_err, dsvert_pois_err),
  stringsAsFactors = FALSE
)

knitr::kable(comparison,
             col.names = c("Family", "Coefficients", "Max |dsVert - local|"),
             align = c("l", "c", "r"),
             caption = "GLM coefficient accuracy (dsVert vs. centralized R glm())")
```

| Family                | Coefficients | Max \|dsVert - local\| |
|:----------------------|:------------:|-----------------------:|
| Gaussian (y=glu)      |      8       |               0.006230 |
| Binomial (y=diabetes) |      8       |               0.000402 |
| Poisson (y=npreg)     |      8       |               0.000184 |

GLM coefficient accuracy (dsVert vs. centralized R glm()) {.table}

### E) Overall validation summary

``` r

summary_df <- data.frame(
  Protocol = c("PSI Alignment", "Correlation", "PCA",
               "GLM Gaussian", "GLM Binomial", "GLM Poisson"),
  Metric   = c("Matched cohort", "Max |r| error", "Max eigenvalue error",
               "Max coef error", "Max coef error", "Max coef error"),
  Value    = c("200/200 (100%)",
               sprintf("%.2e", dsvert_cor_max_error),
               sprintf("%.2e", dsvert_eig_max_error),
               sprintf("%.2e", dsvert_gauss_err),
               sprintf("%.2e", dsvert_binom_err),
               sprintf("%.2e", dsvert_pois_err)),
  stringsAsFactors = FALSE
)

knitr::kable(summary_df,
             col.names = c("Protocol", "Metric", "Value"),
             align = c("l", "l", "r"),
             caption = "Validation summary: dsVert vs. centralized R (Pima dataset, K=3, n=200)")
```

| Protocol      | Metric               |          Value |
|:--------------|:---------------------|---------------:|
| PSI Alignment | Matched cohort       | 200/200 (100%) |
| Correlation   | Max \|r\| error      |       7.37e-06 |
| PCA           | Max eigenvalue error |       3.18e-06 |
| GLM Gaussian  | Max coef error       |       6.23e-03 |
| GLM Binomial  | Max coef error       |       4.02e-04 |
| GLM Poisson   | Max coef error       |       1.84e-04 |

Validation summary: dsVert vs. centralized R (Pima dataset, K=3, n=200)
{.table}

All protocols produce results that are numerically close to centralized
computation. The small errors in correlation and PCA ($`< 10^{-5}`$)
arise from CKKS approximate arithmetic. The slightly larger GLM errors
($`< 10^{-2}`$) reflect both CKKS approximation in encrypted gradient
computation and the ridge regularization ($`\lambda = 10^{-4}`$) used
for BCD convergence.

## 5. Reproducing these results

To reproduce, you need three Opal/Rock servers with dsVert installed:

1.  Prepare the CSV files and upload to each Opal server
2.  Install dsVert on each Rock server
3.  Install dsVertClient locally
4.  Run `demo/run_pima_demo.R`

See the [Getting Started
guide](https://isglobal-brge.github.io/dsVertClient/) for detailed setup
instructions.
