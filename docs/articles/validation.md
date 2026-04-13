# dsVert Validation: Federated Vertical vs Centralised

This vignette compares dsVert’s privacy-preserving federated vertical
analysis against standard centralised R, using two public clinical
datasets hosted on `opal-demo.obiba.org`.

The data partitions were created following the reproducible procedure in
[data-preparation.Rmd](https://isglobal-brge.github.io/dsVertClient/articles/data-preparation.md).
A pruned R script of this vignette is available for download:
[validation.R](https://github.com/isglobal-brge/dsVertClient/blob/main/inst/scripts/validation.R)

``` r

library(DSI)
library(DSOpal)
library(dsVertClient)
library(MASS)
library(survival)
```

## 1. Centralised ground truth

We reconstruct the exact same datasets locally to serve as ground truth.
The analyst can run these locally without any server access.

### birthwt (K=2)

``` r

data(birthwt)
df_bw <- birthwt
df_bw$patient_id <- sprintf("BW%03d", seq_len(nrow(df_bw)))

set.seed(123)
core <- sample(df_bw$patient_id, 160)
remaining <- setdiff(df_bw$patient_id, core)
ids1 <- sample(c(core, sample(remaining, 15)))
ids2 <- sample(c(core, sample(remaining, 10)))

bw_a <- df_bw[df_bw$patient_id %in% ids1, c("patient_id", "age", "lwt", "race", "smoke")]
bw_b <- df_bw[df_bw$patient_id %in% ids2, c("patient_id", "ptl", "ht", "ui", "ftv", "bwt", "low")]
bw_merged <- merge(bw_a, bw_b, by = "patient_id")

nrow(bw_merged)
```

    ## [1] 167

Centralised GLMs:

``` r

bw_formula <- ~ age + lwt + race + smoke + ptl + ht + ui

gt_bw_gaussian <- glm(update(bw_formula, bwt ~ .), data = bw_merged, family = gaussian())
gt_bw_binomial <- glm(update(bw_formula, low ~ .), data = bw_merged, family = binomial())
gt_bw_poisson  <- glm(update(bw_formula, ftv ~ .), data = bw_merged, family = poisson())
```

### CGD (K=3)

``` r

data(cgd)
cgd_one <- cgd[cgd$enum == 1, ]
cgd_one$n_infections <- as.integer(tapply(cgd$status, cgd$id, sum)[as.character(cgd_one$id)])
cgd_one$patient_id <- sprintf("CGD%03d", seq_len(nrow(cgd_one)))
cgd_one$treat_num <- as.integer(cgd_one$treat == "rIFN-g")
cgd_one$sex_num <- as.integer(cgd_one$sex == "male")

nrow(cgd_one)
```

    ## [1] 128

Centralised GLMs:

``` r

gt_cgd_gaussian <- glm(weight ~ age + sex_num + height + treat_num + steroids, data = cgd_one, family = gaussian())
gt_cgd_binomial <- glm(status ~ age + sex_num + height + weight + treat_num + steroids, data = cgd_one, family = binomial())
gt_cgd_poisson  <- glm(n_infections ~ age + sex_num + height + weight + treat_num + steroids, data = cgd_one, family = poisson())
```

## 2. Federated vertical analysis (dsVert)

We connect to `opal-demo.obiba.org`. Each dataset partition is stored as
a separate table; opening K connections to K different tables simulates
K independent hospitals.

### birthwt (K=2)

``` r

opal_url <- "https://opal-demo.obiba.org"

builder <- DSI::newDSLoginBuilder()
builder$append(server = "hospitalA", url = opal_url,
               user = "administrator", password = "password",
               table = "vert_demo.bw_hospital_a")
builder$append(server = "hospitalB", url = opal_url,
               user = "administrator", password = "password",
               table = "vert_demo.bw_hospital_b")
bw_conns <- DSI::datashield.login(logins = builder$build())
```

    ## 
    ## Logging into the collaborating servers

``` r

DSI::datashield.assign.table(bw_conns, "D",
  list(hospitalA = "vert_demo.bw_hospital_a",
       hospitalB = "vert_demo.bw_hospital_b"))
```

Align records via private set intersection:

``` r

bw_psi <- ds.psiAlign("D", "patient_id", "DA", verbose = FALSE, datasources = bw_conns)
bw_psi$n_common
```

    ## [1] 167

Fit all three GLM families using the formula interface:

``` r

dv_bw_gaussian <- ds.vertGLM(bwt ~ age + lwt + race + smoke + ptl + ht + ui,
                              data = "DA", family = "gaussian",
                              verbose = FALSE, datasources = bw_conns)

dv_bw_binomial <- ds.vertGLM(low ~ age + lwt + race + smoke + ptl + ht + ui,
                              data = "DA", family = "binomial",
                              verbose = FALSE, datasources = bw_conns)

dv_bw_poisson  <- ds.vertGLM(ftv ~ age + lwt + race + smoke + ptl + ht + ui,
                              data = "DA", family = "poisson",
                              verbose = FALSE, datasources = bw_conns)
```

Correlation and PCA:

``` r

bw_vars <- c("age", "lwt", "race", "smoke", "ptl", "ht", "ui")
dv_bw_cor <- ds.vertCor("DA", bw_vars, verbose = FALSE, datasources = bw_conns)
dv_bw_pca <- ds.vertPCA(cor_result = dv_bw_cor)
```

    ## Using provided correlation matrix (skipping Ring63 protocol)...

    ## Performing PCA via eigen decomposition...

    ## PCA complete: 7 components extracted.

``` r

DSI::datashield.logout(bw_conns)
```

### CGD (K=3)

``` r

builder <- DSI::newDSLoginBuilder()
builder$append(server = "demographics", url = opal_url,
               user = "administrator", password = "password",
               table = "vert_demo.cgd_demographics")
builder$append(server = "treatment", url = opal_url,
               user = "administrator", password = "password",
               table = "vert_demo.cgd_treatment")
builder$append(server = "outcomes", url = opal_url,
               user = "administrator", password = "password",
               table = "vert_demo.cgd_outcomes")
cgd_conns <- DSI::datashield.login(logins = builder$build())
```

    ## 
    ## Logging into the collaborating servers

``` r

DSI::datashield.assign.table(cgd_conns, "D",
  list(demographics = "vert_demo.cgd_demographics",
       treatment    = "vert_demo.cgd_treatment",
       outcomes     = "vert_demo.cgd_outcomes"))
```

``` r

cgd_psi <- ds.psiAlign("D", "patient_id", "DA", verbose = FALSE, datasources = cgd_conns)
cgd_psi$n_common
```

    ## [1] 128

``` r

dv_cgd_gaussian <- ds.vertGLM(weight ~ age + sex_num + height + treat_num + steroids,
                               data = "DA", family = "gaussian",
                               verbose = FALSE, datasources = cgd_conns)

dv_cgd_binomial <- ds.vertGLM(status ~ age + sex_num + height + weight + treat_num + steroids,
                               data = "DA", family = "binomial",
                               verbose = FALSE, datasources = cgd_conns)

dv_cgd_poisson  <- ds.vertGLM(n_infections ~ age + sex_num + height + weight + treat_num + steroids,
                               data = "DA", family = "poisson",
                               verbose = FALSE, datasources = cgd_conns)
```

``` r

cgd_vars <- c("age", "sex_num", "height", "weight", "treat_num", "steroids")
dv_cgd_cor <- ds.vertCor("DA", cgd_vars, verbose = FALSE, datasources = cgd_conns)
dv_cgd_pca <- ds.vertPCA(cor_result = dv_cgd_cor)
```

    ## Using provided correlation matrix (skipping Ring63 protocol)...

    ## Performing PCA via eigen decomposition...

    ## PCA complete: 6 components extracted.

``` r

DSI::datashield.logout(cgd_conns)
```

## 3. Comparison: centralised vs federated vertical

### Coefficient comparison

``` r

compare_glm <- function(gt, dv, label) {
  nms <- names(dv$coefficients)
  data.frame(
    variable    = nms,
    centralised = round(coef(gt)[nms], 4),
    dsVert      = round(dv$coefficients[nms], 4),
    delta       = round(abs(coef(gt)[nms] - dv$coefficients[nms]), 4),
    SE_central  = round(summary(gt)$coefficients[nms, "Std. Error"], 4),
    SE_dsVert   = round(dv$std_errors[nms], 4),
    row.names   = NULL
  )
}
```

**birthwt Gaussian (birth weight):**

``` r

compare_glm(gt_bw_gaussian, dv_bw_gaussian, "birthwt Gaussian")
```

    ##      variable centralised    dsVert  delta SE_central SE_dsVert
    ## 1 (Intercept)   3128.6096 3128.4739 0.1357   354.6522  358.9306
    ## 2         age     -0.7504   -0.7490 0.0014     9.8192   11.2597
    ## 3         lwt      3.8212    3.8214 0.0002     1.7388    1.8760
    ## 4        race   -212.5562 -212.5182 0.0380    60.0995   64.7204
    ## 5       smoke   -387.3030 -387.2934 0.0096   114.7653  123.1826
    ## 6         ptl    -31.9994  -32.0104 0.0110   104.6260  119.5143
    ## 7          ht   -573.8107 -573.8498 0.0391   229.6754  276.5613
    ## 8          ui   -510.4663 -510.4127 0.0535   145.9039  171.8521

**birthwt Binomial (low birth weight):**

``` r

compare_glm(gt_bw_binomial, dv_bw_binomial, "birthwt Binomial")
```

    ##      variable centralised  dsVert delta SE_central SE_dsVert
    ## 1 (Intercept)      0.0487  0.0493 6e-04     1.3796    1.1495
    ## 2         age     -0.0258 -0.0259 0e+00     0.0384    0.0342
    ## 3         lwt     -0.0169 -0.0169 0e+00     0.0073    0.0059
    ## 4        race      0.5564  0.5562 1e-04     0.2362    0.1934
    ## 5       smoke      1.1767  1.1767 0e+00     0.4459    0.4275
    ## 6         ptl      0.4423  0.4430 7e-04     0.3583    0.4484
    ## 7          ht      1.8431  1.8435 4e-04     0.8114    0.9664
    ## 8          ui      0.6482  0.6488 6e-04     0.4913    0.5932

**birthwt Poisson (physician visits):**

``` r

compare_glm(gt_bw_poisson, dv_bw_poisson, "birthwt Poisson")
```

    ##      variable centralised  dsVert delta SE_central SE_dsVert
    ## 1 (Intercept)     -1.6361 -1.6366 5e-04     0.5957    0.6158
    ## 2         age      0.0433  0.0433 0e+00     0.0156    0.0186
    ## 3         lwt      0.0040  0.0040 0e+00     0.0028    0.0034
    ## 4        race     -0.0571 -0.0570 0e+00     0.1043    0.1098
    ## 5       smoke     -0.0714 -0.0715 1e-04     0.1965    0.1741
    ## 6         ptl     -0.0290 -0.0291 0e+00     0.1951    0.1707
    ## 7          ht     -0.4909 -0.4911 2e-04     0.4686    0.4259
    ## 8          ui     -0.1520 -0.1521 0e+00     0.2790    0.2609

**CGD Gaussian (body weight):**

``` r

compare_glm(gt_cgd_gaussian, dv_cgd_gaussian, "CGD Gaussian")
```

    ##      variable centralised   dsVert  delta SE_central SE_dsVert
    ## 1 (Intercept)    -35.4227 -35.4170 0.0057     4.6350   12.3294
    ## 2         age      0.9035   0.9036 0.0001     0.1319    0.3508
    ## 3     sex_num      5.5956   5.5912 0.0044     2.0177    5.3675
    ## 4      height      0.4242   0.4242 0.0000     0.0410    0.1091
    ## 5   treat_num     -2.4607  -2.4611 0.0004     1.4819    3.9421
    ## 6    steroids     -0.0591  -0.0607 0.0017     5.0197   13.3531

**CGD Binomial (serious infection):**

``` r

compare_glm(gt_cgd_binomial, dv_cgd_binomial, "CGD Binomial")
```

    ##      variable centralised  dsVert delta SE_central SE_dsVert
    ## 1 (Intercept)     -0.4702 -0.4696 6e-04     1.5501    1.5476
    ## 2         age     -0.0826 -0.0823 2e-04     0.0507    0.0507
    ## 3     sex_num      0.2083  0.2087 4e-04     0.5715    0.5699
    ## 4      height      0.0042  0.0042 0e+00     0.0155    0.0155
    ## 5      weight      0.0173  0.0172 1e-04     0.0250    0.0250
    ## 6   treat_num     -1.1495 -1.1499 4e-04     0.4098    0.4104
    ## 7    steroids      1.8046  1.8046 1e-04     1.3111    1.3111

**CGD Poisson (infection count):**

``` r

compare_glm(gt_cgd_poisson, dv_cgd_poisson, "CGD Poisson")
```

    ##      variable centralised  dsVert  delta SE_central SE_dsVert
    ## 1 (Intercept)     -0.8343 -0.7849 0.0495     0.9897    0.9671
    ## 2         age     -0.0906 -0.0880 0.0026     0.0340    0.0322
    ## 3     sex_num      0.3062  0.2951 0.0111     0.3559    0.3428
    ## 4      height      0.0077  0.0072 0.0005     0.0101    0.0099
    ## 5      weight      0.0137  0.0137 0.0001     0.0146    0.0146
    ## 6   treat_num     -1.0864 -1.0757 0.0108     0.2688    0.2524
    ## 7    steroids      1.4986  1.4738 0.0248     0.5824    0.5647

### Deviance comparison

``` r

deviance_table <- data.frame(
  dataset = rep(c("birthwt", "CGD"), each = 3),
  K       = rep(c(2, 3), each = 3),
  family  = rep(c("Gaussian", "Binomial", "Poisson"), 2),
  centralised = c(
    deviance(gt_bw_gaussian), deviance(gt_bw_binomial), deviance(gt_bw_poisson),
    deviance(gt_cgd_gaussian), deviance(gt_cgd_binomial), deviance(gt_cgd_poisson)),
  dsVert = c(
    dv_bw_gaussian$deviance, dv_bw_binomial$deviance, dv_bw_poisson$deviance,
    dv_cgd_gaussian$deviance, dv_cgd_binomial$deviance, dv_cgd_poisson$deviance)
)
deviance_table$error_pct <- round(
  100 * abs(deviance_table$dsVert - deviance_table$centralised) / deviance_table$centralised, 2)
deviance_table$centralised <- round(deviance_table$centralised, 2)
deviance_table$dsVert <- round(deviance_table$dsVert, 2)
deviance_table
```

    ##   dataset K   family centralised      dsVert error_pct
    ## 1 birthwt 2 Gaussian 67414661.46 67423291.54      0.01
    ## 2 birthwt 2 Binomial      174.86      175.06      0.12
    ## 3 birthwt 2  Poisson      216.31      216.71      0.19
    ## 4     CGD 3 Gaussian     8532.95     8539.31      0.07
    ## 5     CGD 3 Binomial      151.06      151.23      0.11
    ## 6     CGD 3  Poisson      161.29      160.95      0.21

### Correlation comparison

**birthwt:**

``` r

gt_bw_cor <- cor(bw_merged[, bw_vars])
max(abs(dv_bw_cor$correlation - gt_bw_cor))
```

    ## [1] 6.101472e-06

**CGD:**

``` r

gt_cgd_cor <- cor(cgd_one[, cgd_vars])
max(abs(dv_cgd_cor$correlation - gt_cgd_cor))
```

    ## [1] 4.196151e-06

### PCA comparison

**birthwt eigenvalues:**

``` r

gt_bw_pca <- eigen(gt_bw_cor)
data.frame(
  component   = paste0("PC", 1:length(gt_bw_pca$values)),
  centralised = round(gt_bw_pca$values, 4),
  dsVert      = round(dv_bw_pca$eigenvalues, 4),
  delta       = round(abs(gt_bw_pca$values - dv_bw_pca$eigenvalues), 6)
)
```

    ##     component centralised dsVert delta
    ## PC1       PC1      1.5642 1.5642 1e-06
    ## PC2       PC2      1.4202 1.4202 0e+00
    ## PC3       PC3      1.0477 1.0477 3e-06
    ## PC4       PC4      0.9905 0.9905 1e-06
    ## PC5       PC5      0.7704 0.7704 1e-06
    ## PC6       PC6      0.6582 0.6582 5e-06
    ## PC7       PC7      0.5488 0.5488 0e+00

**CGD eigenvalues:**

``` r

gt_cgd_pca <- eigen(gt_cgd_cor)
data.frame(
  component   = paste0("PC", 1:length(gt_cgd_pca$values)),
  centralised = round(gt_cgd_pca$values, 4),
  dsVert      = round(dv_cgd_pca$eigenvalues, 4),
  delta       = round(abs(gt_cgd_pca$values - dv_cgd_pca$eigenvalues), 6)
)
```

    ##     component centralised dsVert delta
    ## PC1       PC1      2.7298 2.7298 0e+00
    ## PC2       PC2      1.1988 1.1988 1e-06
    ## PC3       PC3      1.0081 1.0081 0e+00
    ## PC4       PC4      0.7909 0.7909 2e-06
    ## PC5       PC5      0.1851 0.1851 1e-06
    ## PC6       PC6      0.0872 0.0872 0e+00

## 4. Privacy summary

Throughout all analyses above, the analyst (this R session) received
only:

- **p-dimensional gradient sums** (aggregates over all n patients, per
  iteration)
- **p x p correlation scalars** (pairwise sums over n)
- **1 deviance scalar** (sum over n)
- **p standard errors** (derived from aggregate Hessian)

No individual patient values were disclosed.
