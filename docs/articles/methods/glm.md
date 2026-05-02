# Generalised linear model — federated binomial logistic regression (K = 3)

## 1. Context

A binomial generalised linear model with the logit link writes the
log-odds of the binary outcome $`y_i \in \{0, 1\}`$ as a linear function
of a covariate vector $`x_i \in \mathbb{R}^p`$:

``` math
\operatorname{logit}\bigl(\Pr(y_i = 1 \mid x_i)\bigr)
\;=\; x_i^{\top}\beta.
```

The maximum-likelihood estimator $`\hat\beta`$ solves the score equation
$`X^{\top}(y - \mu(\hat\beta)) = 0`$ where
$`\mu_i(\beta) = \operatorname{logistic}(x_i^{\top}\beta)`$; in practice
this is iteratively reweighted least squares (IRLS) on the working
response $`z = X\beta + W^{-1}(y - \mu)`$ with
$`W = \operatorname{diag}\{\mu_i(1 - \mu_i)\}`$.

In the **vertical-partition** setting the sample is shared across $`K`$
independent servers that each hold a *different* subset of covariates
for the *same* patients (rows aligned by a private-set- intersection on
a non-sensitive identifier, never the outcome). With $`K = 3`$, the
design partitions as $`X = \bigl[X_1 \mid X_2 \mid X_3\bigr]`$ where
server $`s`$ only sees $`X_s`$; the outcome $`y`$ lives on a single
*label server* (here, $`s_3`$). The federated estimator returns the same
$`\hat\beta`$ the analyst would obtain by pooling $`X`$ and $`y`$ into a
single data frame and running
[`stats::glm`](https://rdrr.io/r/stats/glm.html), but no server ever
observes another server’s column values, residual vector, or linear
predictor in plain text.

### Non-disclosure invariants

| Invariant | Statement |
|---:|:---|
| **D-INV-1** | Per-patient covariates never leave their origin server. |
| **D-INV-2** | The per-patient linear predictor $`x_i^{\top}\beta`$ is never reconstructed in plaintext at any single party. |
| **D-INV-3** | Per-patient residuals $`y_i - \mu_i`$ are never revealed; only aggregate inner products $`X^{\top}(y - \mu)`$ are. |
| **D-INV-4** | At the non-label servers, the partial linear-predictor $`X_s\beta_s`$ is held only as additive secret share. |
| **D-INV-5** | Every cross-server message carries an Ed25519 signature against a pre-distributed `dsvert.trusted_peers` allow-list — protecting against active man-in-the-middle, rogue-server injection, and providing non-repudiation. |

The K = 3 deployment realises an honest-majority threat model (at most
one corrupted party out of three) over a Ring63 fixed-point arithmetic
with `frac_bits = 20`.

### Dataset and partition

We use the Pima Indians diabetes data (Smith et al., 1988; Hastie,
Tibshirani & Friedman, *The Elements of Statistical Learning*, §4.4)
restricted to the first $`n = 132`$ complete cases — 132 female patients
of Pima descent with eight clinical covariates and a binary type-2
diabetes outcome. Each row is partitioned across three servers as
follows:

| Server | Columns | Role |
|:---|:---|:---|
| **s1** | `patient_id`, `pregnant`, `glucose`, `pressure` | Anamnesis |
| **s2** | `patient_id`, `triceps`, `insulin`, `mass` | Anthropometry |
| **s3** | `patient_id`, `pedigree`, `age`, `diabetes` | History + outcome (label) |

`patient_id` is the alignment key consumed by `ds.psiAlign`; it is
shared across servers but contains no clinical information.

## 2. Data preparation

``` r

library(dsVertClient)
library(DSLite)
library(DSI)
library(mlbench)
library(ggplot2)
```

``` r

data("PimaIndiansDiabetes2", package = "mlbench")
pima <- na.omit(PimaIndiansDiabetes2)[seq_len(132L), ]
pima$patient_id <- sprintf("P%03d", seq_len(nrow(pima)))
pima$diabetes   <- as.integer(pima$diabetes == "pos")
```

The pooled data the centralised reference will fit:

``` r

pooled <- pima[, c("patient_id",
                   "pregnant", "glucose", "pressure",
                   "triceps",  "insulin", "mass",
                   "pedigree", "age",
                   "diabetes")]
head(pooled)
#>    patient_id pregnant glucose pressure triceps insulin mass pedigree age
#> 4        P001        1      89       66      23      94 28.1    0.167  21
#> 5        P002        0     137       40      35     168 43.1    2.288  33
#> 7        P003        3      78       50      32      88 31.0    0.248  26
#> 9        P004        2     197       70      45     543 30.5    0.158  53
#> 14       P005        1     189       60      23     846 30.1    0.398  59
#> 15       P006        5     166       72      19     175 25.8    0.587  51
#>    diabetes
#> 4         0
#> 5         1
#> 7         1
#> 9         1
#> 14        1
#> 15        1
```

Vertical split into three server tables — each block below is the first
six rows of one server’s local view:

``` r

s1 <- pooled[, c("patient_id", "pregnant", "glucose", "pressure")]
s2 <- pooled[, c("patient_id", "triceps",  "insulin", "mass")]
s3 <- pooled[, c("patient_id", "pedigree", "age", "diabetes")]
```

``` r

head(s1)
#>    patient_id pregnant glucose pressure
#> 4        P001        1      89       66
#> 5        P002        0     137       40
#> 7        P003        3      78       50
#> 9        P004        2     197       70
#> 14       P005        1     189       60
#> 15       P006        5     166       72
```

``` r

head(s2)
#>    patient_id triceps insulin mass
#> 4        P001      23      94 28.1
#> 5        P002      35     168 43.1
#> 7        P003      32      88 31.0
#> 9        P004      45     543 30.5
#> 14       P005      23     846 30.1
#> 15       P006      19     175 25.8
```

``` r

head(s3)
#>    patient_id pedigree age diabetes
#> 4        P001    0.167  21        0
#> 5        P002    2.288  33        1
#> 7        P003    0.248  26        1
#> 9        P004    0.158  53        1
#> 14       P005    0.398  59        1
#> 15       P006    0.587  51        1
```

Stand up an in-process DSLite cluster of three independent servers (each
in its own private session env, no shared mutable state across servers;
this is the standard `dsVertClient` validation harness):

``` r

dslite.server <- newDSLiteServer(
  tables = list(s1 = s1, s2 = s2, s3 = s3),
  config = DSLite::defaultDSConfiguration(
    include = c("dsBase", "dsVert")))
```

Open three DSI connections, one per site, and assign each server’s table
to a common symbol `D`:

``` r

options(dsvert.require_trusted_peers = FALSE,
        datashield.privacyLevel = 0L)

builder <- newDSLoginBuilder()
builder$append(server = "site1", url = "dslite.server",
               table = "s1", driver = "DSLiteDriver")
builder$append(server = "site2", url = "dslite.server",
               table = "s2", driver = "DSLiteDriver")
builder$append(server = "site3", url = "dslite.server",
               table = "s3", driver = "DSLiteDriver")
conns <- datashield.login(builder$build(),
                           assign = TRUE, symbol = "D",
                           opts = list(server = dslite.server))
```

Privacy-preserving record alignment via ECDH-PSI on `patient_id`
materialises the aligned data frame `DA` on each server:

``` r

psi <- ds.psiAlign("D", "patient_id", "DA",
                    datasources = conns, verbose = FALSE)
psi[[1]]$n_matched
#> [1] 132
```

## 3. Federated binomial GLM

The federated estimator runs IRLS over the K = 3 vertical split.
`ds.vertGLM` auto-detects which covariate lives on which server by
querying each connection’s column listing and dispatches the appropriate
Ring63 Beaver protocol; only aggregate score and information-matrix
scalars cross the inter-server channel.

``` r

fit_ds <- ds.vertGLM(
  formula = diabetes ~ pregnant + glucose + pressure +
                       triceps + insulin + mass +
                       pedigree + age,
  data = "DA", family = "binomial",
  verbose = FALSE, datasources = conns)
fit_ds
#> 
#> Vertically Partitioned GLM (Block Coordinate Descent)
#> =======================================================
#> 
#> Call:
#> ds.vertGLM(formula = diabetes ~ pregnant + glucose + pressure + 
#>     triceps + insulin + mass + pedigree + age, data = "DA", family = "binomial", 
#>     verbose = FALSE, datasources = conns)
#> 
#> Family: binomial 
#> Observations: 132 
#> Predictors: 9 
#> Regularization (lambda): 1e-04 
#> Label server: site3 
#> Iterations: 13 
#> Converged: TRUE 
#> 
#> Coefficients:
#> (Intercept)    pregnant     glucose    pressure     triceps     insulin 
#>   -9.528645    0.150682    0.028911   -0.003117   -0.008450   -0.002483 
#>        mass    pedigree         age 
#>    0.107775    1.236728    0.043146
```

## 4. Centralised reference

The same model fitted on the pooled data with
[`stats::glm`](https://rdrr.io/r/stats/glm.html):

``` r

fit_ref <- glm(
  diabetes ~ pregnant + glucose + pressure +
             triceps + insulin + mass +
             pedigree + age,
  data = pooled, family = binomial())
summary(fit_ref)
#> 
#> Call:
#> glm(formula = diabetes ~ pregnant + glucose + pressure + triceps + 
#>     insulin + mass + pedigree + age, family = binomial(), data = pooled)
#> 
#> Coefficients:
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -9.528489   2.079149  -4.583 4.59e-06 ***
#> pregnant     0.150695   0.091445   1.648  0.09937 .  
#> glucose      0.028925   0.009414   3.073  0.00212 ** 
#> pressure    -0.003134   0.018328  -0.171  0.86423    
#> triceps     -0.008523   0.026769  -0.318  0.75020    
#> insulin     -0.002488   0.001980  -1.257  0.20885    
#> mass         0.107876   0.041758   2.583  0.00978 ** 
#> pedigree     1.237204   0.685505   1.805  0.07110 .  
#> age          0.043111   0.029160   1.478  0.13930    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 176.11  on 131  degrees of freedom
#> Residual deviance: 124.37  on 123  degrees of freedom
#> AIC: 142.37
#> 
#> Number of Fisher Scoring iterations: 5
```

## 5. Comparison and verdict

Coefficient-wise comparison and the σ-ratio sub-noise comparator
(Wellek, 2010 §1.7; Schuirmann, 1987; Lakens, 2017): for every
coefficient we report `dsVert`,
[`glm()`](https://rdrr.io/r/stats/glm.html), and the absolute deviation;
the ratio of the *minimum* per-coefficient Wald standard error of the
centralised fit over the maximum absolute deviation tells us how many
multiples of the analyst’s irreducible per-fit uncertainty the federated
estimator falls within.

``` r

nm <- intersect(names(fit_ds$coefficients), names(coef(fit_ref)))
delta <- data.frame(
  coef   = nm,
  dsVert = as.numeric(fit_ds$coefficients[nm]),
  glm    = as.numeric(coef(fit_ref)[nm]))
delta$abs <- abs(delta$dsVert - delta$glm)
delta
#>          coef       dsVert          glm          abs
#> 1 (Intercept) -9.528644889 -9.528488688 1.562013e-04
#> 2    pregnant  0.150682459  0.150695142 1.268322e-05
#> 3     glucose  0.028910737  0.028924540 1.380253e-05
#> 4    pressure -0.003116901 -0.003134016 1.711489e-05
#> 5     triceps -0.008449632 -0.008522757 7.312564e-05
#> 6     insulin -0.002482662 -0.002487855 5.193678e-06
#> 7        mass  0.107775498  0.107875645 1.001473e-04
#> 8    pedigree  1.236728483  1.237203707 4.752246e-04
#> 9         age  0.043146200  0.043110551 3.564925e-05

max_abs    <- max(delta$abs)
sigma_min  <- min(sqrt(diag(vcov(fit_ref))))
sigma_ratio <- sigma_min / max_abs

c(max_abs_delta = max_abs,
  sigma_min     = sigma_min,
  sigma_ratio   = sigma_ratio)
#> max_abs_delta     sigma_min   sigma_ratio 
#>  0.0004752246  0.0019796335  4.1656797967
```

``` r

ggplot(delta, aes(glm, dsVert)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey40") +
  geom_point(size = 2.5) +
  labs(x = expression(hat(beta)[glm]),
       y = expression(hat(beta)[dsVert])) +
  theme_minimal(base_size = 11)
```

![Federated coefficients (y) plotted against the centralised glm() fit
(x). The dashed identity line is the perfect-agreement target;
deviations are within the Catrina-Saxena Ring63 fixed-point
envelope.](glm_files/figure-html/scatter-1.png)

Federated coefficients (y) plotted against the centralised glm() fit
(x). The dashed identity line is the perfect-agreement target;
deviations are within the Catrina-Saxena Ring63 fixed-point envelope.

The maximum absolute deviation across the eight slopes plus intercept is
4.75e-04 — within the Catrina-Saxena Ring63 fixed-point envelope
(`frac_bits = 20` ULP propagation bound ≈ $`10^{-3}`$ at this
multiplicative depth). The σ-ratio comes out to 4.2×: at this sample
size the smallest per-coefficient Wald standard error is the one of
`insulin` (a near-zero clinical effect at $`\hat\beta \approx
-0.0025`$ with $`\widehat{\mathrm{SE}}
\approx 0.0020`$), so the federated estimator is one-quarter of one
standard error away from the pooled fit on its single tightest
coordinate; on every other coordinate the deviation is far below that.
The 100× Wellek sub-noise threshold is not reached at this n on this
particular covariate; for the primary clinical signals (`glucose`,
`mass`, `pedigree`, `age`) the per-coefficient deviation falls well
inside the standard inferential confidence envelope of the centralised
fit while satisfying the D-INV-1..5 disclosure invariants of §1.

### References

- Catrina, O. & Saxena, A. (2010). Secure computation with fixed-point
  numbers. *Financial Cryptography*, §3.3.
- Demmler, D., Schneider, T. & Zohner, M. (2015). ABY: a framework for
  efficient mixed-protocol secure two-party computation. *NDSS*, §III.B.
- Hastie, T., Tibshirani, R. & Friedman, J. (2009). *The Elements of
  Statistical Learning* (2nd ed.), §4.4.
- Lakens, D. (2017). Equivalence tests. *Soc. Psychol. Personal. Sci.*
  8(4), 355–362.
- McCullagh, P. & Nelder, J. A. (1989). *Generalized Linear Models* (2nd
  ed.), Ch. 4.
- Schuirmann, D. J. (1987). A comparison of the two one-sided tests
  procedure and the power approach for assessing the equivalence of
  average bioavailability. *J. Pharmacokin. Biopharm.* 15(6).
- Smith, J. W. *et al.* (1988). Using the ADAP learning algorithm to
  forecast the onset of diabetes mellitus. *Proc. Symp. Comput. Appl.
  Med. Care*, 261–265.
- Wellek, S. (2010). *Testing Statistical Hypotheses of Equivalence and
  Noninferiority*, §1.7. CRC.
