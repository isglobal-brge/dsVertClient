Loading required package: progress
Loading required package: R6
Loading required package: opalr
Loading required package: httr
# dsVert Validation Results

Server: opal-demo.obiba.org
Date: 2026-04-13 16:47 UTC 

## birthwt (K=2)


Logging into the collaborating servers
PSI: n = 167 

### Gaussian 

| Variable | Centralised | dsVert | |Delta| |
|----------|------------|--------|--------|
| (Intercept) | 3128.6096 (SE=354.6522) | 3128.4750 (SE=358.9308) | 0.1346 |
| age | -0.7504 (SE=9.8192) | -0.7490 (SE=11.2597) | 0.0014 |
| lwt | 3.8212 (SE=1.7388) | 3.8214 (SE=1.8760) | 0.0002 |
| race | -212.5562 (SE=60.0995) | -212.5181 (SE=64.7205) | 0.0381 |
| smoke | -387.3030 (SE=114.7653) | -387.2937 (SE=123.1826) | 0.0093 |
| ptl | -31.9994 (SE=104.6260) | -32.0110 (SE=119.5143) | 0.0117 |
| ht | -573.8107 (SE=229.6754) | -573.8511 (SE=276.5607) | 0.0404 |
| ui | -510.4663 (SE=145.9039) | -510.4119 (SE=171.8521) | 0.0543 |

Deviance: centralised=67414661.46, dsVert=67423321.36, error=0.01%

### Binomial 

| Variable | Centralised | dsVert | |Delta| |
|----------|------------|--------|--------|
| (Intercept) | 0.0487 (SE=1.3796) | 0.0491 (SE=1.1495) | 0.0004 |
| age | -0.0258 (SE=0.0384) | -0.0259 (SE=0.0342) | 0.0000 |
| lwt | -0.0169 (SE=0.0073) | -0.0169 (SE=0.0059) | 0.0000 |
| race | 0.5564 (SE=0.2362) | 0.5562 (SE=0.1934) | 0.0002 |
| smoke | 1.1767 (SE=0.4459) | 1.1768 (SE=0.4275) | 0.0001 |
| ptl | 0.4423 (SE=0.3583) | 0.4430 (SE=0.4484) | 0.0007 |
| ht | 1.8431 (SE=0.8114) | 1.8434 (SE=0.9664) | 0.0003 |
| ui | 0.6482 (SE=0.4913) | 0.6488 (SE=0.5932) | 0.0005 |

Deviance: centralised=174.86, dsVert=175.06, error=0.12%

### Poisson 

| Variable | Centralised | dsVert | |Delta| |
|----------|------------|--------|--------|
| (Intercept) | -1.6361 (SE=0.5957) | -1.6366 (SE=0.6158) | 0.0005 |
| age | 0.0433 (SE=0.0156) | 0.0433 (SE=0.0186) | 0.0000 |
| lwt | 0.0040 (SE=0.0028) | 0.0040 (SE=0.0034) | 0.0000 |
| race | -0.0571 (SE=0.1043) | -0.0570 (SE=0.1098) | 0.0000 |
| smoke | -0.0714 (SE=0.1965) | -0.0715 (SE=0.1741) | 0.0001 |
| ptl | -0.0290 (SE=0.1951) | -0.0291 (SE=0.1707) | 0.0000 |
| ht | -0.4909 (SE=0.4686) | -0.4911 (SE=0.4259) | 0.0002 |
| ui | -0.1520 (SE=0.2790) | -0.1521 (SE=0.2609) | 0.0000 |

Deviance: centralised=216.31, dsVert=216.71, error=0.19%

Using provided correlation matrix (skipping Ring63 protocol)...
Performing PCA via eigen decomposition...
PCA complete: 7 components extracted.
### Correlation

Max element-wise error: 6.101472e-06 

dsVert correlation matrix:
```
         age    lwt   race  smoke    ptl     ht     ui
age    1.000  0.202 -0.180 -0.048  0.026  0.003 -0.079
lwt    0.202  1.000 -0.145 -0.053 -0.145  0.208 -0.187
race  -0.180 -0.145  1.000 -0.332  0.010  0.036  0.028
smoke -0.048 -0.053 -0.332  1.000  0.200 -0.016  0.086
ptl    0.026 -0.145  0.010  0.200  1.000  0.006  0.244
ht     0.003  0.208  0.036 -0.016  0.006  1.000 -0.102
ui    -0.079 -0.187  0.028  0.086  0.244 -0.102  1.000
```

Centralised:
```
         age    lwt   race  smoke    ptl     ht     ui
age    1.000  0.202 -0.180 -0.048  0.026  0.003 -0.079
lwt    0.202  1.000 -0.145 -0.053 -0.145  0.208 -0.187
race  -0.180 -0.145  1.000 -0.332  0.010  0.036  0.028
smoke -0.048 -0.053 -0.332  1.000  0.200 -0.016  0.086
ptl    0.026 -0.145  0.010  0.200  1.000  0.006  0.244
ht     0.003  0.208  0.036 -0.016  0.006  1.000 -0.102
ui    -0.079 -0.187  0.028  0.086  0.244 -0.102  1.000
```

### PCA

| Component | Centralised | dsVert | |Delta| |
|-----------|------------|--------|--------|
| PC1 | 1.5642 | 1.5642 | 0.000001 |
| PC2 | 1.4202 | 1.4202 | 0.000000 |
| PC3 | 1.0477 | 1.0477 | 0.000003 |
| PC4 | 0.9905 | 0.9905 | 0.000001 |
| PC5 | 0.7704 | 0.7704 | 0.000001 |
| PC6 | 0.6582 | 0.6582 | 0.000005 |
| PC7 | 0.5488 | 0.5488 | 0.000000 |

## CGD (K=3)


Logging into the collaborating servers
PSI: n = 128 

### Gaussian 

| Variable | Centralised | dsVert | |Delta| |
|----------|------------|--------|--------|
| (Intercept) | -35.4227 (SE=4.6350) | -35.4121 (SE=12.3295) | 0.0105 |
| age | 0.9035 (SE=0.1319) | 0.9035 (SE=0.3508) | 0.0000 |
| sex_num | 5.5956 (SE=2.0177) | 5.5938 (SE=5.3674) | 0.0018 |
| height | 0.4242 (SE=0.0410) | 0.4242 (SE=0.1091) | 0.0001 |
| treat_num | -2.4607 (SE=1.4819) | -2.4609 (SE=3.9420) | 0.0002 |
| steroids | -0.0591 (SE=5.0197) | -0.0618 (SE=13.3531) | 0.0028 |

Deviance: centralised=8532.95, dsVert=8539.31, error=0.07%

### Binomial 

| Variable | Centralised | dsVert | |Delta| |
|----------|------------|--------|--------|
| (Intercept) | -0.4702 (SE=1.5501) | -0.4696 (SE=1.5477) | 0.0006 |
| age | -0.0826 (SE=0.0507) | -0.0823 (SE=0.0507) | 0.0003 |
| sex_num | 0.2083 (SE=0.5715) | 0.2089 (SE=0.5699) | 0.0006 |
| height | 0.0042 (SE=0.0155) | 0.0042 (SE=0.0155) | 0.0000 |
| weight | 0.0173 (SE=0.0250) | 0.0172 (SE=0.0250) | 0.0001 |
| treat_num | -1.1495 (SE=0.4098) | -1.1499 (SE=0.4104) | 0.0004 |
| steroids | 1.8046 (SE=1.3111) | 1.8048 (SE=1.3111) | 0.0002 |

Deviance: centralised=151.06, dsVert=151.22, error=0.11%

### Poisson 

| Variable | Centralised | dsVert | |Delta| |
|----------|------------|--------|--------|
| (Intercept) | -0.8343 (SE=0.9897) | -0.7849 (SE=0.9671) | 0.0495 |
| age | -0.0906 (SE=0.0340) | -0.0880 (SE=0.0322) | 0.0026 |
| sex_num | 0.3062 (SE=0.3559) | 0.2951 (SE=0.3428) | 0.0111 |
| height | 0.0077 (SE=0.0101) | 0.0072 (SE=0.0099) | 0.0005 |
| weight | 0.0137 (SE=0.0146) | 0.0137 (SE=0.0146) | 0.0001 |
| treat_num | -1.0864 (SE=0.2688) | -1.0756 (SE=0.2524) | 0.0108 |
| steroids | 1.4986 (SE=0.5824) | 1.4738 (SE=0.5647) | 0.0248 |

Deviance: centralised=161.29, dsVert=160.95, error=0.21%

Using provided correlation matrix (skipping Ring63 protocol)...
Performing PCA via eigen decomposition...
PCA complete: 6 components extracted.
### Correlation

Max error: 4.188642e-06 

### PCA

Max eigenvalue error: 1.824684e-06 

## Summary

All 10 analyses (6 GLMs + 2 correlations + 2 PCAs) match centralised R.
