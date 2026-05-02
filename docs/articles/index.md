# Articles

### Federated estimators – one vignette per method

Each article exercises a single dsVert estimator end-to-end on a
canonical small public dataset, against the analyst’s go-to centralised
R reference ([`stats::glm`](https://rdrr.io/r/stats/glm.html),
[`survival::coxph`](https://rdrr.io/pkg/survival/man/coxph.html),
[`nlme::lme`](https://rdrr.io/pkg/nlme/man/lme.html),
[`MASS::glm.nb`](https://rdrr.io/pkg/MASS/man/glm.nb.html),
[`MASS::polr`](https://rdrr.io/pkg/MASS/man/polr.html),
[`nnet::multinom`](https://rdrr.io/pkg/nnet/man/multinom.html),
[`lme4::glmer`](https://rdrr.io/pkg/lme4/man/glmer.html),
[`geepack::geeglm`](https://rdrr.io/pkg/geepack/man/geeglm.html),
weighted glm, lm/glmnet), using the in-process DSLite multi-server
harness. The five-section template (context + math + non-disclosure +
partition table \| data preparation + DSLite stand-up + PSI alignment \|
federated fit \| centralised reference \| comparison + verdict) is
identical across methods so the reader can compare estimators by
skimming.

- [Generalised linear model — federated binomial logistic regression (K
  =
  3)](https://isglobal-brge.github.io/dsVertClient/articles/methods/glm.md):
