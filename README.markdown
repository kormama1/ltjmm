# `ltjmm`: Latent Time Joint Mixed Effects Models

(cloned from https://bitbucket.org/mdonohue/ltjmm/src/master/)

Diseases that progress over long periods of time are often studied by observing cohorts at different stages of disease for shorter periods of time. We apply MCMC sampling to estimate Latent Time Joint Mixed Effect Models (LTJMM) from short-term observations with unknown relative observation times. The LTJMM is described in [Li, et al. (2017)](https://doi.org/10.1177/0962280217737566). 

## To install, from an R prompt:

```r
install.packages("devtools")
devtools::install_bitbucket("mdonohue/ltjmm")
```

## Model details

The Stan code for LTJMMs is located in the `src/stan_files` subdirectory of this repository. All of the models implemented in this package are restricted to have a random intercept and slope for each outcome, the same fixed effect covariates for each outcome, and Gaussian residuals with identity link function.

The `ltjmm::ltjmm_stan` function can fit models with or without a latent time parameter (`lt=TRUE` or `lt=FALSE`). With `lt=FALSE` the model is a joint (or multivariate) mixed effect model. One can assume all of the random effects are from one multivariate Gaussian distribution (`random_effects='multivariate'`) or each random effect is from separate univariate Gaussian distributions (`random_effects='univariate'`).

Our fork of the [`rstanarm`](https://github.com/mcdonohue/rstanarm) repository includes a function `stan_ltjmm`, which allows different random effects, fixed effects, and exponential family links/distributions for up to 20 outcomes.

## Reference

* Li, D., Iddi, S., Thompson, W. K., Donohue, M. C., for ADNI. (2017). Bayesian latent time joint mixed effect models for multicohort longitudinal data. *Statistical methods in medical research*. [https://doi.org/10.1177/0962280217737566](https://doi.org/10.1177/0962280217737566)

