# **UPG**: Efficient Bayesian Regression Models for Binary and Categorical Outcomes

<!-- badges: start -->
[![CRAN](http://www.r-pkg.org/badges/version/UPG)](https://cran.r-project.org/package=UPG)
[![month](http://cranlogs.r-pkg.org/badges/UPG)](https://www.r-pkg.org/pkg/UPG)
[![total](http://cranlogs.r-pkg.org/badges/grand-total/UPG)](https://www.r-pkg.org/pkg/UPG)
<!-- badges: end -->

**UPG** offers efficient Bayesian implementations of regression models for binary and categorical data. The package can be used to estimate Bayesian versions of probit, logit, multinomial logit and binomial logit models. In this context, the Bayesian paradigm is especially useful for uncertainty quantification and solving issues related to rare events and (quasi-)perfect separation. In fact, **UPG** allows for highly efficient posterior sampling in cases with imbalanced data as the implemented algorithms are based on the marginal data augmentation schemes developed in Frühwirth-Schnatter, Zens, and Wagner (2020). Several functions are available for tabulating and visualizing results as well as for predictive exercises. 

## Installation

**UPG** is available on CRAN and can be installed as follows:

``` r
install.packages("UPG")
```

## Usage

The core function for estimating models is `UPG()`. Given a suitable outcome vector `y` and a suitable design matrix `X`, the four implemented models can be estimated using

- `UPG(y, X, model = "probit")` for probit models
- `UPG(y, X, model = "logit")` for binary logit models
- `UPG(y, X, model = "mnl")` for multinomial logit models
- `UPG(y, X, Ni, model = "binomial")` for binomial logit models

where binomial logit models require the number of trials `Ni` as additional input. 

The estimation output can be analyzed using a variety of tools implemented in **UPG**. To tabulate and visualize the results, `summary()` and `plot()` are available. Predictions can be obtained using `predict()`. Extracting coefficients can be done using `coef()` and `logLik()` returns the log-likelihood of the model. Finally, the user has access to a number of MCMC diagnostics via `UPG.Diag()`.

More details and applied examples may be found in the package vignette.

## References

Frühwirth-Schnatter, S., Zens, G., & Wagner, H. (2020). Ultimate Pólya Gamma Samplers - Efficient MCMC for possibly imbalanced binary and categorical data. arXiv preprint arXiv:2011.06898.
