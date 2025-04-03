mgcv.taps: Implementation of TAPS (Test for Arbitrary Parametric
Structure) by using R package mgcv.
================

# mgcv.taps

`mgcv.taps` is an extension of the `mgcv` package that implements the
**Test for Arbitrary Parametric Structure (TAPS)**. When considering
nonparametric approach like a generalized additive model (GAM), a common
starting point is to evaluate whether a simpler parametric model, such
as a generalized linear model (GLM), can adequately capture the data’s
variation. TAPS directly addresses this question by providing inference
tools to assess whether a parametric structure is sufficient to describe
the target function within a GAM framework.

![](README_files/figure-gfm/TAPS_Introduction.png)

## Installation

You can install the development version of `mgcv.taps` from GitHub:

``` r
devtools::install_github("harryyiheyang/mgcv.taps")
```

## Overview

This package supports smooth function construction, estimation, and
hypothesis testing under arbitrary parametric constraints, by extending
the `mgcv` smoothing interface.

Two types of smoothers are currently implemented:

### 1. Univariate Structured Smooth (`bs="AMatern"`)

This smoother models a function as:

    f(x) = A(x)%*%alpha + B(x)%*%beta

- `A(x)` is a user-defined basis matrix representing a structured
  parametric component (e.g., `cbind(1,x,x^2)` for quadratic functions,
  or `cbind(1,x,pmax(0,x-nu))` for piecewise structures).
- `B(x)` is a set of smooth basis functions adaptively constructed using
  a Matern kernel and constrained to be orthogonal to `A(x)`:

<!-- -->

    t(A)%*%B=0

The statistical method used to enforce `t(A) %*% B = 0` share insights
with the boundary condition applied in thin plate splines ([Wood,
2003](https://academic.oup.com/jrsssb/article/65/1/95/7110632)).

### 2. Bivariate Structured Smooth (`bs="A2Matern"`)

This extension handles bivariate smooths with similar structure:

    f(x1,x2) = A(x1,x2)%*%alpha + B(x1,x2)%*%beta

- `A(x1,x2)` can represent structured interaction effects (e.g.,
  `cbind(1,x1,x2,x1*x2)`).
- `B(x1,x2)` is constructed using Matern kernels and is orthogonal to
  `A(x1,x2)`.

Additional smoother types may be supported in future versions.

## Estimation

Estimation is performed using `mgcv::gam()`. The `xt` argument should be
a list:

- `getA`: a user-specified function returning the parametric basis
  matrix `A`.
- `para`: optional parameters used in constructing `A`.

### Univariate Example

    gam(y~s(x,bs="AMatern",k=10,xt=list(getA=function(x,para) cbind(1,x))))

### Bivariate Example

    gam(y~s(x1,x2,bs="A2Matern",k=10,xt=list(getA=function(x1,x2,para) cbind(1,x1,x2,x1*x2))))

## Hypothesis Testing

The package provides a wrapper for Wald tests on specific components of
the model using `taps_wald_test`, which actually wraps an unobserved
function `mgcv::testStat` ([Wood,
2013](https://academic.oup.com/biomet/article/100/1/221/192816)):

    taps_wald_test(fit,test.component=1)

Alternatively, one can also perform a score test variance component
using `taps_score_test` ([Zhang and Lin,
2003](https://doi.org/10.1093/biostatistics/4.1.57)).

    taps_score_test(fit,test.component=1)

Both can be used to formally test the parametric structure encoded in
`A`. In these two test functions, `test.component` refers to the index
of smooth term to be tested.

## Family of Outcome

Currently, one can use the two novel smoothers, `AMatern` and
`A2Matern`, to estimate models under any family of outcome supported by
`mgcv`. For hypothesis testing, `taps_wald_test` also supports all
outcome families. However, `taps_score_test` is currently limited to
exponential family distributions.

## Examples

The example utilizes the `gov_transfers` dataset from the R package
`causaldata`, which includes data from [Manacorda et
al. (2011)](https://www.aeaweb.org/articles?id=10.1257/app.3.3.1),
regarding a government transfer program allocated based on an income
threshold. The dataset is pre-filtered to include only households near
this income cutoff. For additional details, refer to `?gov_transfer`.

The codes below fit the effect of income on support of the government
yielded by the standard GAM:

``` r
library(mgcv)
library(mgcv.taps)
library(causaldata)
options(bitmapType="cairo")
data("gov_transfers")
gov_transfers$support_success=gov_transfers$Support*2
gov_transfers$support_failure=2-gov_transfers$support_success
fit0=gam(cbind(support_success,support_failure)~s(Income_Centered,bs="gp")+Education+Age,data=gov_transfers,family=binomial(link="probit"))
plot(fit0,main="Standard GAM Fit of Income Effect")
```

![](README_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

We then investigate whether there is a discontinuity at `income=0`,
using the model `A=cbind(1, income, pmax(0, income))`:

``` r
fit1=gam(cbind(support_success,support_failure)~s(Income_Centered,bs="AMatern",xt=list(getA=linearity_discontinuity,para=0))+Education+Age,data=gov_transfers,method="REML",family=binomial(link="probit"))
plot(fit1,main="GAM Fit of Income Effect with one Breakpoint")
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
taps_score_test(fit1)
```

    ##           smooth.term smooth.df smooth.stat smooth.pvalue
    ##                <char>     <num>       <num>         <num>
    ## 1: s(Income_Centered)  2.523577    2.931358     0.3194111

``` r
taps_wald_test(fit1)
```

    ##            mixed.term fix.df fix.chisq   fix.pvalue fix.indices smooth.df
    ##                <char>  <num>     <num>        <num>      <char>     <num>
    ## 1: s(Income_Centered)      3  72.20276 1.440396e-15    reported  4.115946
    ##    smooth.chisq smooth.pvalue
    ##           <num>         <num>
    ## 1:     12.05502    0.01814953

Based on the score test, there is no evidence to suggest that a more
complex model beyond a linearity discontinuity with a breakpoint at
`income = 0` is necessary. However, the Wald test indicates a slight
non-linear effect in addition to the linearity discontinuity.

We conclude by examining the effect pattern of income on support using a
parametric GLM. For simplicity and clarity in presentation, we developed
a fixed-effect smooth term, `s(Income_Centered, fx=TRUE, bs="Fixed")`,
which exactly models the effect of income as a parametric function
defined by linearity_discontinuity. Refer to `?gam`for more details on
fixed-effect smooth terms (`fx=TRUE`). Below are our results:

``` r
fit2 = gam(cbind(support_success, support_failure) ~ s(Income_Centered, fx=TRUE, bs="Fixed", xt=list(getA=linearity_discontinuity, para=0)) + Education + Age, data = gov_transfers, family = binomial(link="probit"))
plot(fit2, main = "Exact GLM Fit of Income Effect with One Breakpoint")
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## License

MIT

## Maintainer

Yihe Yang Email: <yxy1234@case.edu>
