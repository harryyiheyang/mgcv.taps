# Distribution example review

All simulations were run with `n = 1000`, fixed replicate-specific seeds,
REML fitting, and the Davies form of `taps_score_test()`. The nine non-Cox
families use 100 NULL and 100 alternative replicates. Cox uses 500 replicates
per scenario after its data-generating process and censoring were reviewed.
No Wald test was run.

The covariate construction, nuisance functions, response-family parameters,
smooth bases, and fit settings follow the previous examples. Binary retains the
piecewise-linearity/discontinuity null and uses the original smoothing width
`0.05` under the alternative. Binomial retains the no-`x0`-effect NULL and uses
four trials per observation. `sbeta` retains the existing Beta-regression
semantics: `betar(link = "logit")` with generated precision `phi = 0.5`.

## Rejection rates and power

| Family | NULL 0.01 | NULL 0.05 | NULL 0.10 | Power 0.01 | Power 0.05 | Power 0.10 |
|---|---:|---:|---:|---:|---:|---:|
| Gaussian | 0.00 | 0.04 | 0.10 | 1.00 | 1.00 | 1.00 |
| Binary | 0.01 | 0.06 | 0.09 | 0.16 | 0.34 | 0.46 |
| Binomial | 0.02 | 0.09 | 0.11 | 0.72 | 0.84 | 0.88 |
| NB | 0.00 | 0.03 | 0.07 | 0.36 | 0.61 | 0.68 |
| scat | 0.00 | 0.02 | 0.11 | 0.11 | 0.31 | 0.42 |
| ocat | 0.01 | 0.05 | 0.07 | 0.30 | 0.56 | 0.74 |
| Cox | 0.006 | 0.050 | 0.116 | 0.362 | 0.584 | 0.706 |
| sbeta | 0.00 | 0.04 | 0.04 | 0.10 | 0.25 | 0.37 |
| Tweedie | 0.01 | 0.05 | 0.10 | 0.27 | 0.45 | 0.55 |
| ZIP | 0.02 | 0.07 | 0.13 | 0.18 | 0.32 | 0.45 |

## Median p-values and QQ deviations

`QQ max` is `max(abs(sort(p) - ppoints(length(p))))`; `QQ RMSE` is the
corresponding root mean squared deviation. The KS p-value is reported only for
the NULL sample and tests uniformity.

| Family | NULL median | Alt median | NULL QQ max | NULL QQ RMSE | Alt QQ max | Alt QQ RMSE | NULL KS p |
|---|---:|---:|---:|---:|---:|---:|---:|
| Gaussian | 0.5002 | 0.000002 | 0.0599 | 0.0280 | 0.9888 | 0.5769 | 0.7684 |
| Binary | 0.5592 | 0.1103 | 0.1010 | 0.0506 | 0.4095 | 0.2959 | 0.1964 |
| Binomial | 0.5107 | 0.0011 | 0.0689 | 0.0273 | 0.7945 | 0.5220 | 0.6198 |
| NB | 0.5503 | 0.0251 | 0.0810 | 0.0437 | 0.6072 | 0.4157 | 0.4265 |
| scat | 0.5272 | 0.1712 | 0.0880 | 0.0415 | 0.3868 | 0.2705 | 0.3316 |
| ocat | 0.5076 | 0.0307 | 0.0409 | 0.0189 | 0.6498 | 0.4405 | 0.9779 |
| Cox | 0.5245 | 0.0286 | 0.0296 | 0.0123 | 0.6218 | 0.4296 | 0.7247 |
| sbeta | 0.5650 | 0.2063 | 0.1573 | 0.0897 | 0.3470 | 0.2319 | 0.0091 |
| Tweedie | 0.4751 | 0.0696 | 0.0772 | 0.0247 | 0.4777 | 0.3353 | 0.4832 |
| ZIP | 0.4360 | 0.1203 | 0.0726 | 0.0334 | 0.4797 | 0.3159 | 0.5580 |

## Review findings

- The sbeta NULL p-values show the clearest non-uniformity. They are
  conservative in the lower tail: rejection is 0.04 at both 0.05 and 0.10,
  with KS p = 0.0091.
- The revised Cox generator is a correctly parameterized Weibull PH model,
  without observation-level frailty. Independent exponential censoring with
  rate `0.704` produced mean censoring fractions 0.511 under NULL and 0.487
  under the alternative. With 500 NULL replicates, rejection rates are 0.006,
  0.050, and 0.116 at thresholds 0.01, 0.05, and 0.10; KS p = 0.725 and NULL
  QQ max = 0.030. The earlier Cox calibration concern was primarily caused by
  the misspecified simulation setting.
- The 500 Cox NULL fits accumulated 19 fitting warnings, and the 500
  alternative fits accumulated 10. The alternative warnings shown in the log
  are `gam.fit5 step failed` messages with small maximum relative gradients
  (`2.0e-8` to `4.2e-7`). No Cox fit or score test stopped, and all 1,000
  p-values are finite.
- Binomial and ZIP have 0.05 rejection rates of 0.09 and 0.07 respectively;
  with only 100 replicates, these estimates have substantial Monte Carlo
  uncertainty and should be monitored in a larger calibration run.
- Every R invocation emitted locale startup warnings because the current
  Windows R installation could not set the requested `C.UTF-8` categories.
  Tweedie additionally reported that its installed package was built under R
  4.5.2 while the run used R 4.5.1.
- No family failed, and no replicate was skipped or replaced.

The machine-readable results are in `review_summary.csv`, and achieved Cox
censoring is in `cox_censoring_summary.csv`. All ten RDS files were validated
as two-element named lists containing only p-values: 500 NULL and 500
alternative values for Cox, and 100 per scenario for each other family.
