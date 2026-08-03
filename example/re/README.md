# RE (clustered) distribution examples

Same families, nuisance functions, and test settings as the scripts in
`example/`, with three design changes that mimic clustered data such as the
spatial RD applications (running variable constant within a census tract or
survey segment):

1. `n = 1000` observations come from `ncl = 50` clusters of `mcl = 20`.
2. The tested covariate (the `AMatern` term) is drawn at the cluster level,
   so it is constant within cluster. The three nuisance covariates remain
   individual-level (correlated MVN copula among themselves).
3. A cluster random intercept `u ~ N(0, sigma_u^2)` with `sigma_u = 0.5` is
   added on the link scale in both the NULL and alternative DGPs.

Each replicate fits two models and records both score-test p-values:

- `*_iid`: the original model formula (no cluster term), i.e., working
  independence. Under the NULL this is expected to over-reject because the
  cluster noise is confounded with smooth deviations in the tested term.
- `*_re`: the same formula plus `s(cl, bs = "re")`. This is the model-based
  correction; the NULL rejection rates should return to nominal.

100 replicates per scenario (proof-of-concept scale). Seeds are the original
family seed blocks prefixed with 5 (Cox 59, ziP 60). Run each
`*_re_example.R` from the repository root, then `summary_re.R` to produce
`re_summary.csv`.
