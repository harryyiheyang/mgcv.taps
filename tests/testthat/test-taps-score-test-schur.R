getA_linear_test <- function(x, para) cbind(1, x)

simulate_schur_test_data <- function(seed = 1L, n_id = 12L, n_time = 12L) {
  set.seed(seed)
  id <- factor(rep(seq_len(n_id), each = n_time))
  x <- rep(seq(0, 1, length.out = n_time), n_id)
  b0 <- rnorm(n_id, sd = 0.5)
  b1 <- rnorm(n_id, sd = 0.35)
  y <- 1 + 0.7 * x + b0[id] + b1[id] * x + rnorm(length(x), sd = 0.5)
  data.frame(y = y, x = x, id = id)
}

test_that("empty Schur projection reproduces the ordinary Davies test", {
  dat <- simulate_schur_test_data()
  expect_warning(fit <- mgcv::gam(
    y ~ s(
      x, bs = "AMatern", k = 7,
      xt = list(getA = getA_linear_test, para = NULL)
    ) + s(x, id, bs = "fs", k = 5, xt = "cr"),
    data = dat,
    method = "REML"
  ), "repeated 1-d smooths")

  ordinary <- taps_score_test(fit, test.component = 1, method = "davies")
  schur0 <- taps_score_test_schur(
    fit, test.component = 1, schur.component = integer(0),
    method = "davies"
  )

  expect_equal(schur0$smooth.pvalue, ordinary$smooth.pvalue, tolerance = 1e-7)
  expect_equal(schur0$schur.fraction, 0)
  expect_equal(schur0$information, schur0$efficient.information)
  expect_equal(schur0$schur.coefficients, "")
  expect_equal(schur0$schur.couplings, "")
  expect_equal(schur0$schur.contributions, "")
})

test_that("Schur projection removes information along a nuisance smooth", {
  dat <- simulate_schur_test_data(seed = 2L)
  expect_warning(fit <- mgcv::gam(
    y ~ s(
      x, bs = "AMatern", k = 7,
      xt = list(getA = getA_linear_test, para = NULL)
    ) + s(x, id, bs = "fs", k = 5, xt = "cr"),
    data = dat,
    method = "REML"
  ), "repeated 1-d smooths")

  ans <- taps_score_test_schur(
    fit, test.component = 1, schur.component = 2,
    method = "davies"
  )

  expect_true(is.finite(ans$smooth.pvalue))
  expect_gte(ans$smooth.pvalue, 0)
  expect_lte(ans$smooth.pvalue, 1)
  expect_gt(ans$schur.fraction, 0)
  expect_lt(ans$efficient.information, ans$information)
  expect_equal(ans$schur.terms, "s(x,id)")
})

test_that("component indexes are validated", {
  dat <- simulate_schur_test_data(seed = 3L, n_id = 8L, n_time = 10L)
  expect_warning(fit <- mgcv::gam(
    y ~ s(
      x, bs = "AMatern", k = 6,
      xt = list(getA = getA_linear_test, para = NULL)
    ) + s(x, id, bs = "fs", k = 4, xt = "cr"),
    data = dat,
    method = "REML"
  ), "repeated 1-d smooths")

  expect_error(
    taps_score_test_schur(fit, 1, c(2, 2)),
    "must not contain duplicated"
  )
  expect_error(
    taps_score_test_schur(fit, 1, 1),
    "cannot also be"
  )
  expect_error(
    taps_score_test_schur(fit, 1, 3),
    "valid smooth-term indexes"
  )
})

test_that("P-spline factor smooths can be Schur components", {
  dat <- simulate_schur_test_data(seed = 5L, n_id = 10L, n_time = 12L)
  expect_warning(fit <- mgcv::gam(
    y ~ s(
      x, bs = "AMatern", k = 7,
      xt = list(getA = getA_linear_test, para = NULL)
    ) + s(x, id, bs = "fs", k = 5, xt = "ps"),
    data = dat,
    method = "REML"
  ), "repeated 1-d smooths")

  ans <- taps_score_test_schur(
    fit, test.component = 1, schur.component = 2,
    method = "davies"
  )

  expect_true(is.finite(ans$smooth.pvalue))
  expect_gt(ans$schur.fraction, 0)
  expect_lt(ans$efficient.information, ans$information)
})

test_that("a vector of nuisance smooths is projected jointly", {
  dat <- simulate_schur_test_data(seed = 6L, n_id = 12L, n_time = 12L)
  set.seed(6)
  dat$z <- runif(nrow(dat))
  expect_warning(fit <- mgcv::gam(
    y ~ s(
      x, bs = "AMatern", k = 7,
      xt = list(getA = getA_linear_test, para = NULL)
    ) + s(x, id, bs = "fs", k = 5, xt = "ps") +
      s(z, bs = "cr", k = 5),
    data = dat,
    method = "REML"
  ), "repeated 1-d smooths")

  ans <- taps_score_test_schur(
    fit, test.component = 1, schur.component = c(2, 3),
    method = "davies"
  )

  expect_equal(ans$schur.rank, 2)
  expect_match(ans$schur.terms, "s\\(x,id\\)")
  expect_match(ans$schur.terms, "s\\(z\\)")
  expect_match(ans$schur.coefficients, "s\\(x,id\\)=")
  expect_match(ans$schur.coefficients, "s\\(z\\)=")
  expect_match(ans$schur.couplings, "s\\(x,id\\)=")
  expect_match(ans$schur.couplings, "s\\(z\\)=")
  expect_match(ans$schur.contributions, "s\\(x,id\\)=")
  expect_match(ans$schur.contributions, "s\\(z\\)=")
  expect_length(strsplit(ans$schur.couplings, ", ", fixed = TRUE)[[1L]], 2L)
  contribution_text <- strsplit(
    ans$schur.contributions, ", ", fixed = TRUE
  )[[1L]]
  contribution_value <- as.numeric(sub("^.*=", "", contribution_text))
  expect_equal(
    sum(contribution_value), ans$schur.fraction,
    tolerance = 1e-7
  )
  expect_gt(ans$schur.fraction, 0)
  expect_lt(ans$efficient.information, ans$information)
})

test_that("two random-effect couplings are labeled in one result cell", {
  dat <- simulate_schur_test_data(seed = 8L, n_id = 12L, n_time = 12L)
  dat$site <- factor(rep(rep(1:3, each = 4), each = 12))
  fit <- mgcv::gam(
    y ~ s(
      x, bs = "AMatern", k = 7,
      xt = list(getA = getA_linear_test, para = NULL)
    ) + s(id, bs = "re") + s(site, bs = "re"),
    data = dat,
    method = "REML"
  )

  ans <- taps_score_test_schur(
    fit, test.component = 1, schur.component = c(2, 3),
    method = "davies"
  )

  expect_equal(ans$schur.terms, "s(id), s(site)")
  expect_match(ans$schur.couplings, "^s\\(id\\)=")
  expect_match(ans$schur.couplings, ", s\\(site\\)=")
  expect_length(
    strsplit(ans$schur.couplings, ", ", fixed = TRUE)[[1L]], 2L
  )
})

test_that("gamm4 Schur branch handles finite random effects", {
  skip_if_not_installed("gamm4", minimum_version = "0.3.0")
  skip_if_not_installed("lme4", minimum_version = "2.0.6")

  dat <- simulate_schur_test_data(seed = 4L, n_id = 10L, n_time = 10L)
  expect_warning(fit <- gamm4::gamm4(
    y ~ s(
      x, bs = "AMatern", k = 6,
      xt = list(getA = getA_linear_test, para = NULL)
    ) + s(x, id, bs = "fs", k = 4, xt = "cr"),
    random = ~(1 | id),
    data = dat,
    REML = TRUE
  ), "repeated 1-d smooths")

  expect_warning(
    ordinary <- taps_score_test(fit, test.component = 1, method = "davies"),
    "Consider playing with"
  )
  expect_warning(schur0 <- taps_score_test_schur(
    fit, test.component = 1, schur.component = integer(0),
    method = "davies"
  ), "Consider playing with")
  ans <- taps_score_test_schur(
    fit, test.component = 1, schur.component = 2,
    method = "davies"
  )

  expect_equal(schur0$smooth.pvalue, ordinary$smooth.pvalue, tolerance = 1e-7)
  expect_gt(ans$schur.fraction, 0)
  expect_lt(ans$efficient.information, ans$information)
})

test_that("only Davies calibration is accepted", {
  dat <- simulate_schur_test_data(seed = 7L, n_id = 8L, n_time = 10L)
  expect_warning(fit <- mgcv::gam(
    y ~ s(
      x, bs = "AMatern", k = 6,
      xt = list(getA = getA_linear_test, para = NULL)
    ) + s(x, id, bs = "fs", k = 4, xt = "ps"),
    data = dat,
    method = "REML"
  ), "repeated 1-d smooths")

  expect_error(
    taps_score_test_schur(fit, 1, 2, method = "normal"),
    "method must be 'davies'"
  )
})
