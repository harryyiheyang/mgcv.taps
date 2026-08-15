library(mgcv)
library(mgcv.taps)

set.seed(20260814)
n_id <- 14L
n_each <- 7L
id <- factor(rep(seq_len(n_id), each = n_each))
time <- rep(seq(0, 1, length.out = n_each), n_id)
x <- runif(length(id))
u0 <- rnorm(n_id, sd = 0.35)
u1 <- rnorm(n_id, sd = 0.2)
y <- 1 + 0.4 * x + sin(2 * pi * x) + u0[id] + u1[id] * time +
  rnorm(length(id), sd = 0.45)
dat <- data.frame(y = y, x = x, time = time, id = id)

fit_gam <- gam(
  y ~ s(x, k = 6L) + s(id, bs = "re"),
  data = dat, method = "REML"
)
gam_result <- taps_score_test_gamm(fit_gam, method = "liu")
gam_direct <- mgcv.taps:::taps_score_test_re(fit_gam, method = "liu")
if (attr(gam_result, "backend") != "mgcv" ||
    !isTRUE(all.equal(gam_result$smooth.pvalue,
                      gam_direct$smooth.pvalue, tolerance = 1e-12))) {
  stop("taps_score_test_gamm did not dispatch the native gam re model correctly.")
}
gam_default <- taps_score_test(fit_gam, method = "liu")
if (!is.finite(gam_default$smooth.pvalue)) {
  stop("taps_score_test did not retain ordinary mgcv re-smooth handling.")
}

re_component_error <- FALSE
tryCatch(
  taps_score_test_gamm(fit_gam, test.component = 2L, method = "liu"),
  error = function(e) {
    re_component_error <<- grepl("mean-structure smooth", conditionMessage(e),
                                 fixed = TRUE)
  }
)
fit_no_re <- gam(y ~ s(x, k = 6L), data = dat, method = "REML")
no_re_error <- FALSE
tryCatch(
  taps_score_test_gamm(fit_no_re, method = "liu"),
  error = function(e) {
    no_re_error <<- grepl("must contain an re or fs", conditionMessage(e),
                          fixed = TRUE)
  }
)
if (!re_component_error || !no_re_error) {
  stop("The unified GAMM interface did not enforce its mean/re structure contract.")
}

fit_bam <- bam(
  y ~ s(x, k = 6L) + s(time, id, bs = "fs", k = 4L, m = 1L),
  data = dat, method = "fREML", discrete = FALSE
)
bam_result <- taps_score_test_gamm(fit_bam, method = "liu")
if (attr(bam_result, "backend") != "mgcv" ||
    !is.finite(bam_result$smooth.pvalue)) {
  stop("taps_score_test_gamm did not dispatch the native bam fs model correctly.")
}

fit_bam_discrete <- bam(
  y ~ s(x, k = 6L) + s(id, bs = "re"),
  data = dat, method = "fREML", discrete = TRUE
)
bam_discrete_result <- taps_score_test_gamm(
  fit_bam_discrete, method = "liu"
)
if (attr(bam_discrete_result, "backend") != "mgcv" ||
    !is.finite(bam_discrete_result$smooth.pvalue)) {
  stop("taps_score_test_gamm did not support native discrete bam re dispatch.")
}

if (requireNamespace("gamm4", quietly = TRUE) &&
    requireNamespace("lme4", quietly = TRUE)) {
  fit_gamm4 <- gamm4::gamm4(
    y ~ s(x, k = 6L), random = ~(1 | id), data = dat
  )
  gamm4_result <- taps_score_test_gamm(fit_gamm4, method = "liu")
  gamm4_direct <- mgcv.taps:::taps_score_test_gamm4_re(
    fit_gamm4, method = "liu"
  )
  if (attr(gamm4_result, "backend") != "gamm4" ||
      !isTRUE(all.equal(gamm4_result$smooth.pvalue,
                        gamm4_direct$smooth.pvalue,
                        tolerance = 1e-12))) {
    stop("taps_score_test_gamm did not dispatch the gamm4 model correctly.")
  }
  generic_gamm4_error <- FALSE
  tryCatch(
    taps_score_test(fit_gamm4, method = "liu"),
    error = function(e) {
      generic_gamm4_error <<- grepl(
        "use taps_score_test_gamm", conditionMessage(e), fixed = TRUE
      )
    }
  )
  if (!generic_gamm4_error) {
    stop("taps_score_test did not reject a gamm4 fit.")
  }
}
