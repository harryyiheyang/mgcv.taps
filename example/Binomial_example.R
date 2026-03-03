require(mgcv)
require(MASS)
library(devtools)
devtools::document()
gen_data_null <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X  <- mvrnorm(n, rep(0, 4), matrix(0.25, 4, 4) + 0.75 * diag(4))
  x1 <- qbeta(pnorm(X[, 1]), 1.5, 1.5)
  x2 <- qbeta(pnorm(X[, 2]), 1.5, 1.5)
  x3 <- qbeta(pnorm(X[, 3]), 1.5, 1.5)
  x4 <- qbeta(pnorm(X[, 4]), 1.5, 1.5)
  t2  <- 2 * pi * x2
  f2  <- 0.4*sin(t2) + 0.8*cos(t2) + 1.2*sin(t2)^2 + 1.6*cos(t2)^3 + 2*sin(t2)^3
  t3  <- 2 * (x3 - 0.5)
  f3  <- 3 * sin(3 * t3) + 6 * exp(-36 * t3^2)
  # x1 (s(x0)) is the NULL component being tested — no f1 term
  eta <- f2 + f3
  list(x1=x1, x2=x2, x3=x3, x4=x4, eta=eta)
}

n    <- 1000
nsim <- 300

# Two columns: [no PG, PG]
PV_binom <- matrix(NA, nsim, 4)

for (i in 1:nsim) {

  # ---- Binomial ----
  d       <- gen_data_null(n)
  n_trial <- 1
  prob    <- binomial()$linkinv(d$eta - 1)
  y       <- rbinom(n, size = n_trial, prob = prob)
  dat     <- data.frame(x0=d$x1, x1=d$x2, x2=d$x3, x3=d$x4,
                        y=y, n_trial=n_trial)
  b <- gam(cbind(y,n_trial-y) ~ s(x0, bs="AMatern") + s(x1, bs="cr") +
             s(x2, bs="cr") + s(x3, bs="cr"),
           data=dat, family=binomial(), method="REML")

  PV_binom[i, 1] <- taps_score_test(b, test.component=1,
                                    method="davies")$smooth.pvalue

  # ---- Binomial ----
  d       <- gen_data_null(n)
  n_trial <- 2
  prob    <- binomial()$linkinv(d$eta - 1)
  y       <- rbinom(n, size = n_trial, prob = prob)
  dat     <- data.frame(x0=d$x1, x1=d$x2, x2=d$x3, x3=d$x4,
                        y=y, n_trial=n_trial)
  b <- gam(cbind(y,n_trial-y) ~ s(x0, bs="AMatern") + s(x1, bs="cr") +
             s(x2, bs="cr") + s(x3, bs="cr"),
           data=dat, family=binomial(), method="REML")

  PV_binom[i, 2] <- taps_score_test(b, test.component=1,
                                    method="davies")$smooth.pvalue

  # ---- Binomial ----
  d       <- gen_data_null(n)
  n_trial <- 5
  prob    <- binomial()$linkinv(d$eta - 1)
  y       <- rbinom(n, size = n_trial, prob = prob)
  dat     <- data.frame(x0=d$x1, x1=d$x2, x2=d$x3, x3=d$x4,
                        y=y, n_trial=n_trial)
  b <- gam(cbind(y,n_trial-y) ~ s(x0, bs="AMatern") + s(x1, bs="cr") +
             s(x2, bs="cr") + s(x3, bs="cr"),
           data=dat, family=binomial(), method="REML")

  PV_binom[i, 3] <- taps_score_test(b, test.component=1,
                                    method="davies")$smooth.pvalue

  # ---- Binomial ----
  d       <- gen_data_null(n)
  n_trial <- 10
  prob    <- binomial()$linkinv(d$eta - 1)
  y       <- rbinom(n, size = n_trial, prob = prob)
  dat     <- data.frame(x0=d$x1, x1=d$x2, x2=d$x3, x3=d$x4,
                        y=y, n_trial=n_trial)
  b <- gam(cbind(y,n_trial-y) ~ s(x0, bs="AMatern") + s(x1, bs="cr") +
             s(x2, bs="cr") + s(x3, bs="cr"),
           data=dat, family=binomial(), method="REML")

  PV_binom[i, 4] <- taps_score_test(b, test.component=1,
                                    method="davies")$smooth.pvalue

  if (i %% 50 == 0) cat("iteration", i, "\n")
}

# ---- Plot ----
library(ggplot2)

make_qq_df <- function(pv, model_name) {
  pv_clean <- na.omit(pv)
  N  <- length(pv_clean)
  df <- data.frame(
    observed    = sort(pv_clean),
    expected    = (1:N) / (N + 1),
    model       = model_name
  )
  df$ci_lower     <- qbeta(0.025, 1:N, N - 1:N + 1)
  df$ci_upper     <- qbeta(0.975, 1:N, N - 1:N + 1)
  df$log_obs      <- -log10(df$observed)
  df$log_exp      <- -log10(df$expected)
  df$log_ci_lower <- -log10(df$ci_upper)
  df$log_ci_upper <- -log10(df$ci_lower)
  df
}

df_all <- rbind(
  make_qq_df(PV_binom[, 1], "Binomial | 1"),
  make_qq_df(PV_binom[, 2], "Binomial | 2"),
  make_qq_df(PV_binom[, 3], "Binomial | 5"),
  make_qq_df(PV_binom[, 4], "Binomial | 8")
)

ggplot(df_all, aes(x = log_exp, y = log_obs)) +
  geom_ribbon(aes(ymin = log_ci_lower, ymax = log_ci_upper),
              fill = "grey70", alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1,
              color = "red", linetype = "dashed", linewidth = 1) +
  geom_point(shape = 21, fill = "black", color = "white",
             size = 2, alpha = 0.8) +
  facet_wrap(~ model, scales = "free", nrow = 2) +
  labs(
    title = "Null P-value calibration: PG vs No PG",
    x = expression(Expected~~-log[10](italic(p))),
    y = expression(Observed~~-log[10](italic(p)))
  ) +
  theme_bw() +
  theme(
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "grey95"),
    strip.text        = element_text(face = "bold")
  )
