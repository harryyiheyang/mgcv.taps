require(mgcv)
require(MASS)
require(tweedie)

gen_data <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X  <- mvrnorm(n, rep(0, 4), matrix(0.25, 4, 4) + 0.75 * diag(4))
  x1 <- qbeta(pnorm(X[, 1]), 1.5, 1.5)
  x2 <- qbeta(pnorm(X[, 2]), 1.5, 1.5)
  x3 <- qbeta(pnorm(X[, 3]), 1.5, 1.5)
  x4 <- qbeta(pnorm(X[, 4]), 1.5, 1.5)
  t2  <- 2 * pi * x2
  f1  <- x1 * 2
  f2  <- 0.4*sin(t2)+0.8*cos(t2)+1.2*sin(t2)^2+1.6*cos(t2)^3+2*sin(t2)^3
  t3  <- 2 * (x3 - 0.5)
  f3  <- 3 * sin(3 * t3) + 6 * exp(-36 * t3^2)
  eta <- f1 + f2 + f3
  list(x1=x1, x2=x2, x3=x3, x4=x4, eta=eta)
}

n    <- 1000
nsim <- 300

# storage: rows = simulations, cols = [old, new]
PV_nb   <- matrix(NA, nsim, 1)
PV_tw   <- matrix(NA, nsim, 1)
PV_sc   <- matrix(NA, nsim, 1)
PV_br   <- matrix(NA, nsim, 1)

for (i in 1:nsim) {

  # ---- nb ----
  d   <- gen_data(n)
  y   <- rnbinom(n, mu = exp(d$eta/2-1), size = 1.5)
  dat <- data.frame(x0=d$x1, x1=d$x2, x2=d$x3, x3=d$x4, y=y)
  b   <- gam(y ~ s(x0,bs="AMatern") + s(x1,bs="cr") +
     s(x2,bs="cr")+ s(x3,bs="cr"),
   data=dat, family=nb(), method="REML")
  PV_nb[i, 1] <- taps_score_test(b,test.component=1, method="davies")$smooth.pvalue

  # ---- tw ----
  d   <- gen_data(n)
  y   <- rtweedie(n, mu = exp(d$eta/2-1), phi = 3, power = 1.5)
  dat <- data.frame(x0=d$x1, x1=d$x2, x2=d$x3, x3=d$x4, y=y)
  b   <- gam(y ~ s(x0,bs="AMatern") + s(x1,bs="cr") +
     s(x2,bs="cr")      + s(x3,bs="cr"),
   data=dat, family=tw(), method="REML")
  PV_tw[i, 1] <- taps_score_test(b,test.component=1, method="davies")$smooth.pvalue

  # ---- scat ----
  d   <- gen_data(n)
  y   <- d$eta + rt(n, df = 5)*2.5
  dat <- data.frame(x0=d$x1, x1=d$x2, x2=d$x3, x3=d$x4, y=y)
  b   <- gam(y ~ s(x0,bs="AMatern") + s(x1,bs="cr") +
     s(x2,bs="cr")      + s(x3,bs="cr"),
   data=dat, family=scat(), method="REML")
  PV_sc[i, 1] <- taps_score_test(b, test.component=1, method="davies")$smooth.pvalue

  # ---- betar ----
  d   <- gen_data(n)
  mu  <- binomial()$linkinv(d$eta / 2 - 2)
  phi <- 0.5
  y   <- rbeta(n, mu * phi, phi - mu * phi)
  dat <- data.frame(x0=d$x1, x1=d$x2, x2=d$x3, x3=d$x4, y=y)
  b   <- gam(y ~ s(x0,bs="AMatern") + s(x1,bs="cr") +
     s(x2,bs="cr")      + s(x3,bs="cr"),
   data=dat, family=betar(link="logit"), method="REML")
  PV_br[i, 1] <- taps_score_test(b,test.component=1, method="davies")$smooth.pvalue

  if (i %% 50 == 0) cat("iteration", i, "\n")
}

library(ggplot2)

make_qq_df <- function(pv, model_name) {
  pv_clean <- na.omit(pv)
  N <- length(pv_clean)
  df <- data.frame(
    observed = sort(pv_clean),
    expected = (1:N) / (N + 1),
    model = model_name
  )
  df$ci_lower <- qbeta(0.025, 1:N, N - 1:N + 1)
  df$ci_upper <- qbeta(0.975, 1:N, N - 1:N + 1)
  df$log_obs <- -log10(df$observed)
  df$log_exp <- -log10(df$expected)
  df$log_ci_lower <- -log10(df$ci_upper)
  df$log_ci_upper <- -log10(df$ci_lower)
  return(df)
}

df_all <- rbind(
  make_qq_df(PV_nb[, 1], "nb | taps_score_test"),
  make_qq_df(PV_tw[, 1], "tw | taps_score_test"),
  make_qq_df(PV_sc[, 1], "scat | taps_score_test"),
  make_qq_df(PV_br[, 1], "betar | taps_score_test")
)

ggplot(df_all, aes(x = log_exp, y = log_obs)) +
  geom_ribbon(aes(ymin = log_ci_lower, ymax = log_ci_upper), fill = "grey70", alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  geom_point(shape = 21, fill = "black", color = "white", size = 2, alpha = 0.8) +
  facet_wrap(~ model, scales = "free") +
  labs(
    x = expression(Expected~~-log[10](italic(p))),
    y = expression(Observed~~-log[10](italic(p)))
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold")
  )
