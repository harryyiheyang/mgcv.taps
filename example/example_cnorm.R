library(mgcv)
library(ggplot2)
library(devtools)
document()

n <- 1000
n_sim <- 300
p_values <- matrix(NA, nrow = n_sim, ncol = 3)
colnames(p_values) <- c("right", "left", "interval")

for (i in 1:n_sim) {
  X <- MASS::mvrnorm(n, rep(0, 4), matrix(0.25, 4, 4) + 0.75 * diag(4))
  x1 <- qbeta(pnorm(X[, 1]), 1.5, 1.5)
  x2 <- qbeta(pnorm(X[, 2]), 1.5, 1.5)
  x3 <- qbeta(pnorm(X[, 3]), 1.5, 1.5)
  x4 <- qbeta(pnorm(X[, 4]), 1.5, 1.5)
  t2 <- 2 * pi * x2
  t3 <- 2 * (x3 - 0.5)
  eta <- 2 * x1 +
    0.4 * sin(t2) + 0.8 * cos(t2) + 1.2 * sin(t2)^2 + 1.6 * cos(t2)^3 + 2 * sin(t2)^3 +
    3 * sin(3 * t3) + 6 * exp(-36 * t3^2)
  y_true <- eta + rnorm(n, sd = 1.5)
  mu <- mean(y_true)

  # 右删失数据集
  ct_r <- rnorm(n, mean = mu + 1, sd = 2)
  dat_r <- data.frame(
    x1 = x1, x2 = x2, x3 = x3, x4 = x4,
    y1 = pmin(y_true, ct_r),
    y2 = ifelse(y_true > ct_r, Inf, pmin(y_true, ct_r))
  )

  # 左删失数据集
  ct_l <- rnorm(n, mean = mu - 1, sd = 2)
  dat_l <- data.frame(
    x1 = x1, x2 = x2, x3 = x3, x4 = x4,
    y1 = pmax(y_true, ct_l),
    y2 = ifelse(y_true < ct_l, -Inf, pmax(y_true, ct_l))
  )

  # 区间删失数据集
  width <- runif(n, 0.5, 1.5)
  dat_i <- data.frame(
    x1 = x1, x2 = x2, x3 = x3, x4 = x4,
    y1 = y_true - width / 2,
    y2 = y_true + width / 2
  )

  run_test <- function(dat) {
    tryCatch({
      fit <- gam(
        cbind(y1, y2) ~ s(x1, bs = "AMatern", k = 10) +
          s(x2, bs = "cr", k = 10) + s(x3, bs = "cr", k = 15) + s(x4, bs = "cr", k = 10),
        data = dat, method = "REML", family = cnorm()
      )
      taps_score_test(fit, test.component = 1)$smooth.pvalue
    }, error = function(e) { message("Error sim ", i, ": ", e$message); NA })
  }

  p_values[i, "right"]    <- run_test(dat_r)
  p_values[i, "left"]     <- run_test(dat_l)
  p_values[i, "interval"] <- run_test(dat_i)

  if (i %% 20 == 0) cat("Completed", i, "simulations\n")
}

# 画图
make_qq_df <- function(pv) {
  pv <- na.omit(pv)
  n  <- length(pv)
  data.frame(
    log_obs      = -log10(sort(pv)),
    log_exp      = -log10((1:n) / (n + 1)),
    log_ci_lower = -log10(qbeta(0.975, 1:n, n - 1:n + 1)),
    log_ci_upper = -log10(qbeta(0.025, 1:n, n - 1:n + 1))
  )
}

df_all <- do.call(rbind, lapply(c("right", "left", "interval"), function(nm) {
  cbind(make_qq_df(p_values[, nm]), censoring = nm)
}))

ggplot(df_all, aes(x = log_exp, y = log_obs)) +
  geom_ribbon(aes(ymin = log_ci_lower, ymax = log_ci_upper), fill = "grey70", alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  geom_point(shape = 21, fill = "black", color = "white", size = 2, alpha = 0.8) +
  facet_wrap(~ censoring, nrow = 1) +
  labs(x = expression(Expected~~-log[10](italic(p))),
       y = expression(Observed~~-log[10](italic(p))),
       title = "Q-Q Plot for cnorm (null: x1 linear)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))
