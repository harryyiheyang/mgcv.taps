library(mgcv)
library(mgcv.taps)
devtools::document()

n <- 1000
n_sim <- 300
p_values <- numeric(n_sim)

for (i in 1:n_sim) {
  # 生成协变量
  X <- MASS::mvrnorm(n, rep(0, 4), matrix(0.25, 4, 4) + 0.75 * diag(4))
  x1 <- qbeta(pnorm(X[, 1]), 1.5, 1.5)
  x2 <- qbeta(pnorm(X[, 2]), 1.5, 1.5)
  x3 <- qbeta(pnorm(X[, 3]), 1.5, 1.5)
  x4 <- qbeta(pnorm(X[, 4]), 1.5, 1.5)

  # 零假设：x1是线性效应（没有非线性效应）
  f1 <- 2 * x1  # 纯线性

  # 其他非线性效应
  t2 <- 2 * pi * x2
  f2 <- 0.4 * sin(t2) + 0.8 * cos(t2) + 1.2 * sin(t2)^2

  t3 <- 2 * (x3 - 0.5)
  f3 <- 3 * sin(3 * t3) + 6 * exp(-36 * t3^2)

  # Linear predictor
  eta <- f1 + f2 + f3

  # 生成有序分类响应（J=5个类别）
  J <- 5
  cutpoints <- c(-2, -0.5, 0.5, 2)  # K=4个cutpoints

  probs <- matrix(0, n, J)
  for (j in 1:J) {
    if (j == 1) {
      probs[, j] <- plogis(cutpoints[1] - eta)
    } else if (j == J) {
      probs[, j] <- 1 - plogis(cutpoints[J-1] - eta)
    } else {
      probs[, j] <- plogis(cutpoints[j] - eta) - plogis(cutpoints[j-1] - eta)
    }
  }

  y <- apply(probs, 1, function(p) sample(1:J, 1, prob = p))

  data <- data.frame(
    y = as.integer(y),
    x1 = x1,
    x2 = x2,
    x3 = x3,
    x4 = x4
  )

  # 拟合模型
  fit <- mgcv::gam(
    y ~ s(x1, bs = "AMatern", k = 10) + s(x2, bs = "cr", k = 10) +
      s(x3, bs = "cr", k = 15) + s(x4, bs = "cr", k = 10),
    data = data,
    family = ocat(R = J),
    method = "REML"
  )

  # 测试x1的smooth term（在零假设下应该不显著）
  test <- taps_score_test(fit, test.component = 1)
  p_values[i] <- test$smooth.pvalue

  if (i %% 20 == 0) cat("Completed", i, "simulations\n")
}


library(ggplot2)
PV_clean <- na.omit(p_values)
N <- length(PV_clean)
df_qq <- data.frame(
  observed = sort(PV_clean),
  expected = (1:N) / (N + 1)
)
df_qq$ci_lower <- qbeta(0.025, 1:N, N - 1:N + 1)
df_qq$ci_upper <- qbeta(0.975, 1:N, N - 1:N + 1)
df_qq$log_obs <- -log10(df_qq$observed)
df_qq$log_exp <- -log10(df_qq$expected)
df_qq$log_ci_lower <- -log10(df_qq$ci_upper)
df_qq$log_ci_upper <- -log10(df_qq$ci_lower)

ggplot(df_qq, aes(x = log_exp, y = log_obs)) +
  geom_ribbon(aes(ymin = log_ci_lower, ymax = log_ci_upper),
              fill = "grey70", alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  geom_point(shape = 21, fill = "black", color = "white", size = 2, alpha = 0.8) +
  labs(
    x = expression(Expected~~-log[10](italic(p))),
    y = expression(Observed~~-log[10](italic(p))),
    title = "Q-Q Plot of P-values with 95% Confidence Interval"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

