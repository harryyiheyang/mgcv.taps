require(mgcv)
require(MASS)
devtools::document()
n=1000
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

rzip <- function(gamma,theta= c(-2,.3)) {
  ## generate zero inflated Poisson random variables, where
  ## lambda = exp(gamma), eta = theta[1] + exp(theta[2])*gamma
  ## and 1-p = exp(-exp(eta)).
  y <- gamma; n <- length(y)
  lambda <- exp(gamma)
  eta <- theta[1] + exp(theta[2])*gamma
  p <- 1- exp(-exp(eta))
  ind <- p > runif(n)
  y[!ind] <- 0
  np <- sum(ind)
  ## generate from zero truncated Poisson, given presence...
  y[ind] <- qpois(runif(np,dpois(0,lambda[ind]),1),lambda[ind])
  y
}
# Function to manually verify mgcv's ziP pseudo-data against the math from the PDF
check_zip_pseudo_data <- function(fit) {
  y <- as.numeric(fit$y)
  gamma <- as.numeric(fit$linear.predictors) # This is the single linear predictor
  theta <- fit$family$getTheta(TRUE)         # theta[1] = theta_1, theta[2] = theta_2

  # From PDF: eta = theta_1 + e^{theta_2} * gamma
  eta <- theta[1] + exp(theta[2]) * gamma
  eta_gamma <- exp(theta[2]) # Derivative of eta w.r.t gamma

  l_gamma <- numeric(length(y))
  l_gammagamma <- numeric(length(y))

  # --- 1. Derivatives for y == 0 ---
  idx0 <- (y == 0)
  if (any(idx0)) {
    eta0 <- eta[idx0]
    # PDF limit: l = -exp(eta)
    l_eta <- -exp(eta0)
    l_etaeta <- -exp(eta0)

    l_gamma[idx0] <- l_eta * eta_gamma
    l_gammagamma[idx0] <- l_etaeta * (eta_gamma^2)
  }

  # --- 2. Derivatives for y > 0 ---
  idx1 <- (y > 0)
  if (any(idx1)) {
    y1 <- y[idx1]
    eta1 <- eta[idx1]
    gam1 <- gamma[idx1]

    # Partial derivatives w.r.t eta (from PDF)
    l_eta <- exp(eta1) / (exp(exp(eta1)) - 1)
    l_etaeta <- (1 - exp(eta1)) * l_eta - l_eta^2

    # Partial derivatives w.r.t gamma (from PDF)
    alpha <- exp(gam1) / (1 - exp(-exp(gam1)))
    l_gam_direct <- y1 - alpha
    l_gamgam_direct <- alpha^2 - (exp(gam1) + 1) * alpha

    # Total derivatives applying the chain rule
    l_gamma[idx1] <- l_gam_direct + l_eta * eta_gamma
    l_gammagamma[idx1] <- l_gamgam_direct + l_etaeta * (eta_gamma^2)
  }

  # --- 3. The Working Weights and Pseudo-Response ---
  # PIRLS step defines weights as the negative second derivative of the log-likelihood
  w_manual <- -l_gammagamma

  # Pseudo-response z = gamma + (first derivative) / (working weight)
  # which simplifies to z = gamma - l_gamma / l_gammagamma
  z_manual <- gamma + l_gamma / w_manual

  # --- 4. Extract directly from the mgcv gam object ---
  w_mgcv <- fit$weights
  z_mgcv <- fit$linear.predictors + fit$residuals

  # Compare them!
  data.frame(
    Max_Weight_Difference = max(abs(w_manual - w_mgcv), na.rm = TRUE),
    Max_PseudoRes_Difference = max(abs(z_manual - z_mgcv), na.rm = TRUE)
  )
}

PV=c(1:200)
for(i in 1:200){
d   <- gen_data(n)
y   <- rzip(d$eta/3-1)
dat <- data.frame(x0=d$x1, x1=d$x2, x2=d$x3, x3=d$x4, y=y)
b   <- gam(y ~ s(x0,bs="AMatern") + s(x1,bs="cr") +
             s(x2,bs="cr")      + s(x3,bs="cr"),
           data=dat, family=ziP(), method="REML")
PV[i]=taps_score_test_extended(b,test.component=1, method="davies")$smooth.pvalue
}

library(ggplot2)
PV_clean <- na.omit(PV)
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
