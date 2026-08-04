library(devtools)
load_all(".")
library(mgcv)
library(gamm4)
library(jmcm)

data(aids)
dat <- aids
dat$id <- factor(dat$id)
dat$sqrt_cd4 <- sqrt(dat$cd4)
dat$time_s <- as.numeric(scale(dat$time))

tv <- c("packs", "drugs", "sex", "cesd")
within_summary <- data.frame(
  variable = tv,
  varying_subjects = NA_integer_,
  total_subjects = length(unique(dat$id)),
  median_unique_values = NA_real_
)

for (i in seq_along(tv)) {
  v <- tv[i]
  n_unique <- tapply(dat[[v]], dat$id, function(x) length(unique(x)))
  within_summary$varying_subjects[i] <- sum(n_unique > 1)
  within_summary$median_unique_values[i] <- median(n_unique)
  b <- ave(dat[[v]], dat$id, FUN = mean)
  w <- dat[[v]] - b
  dat[[paste0(v, "_b")]] <- as.numeric(scale(b))
  dat[[paste0(v, "_w")]] <- as.numeric(scale(w))
}

getA_t <- function(x, para) cbind(1, x)

m_gamm4_base <- gamm4(
  sqrt_cd4 ~ s(
    time_s, bs = "AMatern", k = 12,
    xt = list(getA = getA_t, para = NULL)
  ) + age + packs + drugs + sex + cesd,
  random = ~(1 + time_s | id),
  data = dat,
  REML = TRUE
)

m_gam_base <- gam(
  sqrt_cd4 ~ s(
    time_s, bs = "AMatern", k = 12,
    xt = list(getA = getA_t, para = NULL)
  ) + age + packs + drugs + sex + cesd +
    s(id, bs = "re") + s(time_s, id, bs = "re"),
  data = dat,
  method = "REML"
)

m_gamm4_tv <- gamm4(
  sqrt_cd4 ~ s(
    time_s, bs = "AMatern", k = 12,
    xt = list(getA = getA_t, para = NULL)
  ) + age + packs_b + packs_w + drugs_b + drugs_w +
    sex_b + sex_w + cesd_b + cesd_w,
  random = ~(1 + time_s + packs_w + drugs_w + sex_w + cesd_w || id),
  data = dat,
  REML = TRUE
)

score <- rbind(
  data.frame(model = "gamm4_base", taps_score_test(m_gamm4_base)),
  data.frame(model = "gam_base", taps_score_test(m_gam_base)),
  data.frame(model = "gamm4_time_varying", taps_score_test(m_gamm4_tv))
)

vc_gamm4_base <- as.data.frame(VarCorr(m_gamm4_base$mer))
vc_gamm4_base$model <- "gamm4_base"
vc_gamm4_tv <- as.data.frame(VarCorr(m_gamm4_tv$mer))
vc_gamm4_tv$model <- "gamm4_time_varying"
vc_gamm4 <- rbind(vc_gamm4_base, vc_gamm4_tv)
vc_gamm4 <- vc_gamm4[, c("model", setdiff(names(vc_gamm4), "model"))]

vc_gam_base <- mgcv::gam.vcomp(m_gam_base)$vc
vc_gam_base <- as.data.frame(vc_gam_base)
vc_gam_base$component <- rownames(vc_gam_base)
vc_gam_base$model <- "gam_base"
vc_gam_base <- vc_gam_base[, c("model", "component", setdiff(names(vc_gam_base), c("model", "component")))]

fixed_gamm4_base <- data.frame(
  model = "gamm4_base",
  term = names(fixef(m_gamm4_base$mer)),
  estimate = as.numeric(fixef(m_gamm4_base$mer))
)
fixed_gamm4_tv <- data.frame(
  model = "gamm4_time_varying",
  term = names(fixef(m_gamm4_tv$mer)),
  estimate = as.numeric(fixef(m_gamm4_tv$mer))
)
fixed_gam_base <- as.data.frame(summary(m_gam_base)$p.table)
fixed_gam_base$term <- rownames(fixed_gam_base)
fixed_gam_base$model <- "gam_base"
fixed_gam_base <- fixed_gam_base[, c("model", "term", setdiff(names(fixed_gam_base), c("model", "term")))]

write.csv(within_summary, "example/gamm4_aids_within_variation.csv", row.names = FALSE)
write.csv(score, "example/gamm4_aids_score_tests.csv", row.names = FALSE)
write.csv(vc_gamm4, "example/gamm4_aids_gamm4_variance_components.csv", row.names = FALSE)
write.csv(vc_gam_base, "example/gamm4_aids_gam_variance_components.csv", row.names = FALSE)
write.csv(rbind(fixed_gamm4_base, fixed_gamm4_tv), "example/gamm4_aids_gamm4_fixed_effects.csv", row.names = FALSE)
write.csv(fixed_gam_base, "example/gamm4_aids_gam_fixed_effects.csv", row.names = FALSE)
