library(mgcv.taps)

set.seed(8401)
ng <- 5L
k <- 2L
n_i <- c(7L, 8L, 6L, 9L, 7L)
id <- rep(seq_len(ng), n_i)
n <- length(id)
x <- runif(n, -1, 1)
t <- runif(n)
X <- cbind(1, x, x^2)
B <- cbind(1, t)
G <- matrix(c(0.7, 0.12, 0.12, 0.35), 2L, 2L)
R_diag <- runif(n, 0.7, 1.4)
S <- diag(c(0, 0.4, 1.1))
response <- rnorm(n)
beta <- rnorm(ncol(X))
empty_probes <- matrix(numeric(), 0L, 0L)
cpp <- mgcv.taps:::gammfast_variance_quadratic(
  response, X, beta, B, id, G, R_diag, S, empty_probes,
  exact = TRUE, n_threads = 2L
)
cpp_diagonal <- mgcv.taps:::gammfast_variance_quadratic(
  response, X, beta, B, id, diag(diag(G)), R_diag, S, empty_probes,
  exact = TRUE, n_threads = 2L
)
random_cpp <- mgcv.taps:::gammfast_variance_quadratic(
  response, X, beta, B, id, G, R_diag, S,
  sqrt(ng * k) * diag(ng * k), exact = FALSE, n_threads = 2L
)

se <- eigen(S, symmetric = TRUE)
positive <- se$values > 1e-8
zero <- !positive
Us <- X %*% se$vectors[, positive, drop = FALSE]
Us <- sweep(Us, 2L, sqrt(se$values[positive]), "/")
Xf <- X %*% se$vectors[, zero, drop = FALSE]
A0 <- diag(R_diag) + tcrossprod(Us)
A0_inv <- solve(A0)
P0 <- A0_inv - A0_inv %*% Xf %*%
  solve(crossprod(Xf, A0_inv %*% Xf)) %*% t(Xf) %*% A0_inv
P0 <- (P0 + t(P0)) / 2

Lg <- chol(G)
L <- matrix(0, n, ng * k)
for (i in seq_len(n)) {
  jj <- (id[i] - 1L) * k + seq_len(k)
  L[i, jj] <- B[i, ] %*% t(Lg)
}
residual <- response - Xf %*% crossprod(se$vectors[, zero, drop = FALSE], beta)
A <- crossprod(L, P0 %*% L)
lambda <- eigen((A + t(A)) / 2, symmetric = TRUE, only.values = TRUE)$values
statistic <- drop(crossprod(residual, P0 %*% L %*% t(L) %*% P0 %*% residual))

stopifnot(max(abs(sort(cpp$lambda) - sort(lambda))) < 1e-9)
stopifnot(abs(cpp$statistic - statistic) < 1e-9)
stopifnot(abs(cpp$statistic - cpp_diagonal$statistic) > 1e-6)
stopifnot(abs(random_cpp$statistic - statistic) < 1e-9)
stopifnot(max(abs(random_cpp$moments - vapply(
  seq_len(4L), function(r) sum(lambda^r), numeric(1)
))) < 1e-8)

set.seed(20260811)
n_i <- rep(c(3L, 5L, 7L, 4L), 12L)
id_fit <- factor(rep(seq_along(n_i), n_i))
time_fit <- unlist(lapply(n_i, function(m) seq(0, 1, length.out = m)))
x_fit <- runif(length(id_fit))
u_fit <- rnorm(nlevels(id_fit), sd = 0.2)
fs_fit <- matrix(rnorm(nlevels(id_fit) * 3L, sd = 0.12),
                 nlevels(id_fit), 3L)
C_fit <- cbind(
  sqrt(2) * cos(pi * time_fit),
  sqrt(2) * cos(2 * pi * time_fit),
  sqrt(2) * cos(3 * pi * time_fit)
)
y_fit <- 1 + 0.4 * x_fit + 0.5 * sin(2 * pi * x_fit) +
  u_fit[id_fit] + rowSums(C_fit * fs_fit[id_fit, ]) +
  rnorm(length(id_fit), sd = 0.5)
d_fit <- data.frame(y = y_fit, x = x_fit, time = time_fit, id = id_fit)
fit <- gammfast(
  y ~ s(x, k = 5L) + s(id, bs = "re") + fs(time, id, k = 3L),
  data = d_fit, inner.max = 5L, nthreads = 2L,
  control = list(max.outer = 1000L, objective.tol = 1e-6)
)
random_fit <- mgcv.taps:::gammfast_random_info(fit)
if (!isTRUE(fit$converged)) stop("Gaussian gammfast fixture did not converge.")
variance_test <- gammfast_variance_test(fit, method = "liu", n_threads = 2L)
fit_diagonal <- fit
fit_diagonal$G <- diag(diag(fit$G))
variance_test_diagonal <- gammfast_variance_test(
  fit_diagonal, method = "liu", n_threads = 2L
)
mean_test <- taps_score_test_gamm(
  fit, method = "liu", n_threads = 2L
)
generic_error <- FALSE
tryCatch(
  taps_score_test(fit, method = "liu", n_threads = 2L),
  error = function(e) {
    generic_error <<- grepl(
      "use taps_score_test_gamm", conditionMessage(e), fixed = TRUE
    )
  }
)
set.seed(8403)
rng_before <- .Random.seed
variance_random <- gammfast_variance_test(
  fit, method = "auto", spectrum = "randomized", n_probe = 64L,
  seed = 8404L, n_threads = 2L
)
fit_summary <- summary(
  fit, method = "liu", spectrum = "exact", n_threads = 2L
)
stopifnot(
  nrow(variance_test) == 1L,
  is.finite(variance_test$statistic),
  variance_test$p.value >= 0,
  variance_test$p.value <= 1,
  abs(variance_test$statistic - variance_test_diagonal$statistic) > 1e-6,
  variance_test$requested.method == "liu",
  variance_test$method == "liu",
  !variance_test$fallback,
  is.na(variance_test$fallback.from),
  variance_test$requested.spectrum == "auto",
  !variance_test$null.refit,
  !variance_test$full.random.design,
  attr(mean_test, "backend") == "gammfast",
  identical(attr(mean_test, "full.random.design"), FALSE),
  generic_error,
  variance_random$method == "randomized-liu",
  variance_random$requested.method == "auto",
  variance_random$requested.spectrum == "randomized",
  !variance_random$fallback,
  variance_random$spectrum == "randomized",
  variance_random$n.probe == 64L,
  is.finite(variance_random$p.value),
  abs(variance_random$statistic - variance_test$statistic) < 1e-8,
  inherits(fit_summary, "summary.gammfast"),
  inherits(fit_summary$mean.fit, "summary.gam"),
  identical(fit_summary$mean.inference.source$full.random.design, FALSE),
  is.matrix(fit_summary$mean.tables$parametric),
  is.matrix(fit_summary$mean.tables$smooth),
  fit_summary$n == nrow(d_fit),
  fit_summary$n.subject == nlevels(id_fit),
  fit_summary$basis.dimension == 4L,
  abs(fit_summary$variance.test$statistic - variance_test$statistic) < 1e-8,
  !fit_summary$variance.test$null.refit,
  identical(.Random.seed, rng_before)
)

fallback_test <- mgcv.taps:::gammfast_quadratic_pvalue(
  statistic = 1e10, lambda = c(1, 0.5), method = "davies"
)
if (!fallback_test$fallback ||
    fallback_test$method != "liu-fallback" ||
    fallback_test$fallback.from != "davies" ||
    fallback_test$p.value < 0 || fallback_test$p.value > 1) {
  stop("Davies numerical underflow did not retain Liu fallback provenance.")
}

summary_fallback <- fit_summary
summary_fallback$variance.test$requested.method <- "davies"
summary_fallback$variance.test$method <- "liu-fallback"
summary_fallback$variance.test$fallback <- TRUE
summary_fallback$variance.test$fallback.from <- "davies"
printed_summary <- capture.output(print(summary_fallback))
if (!any(grepl("Mean-structure fit", printed_summary, fixed = TRUE)) ||
    !any(grepl("Post-estimation variance-component test", printed_summary,
               fixed = TRUE)) ||
    !any(grepl("p.value", printed_summary, fixed = TRUE)) ||
    any(grepl("Liu|Davies|fallback|spectrum", printed_summary,
              ignore.case = TRUE))) {
  stop("summary.gammfast print must show the variance-test p-value without method provenance.")
}

ids_fit <- random_fit$id.levels
selected_ids <- c(ids_fit[2L], ids_fit[1L], ids_fit[2L])
blup_old <- gammblup(fit, ids = selected_ids)
if (!identical(blup_old$subjects$individual, selected_ids) ||
    nrow(blup_old$subjects) != 3L ||
    !identical(blup_old$full.random.design, FALSE) ||
    !identical(blup_old$blup[1L, ], blup_old$blup[3L, ]) ||
    !isTRUE(all.equal(
      unname(blup_old$rows$random.link[
        blup_old$rows$request.index == 1L
      ]),
      unname(fit$random.fitted[random_fit$id.index == 2L]),
      tolerance = 1e-12
    ))) {
  stop("gammblup did not preserve fitted-ID selection, order, or duplicates.")
}

known_rows <- d_fit$id %in% selected_ids[1:2]
new_known <- d_fit[known_rows, , drop = FALSE]
new_known <- new_known[rev(seq_len(nrow(new_known))), , drop = FALSE]
blup_known <- gammblup(fit, ids = selected_ids, newdata = new_known)
if (!identical(blup_known$subjects$individual, selected_ids) ||
    !identical(
      blup_known$rows$data.row[blup_known$rows$request.index == 1L],
      which(as.character(new_known$id) == selected_ids[1L])
    )) {
  stop("gammblup did not preserve newdata row order.")
}

new_zero <- data.frame(
  time = c(-0.2, 0.3, 1.2),
  id = c("new-zero", "new-zero", "new-zero")
)
blup_zero <- gammblup(
  fit, ids = "new-zero", newdata = new_zero, new.level = "zero"
)
if (any(blup_zero$blup != 0) ||
    !all(blup_zero$rows$time.clamped[c(1L, 3L)]) ||
    blup_zero$subjects$source != "zero") {
  stop("gammblup new.level = 'zero' handling is incorrect.")
}

new_estimate <- data.frame(
  y = c(1.1, 1.4, 1.0, 1.5),
  x = c(0.1, 0.35, 0.7, 0.9),
  time = c(0, 0.3, 0.7, 1),
  id = rep("new-estimate", 4L)
)
blup_estimate <- gammblup(
  fit, ids = "new-estimate", newdata = new_estimate,
  new.level = "estimate", n_threads = 2L
)
g_new <- fit$global
g_new$coefficients <- fit$coefficients
eta_new <- as.numeric(predict(g_new, newdata = new_estimate, type = "link"))
B_new <- mgcv.taps:::gammfast_random_design(random_fit, new_estimate)$B
D_new <- solve(solve(fit$G) + crossprod(B_new) / fit$sig2)
u_new <- drop(D_new %*% crossprod(B_new, new_estimate$y - eta_new) /
                fit$sig2)
if (!isTRUE(all.equal(unname(drop(blup_estimate$blup)), unname(u_new),
                      tolerance = 1e-10)) ||
    blup_estimate$subjects$source != "estimated") {
  stop("gammblup new-ID conditional BLUP does not match the dense individual solve.")
}

missing_id_error <- FALSE
tryCatch(
  gammblup(fit),
  error = function(e) {
    missing_id_error <<- grepl("ids is required", conditionMessage(e),
                               fixed = TRUE)
  }
)
missing_variable_error <- FALSE
new_missing <- new_estimate[c("y", "time", "id")]
new_missing$id <- "new-missing"
tryCatch(
  gammblup(
    fit, ids = "new-missing",
    newdata = new_missing,
    new.level = "estimate"
  ),
  error = function(e) {
    missing_variable_error <<- grepl("missing variables", conditionMessage(e),
                                     fixed = TRUE)
  }
)
if (!missing_id_error || !missing_variable_error) {
  stop("gammblup did not reject missing ID selection or mean variables clearly.")
}

phi <- fit$sig2
X_fit <- predict(fit$global, newdata = fit$global$model,
                 type = "lpmatrix") / sqrt(phi)
B_fit <- random_fit$B / sqrt(phi)
p_fit <- ncol(X_fit)
q_fit <- nrow(fit$random.effects) * ncol(B_fit)
Z_fit <- matrix(0, nrow(X_fit), q_fit)
for (i in seq_len(nrow(X_fit))) {
  jj <- (random_fit$id.index[i] - 1L) * ncol(B_fit) + seq_len(ncol(B_fit))
  Z_fit[i, jj] <- B_fit[i, ]
}
S_fit <- matrix(0, p_fit, p_fit)
for (j in seq_along(fit$global$S)) {
  ii <- fit$global$off[j] + seq_len(nrow(fit$global$S[[j]])) - 1L
  S_fit[ii, ii] <- S_fit[ii, ii] +
    fit$global$sp[j] * fit$global$S[[j]] / phi
}
P_fit <- matrix(0, p_fit + q_fit, p_fit + q_fit)
P_fit[seq_len(p_fit), seq_len(p_fit)] <- S_fit
G_fit_inv <- solve(fit$G)
for (i in seq_len(nrow(fit$random.effects))) {
  jj <- p_fit + (i - 1L) * ncol(B_fit) + seq_len(ncol(B_fit))
  P_fit[jj, jj] <- G_fit_inv
}
D_fit <- cbind(X_fit, Z_fit)
y_std <- (fit$y - fit$offset) / sqrt(phi)
theta_fit <- solve(crossprod(D_fit) + P_fit, crossprod(D_fit, y_std))
u_joint <- matrix(
  theta_fit[p_fit + seq_len(q_fit)],
  nrow = ncol(B_fit), ncol = nrow(fit$random.effects)
)
u_joint <- t(u_joint)
stopifnot(max(abs(u_joint - fit$random.effects)) < 2e-5)
