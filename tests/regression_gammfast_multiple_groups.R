library(mgcv)
library(mgcv.taps)

set.seed(20260815)
n_id <- 30L
n_each <- 8L
id <- factor(rep(seq_len(n_id), each = n_each))
id2 <- factor(rep(rep(seq_len(5L), each = 6L), each = n_each))
t <- rep(seq(0, 1, length.out = n_each), n_id)
z <- runif(length(id))
w <- runif(length(id))
u0 <- rnorm(n_id, sd = 0.35)
u1 <- rnorm(n_id, sd = 0.20)
a <- matrix(rnorm(n_id * 4L, sd = 0.12), n_id, 4L)
C1 <- cbind(sqrt(2) * cos(pi * t), sqrt(2) * cos(2 * pi * t))
C2 <- cbind(sqrt(2) * cos(pi * w), sqrt(2) * cos(2 * pi * w))
y <- 1 + sin(2 * pi * z) + u0[id] + u1[id] * z +
  rowSums(C1 * a[id, 1:2]) + rowSums(C2 * a[id, 3:4]) +
  rnorm(length(id), sd = 0.4)
dat <- data.frame(y = y, z = z, w = w, t = t, id = id, id2 = id2)

parsed_bare <- mgcv.taps:::gammfast_parse_formula(
  y ~ s(z, k = 6L) + w + s(id, bs = "re"), dat
)
if (!"w" %in% attr(terms(parsed_bare$formula), "term.labels")) {
  stop("gammfast dropped a bare linear fixed-effect term.")
}
parsed_dt <- mgcv.taps:::gammfast_parse_formula(
  y ~ s(z, k = 6L) + w + s(id, bs = "re"),
  data.table::as.data.table(dat)
)
if (!identical(parsed_dt$id, "id")) {
  stop("gammfast did not parse a data.table input correctly.")
}

fit <- gammfast(
  y ~ s(z, k = 6L) + s(id, bs = "re") +
    s(z, id, bs = "re") + fs(t, id, k = 2L) + fs(w, id, k = 2L),
  data = dat, inner.max = 5L, nthreads = 2L,
  control = list(max.outer = 600L, objective.tol = 1e-6)
)
if (!isTRUE(fit$converged)) {
  stop("The multiple-group gammfast fixture did not converge.")
}

random <- mgcv.taps:::gammfast_random_info(fit)
if (length(random$groups) != 4L || ncol(random$B) != 6L ||
    !identical(vapply(random$groups, `[[`, character(1), "type"),
               c("random_intercept", "random_slope", "fs_cosine", "fs_cosine")) ||
    any(grepl("cos0", random$column.names, fixed = TRUE))) {
  stop("The multiple random-group structure was not constructed correctly.")
}

off_block <- fit$G
for (jj in random$group.index) off_block[jj, jj] <- 0
minimum_eigenvalue <- min(unlist(lapply(random$group.index, function(jj) {
  eigen(fit$G[jj, jj, drop = FALSE], symmetric = TRUE,
        only.values = TRUE)$values
})))
if (max(abs(off_block)) != 0 || minimum_eigenvalue <= 0) {
  stop("gammfast random groups are not independent positive-definite blocks.")
}

selected <- gammblup(fit, ids = c("3", "1", "3"))
if (!identical(selected$requested.ids, c("3", "1", "3")) ||
    nrow(selected$blup) != 3L || ncol(selected$blup) != 6L ||
    length(selected$groups) != 4L) {
  stop("gammblup did not preserve multiple-group ID selection and order.")
}

expect_error_text <- function(expr, pattern) {
  matched <- FALSE
  tryCatch(
    force(expr),
    error = function(e) {
      matched <<- grepl(pattern, conditionMessage(e), fixed = TRUE)
    }
  )
  if (!matched) stop("Expected gammfast error was not produced: ", pattern)
}

expect_error_text(
  gammfast(y ~ s(z, k = 6L) + fs(t, id, k = 2L), data = dat),
  "explicit random intercept"
)
expect_error_text(
  gammfast(
    y ~ s(z, k = 6L) + s(id, bs = "re") +
      s(t, id, bs = "fs", k = 2L),
    data = dat
  ),
  "does not support mgcv s(..., bs = 'fs')"
)
expect_error_text(
  gammfast(
    y ~ s(z, k = 6L) + s(id, bs = "re") + fs(t, id2, k = 2L),
    data = dat
  ),
  "share the same factor ID"
)
expect_error_text(
  gammfast(
    y ~ s(z, k = 6L) + s(id, bs = "re"), data = dat,
    discrete = TRUE
  ),
  "does not support discrete = TRUE"
)

cat("gammfast multiple-group regression checks passed.\n")
