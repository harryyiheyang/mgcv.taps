#' Summarize a gammfast fit
#'
#' Prints an mgcv-style summary for the global mean structure, followed by the
#' fitted random-trajectory covariance and a separate post-estimation
#' variance-component test. Mean inference is reconstructed from the final
#' `gammfast` coefficients, `Vp`, `Ve`, and smoothing parameters in the
#' mean-only mgcv shell. No subject-level random-effect design matrix is made.
#'
#' @param object A fitted `gammfast` object.
#' @param variance.test Whether to calculate the post-estimation variance test.
#' @param method,spectrum,q_threshold,n_probe,seed,max_eps,max_iter,n_threads
#'   Passed to [gammfast_variance_test()].
#' @param ... Unused.
#'
#' @return An object of class `summary.gammfast`. `mean.fit` contains the
#'   mgcv-style mean summary and `variance.test` retains the unmodified
#'   provenance table returned by [gammfast_variance_test()], including the
#'   requested and actual method, fallback flag, spectrum, conditional flag,
#'   and `null.refit` flag.
#' @export
summary.gammfast <- function(object, variance.test = TRUE,
                             method = c("auto", "davies", "liu"),
                             spectrum = c("auto", "exact", "randomized"),
                             q_threshold = 1000L,
                             n_probe = 100L, seed = 20260810L,
                             max_eps = 1e-8, max_iter = 1e5,
                             n_threads = 1L, ...) {
  if (!inherits(object, "gammfast")) {
    stop("object must be a 'gammfast' fit.")
  }
  if (length(variance.test) != 1L || is.na(variance.test)) {
    stop("variance.test must be TRUE or FALSE.")
  }
  method <- match.arg(method)
  spectrum <- match.arg(spectrum)

  mean_fit <- gammfast_mean_summary(object)
  random <- gammfast_random_info(object)
  random_eigen <- lapply(random$group.index, function(jj) {
    eigen(object$G[jj, jj, drop = FALSE], symmetric = TRUE,
          only.values = TRUE)$values
  })
  variance <- NULL
  if (isTRUE(variance.test)) {
    variance <- gammfast_variance_test(
      object,
      method = method,
      spectrum = spectrum,
      q_threshold = q_threshold,
      n_probe = n_probe,
      seed = seed,
      max_eps = max_eps,
      max_iter = max_iter,
      n_threads = n_threads
    )
  }

  ans <- list(
    call = object$call,
    formula = object$formula,
    global.formula = object$global.formula,
    family = object$family$family,
    n = length(object$y),
    n.subject = nrow(object$random.effects),
    basis.dimension = ncol(object$random.effects),
    sig2 = object$sig2,
    smoothing.parameters = object$sp,
    mean.fit = mean_fit,
    mean.tables = list(
      parametric = mean_fit$p.table,
      smooth = mean_fit$s.table
    ),
    mean.inference.source = mean_fit$inference.source,
    random.covariance = object$G,
    random.groups = data.frame(
      group = seq_along(random$groups),
      type = vapply(random$groups, `[[`, character(1), "type"),
      label = vapply(random$groups, `[[`, character(1), "label"),
      dimension = lengths(random$group.index),
      variance.trace = vapply(random_eigen, sum, numeric(1)),
      minimum.eigenvalue = vapply(random_eigen, min, numeric(1)),
      maximum.eigenvalue = vapply(random_eigen, max, numeric(1)),
      stringsAsFactors = FALSE
    ),
    variance.test = variance,
    converged = object$converged,
    outer = object$outer
  )
  class(ans) <- "summary.gammfast"
  ans
}

gammfast_mean_summary <- function(object) {
  g <- object$global
  beta <- as.numeric(object$coefficients)
  names(beta) <- names(object$coefficients)
  Vp <- as.matrix(object$Vp)
  Ve <- as.matrix(object$Ve)
  p <- length(beta)
  if (nrow(Vp) != p || ncol(Vp) != p ||
      nrow(Ve) != p || ncol(Ve) != p ||
      any(!is.finite(Vp)) || any(!is.finite(Ve))) {
    stop("The gammfast coefficient covariance matrices are invalid.")
  }

  F <- matrixMultiply(Ve, matrixGeneralizedInverse(Vp))
  edf <- diag(F)
  edf1 <- 2 * edf - rowSums(t(F) * F)
  names(edf) <- names(edf1) <- names(beta)
  residual_df <- length(object$y) - sum(edf)
  if (!is.finite(residual_df) || residual_df <= 0) {
    stop("The gammfast mean summary has non-positive residual degrees of freedom.")
  }

  family <- object$family
  family_name <- tolower(family$family[1L])
  family_scale <- family$scale
  estimated_dispersion <- identical(family_name, "gaussian") ||
    (is.numeric(family_scale) && length(family_scale) == 1L &&
       is.finite(family_scale) && family_scale < 0)
  se <- sqrt(pmax(0, diag(Vp)))
  nsdf <- min(p, sum(g$nsdf))
  if (nsdf > 0L) {
    ii <- seq_len(nsdf)
    p_coeff <- beta[ii]
    p_se <- se[ii]
    p_stat <- p_coeff / p_se
    if (estimated_dispersion) {
      p_pv <- 2 * stats::pt(abs(p_stat), df = residual_df,
                            lower.tail = FALSE)
      p_table <- cbind(p_coeff, p_se, p_stat, p_pv)
      colnames(p_table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    } else {
      p_pv <- 2 * stats::pnorm(abs(p_stat), lower.tail = FALSE)
      p_table <- cbind(p_coeff, p_se, p_stat, p_pv)
      colnames(p_table) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    }
    rownames(p_table) <- names(p_coeff)
  } else {
    p_coeff <- p_stat <- p_pv <- numeric()
    p_table <- NULL
  }

  smooth_terms <- g$smooth
  m <- length(smooth_terms)
  s_table <- NULL
  chi_sq <- s_pv <- smooth_edf <- numeric()
  if (m > 0L) {
    X <- as.matrix(g$X)
    if (ncol(X) != p || nrow(X) != length(object$y)) {
      stop("The gammfast global mean design is incompatible with its coefficients.")
    }
    test_stat <- utils::getFromNamespace("testStat", "mgcv")
    smooth_edf <- ref_df <- chi_sq <- s_pv <- numeric(m)
    smooth_names <- character(m)
    for (i in seq_len(m)) {
      s <- smooth_terms[[i]]
      ii <- s$first.para:s$last.para
      smooth_edf[i] <- sum(edf[ii])
      ref_edf <- min(ncol(X[, ii, drop = FALSE]), sum(edf1[ii]))
      result <- test_stat(
        beta[ii], X[, ii, drop = FALSE], Vp[ii, ii, drop = FALSE],
        ref_edf, type = 0,
        res.df = if (estimated_dispersion) residual_df else -1
      )
      ref_df[i] <- result$rank
      chi_sq[i] <- result$stat
      s_pv[i] <- result$pval
      smooth_names[i] <- s$label
    }
    if (estimated_dispersion) {
      s_table <- cbind(smooth_edf, ref_df, chi_sq / ref_df, s_pv)
      colnames(s_table) <- c("edf", "Ref.df", "F", "p-value")
    } else {
      s_table <- cbind(smooth_edf, ref_df, chi_sq, s_pv)
      colnames(s_table) <- c("edf", "Ref.df", "Chi.sq", "p-value")
    }
    rownames(s_table) <- smooth_names
    names(chi_sq) <- smooth_names
  }

  prior_weights <- object$prior.weights
  if (is.null(prior_weights)) prior_weights <- rep(1, length(object$y))
  prior_weights <- as.numeric(prior_weights)
  eta_mean <- as.numeric(object$global.fitted)
  if (length(eta_mean) != length(object$y)) {
    stop("The gammfast global fitted values have an invalid length.")
  }
  mu_mean <- family$linkinv(eta_mean)
  mean_y <- stats::weighted.mean(as.numeric(object$y), prior_weights)
  null_mu <- rep(mean_y, length(object$y))
  deviance <- sum(family$dev.resids(object$y, mu_mean, prior_weights))
  null_deviance <- sum(family$dev.resids(
    object$y, null_mu, prior_weights
  ))
  dev_expl <- if (is.finite(null_deviance) && null_deviance > 0) {
    (null_deviance - deviance) / null_deviance
  } else {
    numeric()
  }
  r_sq <- if (inherits(family, "general.family") || !is.null(family$no.r.sq)) {
    NULL
  } else {
    w <- sqrt(prior_weights)
    n <- length(object$y)
    1 - stats::var(w * (as.numeric(object$y) - mu_mean)) * (n - 1) /
      (stats::var(w * (as.numeric(object$y) - mean_y)) * residual_df)
  }

  out <- list(
    p.coeff = p_coeff,
    se = se,
    p.t = p_stat,
    p.pv = p_pv,
    residual.df = residual_df,
    m = m,
    chi.sq = chi_sq,
    s.pv = s_pv,
    scale = object$sig2,
    r.sq = r_sq,
    family = family,
    formula = object$global.formula,
    n = length(object$y),
    dev.expl = dev_expl,
    edf = smooth_edf,
    dispersion = object$sig2,
    pTerms.pv = numeric(),
    pTerms.chi.sq = numeric(),
    pTerms.df = numeric(),
    cov.unscaled = Vp / object$sig2,
    cov.scaled = Vp,
    p.table = p_table,
    pTerms.table = NULL,
    s.table = s_table,
    method = "gammfast fREML",
    sp.criterion = object$objective,
    rank = qr(F)$rank,
    np = p,
    inference.source = list(
      coefficients = "gammfast$coefficients",
      covariance = "gammfast$Vp and gammfast$Ve",
      smoothing.parameters = "gammfast$sp",
      shell = "global mean-only mgcv setup",
      full.random.design = FALSE
    )
  )
  class(out) <- "summary.gam"
  out
}

#' @export
print.summary.gammfast <- function(x, digits = max(3L, getOption("digits") - 3L),
                                   ...) {
  cat("Mean-structure fit (mgcv style; global mean-only shell):\n\n")
  print(x$mean.fit, digits = digits, ...)
  cat(
    "Mean inference source: final gammfast coefficients/Vp/Ve/sp; ",
    "no subject-level random-effect design matrix.\n",
    sep = ""
  )
  cat("\nRandom-trajectory covariance:\n")
  print(x$random.groups, row.names = FALSE)
  print(x$random.covariance, digits = digits)
  if (!is.null(x$variance.test)) {
    shown <- data.frame(
      component = x$variance.test$component,
      statistic = x$variance.test$statistic,
      p.value = x$variance.test$p.value,
      check.names = FALSE
    )
    cat("\nPost-estimation variance-component test:\n")
    print(shown, row.names = FALSE, digits = digits)
  }
  invisible(x)
}
