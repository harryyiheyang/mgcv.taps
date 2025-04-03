#' Perform Wald Tests on mgcv Smooth Terms (Null and Penalized Components)
#'
#' This function takes a fitted mgcv `gam` or `bam` object,
#' then performs separate Wald tests on the null (unpenalized/fixed) and penalized
#' (smooth) components of each smooth term.
#' It returns a clear summary as a `data.table`.
#'
#' @param fit A fitted `gam` or `bam` model from the `mgcv` package.
#' @param test.component The index of smooth component to be tested.
#' @param null.tol The tolerance of row norm to detect indices of null space. Default to 1e-10.
#' @return A `data.table` summarizing Wald statistics and p-values separately
#'         for null (unpenalized) and penalized (smooth) components of each term.
#'
#' @importFrom mgcv vcov.gam
#' @importFrom data.table data.table rbindlist
#' @importFrom ordgam testStat
#' @importFrom stats pchisq
#' @import CppMatrix
#'
#' @examples
#' \dontrun{
#' fit <- gam(y ~ s(x0) + s(x1), data=dat, method="REML")
#' taps_wald_test(fit)
#' }
#'
#' @export
taps_wald_test <- function(fit,test.component=1,null.tol=1e-10) {
if (!inherits(fit, "gam")) stop("fit must be a 'gam' or 'bam' object.")

# Extract coefficients and covariance matrix from the fitted model
beta <- fit$coefficients
Vb   <- vcov.gam(fit)

# Generate the global design matrix (lpmatrix)
X <- predict(fit, fit$model, type = "lpmatrix")

# Residual degrees of freedom (global)
rdf <- fit$df.residual

# Smooth terms list
smooth_terms <- fit$smooth

# Prepare a list to store results for each smooth term
out_list <- vector("list", length(smooth_terms))

s <- smooth_terms[[test.component]]
term_name <- s$label
indices   <- s$first.para:s$last.para
reported_null_dim  <- s$null.space.dim  # The reported null space dimension

if(reported_null_dim==0) stop("No null space detected. Using summary() for p-value of this component.")
# Dynamically determine null space from penalty matrix S
S_matrix <- s$S[[1]]  # Extract first penalty matrix for the term

# Take intersection of zero rows and columns, then map back to indices
if(is.null(s$getA)==0){
detected_null_indices <- indices[1:reported_null_dim]
detect_method="reported"
}else{
col_norms <- apply(S_matrix, 2, function(x) sqrt(sum(x^2)))
detected_null_indices=which(col_norms<null.tol)
detect_method="detected"
}

# Ensure consistency with reported null space
if (length(detected_null_indices) != reported_null_dim) {
warning(sprintf("Mismatch in null space detection for term '%s': S reports %d null dimensions, detected %d from S.",
              term_name, reported_null_dim, length(detected_null_indices)))
reported_null_dim <- length(detected_null_indices)  # Override with detected value
}

# Smooth indices: remove detected null indices from full set
smooth_indices <- setdiff(indices, detected_null_indices)

# 1) Wald test for the unpenalized (fixed) component
fix_stat <- fix_df <- fix_p <- NA
if (length(detected_null_indices) > 0) {
est_null <- beta[detected_null_indices]
V_null   <- Vb[detected_null_indices, detected_null_indices, drop = FALSE]
tmp_fixed <- mgcv_wald(est_null, V_null, df_fixed = TRUE)

fix_stat <- tmp_fixed["statistic"]
fix_df   <- tmp_fixed["df"]
fix_p    <- tmp_fixed["p.value"]
} else {
fix_stat <- NA
fix_df   <- 0
fix_p    <- NA
}

# 2) Wald-like (or pseudo-likelihood ratio) test for the penalized (smooth) component
sm_stat <- sm_df <- sm_p <- NA
if (length(smooth_indices) > 0) {
p_smooth <- beta[smooth_indices]
V_smooth <- Vb[smooth_indices, smooth_indices, drop = FALSE]

# Extract the relevant columns from X
Xt <- X[, smooth_indices, drop = FALSE]

# Effective degrees of freedom for this smooth
edf_for_this_smooth <- sum(fit$edf[indices])
if (!is.null(fit$edf1)) {
edf_for_this_smooth <- sum(fit$edf1[indices])
}
edf_penalized <- edf_for_this_smooth - reported_null_dim
the_rank <- min(ncol(Xt), edf_penalized)

tmp_smooth <- testStat_wrapper(p_smooth, Xt, V_smooth, rank = the_rank, rdf = rdf)
sm_stat <- tmp_smooth["statistic"]
sm_df   <- tmp_smooth["df"]
sm_p    <- tmp_smooth["p.value"]
}

# Combine the results for this smooth term
out_list <- data.table(
mixed.term           = term_name,
`fix.df`       = fix_df,
`fix.chisq`    = fix_stat,
`fix.pvalue`   = fix_p,
`fix.indices`  = detect_method,
`smooth.df`    = sm_df,
`smooth.chisq` = sm_stat,
`smooth.pvalue`= sm_p
)

return(out_list)
}
