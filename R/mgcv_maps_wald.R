#' Perform Wald Tests on mgcv Smooth Terms (Null and Penalized Components)
#'
#' This function takes a fitted mgcv `gam` or `bam` object,
#' then performs separate Wald tests on the null (unpenalized/fixed) and penalized
#' (smooth) components of each smooth term.
#' It returns a clear summary as a `data.table`.
#'
#' @param fit A fitted `gam` or `bam` model from the `mgcv` package.
#'
#' @return A `data.table` summarizing Wald statistics and p-values separately
#'         for null (unpenalized) and penalized (smooth) components of each term.
#'
#' @importFrom mgcv vcov.gam
#' @importFrom data.table data.table rbindlist
#' @importFrom ordgam testStat
#' @importFrom stats pchisq
#'
#' @examples
#' \dontrun{
#' fit <- gam(y ~ s(x0) + s(x1), data=dat, method=\"REML\")
#' mgcv_maps_wald(fit)
#' }
#'
#' @export
mgcv_maps_wald <- function(fit) {
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

for (i in seq_along(smooth_terms)) {
s <- smooth_terms[[i]]
term_name <- s$label
indices   <- s$first.para:s$last.para
reported_null_dim  <- s$null.space.dim  # The reported null space dimension

# Dynamically determine null space from penalty matrix S
S_matrix <- s$S[[1]]  # Extract first penalty matrix for the term

# Identify rows/columns where S is completely zero (null space)
zero_rows <- which(rowSums(abs(S_matrix)) == 0)
zero_cols <- which(colSums(abs(S_matrix)) == 0)

# Take intersection of zero rows and columns, then map back to indices
detected_null_indices <- indices[intersect(zero_rows, zero_cols)]

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
out_list[[i]] <- data.table(
term           = term_name,
`fix.df`       = fix_df,
`fix.chisq`    = fix_stat,
`fix.pvalue`   = fix_p,
`smooth.df`    = sm_df,
`smooth.chisq` = sm_stat,
`smooth.pvalue`= sm_p
)
}

# Bind the list into a single data.table
out_dt <- rbindlist(out_list)
return(out_dt)
}
