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
#' @importFrom stats pchisq
#'
#' @examples
#' \dontrun{
#' fit <- gam(y ~ s(x0) + s(x1), data=dat, method=\"REML\")
#' mgcv_smooth_wald(fit)
#' }
#'
#' @export
mgcv_smooth_wald <- function(fit) {
if (!inherits(fit, "gam")) stop("fit must be a 'gam' or 'bam' object.")

beta <- fit$coefficient
Vb <- vcov.gam(fit)
smooth_terms <- fit$smooth

results_list <- vector("list", length(smooth_terms))

for (i in seq_along(smooth_terms)) {
s <- smooth_terms[[i]]
term_name <- s$label
indices <- s$first.para:s$last.para
null_dim <- s$null.space.dim

# Split indices
if (null_dim > 0) {
null_indices <- indices[1:null_dim]
smooth_indices <- indices[(null_dim+1):length(indices)]
} else {
null_indices <- integer(0)
smooth_indices <- indices
}

# Wald test for null space
null_results <- mgcv_wald(beta, Vb, null_indices)
smooth_results <- mgcv_wald(beta, Vb, smooth_indices)

results_list[[i]] <- data.table(
`Term` = term_name,
`Fix df` = null_results[2],
`Fix chisq` = null_results[1],
`Fix p` = null_results[3],
`Random df` = smooth_results[2],
`Random chisq` = smooth_results[1],
`Random p` = smooth_results[3]
)
}

results_df <- do.call(rbind, results_list)
rownames(results_df) <- NULL
return(results_df)
}


