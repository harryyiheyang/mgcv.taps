#' Perform Score Tests on mgcv Smooth Terms (Null and Penalized Components)
#'
#' This function takes a fitted mgcv `gam` or `bam` object,
#' then performs Score test on the null (unpenalized/fixed) and penalized
#' (smooth) components of each smooth term.
#' It returns a clear summary as a `data.table`.
#'
#' @param fit A fitted `gam` or `bam` model from the `mgcv` package.
#' @param test.component The index of smooth component to be tested.
#' @param null.tol The tolerance of row norm to detect indices of null space. Default to 1e-10.
#' @param method Method for p-value calculation: `"satterthwaite"` or `"liu"` or `"davies"` or `"imhof"` or `"hall"` or `"wood"` (default).
#' @param max_iter The maximum iteration used in CompQuadForm to compute the p-values. Default to 1e5.
#' @param max_eps The tolerance used in CompQuadForm to compute the p-values. Default to 1e-8.

#' @return A `data.table` summarizing Score statistics and p-values separately
#'         for penalized (smooth) components of each term.
#'
#' @importFrom mgcv vcov.gam
#' @importFrom data.table data.table rbindlist
#' @importFrom ordgam testStat
#' @importFrom Matrix bdiag
#' @import CppMatrix
#' @import CompQuadForm
#'
#' @examples
#' \dontrun{
#' fit <- gam(y ~ s(x0) + s(x1), data=dat, method="REML")
#' taps_score_test(fit)
#' }
#'
#' @export
taps_score_test <- function(fit,test.component=1,null.tol=1e-10,method="wood",max_eps=1e-8,max_iter=1e5) {
if (!inherits(fit, "gam")) stop("fit must be a 'gam' or 'bam' object.")
prior_weights <- fit$prior.weights
if (any(prior_weights != 1)) {
stop("Score test currently only supports models with prior.weights == 1 for all observations.")
}
################################################################################################
test_result=taps_score_test_operator(fit=fit,test.component=test.component,null.tol=null.tol,method=method,max_eps=max_eps,max_iter=max_iter)
return(test_result)
}
