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
#' @param method Method for p-value calculation: `"satterthwaite"` (default) or `"liu"` or `"davies"` or `"imhof"` or `"farebrother"`.
#' @return A `data.table` summarizing Score statistics and p-values separately
#'         for penalized (smooth) components of each term.
#'
#' @importFrom mgcv vcov.gam
#' @importFrom data.table data.table rbindlist
#' @importFrom ordgam testStat
#' @importFrom stats pchisq
#' @importFrom Matrix bdiag
#' @import CppMatrix
#'
#' @examples
#' \dontrun{
#' fit <- gam(y ~ s(x0) + s(x1), data=dat, method="REML")
#' taps_score_test(fit)
#' }
#'
#' @export
taps_score_test <- function(fit,test.component=1,null.tol=1e-10,method="satterthwaite") {
if (!inherits(fit, "gam")) stop("fit must be a 'gam' or 'bam' object.")
prior_weights <- fit$prior.weights
if (any(prior_weights != 1)) {
stop("Score test currently only supports models with prior.weights == 1 for all observations.")
}
################################################################################################
test_result=taps_score_test_operator(fit=fit,test.component=test.component,null.tol=null.tol,method=method)
return(test_result)
}
