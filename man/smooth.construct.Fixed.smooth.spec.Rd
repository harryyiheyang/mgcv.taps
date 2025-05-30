% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smooth.construct.Fixed.smooth.spec.R
\name{smooth.construct.Fixed.smooth.spec}
\alias{smooth.construct.Fixed.smooth.spec}
\title{Construct a Fixed-Effect-Only Smooth Term}
\usage{
\method{smooth.construct}{Fixed.smooth.spec}(object, data, knots)
}
\arguments{
\item{object}{A smooth specification object created by `s()`, containing user-defined smoothing parameters.}

\item{data}{A data frame containing the covariate for the smooth term.}

\item{knots}{A list of knots (not used in this fixed-effect-only implementation).}
}
\value{
A smooth term object of class `"Fixed.smooth"`, `"mgcv.smooth"`, containing:
\itemize{
  \item `X`: The design matrix consisting solely of the fixed effect components from `A`.
  \item `S`: A zero penalty matrix with dimensions matching the number of columns in `A`.
  \item `rank`: The rank of the penalty matrix (set to 0 as it is a zero matrix).
  \item `null.space.dim`: The dimension of the fixed effect (null space), equal to the number of columns in `A`.
  \item `df`: The degrees of freedom, equal to the number of columns in `A`.
}
}
\description{
This function constructs a smooth term based solely on a fixed effect structure defined by `getA`.
The fixed effect structure is user-defined via `xt$getA`, with corresponding parameters in `xt$para`.
The function ensures that only fixed effects are modeled by checking that `fx=TRUE`. If `fx=FALSE`, it
stops with an error message. The penalty matrix `S` is set to a zero matrix with dimensions matching
the number of columns in `A`.
}
