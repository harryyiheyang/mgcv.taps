% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linearity_piecewise_discontinuity.R
\name{linearity_piecewise_discontinuity}
\alias{linearity_piecewise_discontinuity}
\title{Generate a basis matrix for piecewise linear and regression discontinuity effects}
\usage{
linearity_piecewise_discontinuity(x, para)
}
\arguments{
\item{x}{A numeric vector representing the covariate values.}

\item{para}{A list with two numeric vectors:
\itemize{
  \item `para[[1]]`: Numeric vector of changepoints where slope changes occur.
  \item `para[[2]]`: Numeric vector of discontinuity points where jumps occur.
}}
}
\value{
A numeric matrix where columns represent:
  - An intercept term (1)
  - The original covariate `x`
  - Additional columns for each changepoint capturing slope changes
  - Additional columns for each discontinuity capturing jumps and post-breakpoint slope changes
}
\description{
This function constructs a basis matrix that captures both piecewise linear changes at
changepoints and discontinuous jumps at breakpoints.
}
