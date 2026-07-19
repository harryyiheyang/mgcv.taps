// Cox score and observed information using mgcv's Peto approximation for ties.
// The calculation stays in coefficient space; no n by n matrix is formed.

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List cox_peto_suffstat(const arma::mat& X,
                             const arma::vec& eta,
                             const arma::vec& time,
                             const arma::ivec& status) {
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;

  if (eta.n_elem != n || time.n_elem != n || status.n_elem != n) {
    Rcpp::stop("Cox design, eta, time, and status must have matching rows.");
  }
  if (n == 0) Rcpp::stop("Cox sufficient statistics require observations.");

  for (arma::uword i = 0; i < n; ++i) {
    if (status(i) != 0 && status(i) != 1) {
      Rcpp::stop("Cox status must contain only zero and one.");
    }
  }

  const arma::uvec ord = arma::sort_index(time, "descend");
  const arma::vec w = arma::exp(eta - eta.max());

  double S0 = 0.0;
  arma::rowvec S1(p, arma::fill::zeros);
  arma::mat S2(p, p, arma::fill::zeros);
  arma::vec score(p, arma::fill::zeros);
  arma::mat information(p, p, arma::fill::zeros);
  int events = 0;

  arma::uword i = 0;
  while (i < n) {
    arma::uword j = i + 1;
    const double tied_time = time(ord(i));
    while (j < n && time(ord(j)) == tied_time) ++j;

    for (arma::uword r = i; r < j; ++r) {
      const arma::uword idx = ord(r);
      const double wr = w(idx);
      const arma::rowvec xr = X.row(idx);
      S0 += wr;
      S1 += wr * xr;
      S2 += wr * xr.t() * xr;
    }

    int d = 0;
    arma::rowvec event_sum(p, arma::fill::zeros);
    for (arma::uword r = i; r < j; ++r) {
      const arma::uword idx = ord(r);
      if (status(idx) == 1) {
        ++d;
        event_sum += X.row(idx);
      }
    }

    if (d > 0) {
      const arma::rowvec xbar = S1 / S0;
      score += event_sum.t() - static_cast<double>(d) * xbar.t();
      information += static_cast<double>(d) *
        (S2 / S0 - xbar.t() * xbar);
      events += d;
    }

    i = j;
  }

  information = 0.5 * (information + information.t());
  return Rcpp::List::create(
    Rcpp::Named("score") = score,
    Rcpp::Named("information") = information,
    Rcpp::Named("events") = events
  );
}
