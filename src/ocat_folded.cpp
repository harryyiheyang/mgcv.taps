// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

inline double clip(double x, double lo, double hi){
  return x < lo ? lo : (x > hi ? hi : x);
}

inline double logistic(double x){
  return 1.0 / (1.0 + std::exp(-x));
}

// [[Rcpp::export]]
List ocat_folded(
  NumericVector eta,
  IntegerVector y_int,      // 1-indexed categories, values in 1..J
  NumericVector alpha,      // J-1 cut points (not padded), from fit$family$getTheta(TRUE)
  double eps_mu       = 1e-12,
  double eps_w        = 1e-12,
  double clip_eta_min = -20.0,
  double clip_eta_max =  20.0,
  int n_threads = 1
){
  const int n   = eta.size();
  const int Jm1 = alpha.size();   // J - 1
  const int J   = Jm1 + 1;

  if(y_int.size() != n)
    stop("eta and y_int must have same length");

  NumericVector w_star(n, 0.0);
  NumericVector z_star(n, NA_REAL);

  #ifdef _OPENMP
  if(n_threads < 1) n_threads = 1;
  #pragma omp parallel for num_threads(n_threads) schedule(static)
  #endif
  for(int i = 0; i < n; ++i){
    double etai = clip(eta[i], clip_eta_min, clip_eta_max);
    int ci = y_int[i];   // 1-indexed observed category

    // v[j] = P(y <= j), v[0]=0, v[J]=1
    std::vector<double> v(J + 1);
    v[0] = 0.0;
    for(int j = 1; j <= Jm1; ++j)
      v[j] = logistic(alpha[j-1] - etai);
    v[J] = 1.0;

    // mu[j] = P(y == j) = v[j] - v[j-1], j = 1..J
    std::vector<double> mu(J + 1);
    for(int j = 1; j <= J; ++j)
      mu[j] = v[j] - v[j-1];

    double sum_Rmu  = 0.0;
    double sum_R2mu = 0.0;
    double score    = 0.0;
    for(int j = 1; j <= Jm1; ++j){
      double Rij = (v[j] + v[j-1] - 1.0) - v[Jm1];
      double yij = (ci == j) ? 1.0 : 0.0;
      sum_Rmu  += Rij * mu[j];
      sum_R2mu += Rij * Rij * mu[j];
      score    += Rij * (yij - mu[j]);
    }
    double wi = sum_R2mu - sum_Rmu * sum_Rmu;

    if(wi < eps_w){
      w_star[i] = 0.0;
      z_star[i] = NA_REAL;
      continue;
    }
    w_star[i] = wi;
    z_star[i] = etai + score / wi;
  }

  return List::create(
    _["w_star"] = w_star,
    _["z_star"] = z_star
  );
}
