// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;
using namespace std;
using namespace arma;
using namespace R;

//[[Rcpp::export]]
arma::vec weight_check(arma::vec w_log,
                       const double min_lim) {
  const int n = w_log.n_rows;
  bool min_check;
  
  for (int i=0; i<n; ++i) {
    min_check = (w_log(i) < min_lim);
    if (min_check) {
      w_log(i) = min_lim;
    }
  }
  
  return w_log;
}


//[[Rcpp::export]]
arma::vec weight_update_binom(arma::vec data,
                              arma::vec p,
                              int n_trial,
                              const double min_lim) {
  const int N = p.n_rows;
  const int n = data.n_rows;
  
  arma::vec res(N);
  
  for (int i=0; i<N; ++i) {
    double w = 0;
    for (int j=0; j<n; ++j) {
      w = w + R::dbinom(data[j], n_trial, p[i], TRUE);
    }
    
    bool min_check;
    min_check = (w < min_lim);
    if (min_check) {
      w = min_lim;
    }
    
    res[i] = w;
  }
  
  return res;
}



//[[Rcpp::export]]
arma::vec likelihood_particle_beta(arma::mat obs,
                                   arma::mat alpha,
                                   arma::mat beta) {
  const int N = obs.n_rows;
  const int p = obs.n_cols;
  
  arma::vec res(N, arma::fill::zeros);
  
  for (int i=0; i<N; ++i) {
    for (int j=0; j<p; ++j) {
      res(i) = res(i) + R::dbeta(obs(i,j), alpha(i,j), beta(i,j), TRUE);
    }
  }
  
  return res;
}



//[[Rcpp::export]]
arma::vec likelihood_particle_gamma(arma::mat obs,
                                    arma::mat alpha,
                                    arma::mat beta) {
  const int N = obs.n_rows;
  const int p = obs.n_cols;
  
  arma::vec res(N, arma::fill::zeros);
  
  for (int i=0; i<N; ++i) {
    for (int j=0; j<p; ++j) {
      res(i) = res(i) + R::dgamma(obs(i,j), alpha(i,j), beta(i,j), TRUE);
    }
  }
  
  return res;
}


