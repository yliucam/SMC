// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
//#include <math.h>
#include <cmath>
  
using namespace Rcpp;
using namespace std;
using namespace arma;
using namespace R;

//[[Rcpp::export]]
arma::vec dnorm_log_uni(const arma::vec& data,
                        const arma::vec& mu,
                        const arma::vec& sigma) {
  const int n = data.n_rows;
  
  arma::vec res(n, arma::fill::zeros);
  for (int i=0; i<n; ++i) {
    res(i) = R::dnorm(data(i), mu(i), sigma(i), TRUE);
  }
  
  return res;
}

//[[Rcpp::export]]
arma::vec dlnorm_log_uni(const arma::vec& data,
                         const arma::vec& mulog,
                         const arma::vec& sigmalog) {
  const int n = data.n_rows;
  
  arma::vec res(n, arma::fill::zeros);
  for (int i=0; i<n; ++i) {
    res(i) = R::dlnorm(data(i), mulog(i), sigmalog(i), TRUE);
  }
  
  return res;
}

//[[Rcpp::export]]
arma::vec dgamma_log(const arma::vec& data,
                     const arma::vec& alpha,
                     const arma::vec& beta) {
  const int n = data.n_rows;
  
  arma::vec res(n, arma::fill::zeros);
  for (int i=0; i<n; ++i) {
    res(i) = R::dgamma(data(i), alpha(i), beta(i), TRUE);
  }
  
  return res;
}


//[[Rcpp::export]]
arma::vec dnorm_log_help(const arma::vec& data,
                         const arma::vec& mu,
                         const arma::vec& sigma) {
  const int n = data.n_rows;
  const int N = mu.n_rows;
  
  arma::vec res(N, arma::fill::zeros);
  for (int i=0; i<N; ++i) {
    for (int j=0; j<n; ++j) {
      res(i) += R::dnorm(data(j), mu(i), sigma(i), TRUE);
    }
  }
  
  return res;
}

