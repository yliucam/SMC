// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;
using namespace std;
using namespace arma;
using namespace R;

//[[Rcpp::export]]
arma::vec exponential_psd_log_sum(arma::vec data,
                                  arma::vec lambda) {
  const int N = lambda.n_rows;
  const int n = data.n_rows;
  
  arma::vec res(N, arma::fill::zeros);
  
  for (int i=0; i<N; ++i) {
    for (int j=0; j<n; ++j) {
      res(i) = res(i) + R::dexp(data(j), 1/lambda(i), TRUE);
    }
  }
  
  return res;
}


//[[Rcpp::export]]
arma::vec normal_nodes_log_sum(arma::mat data,
                               arma::vec mu,
                               arma::vec sigma) {
  const int N = data.n_rows;
  
  arma::vec res(N, arma::fill::zeros);
  
  for (int i=0; i<N; ++i) {
    res(i) = res(i) + R::dnorm(data(i,0), mu(i), sigma(0), TRUE) + R::dnorm(data(i,1), mu(i), sigma(1), TRUE);
  }
  
  return res;
}


//[[Rcpp::export]]
arma::vec normal_nodes_log_sum_new(arma::mat data,
                               arma::mat mu,
                               arma::mat sigma) {
  const int N = data.n_rows;
  const int d = data.n_cols;
  
  arma::vec res(N, arma::fill::zeros);
  
  for (int i=0; i<N; ++i) {
    for (int j=0; j<d; ++j) {
      res(i) = res(i) + R::dnorm(data(i,j), mu(i,j), sigma(i,j), TRUE);
    }
  }
  
  return res;
}




