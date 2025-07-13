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
arma::vec normal_uni_log(arma::vec data,
                         const double mu,
                         const double sigma) {
  const int N = data.n_rows;
  
  arma::vec res(N, arma::fill::zeros);
  for (int i=0; i<N; ++i) {
    res(i) = R::dnorm(data(i), mu, sigma, TRUE);
  }
  
  return res;
}


//[[Rcpp::export]]
arma::vec trunc_normal_uni_log(arma::vec data,
                               const double mu,
                               const double sigma,
                               const double Z) {
  const int N = data.n_rows;
  
  arma::vec res(N, arma::fill::zeros);
  for (int i=0; i<N; ++i) {
    res(i) = R::dnorm(data(i), mu, sigma, TRUE) - std::log(Z);
  }
  
  return res;
}


//[[Rcpp::export]]
arma::vec unif_log(arma::vec data,
                   const double a,
                   const double b) {
  const int N = data.n_rows;
  
  arma::vec res(N, arma::fill::zeros);
  for (int i=0; i<N; ++i) {
    res(i) = R::dunif(data(i), a, b, TRUE);
  }
  
  return res;
}


//[[Rcpp::export]]
arma::vec pois_log(const int data,
                   arma::vec rate) {
  const int N = rate.n_rows;
  
  arma::vec res(N, arma::fill::zeros);
  for (int i=0; i<N; ++i) {
    res(i) = R::dpois(data, rate(i), TRUE);
  }
  
  return res;
}







