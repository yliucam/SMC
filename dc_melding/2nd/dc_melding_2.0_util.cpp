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
                         const arma::vec mu,
                         const arma::vec sigma) {
  const int N = data.n_rows;
  
  arma::vec res(N, arma::fill::zeros);
  for (int i=0; i<N; ++i) {
    res(i) = R::dnorm(data(i), mu(i), sigma(i), TRUE);
  }
  
  return res;
}


//[[Rcpp::export]]
arma::vec normal_mult_log(arma::mat data,
                         const arma::vec mu,
                         const arma::mat Sigma) {
  const int N = data.n_rows;
  const int d = data.n_cols;
  
  arma::vec res(N, arma::fill::zeros);
  
  double log_determ = -std::log(2 * M_PI) - 1/2.0 * arma::log_det(Sigma).real();
  for (int i=0; i<N; ++i) {
    arma::vec z = data.row(i).st() - mu;
    arma::mat zSz = arma::trans(z) * arma::inv(Sigma) * z;
    res(i) = log_determ - 1/2.0 * zSz(0,0);
  }
  
  return res;
}



//[[Rcpp::export]]
arma::vec normal_uni_log_sum(arma::vec data,
                             arma::vec mu,
                             arma::vec sigma) {
  const int n = data.n_rows;
  const int N = mu.n_rows;
  
  arma::vec res(N, arma::fill::zeros);
  
  for (int i=0; i<N; ++i) {
    for (int j=0; j<n; ++j) {
      res(i) += R::dnorm(data(j), mu(i), sigma(i), TRUE);
    }
  }
  
  return (res);
}









