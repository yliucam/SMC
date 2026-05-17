// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
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
arma::vec gamma_log(arma::vec data,
                    const double alpha,
                    const double beta) {
  const int N = data.n_rows;
  
  arma::vec res(N, arma::fill::zeros);
  for (int i=0; i<N; ++i) {
    res(i) = R::dgamma(data(i), alpha, 1/beta, TRUE);
  }
  
  return res;
}


//[[Rcpp::export]]
arma::vec normal_uni_log_particle(arma::vec data,
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
  
  return res;
}


//[[Rcpp::export]]
arma::vec tstudent_log_particle(arma::vec data,
                                arma::vec mu,
                                arma::vec tau,
                                arma::vec nu) {
  const int n = data.n_rows;
  const int N = mu.n_rows;
  
  arma::vec res(N, arma::fill::zeros);
  for (int i=0; i<N; ++i) {
    for (int j=0; j<n; ++j) {
      res(i) += std::lgamma((nu(i)+1.0)/2.0) - std::lgamma(nu(i)/2.0) +
                                    .5 * (std::log(tau(i)) - std::log(nu(i))) - (nu(i)+1) / 2 *
      std::log(1.0 + tau(i) / nu(i) * std::pow(data(j)-mu(i), 2));
    }
  }
  
  return res;
}

