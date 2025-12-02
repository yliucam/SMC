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
arma::vec pois_log(const arma::vec& data,
                   const arma::vec& rate) {
  const int N = rate.n_rows;
  
  arma::vec res(N, arma::fill::zeros);
  for (int i=0; i<N; ++i) {
    res(i) = R::dpois(data(i), rate(i), TRUE);
  }
  
  return res;
}


//[[Rcpp::export]]
arma::vec binom_log(const arma::vec& data,
                    const arma::vec& size,
                    const arma::vec& prob) {
  const int N = size.n_rows;
  
  arma::vec res(N, arma::fill::zeros);
  for (int i=0; i<N; ++i) {
    res(i) = R::dbinom(data(i), size(i), prob(i), TRUE);
  }
  
  return res;
}


//[[Rcpp::export]]
arma::vec latent_transition_density_log(arma::mat xJ,
                                        arma::mat sur,
                                        arma::mat imm,
                                        arma::vec alpha6,
                                        arma::vec deltaJF,
                                        arma::vec deltaAF,
                                        arma::vec rho) {
  const int N = xJ.n_rows;
  const int tt = xJ.n_cols;
  
  arma::vec res(N, arma::fill::zeros);
  
  // Sum of log priors for alpha6, xJ[1], sur[1] and imm[1]
  for (int i=0; i<N; ++i) {
    res(i) += R::dnorm(alpha6(i), -2, 2, TRUE) - std::log(1); // N(0,2^2)[-10,10] prior for alpha6. Note that pnorm(8,-2,2) - pnorm(-12,-2,2) is extremely close to 1, so we subtract log(1).
    res(i) += 3/50; // DU(0,50) priors for xJ[1], sur[1] and imm[1].
  }
  
  arma::vec rateJ(N, arma::fill::zeros);
  arma::vec rate_imm(N, arma::fill::zeros);
  arma::mat x = xJ + sur + imm;
  
  for (int i=1; i<tt; ++i) {
    for (int j=0; j<N; ++j) {
      rateJ(j) = 0.5 * x(j,i) * rho(j) * deltaJF(j);
      rate_imm(j) = x(j,i) * std::exp(alpha6(j));
      res(j) += R::dpois(xJ(j,i), rateJ(j), TRUE) +
        R::dbinom(sur(j,i), x(j,i), deltaAF(j), TRUE) +
        R::dpois(imm(j,i), rate_imm(j), TRUE);
    }
  }
  
  return res;
}


//[[Rcpp::export]]
arma::vec p_common_likelihood(arma::vec data,
                              arma::mat xJ,
                              arma::mat sur,
                              arma::mat imm) {
  const int N = xJ.n_rows;
  const int tt = xJ.n_cols;
  
  arma::vec res(N, arma::fill::zeros);
  arma::mat x = xJ + sur + imm;
  
  for (int i=0; i<tt; ++i) {
    for (int j=0; j<N; ++j) {
      res(j) += R::dpois(data(i), x(j,i), TRUE);
    }
  }
  
  return res;
}


// [[Rcpp::export]]
arma::vec fecundity_llike(arma::mat data,
                          arma::vec rho) {
  const int n = data.n_rows;
  const int N = rho.n_rows;
  
  arma::vec res(N, arma::fill::zeros);
  
  for (int i=0; i<N; ++i) {
    for (int j=0; j<n; ++j) {
      res(i) += R::dpois(data(j,1), rho(i)*data(j,0), TRUE);
    }
  }
  
  return res;
}


// [[Rcpp::export]]
arma::vec merging_likelihood_cpp(const double alpha6,
                                 arma::vec xJ,
                                 arma::vec sur,
                                 arma::vec imm,
                                 arma::vec alpha0,
                                 arma::vec alpha2,
                                 arma::vec rho,
                                 const int N,
                                 const int tt) {

  arma::vec res(N, arma::fill::zeros);
  arma::vec x = xJ + sur + imm;
  
  
  for (int i=0; i<N; ++i) {
    double deltaJF = std::exp(alpha0(i)) / (1 + std::exp(alpha0(i)));
    double deltaAF = std::exp(alpha0(i) + alpha2(i)) / (1 + std::exp(alpha0(i) + alpha2(i)));
    for (int j=0; j<tt; ++j) {
      if (j == 0) {
        res(i) += std::log(1.0 / 50.0) * 3;
      } else {
        double rate_J = .5 * x(j-1) * rho(i) * deltaJF;
        res(i) += R::dpois(xJ(j), rate_J, TRUE) + 
          R::dbinom(sur(j), x(j-1), deltaAF, TRUE);
      }
    }
  }
  
  return res;
}
  
  
  
  
  
  
  
  
  
  





