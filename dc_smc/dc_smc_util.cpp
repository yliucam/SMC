// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;
using namespace std;
using namespace arma;
using namespace R;

//[[Rcpp::export]]
arma::vec weight_check(arma::vec data,
                    arma::vec mu,
                    const double min_lim) {
  const int N = mu.n_rows;
  const int p = data.n_rows;
  
  arma::vec res(N);
  
  for (int i=0; i<N; ++i) {
    double w = 0;
    for (int j=0; j<p; ++j) {
      w = w + std::log(1 / std::sqrt(2*M_PI)) - std::pow((data[j] - mu[i]), 2.0) / 2;
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



















