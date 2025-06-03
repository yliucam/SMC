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
double normal_log_sum(arma::vec data,
                      const double mu,
                      const double sigma) {
  const int n = data.n_rows;
  
  double res = 0.0;
  
  for (int i=0; i<n; ++i) {
    res += R::dnorm(data(i), mu, sigma, TRUE);
  }
  
  return res;
}






