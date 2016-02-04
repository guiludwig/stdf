#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <cmath>
#include <valarray>
#include <iostream>

using namespace Eigen;
using namespace Rcpp;

//' @export
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
MatrixXd covGauss(const MatrixXd DTR, double theta){
  return((-1*DTR.array().pow(2)/theta).exp().matrix());
}

/*** R
# tests
test <- covGauss(A <- matrix(c(0,sqrt(2),sqrt(2),0), ncol=2), .4)
all.equal(test, exp(-(A^2)/.4))
*/