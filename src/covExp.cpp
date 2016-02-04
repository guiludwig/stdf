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
MatrixXd covExp(const MatrixXd DTR, double theta){
  return((-1*DTR/theta).array().exp().matrix());
}

/*** R
# tests
test <- covExp(A <- matrix(c(0,sqrt(2),sqrt(2),0), ncol=2), .4)
all.equal(test, exp(-abs(A)/.4))
*/