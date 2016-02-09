#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <cmath>
#include <Rmath.h>
#include <valarray>
#include <iostream>
#include "covMat.h"

using namespace Eigen;
using namespace Rcpp;

//' @export
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd covMat(const Eigen::MatrixXd DTR, double theta, double nu){
  int N = DTR.rows();
  double con = (pow(2.0, 1-nu))/tgamma(nu);
  MatrixXd sDTR = sqrt(2*nu)*DTR/theta;
  MatrixXd BsDTR(MatrixXd(N,N).setZero());
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      BsDTR(i,j) = R::bessel_k(sDTR(i,j), nu, 1.0); 
      if(BsDTR(i,j) == R_PosInf){
        // tolerance, assuming
        // > .Machine$double.eps ^ 0.5
        //   [1] 1.490116e-08
        sDTR(i,j) = 0.000000001;
        BsDTR(i,j) = R::bessel_k(0.000000001, nu, 1.0); 
      }
    }
  }
  MatrixXd M(N,N);
  M = con*(sDTR.array().pow(nu).matrix()).cwiseProduct(BsDTR);
  return(M);
}

/*** R
# tests
nu <- 6.5
test <- covMat(A <- matrix(c(0,sqrt(2),sqrt(2),0), ncol=2), .4, nu)
testR <- 2^(1-nu)/gamma(nu)*((sqrt(2*nu)*A/.4)^nu)*besselK((sqrt(2*nu)*A/.4), nu); diag(testR) <- 1
all.equal(test, testR)
*/