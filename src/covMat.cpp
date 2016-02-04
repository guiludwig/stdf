#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <cmath>
#include <Rmath.h>
#include <valarray>
#include <iostream>
#include <covMat.h>

using namespace Eigen;
using namespace Rcpp;

//' @export
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
MatrixXd covMat(const MatrixXd DTR, double theta, double nu){
  int N = DTR.rows();
  double con = (pow(2.0, nu-1))/tgamma(nu);
  MatrixXd sDTR = sqrt(2*nu)*DTR/theta;
  MatrixXd BsDTR(MatrixXd(N,N).setZero());
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      if(i != j) { // else zero
        BsDTR(i,j) = R::bessel_k(sDTR(i,j), nu, 1.0); 
      }
    }
  }
  MatrixXd M(N,N);
  M = con*(sDTR.array().pow(nu).matrix()).cwiseProduct(BsDTR);
  for(int i=0; i<N; i++){
    M(i,i) = 1;
  }
  // BsDTR.setZero()
  return(M);
}

/*** R
# tests
nu <- 6.5
test <- covMat(A <- matrix(c(0,sqrt(2),sqrt(2),0), ncol=2), .4, nu)
testR <- 2^(nu-1)/gamma(nu)*((sqrt(2*nu)*A/.4)^nu)*besselK((sqrt(2*nu)*A/.4), nu); diag(testR) <- 1
all.equal(test, testR)
*/