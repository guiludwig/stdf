#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <cmath>        // std::exp(double)
#include <valarray>     // std::valarray, std::exp(valarray)
#include <iostream>

using namespace Eigen;
using namespace Rcpp;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double logProfileCpp(const VectorXd theta, const MatrixXd DTR, const VectorXd Y, 
                     const MatrixXd XTR, const MatrixXd PhiTime, const VectorXd LambEst) {
  /* 
   theta: J x 1
   DTR: N x N
   Y: N x 1
   XTR: N x 6 (b0, bx, by, spline-basis)
   PhiTime: N x J
   LambEst: (J-1) x 1
  
   removed arg: const MatrixXd&NoiTR,
     pass NoiTR[,1] to it
   removed arg: const VectorXd& PhiEst, 
     pass Phi.est[NoiTR[,2]] to it
     
   More about eigen: http://home.uchicago.edu/~skrainka/pdfs/Talk.Eigen.pdf
   
     */
  int N = Y.size();
  int J = LambEst.size();
  MatrixXd psi; // Dynamic size means: not known at compilation time.
  psi.setZero(N,N);
  for(int j = 0; j < J; j++){ 
    psi += LambEst(j)*((-1*DTR/theta(j)).array().exp().matrix()).cwiseProduct(PhiTime.col(j)*PhiTime.col(j).transpose());
  }
  psi += theta(J)*MatrixXd::Identity(N,N); // theta has J+1 elements
  MatrixXd U = psi.llt().matrixL().transpose(); // same as chol(psi) in R
  // This step finds beta by Generalized Least Squares
  MatrixXd SX = psi.llt().solve(XTR); // A\b by Cholesky's decomposition
  VectorXd beta = ((XTR.adjoint())*SX).ldlt().solve((SX.adjoint())*Y);
  // End GLS
  VectorXd resid = Y - XTR*beta;
  double quadForm = (resid.transpose())*(psi.ldlt().solve(resid));
  double logUdet = U.diagonal().array().sum().log();
  return(quadForm/2 + logUdet + N/2);
}

// Test:
// logProfileCpp(testTheta <- 1:3, testDTR <- matrix(1,3,3), testY <- 1:3, testXTR <- matrix(0,3,3), testPhiTime <- matrix(rep(1:3,2), ncol=2), testLambEst <- 1:3)
// system.time({logProfileCpp(theta0, DTR, NoiTR[,1], XTR, Phi.est[NoiTR[,2],], lamb.est)})
// Compare with
// system.time({logProfile(theta0, DTR, NoiTR, XTR, Phi.est, lamb.est)})