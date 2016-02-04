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

//' @export
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double logProfileCpp(const Eigen::VectorXd theta, const Eigen::MatrixXd DTR, 
                     const Eigen::VectorXd Y, const Eigen::MatrixXd XTR, 
                     const Eigen::VectorXd subsetStatic, const Eigen::MatrixXd PhiTime, 
                     const Eigen::VectorXd LambEst) {
  /* 
   theta: (J + 2) x 1; first J elements have temporal components,
          J+1 has static variance, J+2 has roving variance.
   DTR: N x N
   Y: N x 1
   XTR: N x (3 + spline.df) (b0, bx, by, spline-basis)
   PhiTime: N x J
   LambEst: J x 1
  
   More about eigen: http://home.uchicago.edu/~skrainka/pdfs/Talk.Eigen.pdf
                     http://eigen.tuxfamily.org/dox/AsciiQuickReference.txt
   
     */
  
  int N = Y.size();
  int J = LambEst.size();
  MatrixXd psi(MatrixXd(N,N).setZero()); // Dynamic size means: not known at compilation time.
  for(int j = 0; j < J; j++){ 
    // psi += LambEst(j)*((-1*DTR/theta(j)).array().exp().matrix()).cwiseProduct(PhiTime.col(j)*PhiTime.col(j).adjoint()); // PhiPhit.selfadjointView<Lower>().rankUpdate(PhiTime.col(j))
    psi += LambEst(j)*covExp(DTR, theta(j)).cwiseProduct(PhiTime.col(j)*PhiTime.col(j).adjoint()); // PhiPhit.selfadjointView<Lower>().rankUpdate(PhiTime.col(j))
  }
  // This is weird but: sigma_R I + (sigma_S - sigma_R) 1{static}
  VectorXd RandNoise = theta(J+1)*Eigen::VectorXd::Constant(N,1) + (theta(J) - theta(J+1))*subsetStatic;
  psi += RandNoise.asDiagonal(); // theta has J+2 elements
  Eigen::MatrixXd U = psi.llt().matrixL().adjoint(); // same as chol(psi) in R
  // This step finds beta by Generalized Least Squares
  Eigen::MatrixXd SX = psi.llt().solve(XTR); // A\b by Cholesky's decomposition
  Eigen::VectorXd beta = ((XTR.adjoint())*SX).ldlt().solve((SX.adjoint())*Y);
  // End GLS
  Eigen::VectorXd resid = Y - XTR*beta;
  double quadForm = (resid.adjoint())*(psi.ldlt().solve(resid));
  double logUdet = 2*U.diagonal().array().log().sum(); // = log(|Sigma|)
  return(quadForm + logUdet);
}

// Test:
// system.time({logProfileCpp(testTheta <- 1:3, testDTR <- matrix(1,420,420), testY <- 1:420, testXTR <- matrix(1,420,6), testSubsetStatic <- rep(c(0,1), each = 210), testPhiTime <- matrix(rep(1:420,2), ncol=2), testLambEst <- 1:2)})
// replicate(10, system.time({logProfileCpp(testTheta <- 1:3, testDTR <- matrix(1,420,420), testY <- 1:420, testXTR <- matrix(1,420,6), testSubsetStatic <- rep(c(0,1), each = 210), testPhiTime <- matrix(rep(1:420,2), ncol=2), testLambEst <- 1:2)})[3])
// system.time({logProfileCpp(theta0, DTR, NoiTR[,1], XTR, Phi.est[NoiTR[,2],], lamb.est)})
// Compare with
// system.time({logProfile(theta0, DTR, NoiTR, XTR, Phi.est, lamb.est)})