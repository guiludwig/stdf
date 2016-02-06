#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <cmath>        // std::exp(double)
#include <valarray>     // std::valarray, std::exp(valarray)
#include <iostream>
#include <Rmath.h>
#include "covGauss.h"
#include "covMat.h"
#include "covExp.h"

using namespace Eigen;
using namespace Rcpp;

//' @export
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd dlogProfileCpp(const Eigen::VectorXd theta, const Eigen::MatrixXd DTR, 
const Eigen::VectorXd Y, const Eigen::MatrixXd XTR, 
const Eigen::VectorXd subsetStatic, const Eigen::MatrixXd PhiTime, 
const Eigen::VectorXd LambEst, const double nu) {
  /* 
  theta: J x 1
  DTR: N x N
  Y: N x 1
  XTR: N x 6 (b0, bx, by, spline-basis)
  PhiTime: N x J
  LambEst: (J-1) x 1
    
  More about eigen: http://home.uchicago.edu/~skrainka/pdfs/Talk.Eigen.pdf
  http://eigen.tuxfamily.org/dox/AsciiQuickReference.txt
  
  */
  
  int N = Y.size();
  int J = LambEst.size();
  MatrixXd psi(MatrixXd(N,N).setZero()); // Dynamic size means: not known at compilation time.
  if(nu == 0.5){
    for(int j = 0; j < J; j++){ 
      psi += LambEst(j)*covExp(DTR, theta(j)).cwiseProduct(PhiTime.col(j)*PhiTime.col(j).adjoint()); // PhiPhit.selfadjointView<Lower>().rankUpdate(PhiTime.col(j))
    }
  } else if(nu > 10){
    for(int j = 0; j < J; j++){ 
      psi += LambEst(j)*covGauss(DTR, theta(j)).cwiseProduct(PhiTime.col(j)*PhiTime.col(j).adjoint()); // PhiPhit.selfadjointView<Lower>().rankUpdate(PhiTime.col(j))
    }
  } else {
    for(int j = 0; j < J; j++){ 
      psi += LambEst(j)*covMat(DTR, theta(j), nu).cwiseProduct(PhiTime.col(j)*PhiTime.col(j).adjoint()); // PhiPhit.selfadjointView<Lower>().rankUpdate(PhiTime.col(j))
    }
  }
  VectorXd RandNoise = theta(J+1)*VectorXd::Constant(N,1) + (theta(J) - theta(J+1))*subsetStatic;
  psi += RandNoise.asDiagonal(); // theta has J+2 elements
  MatrixXd psiInv = psi.llt().solve(MatrixXd::Identity(N,N));
  // Eigen::MatrixXd U = psi.llt().matrixL().adjoint(); // same as chol(psi) in R
  // This step finds beta by Generalized Least Squares
  // Eigen::MatrixXd SX = psi.llt().solve(XTR); // A\b by Cholesky's decomposition
  // Eigen::VectorXd beta = ((XTR.adjoint())*SX).ldlt().solve((SX.adjoint())*Y);
  VectorXd beta = ((XTR.adjoint())*psiInv*XTR).ldlt().solve((XTR.adjoint())*psiInv*Y);
  // End GLS
  VectorXd resid = Y - XTR*beta;
  
  // Difference to logProfile() starts here 
  
  // Eigen::VectorXd sigmaRes = psi.ldlt().solve(resid);
  VectorXd sigmaRes = psiInv * resid;
  MatrixXd psij(MatrixXd(N,N).setZero());
  VectorXd dTheta = VectorXd::Constant(J+2,0);
  // psij.setZero(N,N);
  // dTheta.resize(J+2);
  
  if(nu == 0.5){
    for(int j = 0; j < J; j++){
      psij = DTR.cwiseProduct((LambEst(j)*covExp(DTR, theta(j))).cwiseProduct(PhiTime.col(j)*PhiTime.col(j).transpose()))/(theta(j)*theta(j)); // if using c = 1/theta, then -1
      dTheta(j) = -1*(sigmaRes.adjoint())*psij*sigmaRes + (psiInv*psij).trace();
    }
  } else if(nu > 10){
    for(int j = 0; j < J; j++){
      psij = (DTR).cwiseProduct((LambEst(j)*covGauss(DTR, theta(j))).cwiseProduct(PhiTime.col(j)*PhiTime.col(j).transpose()))/(pow(theta(j),4)); // if using c = 1/theta, then -1
      dTheta(j) = -1*(sigmaRes.adjoint())*psij*sigmaRes + (psiInv*psij).trace();
    }
  } else {
    for(int j = 0; j < J; j++){
      // http://functions.wolfram.com/Bessel-TypeFunctions/BesselK/20/ShowAll.html
      MatrixXd sDTR = (sqrt(2*nu)*DTR/theta(j)).array().pow(nu).matrix();
      psij = (pow(2, nu-1)/tgamma(nu))*sDTR.cwiseProduct(DTR.cwiseProduct((LambEst(j)*covMat(DTR, theta(j), nu - 1)).cwiseProduct(PhiTime.col(j)*PhiTime.col(j).transpose()))/(theta(j)*theta(j))); 
      dTheta(j) = -1*(sigmaRes.adjoint())*psij*sigmaRes + (psiInv*psij).trace();
    }
  }
  // need eig equivalent of which
  int staticObs = (int) subsetStatic.sum();
  
  VectorXd diagInverse = psiInv.diagonal();
  double diagS = 0.0;
  double diagR = 0.0;
  
  VectorXd sigmaResStatic = VectorXd::Constant(staticObs,1);
  VectorXd sigmaResRoving = VectorXd::Constant(N - staticObs,1);
  int wS = 0;
  int wR = 0;
  for(int w = 0; w < N; w++){
    if(subsetStatic(w) == 1) {
      sigmaResStatic(wS) = sigmaRes(w);
      diagS += diagInverse(w);
      wS++;
    } else {
      sigmaResRoving(wR) = sigmaRes(w);
      diagR += diagInverse(w);
      wR++;
    }
  }
  
    dTheta(J) = -1*(sigmaResStatic.adjoint())*sigmaResStatic + diagS;
    dTheta(J+1) = -1*(sigmaResRoving.adjoint())*sigmaResRoving + diagR;
  
  return(dTheta);
}

// Test:
// system.time({dlogProfileCpp(testTheta <- 1:3, testDTR <- matrix(1,420,420), testY <- 1:420, testXTR <- matrix(1,420,6), testSubsetStatic <- rep(c(0,1), each = 210), testPhiTime <- matrix(rep(1:420,2), ncol=2), testLambEst <- 1:2)})
// replicate(10, system.time({dlogProfileCpp(testTheta <- 1:3, testDTR <- matrix(1,420,420), testY <- 1:420, testXTR <- matrix(1,420,6), testSubsetStatic <- rep(c(0,1), each = 210), testPhiTime <- matrix(rep(1:420,2), ncol=2), testLambEst <- 1:2)})[3])
// system.time({dlogProfileCpp(theta0, DTR, NoiTR[,1], XTR, Phi.est[NoiTR[,2],], lamb.est)})