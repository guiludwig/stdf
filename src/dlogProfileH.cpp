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
Eigen::VectorXd dlogProfileCppH(const Eigen::VectorXd theta, const Eigen::MatrixXd DTR, 
                                const Eigen::VectorXd Y, const Eigen::MatrixXd XTR, 
                                const Eigen::MatrixXd PhiTime, const Eigen::VectorXd LambEst, 
                                const double nu) {
  /* 
   theta: J x 1
   DTR: N x N
   Y: N x 1
   XTR: N x 6 (b0, bx, by, spline-basis)
   PhiTime: N x J
   LambEst: (J-1) x 1
     
   More about eigen: http://home.uchicago.edu/~skrainka/pdfs/Talk.Eigen.pdf
   
     */
  
  int N = Y.size();
  int J = LambEst.size();
  MatrixXd psi(MatrixXd(N,N).setZero()); // Dynamic size means: not known at compilation time.
  for(int j = 0; j < J; j++){ 
    psi += LambEst(j)*((-1*DTR/theta(j)).array().exp().matrix()).cwiseProduct(PhiTime.col(j)*PhiTime.col(j).transpose());
  }
  psi += theta(J)*MatrixXd::Identity(N,N); // theta has J+1 elements
  // This step finds beta by Generalized Least Squares
  Eigen::MatrixXd SX = psi.llt().solve(XTR); // A\b by Cholesky's decomposition
  Eigen::VectorXd beta = ((XTR.adjoint())*SX).ldlt().solve((SX.adjoint())*Y);
  // End GLS
  Eigen::VectorXd resid = Y - XTR*beta;
  
  // Difference to logProfile() starts here 
  
  Eigen::VectorXd sigmaRes = psi.ldlt().solve(resid);
  Eigen::MatrixXd psij(Eigen::MatrixXd(N,N).setZero());
  Eigen::VectorXd dTheta = Eigen::VectorXd::Constant(J+1,0);
  
  // Eigen::SelfAdjointEigenSolver<MatrixXd> eig;
  
  if(nu == 0.5){
    for(int j = 0; j < J; j++){
      psij = DTR.cwiseProduct((LambEst(j)*covExp(DTR, theta(j))).cwiseProduct(PhiTime.col(j)*PhiTime.col(j).transpose()))/(theta(j)*theta(j)); // if using c = 1/theta, then -1
      dTheta(j) = -1*(sigmaRes.adjoint())*psij*sigmaRes + (psi.ldlt().solve(psij)).trace(); 
    }
  } else if(nu > 10){
    for(int j = 0; j < J; j++){
      psij = (DTR).cwiseProduct((LambEst(j)*covGauss(DTR, theta(j))).cwiseProduct(PhiTime.col(j)*PhiTime.col(j).transpose()))/(pow(theta(j),4)); // if using c = 1/theta, then -1
      dTheta(j) = -1*(sigmaRes.adjoint())*psij*sigmaRes + (psi.ldlt().solve(psij)).trace(); 
    }
  } else {
    for(int j = 0; j < J; j++){
      // http://functions.wolfram.com/Bessel-TypeFunctions/BesselK/20/ShowAll.html
      MatrixXd sDTR = (sqrt(2*nu)*DTR/theta(j)).array().pow(nu).matrix();
      psij = (pow(2, nu-1)/tgamma(nu))*sDTR.cwiseProduct(DTR.cwiseProduct((LambEst(j)*covMat(DTR, theta(j), nu - 1)).cwiseProduct(PhiTime.col(j)*PhiTime.col(j).transpose()))/(theta(j)*theta(j))); 
      dTheta(j) = -1*(sigmaRes.adjoint())*psij*sigmaRes + (psi.ldlt().solve(psij)).trace(); 
    }
  }
//   for(int j = 0; j < J; j++){
//     psij = DTR.cwiseProduct((LambEst(j)*((-1*DTR/theta(j)).array().exp().matrix())).cwiseProduct(PhiTime.col(j)*PhiTime.col(j).transpose()))/(theta(j)*theta(j)); // if using c = 1/theta, then -1
//     dTheta(j) = -1*(sigmaRes.adjoint())*psij*sigmaRes + (psi.ldlt().solve(psij)).trace(); //.diagonal().array().sum();
//   }
  dTheta(J) = -1*(sigmaRes.adjoint())*sigmaRes + (psi.llt().solve(Eigen::MatrixXd::Identity(N, N))).trace(); // eig.compute(psi).eigenvalues().array().pow(-1).sum();
  
  return(dTheta);
}

// Test:
// system.time({dlogProfileCpp(theta0, DTR, NoiTR[,1], XTR, Phi.est[NoiTR[,2],], lamb.est)})
// Compare to
// system.time({dlogProfile(theta0, DTR, NoiTR, XTR, Phi.est, lamb.est)})