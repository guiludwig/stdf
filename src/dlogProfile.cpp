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
Eigen::VectorXd dlogProfileCpp(const Eigen::VectorXd theta, const Eigen::MatrixXd DTR, const Eigen::VectorXd Y, 
                        const Eigen::MatrixXd XTR, const Eigen::MatrixXd PhiTime, const Eigen::VectorXd LambEst) {
  /* 
   theta: J x 1
   DTR: N x N
   Y: N x 1
   XTR: N x 6 (b0, bx, by, spline-basis)
   PhiTime: N x J
   LambEst: (J-1) x 1
  
   removed arg: const Eigen::MatrixXd&NoiTR,
     pass NoiTR[,1] to it
   removed arg: const Eigen::VectorXd& PhiEst, 
     pass Phi.est[NoiTR[,2]] to it
     
     
   More about eigen: http://home.uchicago.edu/~skrainka/pdfs/Talk.Eigen.pdf
   
     */
  
  int N = Y.size();
  int J = LambEst.size();
  Eigen::MatrixXd psi; // Dynamic size means: not known at compilation time.
  psi.setZero(N,N);
  for(int j = 0; j < J; j++){ 
    psi += LambEst(j)*((-1*DTR/theta(j)).array().exp().matrix()).cwiseProduct(PhiTime.col(j)*PhiTime.col(j).transpose());
  }
  psi += theta(J)*Eigen::MatrixXd::Identity(N,N); // theta has J+1 elements
  Eigen::MatrixXd U = psi.llt().matrixL().transpose(); // same as chol(psi) in R
  // This step finds beta by Generalized Least Squares
  Eigen::MatrixXd SX = psi.llt().solve(XTR); // A\b by Cholesky's decomposition
  Eigen::VectorXd beta = ((XTR.adjoint())*SX).ldlt().solve((SX.adjoint())*Y);
  // End GLS
  Eigen::VectorXd resid = Y - XTR*beta;
  
  // Difference to logProfile() starts here 
  
  Eigen::VectorXd sigmaRes = psi.llt().solve(resid);
  Eigen::MatrixXd psij;
  Eigen::VectorXd dTheta;
  dTheta.resize(J+1);
  
  Eigen::SelfAdjointEigenSolver<MatrixXd> eig;
  
  for(int j = 0; j < J; j++){
    psij = DTR.cwiseProduct((LambEst(j)*((-1*DTR/theta(j)).array().exp().matrix())).cwiseProduct(PhiTime.col(j)*PhiTime.col(j).transpose()))/(theta(j)*theta(j)); // if using c = 1/theta, then -1
    dTheta(j) = -1*(sigmaRes.transpose())*psij*sigmaRes + (psi.llt().solve(psij)).diagonal().array().sum();
    dTheta(j) = (psi.llt().solve(psij)).diagonal().array().sum();
  }
  dTheta(J) = -1*(sigmaRes.transpose())*sigmaRes + eig.compute(psi).eigenvalues().array().pow(-1).sum();
  
  return(dTheta);
}

// Test:
// system.time({dlogProfileCpp(theta0, DTR, NoiTR[,1], XTR, Phi.est[NoiTR[,2],], lamb.est)})
// Compare to
// system.time({dlogProfile(theta0, DTR, NoiTR, XTR, Phi.est, lamb.est)})