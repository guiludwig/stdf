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
Eigen::VectorXd dlogProfileCpp(const Eigen::VectorXd theta, const Eigen::MatrixXd DTR, 
                               const Eigen::VectorXd Y, const Eigen::MatrixXd XTR, 
                               const Eigen::VectorXd subsetStatic, const Eigen::MatrixXd PhiTime, 
                               const Eigen::VectorXd LambEst) {
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
                     http://eigen.tuxfamily.org/dox/AsciiQuickReference.txt

     */
  
  int N = Y.size();
  int J = LambEst.size();
  Eigen::MatrixXd psi(MatrixXd(N,N).setZero()); // Dynamic size means: not known at compilation time.
  for(int j = 0; j < J; j++){ 
    psi += LambEst(j)*((-1*DTR/theta(j)).array().exp().matrix()).cwiseProduct(PhiTime.col(j)*PhiTime.col(j).adjoint()); // PhiPhit.selfadjointView<Lower>().rankUpdate(PhiTime.col(j))
  }
  VectorXd RandNoise = theta(J)*Eigen::VectorXd::Constant(N,1) + (theta(J+1) - theta(J))*subsetStatic;
  psi += RandNoise.asDiagonal(); // theta has J+2 elements
  Eigen::MatrixXd U = psi.llt().matrixL().adjoint(); // same as chol(psi) in R
  // This step finds beta by Generalized Least Squares
  Eigen::MatrixXd SX = psi.llt().solve(XTR); // A\b by Cholesky's decomposition
  Eigen::VectorXd beta = ((XTR.adjoint())*SX).ldlt().solve((SX.adjoint())*Y);
  // End GLS
  Eigen::VectorXd resid = Y - XTR*beta;
  
  // Difference to logProfile() starts here 
  
  Eigen::VectorXd sigmaRes = psi.ldlt().solve(resid);
  Eigen::MatrixXd psij = Eigen::MatrixXd::Identity(N, N);
  Eigen::VectorXd dTheta = Eigen::VectorXd::Constant(J+2,0);
  // psij.setZero(N,N);
  // dTheta.resize(J+2);
  
  for(int j = 0; j < J; j++){
    psij = DTR.cwiseProduct((LambEst(j)*((-1*DTR/theta(j)).array().exp().matrix())).cwiseProduct(PhiTime.col(j)*PhiTime.col(j).transpose()))/(theta(j)*theta(j)); // if using c = 1/theta, then -1
    dTheta(j) = -1*(sigmaRes.adjoint())*psij*sigmaRes + (psi.ldlt().solve(psij)).diagonal().array().sum();
  }
  // need eig equivalent of which
  int staticObs = (int) subsetStatic.sum();
  Eigen::VectorXd sigmaResStatic = Eigen::VectorXd::Constant(staticObs,1);
  Eigen::VectorXd sigmaResRoving = Eigen::VectorXd::Constant(N - staticObs,1);
  int wS = 0;
  int wR = 0;
  for(int w = 0; w < N; w++){
    if(subsetStatic(w) == 1) {
      sigmaResStatic(wS) = sigmaRes(w);
      wS++;
    } else {
      sigmaResRoving(wR) = sigmaRes(w);
      wR++;
    }
  }
  
  /*
  /// WRONG! subset of inverse, not inverse of subset
  Eigen::MatrixXd psiStatic = Eigen::MatrixXd::Identity(staticObs, staticObs);
  Eigen::MatrixXd psiRoving = Eigen::MatrixXd::Identity(N - staticObs, N - staticObs);
  wS = 0;
  wR = 0;
  int w2;
  int wScol = 0;
  int wRcol = 0;
  for(int w = 0; w < N; w++){
    if(subsetStatic(w) == 1){
      for(w2 = 0; w2 < N; w2++){
        if(subsetStatic(w2) == 1){
          psiStatic(wS, wScol) = psi(w,w2);
          wScol++;
        }
        wS++;
      }
    } else { // subsetStatic(w) == 0
    for(int w2 = 0; w2 < N; w2++){
      if(subsetStatic(w2) == 0){
        psiStatic(wR, wRcol) = psi(w,w2);
        wRcol++;
      }
      wR++;
    }
    }
  }
  
  // STILL WRONG
  Eigen::MatrixXd whichCells = Eigen::MatrixXd::Constant(N, N, 0.0);
  wS = 0;
  wR = 0;
  int w2;
  int wSrow = 0;
  int wRrow = 0;
  for(int w = 0; w < N; w++){
    if(subsetStatic(w) == 1){
      for(w2 = 0; w2 < N; w2++){
        if(subsetStatic(w2) == 1){
          whichCells(wSrow, wS) = 1.0;
          wSrow++;
        }
        wS++;
      }
    } else { // subsetStatic(w) == 0
    for(int w2 = 0; w2 < N; w2++){
      if(subsetStatic(w2) == 0){
        whichCells(wR, wRrow) = 1.0;
        wRrow++;
      }
      wR++;
    }
    }
  }
  
  */
  
  // Eigen::SelfAdjointEigenSolver<MatrixXd> eig; // specific routine for symmmetric matrices
  Eigen::SelfAdjointEigenSolver<MatrixXd> eig(psi); // specific routine for symmmetric matrices
  
  dTheta(J) = -1*(sigmaResStatic.adjoint())*sigmaResStatic + eig.compute(psi.cwiseProduct(whichCells)).eigenvalues().array().pow(-1).sum();
  dTheta(J+1) = -1*(sigmaResRoving.adjoint())*sigmaResRoving + eig.compute(psi.cwiseProduct(1-whichCells)).eigenvalues().array().pow(-1).sum();
  
  return(dTheta);
}

// Test:
// system.time({dlogProfileCpp(theta0, DTR, NoiTR[,1], XTR, Phi.est[NoiTR[,2],], lamb.est)})
// Compare to
// system.time({dlogProfile(theta0, DTR, NoiTR, XTR, Phi.est, lamb.est)})