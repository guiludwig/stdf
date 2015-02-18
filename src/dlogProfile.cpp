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
  // This is weird but: sigma_R I + (sigma_S - sigma_R) 1{static}
  VectorXd RandNoise = theta(J+1)*Eigen::VectorXd::Constant(N,1) + (theta(J) - theta(J+1))*subsetStatic;
  psi += RandNoise.asDiagonal(); // theta has J+2 elements
  Eigen::MatrixXd psiInv = psi.ldlt().solve(Eigen::MatrixXd::Identity(N,N));
  // Eigen::MatrixXd U = psi.llt().matrixL().adjoint(); // same as chol(psi) in R
  // This step finds beta by Generalized Least Squares
  // Eigen::MatrixXd SX = psi.llt().solve(XTR); // A\b by Cholesky's decomposition
  // Eigen::VectorXd beta = ((XTR.adjoint())*SX).ldlt().solve((SX.adjoint())*Y);
  Eigen::VectorXd beta = ((XTR.adjoint())*psiInv*XTR).ldlt().solve((XTR.adjoint())*psiInv*Y);
  // End GLS
  Eigen::VectorXd resid = Y - XTR*beta;
  
  // Difference to logProfile() starts here 
  
  // Eigen::VectorXd sigmaRes = psi.ldlt().solve(resid);
  Eigen::VectorXd sigmaRes = psiInv * resid;
  Eigen::MatrixXd psij = Eigen::MatrixXd::Identity(N, N);
  Eigen::VectorXd dTheta = Eigen::VectorXd::Constant(J+2,0);
  // psij.setZero(N,N);
  // dTheta.resize(J+2);
  
  for(int j = 0; j < J; j++){
    psij = DTR.cwiseProduct((LambEst(j)*((-1*DTR/theta(j)).array().exp().matrix())).cwiseProduct(PhiTime.col(j)*PhiTime.col(j).transpose()))/(theta(j)*theta(j)); // if using c = 1/theta, then -1
    // dTheta(j) = -1*(sigmaRes.adjoint())*psij*sigmaRes + (psi.ldlt().solve(psij)).diagonal().array().sum();
    dTheta(j) = -1*(sigmaRes.adjoint())*psij*sigmaRes + (psiInv*psij).trace();
  }
  // need eig equivalent of which
  int staticObs = (int) subsetStatic.sum();
  
  Eigen::VectorXd diagInverse = psiInv.diagonal();
  double diagS = 0.0;
  double diagR = 0.0;
  
  Eigen::VectorXd sigmaResStatic = Eigen::VectorXd::Constant(staticObs,1);
  Eigen::VectorXd sigmaResRoving = Eigen::VectorXd::Constant(N - staticObs,1);
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
  
  // Eigen::SelfAdjointEigenSolver<MatrixXd> eig; // specific routine for symmmetric matrices
  // Eigen::SelfAdjointEigenSolver<MatrixXd> eig(psi); // specific routine for symmmetric matrices
  
  // dTheta(J) = -1*(sigmaResStatic.adjoint())*sigmaResStatic + eig.compute(psi.cwiseProduct(whichCells)).eigenvalues().array().pow(-1).sum();
  // dTheta(J+1) = -1*(sigmaResRoving.adjoint())*sigmaResRoving + eig.compute(psi.cwiseProduct(1-whichCells)).eigenvalues().array().pow(-1).sum();
//  if(wS == 0){
//    dTheta(J) = 0;
//    dTheta(J+1) = -1*(sigmaResRoving.adjoint())*sigmaResRoving + diagR;
//  } else if(wR == 0) {
//    dTheta(J) = -1*(sigmaResStatic.adjoint())*sigmaResStatic + diagS;
//    dTheta(J+1) = 0;
//  } else {
    dTheta(J) = -1*(sigmaResStatic.adjoint())*sigmaResStatic + diagS;
    dTheta(J+1) = -1*(sigmaResRoving.adjoint())*sigmaResRoving + diagR;
//  }
  
  return(dTheta);
}

// Test:
// system.time({dlogProfileCpp(testTheta <- 1:3, testDTR <- matrix(1,420,420), testY <- 1:420, testXTR <- matrix(1,420,6), testSubsetStatic <- rep(c(0,1), each = 210), testPhiTime <- matrix(rep(1:420,2), ncol=2), testLambEst <- 1:2)})
// replicate(10, system.time({dlogProfileCpp(testTheta <- 1:3, testDTR <- matrix(1,420,420), testY <- 1:420, testXTR <- matrix(1,420,6), testSubsetStatic <- rep(c(0,1), each = 210), testPhiTime <- matrix(rep(1:420,2), ncol=2), testLambEst <- 1:2)})[3])
// system.time({dlogProfileCpp(theta0, DTR, NoiTR[,1], XTR, Phi.est[NoiTR[,2],], lamb.est)})