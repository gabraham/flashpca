#pragma once

#include <Eigen/QR>
#include <Eigen/SVD>

using namespace Eigen;

const double tol = 1e-7;

#define METHOD_EIGEN 1
#define METHOD_SVD 2

class RandomPCA {
   public:
      MatrixXd M;
      MatrixXd U, V, W, P;
      VectorXd d;

      int stand_method;

      void pca(MatrixXd &X,
	    int method=METHOD_EIGEN,
	    bool transpose=false,
	    unsigned int ndim=10, unsigned int nextra=10,
	    unsigned int maxiter=10);
      void zca_whiten();
};

