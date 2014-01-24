#pragma once

#include <Eigen/QR>
#include <Eigen/SVD>

using namespace Eigen;

#define METHOD_EIGEN 1
#define METHOD_SVD 2

class RandomPCA {
   public:
      MatrixXd M;
      MatrixXd U, V, W, P;
      VectorXd d;

      int stand_method;

      void pca(MatrixXd &X,
	    int method, bool transpose,
	    unsigned int ndim, unsigned int nextra,
	    unsigned int maxiter, double tol);
      void zca_whiten(bool transpose);
};

