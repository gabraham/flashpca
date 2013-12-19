#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <Eigen/SVD>
#include <Eigen/QR>

using namespace Eigen;

const double tol = 1e-6;

class RandomPCA {
   public:
      MatrixXd M;
      MatrixXd U, V, W, P;
      VectorXd d;

      MatrixXd pca(MatrixXd X, bool transpose=false,
	    unsigned int ndim=10, unsigned int nextra=10,
	    unsigned int maxiter=500);
      MatrixXd zca_whiten();
};

