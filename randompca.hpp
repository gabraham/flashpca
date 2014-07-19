#pragma once

#include <Eigen/QR>
#include <Eigen/SVD>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>

using namespace Eigen;

#define METHOD_EIGEN 1
#define METHOD_SVD 2

#define KERNEL_LINEAR 1
#define KERNEL_RBF 2

class RandomPCA {
   public:
      MatrixXd X;
      MatrixXd U, V, W;
      MatrixXd P; // projected X
      VectorXd d;
      double trace;
      VectorXd pve;
      MatrixXd X_meansd;

      int stand_method;
      long seed;
      bool verbose;

      void pca(MatrixXd &X,
	    int method, bool transpose,
	    unsigned int ndim, unsigned int nextra,
	    unsigned int maxiter, double tol, long seed,
	    int kernel, double sigma, bool rbf_center,
	    unsigned int rbf_sample, bool save_kernel,
	    bool do_orth, bool do_loadings);
	    
      void zca_whiten(bool transpose);
};

