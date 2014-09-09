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

#define MODE_PCA 1
#define MODE_CCA 2
#define MODE_SCCA 3

#define SCCA_LOWMEM 1
#define SCCA_HIGHMEM 2

class RandomPCA {
   public:
      MatrixXd M;
      MatrixXd U, V, W, Px, Py;
      VectorXd d;
      double trace;
      VectorXd pve;

      int stand_method;
      long seed;
      bool debug;

      void pca(MatrixXd &X,
	    int method, bool transpose,
	    unsigned int ndim, unsigned int nextra,
	    unsigned int maxiter, double tol, long seed,
	    int kernel, double sigma, bool rbf_center,
	    unsigned int rbf_sample, bool save_kernel,
	    bool do_orth);
      void cca(MatrixXd &X, MatrixXd &Y, double lambda1, double lambda2,
	    long seed);
      void scca(MatrixXd &X, MatrixXd &Y, double lambda1, double lambda2,
	    long seed, unsigned int ndim, int scca_method,
	    unsigned int maxiter, double tol);
      void zca_whiten(bool transpose);
};

