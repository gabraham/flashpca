/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2014 Gad Abraham
 * All rights reserved.
 */

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

#define LOWMEM 1
#define HIGHMEM 2

class RandomPCA {
   public:
      MatrixXd U, V, W, Px, Py;
      VectorXd d;
      double trace;
      VectorXd pve;
      MatrixXd X_meansd, Y_meansd;

      int stand_method;
      long seed;
      bool verbose;
      bool debug;

      void pca(MatrixXd &X,
	    int method, bool transpose,
	    unsigned int ndim, unsigned int nextra,
	    unsigned int maxiter, double tol, long seed,
	    int kernel, double sigma, bool rbf_center,
	    unsigned int rbf_sample, bool save_kernel,
	    bool do_orth, bool do_loadings, int mem,
	    bool divide_n);
      void scca(MatrixXd &X, MatrixXd &Y, double lambda1, double lambda2,
	    long seed, unsigned int ndim, int mem,
	    unsigned int maxiter, double tol);
      void scca(MatrixXd &X, MatrixXd &Y, double lambda1, double lambda2,
	    long seed, unsigned int ndim, int mem,
	    unsigned int maxiter, double tol, MatrixXd &V);
      void zca_whiten(bool transpose);
};

