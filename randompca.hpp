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
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/fisher_f.hpp>

#include "spectra/SymEigsSolver.h"


#include "data.hpp"

using namespace Eigen;

#define METHOD_EIGEN 1
#define METHOD_SVD 2

#define KERNEL_LINEAR 1
#define KERNEL_RBF 2

#define MODE_PCA 1
#define MODE_CCA 2
#define MODE_SCCA 3
#define MODE_UCCA 4

#define MEM_MODE_OFFLINE 1
#define MEM_MODE_ONLINE 2

#define LOWMEM 1
#define HIGHMEM 2

class RandomPCA {
   public:
      MatrixXd U, V, W, Px, Py;
      VectorXd d;
      double trace;
      VectorXd pve;
      MatrixXd X_meansd, Y_meansd;
      ArrayXXd res;

      int stand_method_x, stand_method_y;
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
      void pca_fast(MatrixXd &X, unsigned int block_size,
	    int method, bool transpose,
	    unsigned int ndim, unsigned int nextra,
	    unsigned int maxiter, double tol, long seed,
	    int kernel, double sigma, bool rbf_center,
	    unsigned int rbf_sample, bool save_kernel,
	    bool do_orth, bool do_loadings, int mem,
	    bool divide_n);
      void pca_fast(Data &dat, unsigned int block_size,
	    int method, bool transpose,
	    unsigned int ndim, unsigned int nextra,
	    unsigned int maxiter, double tol, long seed,
	    int kernel, double sigma, bool rbf_center,
	    unsigned int rbf_sample, bool save_kernel,
	    bool do_orth, bool do_loadings, int mem,
	    bool divide_n);
      void pca(Data &dat,
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
      void ucca(MatrixXd &X, MatrixXd &Y);
      void ucca(Data &dat);
};

