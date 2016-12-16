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

#ifdef RENV
#include <Spectra/SymEigsSolver.h> // one or the other will be available
#else
#include <SymEigsSolver.h>
#endif

#include "data.h"

using namespace Eigen;

#define _8GiB 8589934592

#define METHOD_EIGEN 1
#define METHOD_SVD 2

#define KERNEL_LINEAR 1
#define KERNEL_RBF 2

#define MODE_PCA 1
#define MODE_CCA 2
#define MODE_SCCA 3
#define MODE_UCCA 4
#define MODE_CHECK_PCA 5
#define MODE_PREDICT_PCA 6

#define MEM_MODE_OFFLINE 1
#define MEM_MODE_ONLINE 2

#define LOWMEM 1
#define HIGHMEM 2

#define DIVISOR_NONE 0
#define DIVISOR_N1 1
#define DIVISOR_P 2

class RandomPCA {
   public:
      MatrixXd U, V, W, Px, Py;
      VectorXd d;
      MatrixXd V0;
      double trace;
      VectorXd pve;
      MatrixXd X_meansd, Y_meansd;
      ArrayXXd res;
      VectorXd err;
      double mse;
      double rmse;

      int stand_method_x, stand_method_y;
      bool verbose;
      bool debug;
      int divisor;

      void pca_fast(MatrixXd &X, unsigned int block_size,
	    unsigned int ndim,
	    unsigned int maxiter, double tol, long seed,
	    bool do_loadings);
      void pca_fast(Data &dat, unsigned int block_size,
	    unsigned int ndim,
	    unsigned int maxiter, double tol, long seed,
	    bool do_loadings);
      void scca(MatrixXd &X, MatrixXd &Y, double lambda1, double lambda2,
	    long seed, unsigned int ndim, int mem,
	    unsigned int maxiter, double tol);
      void scca(MatrixXd &X, MatrixXd &Y, double lambda1, double lambda2,
	    long seed, unsigned int ndim, int mem,
	    unsigned int maxiter, double tol, MatrixXd &V);
      void scca(Data &dat, double lambda1, double lambda2,
	    long seed, unsigned int ndim, int mem,
	    unsigned int maxiter, double tol,
	    unsigned int block_size);
      void scca(Data &dat, double lambda1, double lambda2,
	    long seed, unsigned int ndim, int mem,
	    unsigned int maxiter, double tol,
	    unsigned int block_size, MatrixXd &V);
      void zca_whiten(bool transpose);
      void ucca(MatrixXd &X, MatrixXd &Y);
      void ucca(Data &dat);
      void check(Data& dat, unsigned int block_size,
	    std::string evec_file, std::string eval_file);
      void check(Data& dat, unsigned int block_size,
	    MatrixXd& evec, VectorXd& eval);
      void check(MatrixXd& X, MatrixXd& evec, VectorXd& eval);
      void project(Data& dat, unsigned int block_size,
	 std::string loadings_file, std::string maf_file,
	 std::string meansd_file);
      void project(Data& dat, unsigned int block_size);
};

