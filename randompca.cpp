/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2014-2016 Gad Abraham
 * All rights reserved.
 */

#include "randompca.h"
#include "util.h"
#include "svdwide.h"
#include "svdtall.h"

MatrixXd make_gaussian(unsigned int rows, unsigned int cols, long seed)
{
   boost::random::mt19937 rng;
   rng.seed(seed);
   boost::random::normal_distribution<double> nrm;
   boost::random::variate_generator<boost::random::mt19937&,
      boost::random::normal_distribution<double> > randn(rng, nrm);

   MatrixXd G(rows, cols);
   std::generate(G.data(), G.data() + G.size(), randn);
   return G;
}

// Compute median of pairwise distances on sample of size n from the matrix X
// We're sampling with replacement
// Based on http://www.machinedlearnings.com/2013/08/cosplay.html
// Todo: use Boost accummulators for quantiles
// http://boost-sandbox.sourceforge.net/libs/accumulators/doc/html/accumulators/user_s_guide/the_statistical_accumulators_library.html#accumulators.user_s_guide.the_statistical_accumulators_library.p_square_quantile
double median_dist(MatrixXd& X, unsigned int n, long seed, bool verbose)
{
   boost::random::mt19937 rng;
   rng.seed(seed);
   boost::random::uniform_real_distribution<> dist(0, 1);
   double prop = (double)X.rows() / n;

   verbose && STDOUT << timestamp() <<
      " Computing median Euclidean distance (" << n << " samples)" <<
      std::endl;

   MatrixXd X2(n, X.cols());
   if(n < X.rows())
   {
      verbose && STDOUT << timestamp() << "Sampling" << std::endl;

      // Sample n rows from X
      for(unsigned int i = 0, k = 0 ; i < X.rows() ; i++)
      {
         if(dist(rng) < prop)
         {
            X2.row(k++) = X.row(i);
            if(k == n)
               break;
         }
      }
   }
   else
      X2.noalias() = X;

   VectorXd norms = X2.array().square().rowwise().sum();
   VectorXd ones = VectorXd::Ones(n);
   MatrixXd R = norms * ones.transpose();
   MatrixXd D = R + R.transpose() - 2 * X2 * X2.transpose();

   unsigned int m = D.size();
   double *d = D.data();
   std::sort(d, d + m);
   double med;

   if(m % 2 == 0)
      med = (d[m / 2 - 1] + d[m / 2]) / 2;
   else
      med = d[m / 2];

   verbose && STDOUT << timestamp() << "Median Euclidean distance: "
      << med << std::endl;

   return med;
}

template <typename Derived>
double var(const MatrixBase<Derived>& x)
{
   const unsigned int n = x.size();
   const double mean = x.mean();
   return (x.array() - mean).square().sum() / (n - 1);
}

template <typename DerivedA, typename DerivedB>
RowVectorXd cov(const MatrixBase<DerivedA>& X, const MatrixBase<DerivedB>& Y)
{
   const unsigned int n = X.rows();
   const RowVectorXd xmean = X.colwise().mean();
   const RowVectorXd ymean = Y.colwise().mean();

   return (X.rowwise() - xmean).transpose() * (Y.rowwise() - ymean) / (n - 1);
}

ArrayXXd wilks(const ArrayXd& r2, unsigned int n, unsigned int k)
{
   ArrayXd lambda = 1 - r2;
   boost::math::fisher_f pf(k, n - k - 1);
   ArrayXd pval(r2.size());
   ArrayXXd res(r2.size(), 3);
   res.col(0) = r2.sqrt();
   res.col(1) = (1 - lambda) / lambda * (n - k - 1) / k;

   for(unsigned int i = 0 ; i < lambda.size() ; i++)
   {
      //double fstat = (1 - lambda(i)) / lambda(i) * (n - k - 1) / k;
      res(i, 2) = cdf(complement(pf, res(i, 1)));
   }

   return res;
}

void RandomPCA::pca_fast(MatrixXd& X, unsigned int block_size,
   unsigned int ndim, unsigned int maxiter,
   double tol, long seed, bool do_loadings)
{
   unsigned int N, p;

   X_meansd = standardise(X, stand_method_x, verbose);
   N = X.rows();
   p = X.cols();

   SVDWide op(X, verbose);
   Spectra::SymEigsSolver<double,
      Spectra::LARGEST_ALGE, SVDWide> eigs(&op, ndim, ndim * 2 + 1);

   eigs.init();
   eigs.compute(maxiter, tol);

   double div = 1;
   if(divisor == DIVISOR_N1)
      div = N - 1;
   else if(divisor == DIVISOR_P)
      div = p;

   if(eigs.info() == Spectra::SUCCESSFUL)
   {
      U = eigs.eigenvectors();
      // Note: _eigenvalues_, not singular values
      d = eigs.eigenvalues().array() / div;
      if(do_loadings)
      {
         VectorXd s = d.array().sqrt().inverse() / sqrt(div);
         V.noalias() = X.transpose() * U * s.asDiagonal();
      }
      trace = X.array().square().sum() / div;
      pve = d / trace;
      Px = U * d.array().sqrt().matrix().asDiagonal();

      verbose && STDOUT << timestamp() << "GRM trace: " << trace << std::endl;
   }
   else
   {
      throw new std::runtime_error(
	 std::string("Spectra eigen-decomposition was not successful")
	    + ", status: " + std::to_string(eigs.info()));
   }
}

void RandomPCA::pca_fast(Data& dat, unsigned int block_size,
   unsigned int ndim, unsigned int maxiter, double tol,
   long seed, bool do_loadings)
{
   unsigned int N = dat.N, p = dat.nsnps;
   SVDWideOnline op(dat, block_size, stand_method_x, verbose);
   Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE,
      SVDWideOnline> eigs(&op, ndim, ndim * 2 + 1);

   eigs.init();
   eigs.compute(maxiter, tol);

   double div = 1;
   if(divisor == DIVISOR_N1)
      div = N - 1;
   else if(divisor == DIVISOR_P)
      div = p;

   if(eigs.info() == Spectra::SUCCESSFUL)
   {
      U = eigs.eigenvectors();
      // Note: _eigenvalues_, not singular values
      d = eigs.eigenvalues().array() / div;
      if(do_loadings)
      {
         V = MatrixXd::Zero(dat.nsnps, U.cols());
         verbose && STDOUT << "Computing loadings" << std::endl;
         VectorXd v(dat.nsnps);
         for(unsigned int j = 0 ; j < U.cols() ; j++)
         {
	    verbose && STDOUT << "loading " << j << std::endl;
            VectorXd u = U.col(j);
            op.crossprod(u.data(), v.data());
            double s = d(j);
            V.col(j) = v * (1.0 / sqrt(s)) / sqrt(div);
         }
      }
      trace = op.trace / div;
      pve = d / trace;
      Px = U * d.array().sqrt().matrix().asDiagonal();
      X_meansd = dat.X_meansd; // TODO: duplication

      verbose && STDOUT << timestamp() << "GRM trace: " << trace << std::endl;
   }
   else
   {
      throw new std::runtime_error(
	 std::string("Spectra eigen-decomposition was not successful")
	    + ", status: " + std::to_string(eigs.info()));
   }
}

double inline sign_scalar(double x)
{
   return (0 < x) - (x < 0);
}

VectorXd inline soft_thresh(VectorXd& a, double b)
{
   VectorXd s = a.unaryExpr(std::ptr_fun(sign_scalar));
   VectorXd d = a.array().abs() - b;
   VectorXd z = (d.array() < 0).select(0, d);
   return s.array() * z.array();
}

VectorXd norm_thresh(VectorXd& x, double lambda)
{
   double s = x.norm();
   if(s > 0)
   {
      x = x.array() / s;
      x = soft_thresh(x, lambda);
      s = x.norm();
      if(s > 0)
	 x = x.array() / s;
   }
   return x;
}

void scca_lowmem(MatrixXd& X, MatrixXd &Y, MatrixXd& U, MatrixXd& V,
   VectorXd& d, double lambda1, double lambda2,
   unsigned int maxiter, double tol, bool verbose)
{
   verbose && STDOUT << timestamp() 
      << "[scca_lowmem] " << std::endl;

   // TODO: X2 and Y2 take up lots of memory
   MatrixXd X2 = MatrixXd::Zero(X.rows() + U.cols(), X.cols());
   MatrixXd Y2 = MatrixXd::Zero(Y.rows() + U.cols(), Y.cols());
   X2.block(0, 0, X.rows(), X.cols()) = X;
   Y2.block(0, 0, Y.rows(), Y.cols()) = Y;
   VectorXd u, v, u_old, v_old;

   for(unsigned int j = 0 ; j < U.cols() ; j++)
   {
      if(j > 0)
      {
	 X2.row(X.rows() + j) = sqrt(d[j - 1]) * U.col(j - 1).transpose();
	 Y2.row(Y.rows() + j) = -sqrt(d[j - 1]) * V.col(j - 1).transpose();
      }

      unsigned int iter = 0;
      for( ; iter < maxiter ; iter++)
      {
	 u_old = u = U.col(j);
	 v_old = v = V.col(j);

	 u = X2.transpose() * (Y2 * v);
	 u = norm_thresh(u, lambda1);
	 U.col(j) = u;

	 v = Y2.transpose() * (X2 * U.col(j));
	 v = norm_thresh(v, lambda2);
	 V.col(j) = v;

	 if(iter > 0
	    && (v_old.array() - v.array()).abs().maxCoeff() < tol
	       && (u_old.array() - u.array()).abs().maxCoeff() < tol)
	 {
	    verbose && STDOUT << timestamp() << "dim " << j << " finished in "
	       << iter << " iterations" << std::endl;

	    break;
	 }
      }

      if(iter >= maxiter)
      {
	 verbose && STDOUT << timestamp()
	    << " SCCA did not converge in " << maxiter << " iterations" <<
	    std::endl;
      }

      long long nzu = (U.col(j).array() != 0).count();
      long long nzv = (V.col(j).array() != 0).count();

      verbose && STDOUT << timestamp() << "U_" << j
	 << " non-zeros: " << nzu << ", V_" << j
	 << " non-zeros: " << nzv << std::endl;

      // Use X and Y, not X2 and Y2
      d[j] = (X * U.col(j)).transpose() * (Y * V.col(j));
      verbose && STDOUT << timestamp() << "d[" << j << "]: "
	 << d[j] << std::endl;
   }
}

void scca_highmem(MatrixXd& X, MatrixXd &Y, MatrixXd& U, MatrixXd& V,
   VectorXd& d, double lambda1, double lambda2,
   unsigned int maxiter, double tol, bool verbose)
{
   verbose && STDOUT << timestamp() 
      << "[scca_highmem] Begin computing X^T Y" << std::endl;

   MatrixXd XY = X.transpose() * Y;

   verbose && STDOUT << timestamp() 
      << "[scaa_highmem] End computing X^T Y" << std::endl;

   MatrixXd XYj;
   VectorXd u, v, u_old, v_old;

   for(unsigned int j = 0 ; j < U.cols() ; j++)
   {
      verbose && STDOUT << timestamp() << "dim " << j << std::endl;

      if(j == 0)
	 XYj = XY;
      else
	 XYj = XYj - d[j - 1] * U.col(j - 1) * V.col(j - 1).transpose();

      unsigned int iter = 0;
      for(; iter < maxiter ; iter++)
      {
	 u_old = u = U.col(j);
	 v_old = v = V.col(j);

	 u = XYj * v;
	 u = norm_thresh(u, lambda1);
	 U.col(j) = u;

	 v = XYj.transpose() * U.col(j);
	 v = norm_thresh(v, lambda2);
	 V.col(j) = v;

	 if(iter > 0
	    && (v_old.array() - v.array()).abs().maxCoeff() < tol
	       && (u_old.array() - u.array()).abs().maxCoeff() < tol)
	 {
	    verbose && STDOUT << timestamp() << "dim " << j << " finished in "
	       << iter << " iterations" << std::endl;
	    break;
	 }
      }

      if(iter >= maxiter)
      {
	 verbose && STDOUT << timestamp()
	    << "SCCA did not converge in " << maxiter << " iterations" <<
	    std::endl;
      }

      long long nzu = (U.col(j).array() != 0).count();
      long long nzv = (V.col(j).array() != 0).count();

      verbose && STDOUT << timestamp() << "U_" << j
	 << " non-zeros: " << nzu << ", V_" << j
	 << " non-zeros: " << nzv << std::endl;

      d[j] = U.col(j).transpose() * XYj * V.col(j);
   }
}

void RandomPCA::scca(MatrixXd &X, MatrixXd &Y, double lambda1, double lambda2,
   long seed, unsigned int ndim, int mem, unsigned int maxiter, double tol)
{
   unsigned int k = Y.cols();
   MatrixXd M = make_gaussian(k, ndim, seed);
   this->scca(X, Y, lambda1, lambda2, seed, ndim, mem, maxiter, tol, M);
}

void RandomPCA::scca(MatrixXd &X, MatrixXd &Y, double lambda1, double lambda2,
   long seed, unsigned int ndim, int mem, unsigned int maxiter, double tol,
   MatrixXd &V0)
{
   X_meansd = standardise(X, stand_method_x);
   Y_meansd = standardise(Y, stand_method_y);

   verbose && STDOUT << timestamp() << "dim(X): " << dim(X) << std::endl;
   verbose && STDOUT << timestamp() << "dim(Y): " << dim(Y) << std::endl;
   verbose && STDOUT << timestamp() << "lambda1: " << lambda1
      << " lambda2: " << lambda2 << std::endl;

   unsigned int p = X.cols();

   //V = make_gaussian(k, ndim, seed);
   this->V0 = V0;
   V = V0;
   U = MatrixXd::Zero(p, ndim);
   d = VectorXd::Zero(ndim);

   if(mem == HIGHMEM)
      scca_highmem(X, Y, U, V, d, lambda1, lambda2, maxiter, tol, verbose);
   else
      scca_lowmem(X, Y, U, V, d, lambda1, lambda2, maxiter, tol, verbose);

   Px = X * U;
   Py = Y * V;
}

void RandomPCA::scca(Data &dat, double lambda1, double lambda2,
   long seed, unsigned int ndim, int mem, unsigned int maxiter, double tol,
   unsigned int block_size)
{
   unsigned int k = dat.Y.cols();
   MatrixXd M = make_gaussian(k, ndim, seed);
   this->scca(dat, lambda1, lambda2, seed, ndim, mem, maxiter, tol, block_size, M);
}

// Block loading of X (genotypes)
// Assumes that data.Y has been set
void RandomPCA::scca(Data &dat, double lambda1, double lambda2,
   long seed, unsigned int ndim, int mem, unsigned int maxiter, double tol,
   unsigned int block_size, MatrixXd &V0)
{
   Y_meansd = standardise(dat.Y, stand_method_y);

   verbose && STDOUT << timestamp() << "dim(Y): " << dim(dat.Y) << std::endl;
   verbose && STDOUT << timestamp() << "lambda1: " << lambda1
      << " lambda2: " << lambda2 << std::endl;

   SVDWideOnline op(dat, block_size, stand_method_x, verbose);

   unsigned int p = dat.nsnps;

   this->V0 = V0;
   V = V0;
   U = MatrixXd::Zero(p, ndim);
   d = VectorXd::Zero(ndim);

   VectorXd uj, vj, uj_old, vj_old;

   for(unsigned int j = 0 ; j < U.cols() ; j++)
   {
      unsigned int iter = 0;
      for( ; iter < maxiter ; iter++)
      {
	 uj_old = uj = U.col(j);
	 vj_old = vj = V.col(j);

	 // u = X.transpose() * (Y * v);
	 MatrixXd Yvj = dat.Y * vj;
	 uj = op.crossprod2(Yvj);
	 
	 // deflate u
	 if(j > 0)
	 {
	    uj -= U.leftCols(j) * d.head(j - 1).asDiagonal()
		  * V.leftCols(j).transpose() * vj;
	 }

	 uj = norm_thresh(uj, lambda1);
	 U.col(j) = uj;

	 // v = Y.transpose() * (X * u);
	 VectorXd Xuj = op.prod3(uj);
	 vj = dat.Y.transpose() * Xuj;

	 // deflate v
	 if(j > 0)
	 {
	    vj -= V.leftCols(j) * d.head(j - 1).asDiagonal()
		  * U.leftCols(j).transpose() * uj;
	 }

	 vj = norm_thresh(vj, lambda2);
	 V.col(j) = vj;

	 if(iter > 0
	    && (vj_old.array() - vj.array()).abs().maxCoeff() < tol
	       && (uj_old.array() - uj.array()).abs().maxCoeff() < tol)
	 {
	    verbose && STDOUT << timestamp() << "dim "
	       << j << " finished in "
	       << iter << " iterations" << std::endl;

	    break;
	 }
      }

      if(iter >= maxiter)
      {
	 verbose && STDOUT << timestamp()
	    << " SCCA did not converge in " << maxiter
	    << " iterations" << std::endl;
      }

      long long nzu = (U.col(j).array() != 0).count();
      long long nzv = (V.col(j).array() != 0).count();

      verbose && STDOUT << timestamp() << "U_" << j
	 << " non-zeros: " << nzu << ", V_" << j
	 << " non-zeros: " << nzv << std::endl;

      VectorXd Xuj = op.prod3(U.col(j));
      d[j] = Xuj.transpose() * (dat.Y * V.col(j));
      verbose && STDOUT << timestamp() << "d[" << j << "]: "
	 << d[j] << std::endl;
   }

   Px = op.prod3(U);
   Py = dat.Y * V;
}

// Single-SNP CCA (like plink.multivariate), offline version (loading all SNPs
// into memory)
void RandomPCA::ucca(MatrixXd &X, MatrixXd &Y)
{
   X_meansd = standardise(X, stand_method_x);
   Y_meansd = standardise(Y, stand_method_y);

   unsigned int n = X.rows();
   unsigned int p = X.cols();
   double varx;
   RowVectorXd covXY;

   JacobiSVD<MatrixXd> svd(Y, ComputeThinU | ComputeThinV);
   ArrayXd d = svd.singularValues();
   MatrixXd V = svd.matrixV();
   ArrayXd s(V.cols());
   ArrayXd r2 = ArrayXd(p);

   for(unsigned int j = 0 ; j < p ; j++)
   {
      varx = var(X.col(j));
      covXY = cov(X.col(j), Y);

      // take absolute value to prevent numerical issues with negative numbers
      // close to zero
      // r2(j) = fabs(1.0 / varx * covXY * covYinv * covXY.transpose());

      // covXY * covYinv * covXY'
      // = covXY * (V * D^(-2) * V') * covXY'
      // = ss'
      // where s = covXY * V * D^(-1)
      //
      // The factor of sqrt(n - 1) comes from the fact that we did SVD of Y,
      // but cov(Y) is 1/(n-1) Y'Y when Y is standardised.
      s = covXY * V * std::sqrt(double(n - 1));
      r2(j) = std::abs((s / d.array()).square().sum() / varx);
   }

   res = wilks(r2, X.rows(), Y.cols());
}

// Single-SNP CCA (like plink.multivariate), online version (loading one SNP
// at a time)
// 
// Assumes data.Y has been set
void RandomPCA::ucca(Data& data)
{
   Y_meansd = standardise(data.Y, stand_method_y);

   unsigned int n = data.N;
   unsigned int p = data.nsnps;
   double varx;
   RowVectorXd covXY;

   verbose && STDOUT << timestamp()
      << "UCCA online mode, N=" << n << " p=" << p << std::endl;

   // QR might be faster here
   JacobiSVD<MatrixXd> svd(data.Y, ComputeThinU | ComputeThinV);
   ArrayXd d = svd.singularValues();
   MatrixXd V = svd.matrixV();
   ArrayXd s(V.cols());
   ArrayXd r2 = ArrayXd(p);

   for(unsigned int j = 0 ; j < p ; j++)
   {
      // No need to explicitly standardise X, since read_snp_block will
      // already standardise it internally, assuming that
      // data.stand_method_x has been set previously
      data.read_snp_block(j, j, false, true);
      varx = var(data.X.col(0));
      covXY = cov(data.X.col(0), data.Y);

      // take absolute value to prevent numerical issues with negative numbers
      // that are very close to zero
      // r2(j) = fabs(1.0 / varx * covXY * covYinv * covXY.transpose());

      // covXY * covYinv * covXY'
      // = covXY * (V * D^(-2) * V') * covXY'
      // = ss'
      // where s = covXY * V * D^(-1)
      //
      // The factor of sqrt(n - 1) comes from the fact that we did SVD of Y,
      // but cov(Y) is 1/(n-1) Y'Y when Y is standardised.
      s = covXY * V * std::sqrt(double(n - 1));
      r2(j) = std::abs((s / d.array()).square().sum() / varx);
   }

   res = wilks(r2, n, data.Y.cols());
}

void RandomPCA::check(Data& dat, unsigned int block_size,
   std::string evec_file, std::string eval_file)
{
   // Read eigenvalues
   // Expects no header, no rownames, one eigenvalue per row
   verbose && STDOUT << timestamp() << "Loading eigenvalue file '"
       << eval_file << "'" << std::endl;
   NamedMatrixWrapper M1 = read_text(eval_file.c_str(), 1, -1, 0);
   MatrixXd ev = M1.X;
   if(ev.rows() == 0)
      throw std::runtime_error("No eigenvalues found in file");
   VectorXd eval = ev.col(0);

   // Read eigenvectors
   // Expects header (colnames), FID and IID cols
   verbose && STDOUT << timestamp() << "Loading eigenvector file '"
       << evec_file << "'" << std::endl;
   NamedMatrixWrapper M2 = read_text(evec_file.c_str(), 3, -1, 1);
   MatrixXd& evec = M2.X;

   if(evec.rows() != dat.N)
      throw std::runtime_error(
	 std::string("Eigenvector dimension doesn't match data dimension")
	    + " (evec.rows = " + std::to_string(evec.rows())
	    + "; dat.N = " + std::to_string(dat.N) + ")");

   if(eval.size() != evec.cols())
      throw std::runtime_error(
	 "Eigenvector dimension doesn't match the number of eigenvalues");

   check(dat, block_size, evec, eval);
}

// Check the eigenvalues/eigenvectors, computing the root mean squared error
// of E = X X' U / div - U D^2, i.e., averaged over the n * k dimensions.
// assumes WIDE
void RandomPCA::check(Data& dat, unsigned int block_size,
   MatrixXd& evec, VectorXd& eval)
{
   SVDWideOnline op(dat, block_size, 1, verbose);

   unsigned int K = std::min(evec.cols(), eval.size());

   // X X' U / div = U D^2
   verbose && STDOUT << timestamp()
      << "Checking mean square error between (X X' U) / div and (U D^2)"
      << " for " << K << " dimensions"
      << std::endl;

   double div = 1;
   if(divisor == DIVISOR_N1)
      div = dat.N - 1;
   else if(divisor == DIVISOR_P)
      div = dat.nsnps;

   MatrixXd XXU = op.perform_op_mat(evec);
   XXU /= div;
   MatrixXd UD2 = evec * eval.asDiagonal();

   RowVectorXd rerr = (XXU - UD2).colwise().squaredNorm();
   err = rerr.transpose();

   for(unsigned int j = 0 ; j < K ; j++)
   {
      verbose && STDOUT << timestamp() << "eval(" << (j + 1)
         << "): " << eval(j) << ", sum squared error: "
         << err(j) << std::endl;
   }

   mse = err.sum() / (dat.N * K);
   rmse = std::sqrt(mse);

   verbose && STDOUT << timestamp() << "Mean squared error: " << mse
      << ", Root mean squared error: " << rmse
      << " (n=" << dat.N << ")" << std::endl;

}

void RandomPCA::check(MatrixXd& X, MatrixXd& evec, VectorXd& eval)
{
   X_meansd = standardise(X, stand_method_x, verbose);
   unsigned int K = std::min(evec.cols(), eval.size());
   unsigned int N = X.rows();
   unsigned int p = X.cols();

   // X X' U / div = U D^2
   verbose && STDOUT << timestamp()
      << "Checking mean square error between (X X' U) and (U D^2)"
      << " for " << K << " dimensions"
      << std::endl;

   double div = 1;
   if(divisor == DIVISOR_N1)
      div = N - 1;
   else if(divisor == DIVISOR_P)
      div = p;

   MatrixXd XXU = X * (X.transpose() * evec) / div;
   MatrixXd UD2 = evec * eval.asDiagonal();

   RowVectorXd rerr = (XXU - UD2).colwise().squaredNorm();
   err = rerr.transpose();

   for(unsigned int j = 0 ; j < K ; j++)
   {
      verbose && STDOUT << timestamp() << "eval(" << (j + 1)
         << "): " << eval(j) << ", sum squared error: "
         << err(j) << std::endl;
   }

   mse = err.sum() / (N * K);
   rmse = std::sqrt(mse);

   verbose && STDOUT << timestamp() << "Mean squared error: " << mse
      << ", Root mean squared error: " << rmse
      << " (n=" << N << ")" << std::endl;
}

MatrixXd maf2meansd(MatrixXd maf)
{
   MatrixXd X_meansd(maf.rows(), 2);
   X_meansd.col(0) = maf * 2.0;
   X_meansd.col(1) = maf.array() * 2.0 * (1.0 - maf.array());
   return X_meansd;
}

void RandomPCA::project(Data& dat, unsigned int block_size,
   std::string loadings_file, std::string maf_file,
   std::string meansd_file)
{
   // Read the loadings
   // TODO: check that SNP ids match
   NamedMatrixWrapper M = read_text(loadings_file.c_str(), 3, -1, 1);
   V = M.X;

   // Read the means+sds or the MAF (and convert MAF to means+sds)
   if(maf_file != "")
   {
      // TODO: missing/non-numeric values?
      verbose && STDOUT << timestamp() << "Reading MAF file "
	 << maf_file << std::endl;
      NamedMatrixWrapper M2 = read_text(maf_file.c_str(), 3, -1, 1);
      dat.X_meansd = maf2meansd(M2.X);
      dat.use_preloaded_maf = true;
   }
   else if(meansd_file != "")
   {
      verbose && STDOUT << timestamp()
	 << " Reading mean/stdev file " << meansd_file << std::endl;
      NamedMatrixWrapper M2 = read_text(meansd_file.c_str(), 3, -1, 1);
      dat.X_meansd = M2.X;
      dat.use_preloaded_maf = true;
   }
   else
   {
      verbose && STDOUT << timestamp()
	 << " Using MAF from the data" << std::endl;
      dat.use_preloaded_maf = false;
   }

   project(dat, block_size);
}

// Project new samples onto existing principal components.
//
// Doesn't do a lot of sanity checking.
//
// Assumes:
// - The loadings matrix V must have been set already
// - dat.X_meansd has beeen set
// - dat.use_preloaded_maf has been set, if needed
void RandomPCA::project(Data& dat, unsigned int block_size)
{
   // Check that the SNPs in the data match

   SVDWideOnline op(dat, block_size, 1, verbose);

   unsigned int k = V.cols();
   Px = MatrixXd::Zero(dat.N, k);

   double div = 1;
   if(divisor == DIVISOR_N1)
      div = dat.N - 1;
   else if(divisor == DIVISOR_P)
      div = V.rows();

   VectorXd pxi(dat.N);
   for(unsigned int i = 0 ; i < k ; i++)
   {
      MatrixXd v = V.col(i);
      op.prod(v.data(), pxi.data());
      Px.col(i) = pxi.array() / sqrt(div); // X V = U D
   }
}

