/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2014-2016 Gad Abraham
 * All rights reserved.
 */

#include "randompca.hpp"
#include "util.hpp"

class SVDWide
{
   private:
      const MatrixXd& mat;
      const unsigned int n;
      bool verbose;
      unsigned int nops;

   public:
      SVDWide(const MatrixXd& mat_, bool verbose_=false):
	 mat(mat_), n(mat_.rows())
      {
	 verbose = verbose_;
	 nops = 1;
      }

      inline unsigned int rows() const { return n; }
      inline unsigned int cols() const { return n; }

      void perform_op(double *x_in, double* y_out)
      {
	 Map<VectorXd> x(x_in, n);
	 Map<VectorXd> y(y_out, n);
	 verbose && STDOUT << timestamp() << "Matrix operation "
	    << nops << std::endl;
	 y.noalias() = mat * (mat.transpose() * x);
	 nops++;
      }
};

class SVDWideOnline
{
   public:
      // Trace of X X'
      double trace;

   private:
      Data& dat;
      const unsigned int n, p;
      unsigned int nblocks;
      unsigned int *start, *stop;
      int stand_method;
      bool verbose;
      unsigned int nops;
      unsigned int block_size;
      bool trace_done;

   public:
      SVDWideOnline(Data& dat_, unsigned int block_size_, int stand_method_,
	 bool verbose_): dat(dat_), n(dat_.N), p(dat_.nsnps)
      {
	 verbose = verbose_;
	 block_size = block_size_;
	 stand_method = stand_method_;
	 nblocks = (unsigned int)ceil((double)p / block_size);
	 verbose && STDOUT << timestamp()
	    << "Using blocksize " << block_size << ", " <<
	    nblocks << " blocks"<< std::endl;
	 start = new unsigned int[nblocks];
	 stop = new unsigned int[nblocks];
	 for(unsigned int i = 0 ; i < nblocks ; i++)
	 {
	    start[i] = i * block_size;
	    stop[i] = start[i] + block_size - 1;
	    stop[i] = stop[i] >= p ? p - 1 : stop[i];
	 }

	 nops = 1;
	 trace = 0;
	 trace_done = false;
      }

      ~SVDWideOnline()
      {
	 delete[] start;
	 delete[] stop;
      }

      inline unsigned int rows() const { return n; }
      inline unsigned int cols() const { return n; }

      // y = X X' * x
      void perform_op(double *x_in, double* y_out)
      {
	 Map<VectorXd> x(x_in, n);
	 Map<VectorXd> y(y_out, n);
	 unsigned int actual_block_size;

	 verbose && STDOUT << timestamp()
	    << "Matrix operation " << nops << std::endl;

	 actual_block_size = stop[0] - start[0] + 1;

	 // If we only have one block, keep it in memory instead of reading it
	 // over again
	 if(nblocks > 1 || nops == 1)
	 {
	    dat.read_snp_block(start[0], stop[0], false, false);
	    verbose && STDOUT << timestamp() << "Reading block " <<
	       0 << " (" << start[0] << ", " << stop[0]
	       << ")"  << std::endl;
	 }

	 y.noalias() = dat.X.leftCols(actual_block_size) *
	    (dat.X.leftCols(actual_block_size).transpose() * x);
	 if(!trace_done)
	    trace = dat.X.leftCols(actual_block_size).array().square().sum();

	 // If there's only one block, this loop doesn't run anyway
	 for(unsigned int k = 1 ; k < nblocks ; k++)
	 {
	    verbose && STDOUT << timestamp() << "Reading block " <<
	       k << " (" << start[k] << ", " << stop[k] << ")"  << std::endl;
	    actual_block_size = stop[k] - start[k] + 1;
	    dat.read_snp_block(start[k], stop[k], false, false);
	    //TODO: Kahan summation better here?
	    y.noalias() = y + dat.X.leftCols(actual_block_size) *
	       (dat.X.leftCols(actual_block_size).transpose() * x);
	    if(!trace_done)
	       trace += dat.X.leftCols(actual_block_size).array().square().sum();
	 }

	 if(!trace_done)
	    trace_done = true;

	 nops++;
      }

      // y = X X' * x
      MatrixXd perform_op_mat(MatrixXd x)
      {
	 unsigned int actual_block_size;

	 verbose && STDOUT << timestamp()
	    << "Matrix operation " << nops << std::endl;

	 actual_block_size = stop[0] - start[0] + 1;

	 // If we only have one block, keep it in memory instead of reading it
	 // over again
	 if(nblocks > 1 || nops == 1)
	 {
	    dat.read_snp_block(start[0], stop[0], false, false);
	    verbose && STDOUT << timestamp() << "Reading block " <<
	       0 << " (" << start[0] << ", " << stop[0]
	       << ")"  << std::endl;
	 }
	 
	 MatrixXd Y(n, x.cols());
	 Y.noalias() = dat.X.leftCols(actual_block_size) *
	    (dat.X.leftCols(actual_block_size).transpose() * x);
	 if(!trace_done)
	    trace = dat.X.leftCols(actual_block_size).array().square().sum();

	 // If there's only one block, this loop doesn't run anyway
	 for(unsigned int k = 1 ; k < nblocks ; k++)
	 {
	    verbose && STDOUT << timestamp() << "Reading block " <<
	       k << " (" << start[k] << ", " << stop[k] << ")"  << std::endl;
	    actual_block_size = stop[k] - start[k] + 1;
	    dat.read_snp_block(start[k], stop[k], false, false);
	    //TODO: Kahan summation better here?
	    Y.noalias() = Y + dat.X.leftCols(actual_block_size) *
	       (dat.X.leftCols(actual_block_size).transpose() * x);
	    if(!trace_done)
	       trace += dat.X.leftCols(actual_block_size).array().square().sum();
	 }

	 if(!trace_done)
	    trace_done = true;

	 nops++;
	 return Y;
      }
      // Like R crossprod(): y = X' * x
      // Note: size of x must be number of samples, size y must number of SNPs
      void crossprod(double *x_in, double *y_out)
      {
	 Map<VectorXd> x(x_in, n);
	 Map<VectorXd> y(y_out, p);
	 unsigned int actual_block_size = stop[0] - start[0] + 1;

	 verbose && STDOUT << timestamp()
	    << "Matrix operation " << nops << std::endl;

	 dat.read_snp_block(start[0], stop[0], false, false);
	 verbose && STDOUT << timestamp() << "Reading block " <<
	    0 << " (" << start[0] << ", " << stop[0]
	    << ")"  << std::endl;

	 y.segment(start[0], actual_block_size) 
	    = dat.X.leftCols(actual_block_size).transpose() * x;

	 for(unsigned int k = 1 ; k < nblocks ; k++)
	 {
	    verbose && STDOUT << timestamp() << "Reading block " <<
	       k << " (" << start[k] << ", " << stop[k] << ")"  << std::endl;
	    actual_block_size = stop[k] - start[k] + 1;
	    dat.read_snp_block(start[k], stop[k], false, false);
	    //TODO: Kahan summation better here?
	    y.segment(start[k], actual_block_size)
	       = dat.X.leftCols(actual_block_size).transpose() * x;
	 }
	 nops++;
      }

      // Like R crossprod(): y = X' * x
      // Note: size of x must be number of samples, size y must number of SNPs
      MatrixXd crossprod2(MatrixXd& x)
      {
	 unsigned int actual_block_size = stop[0] - start[0] + 1;

	 verbose && STDOUT << timestamp()
	    << "Matrix operation " << nops << std::endl;

	 dat.read_snp_block(start[0], stop[0], false, false);
	 verbose && STDOUT << timestamp() << "Reading block " <<
	    0 << " (" << start[0] << ", " << stop[0]
	    << ")"  << std::endl;

	 MatrixXd Y(p, x.cols());
	 Y.middleRows(start[0], actual_block_size) =
	    dat.X.leftCols(actual_block_size).transpose() * x;

	 for(unsigned int k = 1 ; k < nblocks ; k++)
	 {
	    verbose && STDOUT << timestamp() << "Reading block " <<
	       k << " (" << start[k] << ", " << stop[k] << ")"  << std::endl;
	    actual_block_size = stop[k] - start[k] + 1;
	    dat.read_snp_block(start[k], stop[k], false, false);
	    //TODO: Kahan summation better here?
	    Y.middleRows(start[k], actual_block_size) =
	       dat.X.leftCols(actual_block_size).transpose() * x;
	 }
	 nops++;
	 return Y;
      }

      // Like y = X %*% x
      // Note: size of x must be number of SNPs,
      // size y must number of samples
      void prod(double *x_in, double *y_out)
      {
	 Map<VectorXd> x(x_in, p);
	 Map<VectorXd> y(y_out, n);
	 unsigned int actual_block_size = stop[0] - start[0] + 1;

	 verbose && STDOUT << timestamp()
	    << "Matrix operation " << nops << std::endl;

	 dat.read_snp_block(start[0], stop[0], false, false);
	 verbose && STDOUT << timestamp() << "Reading block " <<
	    0 << " (" << start[0] << ", " << stop[0]
	    << ")"  << std::endl;

	 y.noalias() =
	    dat.X.leftCols(actual_block_size)
	       * x.segment(start[0], actual_block_size);

	 for(unsigned int k = 1 ; k < nblocks ; k++)
	 {
	    verbose && STDOUT << timestamp() << "Reading block " <<
	       k << " (" << start[k] << ", " << stop[k] << ")"  << std::endl;
	    actual_block_size = stop[k] - start[k] + 1;
	    dat.read_snp_block(start[k], stop[k], false, false);
	    //TODO: Kahan summation better here?
	    y.noalias() =
	       y + dat.X.leftCols(actual_block_size)
		  * x.segment(start[k], actual_block_size);
	 }
	 nops++;
      }

      // return Y = X * X' * x where x is a matrix (despite x being lower case)
      MatrixXd perform_op_multi(MatrixXd& x)
      {
	 unsigned int actual_block_size;

	 verbose && STDOUT << timestamp()
	    << "Matrix operation " << nops << std::endl;

	 actual_block_size = stop[0] - start[0] + 1;

	 // If we only have one block, keep it in memory instead of reading it
	 // over again
	 if(nblocks > 1 || nops == 1)
	 {
	    dat.read_snp_block(start[0], stop[0], false, false);
	    verbose && STDOUT << timestamp() << "Reading block " <<
	       0 << " (" << start[0] << ", " << stop[0]
	       << ")"  << std::endl;
	 }

	 MatrixXd Y = dat.X.leftCols(actual_block_size) *
	    (dat.X.leftCols(actual_block_size).transpose() * x);
	 if(!trace_done)
	    trace = dat.X.leftCols(actual_block_size).array().square().sum();

	 // If there's only one block, this loop doesn't run anyway
	 for(unsigned int k = 1 ; k < nblocks ; k++)
	 {
	    verbose && STDOUT << timestamp() << "Reading block " <<
	       k << " (" << start[k] << ", " << stop[k] << ")"  << std::endl;
	    actual_block_size = stop[k] - start[k] + 1;
	    dat.read_snp_block(start[k], stop[k], false, false);
	    //TODO: Kahan summation better here?
	    Y.noalias() = Y + dat.X.leftCols(actual_block_size) *
	       (dat.X.leftCols(actual_block_size).transpose() * x);
	    if(!trace_done)
	       trace += dat.X.leftCols(actual_block_size).array().square().sum();
	 }

	 nops++;
	 if(!trace_done)
	    trace_done = true;

	 return Y;
      }

      // Like Y = x' * X, where X is genotypes, x is a matrix
      MatrixXd prod2(MatrixXd& x)
      {
	 unsigned int actual_block_size = stop[0] - start[0] + 1;

	 verbose && STDOUT << timestamp()
	    << "Matrix operation " << nops << std::endl;

	 dat.read_snp_block(start[0], stop[0], false, false);
	 verbose && STDOUT << timestamp() << "Reading block " <<
	    0 << " (" << start[0] << ", " << stop[0]
	    << ")"  << std::endl;

	 MatrixXd Y(x.cols(), p);
	 Y.middleCols(start[0], actual_block_size) =
	    x.transpose() * dat.X.leftCols(actual_block_size);

	 for(unsigned int k = 1 ; k < nblocks ; k++)
	 {
	    verbose && STDOUT << timestamp() << "Reading block " <<
	       k << " (" << start[k] << ", " << stop[k] << ")"  << std::endl;
	    actual_block_size = stop[k] - start[k] + 1;
	    dat.read_snp_block(start[k], stop[k], false, false);
	    //TODO: Kahan summation better here?
	    Y.middleCols(start[k], actual_block_size) =
	       x.transpose() * dat.X.leftCols(actual_block_size);
	 }
	 nops++;
	 return Y;
      }

      // Like Y = X * x, where X is genotypes, x is a matrix
      MatrixXd prod3(MatrixXd& x)
      {
	 unsigned int actual_block_size = stop[0] - start[0] + 1;

	 verbose && STDOUT << timestamp()
	    << "Matrix operation " << nops << std::endl;

	 dat.read_snp_block(start[0], stop[0], false, false);
	 verbose && STDOUT << timestamp() << "Reading block " <<
	    0 << " (" << start[0] << ", " << stop[0]
	    << ")"  << std::endl;

	 MatrixXd Y = dat.X.leftCols(actual_block_size) 
	    * x.middleRows(start[0], actual_block_size);

	 for(unsigned int k = 1 ; k < nblocks ; k++)
	 {
	    verbose && STDOUT << timestamp() << "Reading block " <<
	       k << " (" << start[k] << ", " << stop[k] << ")"  << std::endl;
	    actual_block_size = stop[k] - start[k] + 1;
	    dat.read_snp_block(start[k], stop[k], false, false);
	    //TODO: Kahan summation better here?
	    Y.noalias() =
	       Y + dat.X.leftCols(actual_block_size) 
		  * x.middleRows(start[k], actual_block_size);
	 }
	 nops++;
	 return Y;
      }
};

MatrixXd make_gaussian(unsigned int rows, unsigned int cols, long seed)
{
   boost::random::mt19937 rng;
   rng.seed(seed);
   boost::random::normal_distribution<double> nrm;
   boost::random::variate_generator<boost::random::mt19937&,
      boost::random::normal_distribution<double> > rand(rng, nrm);

   MatrixXd G(rows, cols);
   for(unsigned int i = 0 ; i < rows ; i++)
      for(unsigned int j = 0 ; j < cols ; j++)
	 G(i, j) = rand();
   return G;
}

// normalize each column of X to unit l2 norm
inline void normalize(MatrixXd& X)
{
   unsigned int p = X.cols();
   for(unsigned int j = 0 ; j < p ; j++)
   {
      double s = 1 / sqrt(X.col(j).array().pow(2).sum());
      X.col(j) = X.col(j).array() * s;
   }
}

void pca_small(MatrixXd &B, int method, MatrixXd& U, VectorXd &d, bool verbose)
{
   if(method == METHOD_SVD)
   {
      verbose && STDOUT << timestamp() << "SVD begin" << std::endl;

      JacobiSVD<MatrixXd> svd(B, ComputeThinU | ComputeThinV);
      U = svd.matrixU();
      MatrixXd V = svd.matrixV();
      d = svd.singularValues().array().pow(2);

      verbose && STDOUT << timestamp() << "SVD done" << std::endl;
   }
   else if(method == METHOD_EIGEN)
   {
      verbose && STDOUT << timestamp() << "Eigen-decomposition begin" << std::endl;

      MatrixXd BBT = B * B.transpose();

      verbose && STDOUT << timestamp() << "dim(BBT): " << dim(BBT) << std::endl;

      SelfAdjointEigenSolver<MatrixXd> eig(BBT);

      // The eigenvalues come out sorted in *increasing* order,
      // but we need decreasing order
      VectorXd eval = eig.eigenvalues();
      MatrixXd evec = eig.eigenvectors();
      d.resize(eval.size());
      U.resize(BBT.rows(), BBT.rows());

      unsigned int k = 0;
      for(unsigned int i = d.size() - 1 ; i != -1 ; --i)
      {
	 // we get eigenvalues, which are the squared singular values
	 d(k) = eval(i);
	 U.col(k) = evec.col(i);
	 k++;
      }
   }
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
   const double mean = x.array().sum() / n;
   return (x.array() - mean).array().pow(2).sum() / (n - 1);
}

template <typename DerivedA, typename DerivedB>
RowVectorXd cov(const MatrixBase<DerivedA>& X, const MatrixBase<DerivedB>& Y)
{
   const unsigned int n = X.rows();
   const RowVectorXd xmean = X.colwise().sum() / n;
   const RowVectorXd ymean = Y.colwise().sum() / n;

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

void RandomPCA::pca_fast(MatrixXd& X, unsigned int block_size, int method,
   unsigned int ndim, unsigned int nextra, unsigned int maxiter,
   double tol, long seed, bool do_loadings, int mem)
{
   unsigned int N, p;

   if(stand_method_x != STANDARDIZE_NONE)
      X_meansd = standardize(X, stand_method_x, verbose);
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
   }
   else
   {
      std::cerr 
	 << "Spectra eigen-decomposition was not successful, status: " 
	 << eigs.info() << std::endl;
   }
}

void RandomPCA::pca_fast(Data& dat, unsigned int block_size, int method,
   unsigned int ndim, unsigned int nextra,
   unsigned int maxiter, double tol, long seed, bool do_loadings, int mem)
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
         std::cout << "Computing loadings" << std::endl;
         VectorXd v(dat.nsnps);
         for(unsigned int j = 0 ; j < U.cols() ; j++)
         {
            std::cout << "loading " << j << std::endl;
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
   }
   else
   {
      std::cerr
	 << "Spectra eigen-decomposition was not successful, status: " 
	 << eigs.info() << std::endl;
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
   verbose && STDOUT << timestamp() << "Begin computing X^T Y" << std::endl;

   MatrixXd XY = X.transpose() * Y;

   verbose && STDOUT << timestamp() << "End computing X^T Y" << std::endl;

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
   if(stand_method_x != STANDARDIZE_NONE)
      X_meansd = standardize(X, stand_method_x);

   if(stand_method_y != STANDARDIZE_NONE)
      Y_meansd = standardize(Y, stand_method_y);

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

   //if(mem == HIGHMEM)
   //   scca_highmem(X, Y, U, V, d, lambda1, lambda2, maxiter, tol, verbose);
   //else
      scca_lowmem(X, Y, U, V, d, lambda1, lambda2, maxiter, tol, verbose);

   Px = X * U;
   Py = Y * V;
}

void RandomPCA::scca(Data &dat, MatrixXd &Y, double lambda1, double lambda2,
   long seed, unsigned int ndim, int mem, unsigned int maxiter, double tol)
{
   unsigned int k = Y.cols();
   MatrixXd M = make_gaussian(k, ndim, seed);
   this->scca(dat, Y, lambda1, lambda2, seed, ndim, mem, maxiter, tol, M);
}

// Block loading of X (genotypes)
void RandomPCA::scca(Data &dat, MatrixXd &Y, double lambda1, double lambda2,
   long seed, unsigned int ndim, int mem, unsigned int maxiter, double tol,
   MatrixXd &V0)
{
//   //if(stand_method_x != STANDARDIZE_NONE)
//   //   X_meansd = standardize(X, stand_method_x);
//
//   if(stand_method_y != STANDARDIZE_NONE)
//      Y_meansd = standardize(Y, stand_method_y);
//
//   //verbose && STDOUT << timestamp() << "dim(X): " << dim(X) << std::endl;
//   verbose && STDOUT << timestamp() << "dim(Y): " << dim(Y) << std::endl;
//   verbose && STDOUT << timestamp() << "lambda1: " << lambda1 
//      << " lambda2: " << lambda2 << std::endl;
//
//   unsigned int p = dat.nsnps;
//
//   this->V0 = V0;
//   V = V0;
//   U = MatrixXd::Zero(p, ndim);
//   d = VectorXd::Zero(ndim); 
//
//   VectorXd u, v, u_old, v_old;
//  
//   for(unsigned int j = 0 ; j < U.cols() ; j++)
//   {
//      unsigned int iter = 0;
//      for( ; iter < maxiter ; iter++)
//      {
//	 u_old = u = U.col(j);
//	 v_old = v = V.col(j);
//
//	 if(j == 0)
//	 {
//	    //u = X.transpose() * (Y * v);
//	    MatrixXd Yv = Y * v;
//	    u = crossprod2(Yv);
//	 }
//	 else // deflation
//	 {
//	  TODO: broken
//	    u = (X.array() + sqrt(d[j - 1]) * u.transpose()).transpose() 
//	       * (Y.array() - sqrt(d[j - 1] * v.transpose())  * v);
//	 }
//	 u = norm_thresh(u, lambda1);
//	 U.col(j) = u;
//
//	 if(j == 0)
//	 {
//	    //v = Y.transpose() * (X * u);
//	    VectorXd Xu = prod2(u);
//	    v = Y.transpose() * Xu;
//	 }
//	 else // deflation
//	 {
//	    TODO: broken
//	    // (Y - sqrt(d[j-1)) * v')' * ((X - sqrt(d[j-1]) * u') * u)
//	    v = (Y.array() - sqrt(d[j - 1]) * v.transpose()).transpose() 
//	       * (X.array() + sqrt(d[j - 1] * u.transpose())  * u);
//	 }
//
//	 v = norm_thresh(v, lambda2);
//	 V.col(j) = v;
//
//	 if(iter > 0 
//	    && (v_old.array() - v.array()).abs().maxCoeff() < tol
//	       && (u_old.array() - u.array()).abs().maxCoeff() < tol)
//	 {
//	    verbose && STDOUT << timestamp() << "dim " << j << " finished in "
//	       << iter << " iterations" << std::endl;
//
//	    break;
//	 }
//      }
//
//      if(iter >= maxiter)
//      {
//	 verbose && STDOUT << timestamp()
//	    << " SCCA did not converge in " << maxiter
//	    << " iterations" << std::endl;
//      }
//
//      long long nzu = (U.col(j).array() != 0).count();
//      long long nzv = (V.col(j).array() != 0).count();
//
//      verbose && STDOUT << timestamp() << "U_" << j 
//	 << " non-zeros: " << nzu << ", V_" << j
//	 << " non-zeros: " << nzv << std::endl;
//
//      d[j] = (X * U.col(j)).transpose() * (Y * V.col(j)); 
//      verbose && STDOUT << timestamp() << "d[" << j << "]: "
//	 << d[j] << std::endl;
//   }
//
//   //Px = X * U;
//   //Py = Y * V;
}

// Single-SNP CCA (like plink.multivariate), offline version (loading all SNPs
// into memory)
void RandomPCA::ucca(MatrixXd &X, MatrixXd &Y)
{
   if(stand_method_x != STANDARDIZE_NONE)
      X_meansd = standardize(X, stand_method_x);

   if(stand_method_y != STANDARDIZE_NONE)
      Y_meansd = standardize(Y, stand_method_y);

   unsigned int n = X.rows();
   unsigned int p = X.cols();
   double varx;
   RowVectorXd covXY;

   JacobiSVD<MatrixXd> svd(Y, ComputeThinU | ComputeThinV);
   VectorXd d = svd.singularValues() / sqrt(n - 1);
   VectorXd d2 = d.array().pow(-2);
   MatrixXd V = svd.matrixV();
   MatrixXd covYinv =
      svd.matrixV() * d2.asDiagonal() * svd.matrixV().transpose();
   ArrayXd r2 = ArrayXd(p);

   for(unsigned int j = 0 ; j < p ; j++)
   {
      varx = var(X.col(j));
      covXY = cov(X.col(j), Y);
      
      // take absolute value to prevent numerical issues with negative numbers
      // close to zero
      r2(j) = fabs(1.0 / varx * covXY * covYinv * covXY.transpose());
   }

   res = wilks(r2, X.rows(), Y.cols());
}

// Single-SNP CCA (like plink.multivariate), online version (loading one SNP
// at a time)
void RandomPCA::ucca(Data& data)
{
   if(stand_method_y != STANDARDIZE_NONE)
      Y_meansd = standardize(data.Y, stand_method_y);

   unsigned int n = data.N;
   unsigned int p = data.nsnps;
   double varx;
   RowVectorXd covXY;

   verbose && STDOUT << timestamp()
      << "ucca online mode, N=" << n << " p=" << p << std::endl;

   // QR might be faster here
   JacobiSVD<MatrixXd> svd(data.Y, ComputeThinU | ComputeThinV);
   VectorXd d = svd.singularValues() / sqrt(n - 1);
   VectorXd d2 = d.array().pow(-2);
   MatrixXd V = svd.matrixV();
   MatrixXd covYinv
      = svd.matrixV() * d2.asDiagonal() * svd.matrixV().transpose();
   ArrayXd r2 = ArrayXd(p);

   for(unsigned int j = 0 ; j < p ; j++)
   {
      data.read_snp_block(j, j, false, true);
      if(stand_method_x != STANDARDIZE_NONE)
	 X_meansd = standardize(data.X, stand_method_x);
      varx = var(data.X.col(0));
      covXY = cov(data.X.col(0), data.Y);

      // take absolute value to prevent numerical issues with negative numbers
      // that are very close to zero
      r2(j) = fabs(1.0 / varx * covXY * covYinv * covXY.transpose());
   }

   res = wilks(r2, n, data.Y.cols());
}

// Check the eigenvalues/eigenvectors, computing the root mean squared error
// of E = X X' U / div - U D^2, i.e., averaged over the n * k dimensions.
// assumes WIDE
void RandomPCA::check(Data& dat, unsigned int block_size,
   std::string evec_file, std::string eval_file)
{
   SVDWideOnline op(dat, block_size, 1, verbose);

   NamedMatrixWrapper M1 = read_text(
      eval_file.c_str(), 1, -1, 0);
   MatrixXd ev = M1.X;
   if(ev.rows() == 0)
      throw std::runtime_error("No eigenvalues found in file");

   VectorXd eval = ev.col(0);
   NamedMatrixWrapper M2 = read_text(evec_file.c_str(), 1, -1, 0);
   MatrixXd evec = M2.X;
   MatrixXd out = MatrixXd::Zero(evec.rows(), 1);
   VectorXd uexp;

   if(evec.rows() != dat.N)
      throw std::runtime_error(
	 "Eigenvector dimension doesn't match data dimension");

   if(eval.size() != evec.cols())
      throw std::runtime_error(
	 "Eigenvector dimension doesn't match the number of eigenvalues");

   unsigned int K = fminl(evec.cols(), eval.size());

   // X X' U / div = U D^2
   STDOUT << timestamp()
      << "Checking mean square error between (X X' U) and (U D^2)"
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

   VectorXd err(eval.size());

   for(unsigned int j = 0 ; j < K ; j++)
   {
      err(j) = (XXU.col(j).array() - UD2.col(j).array()).square().sum();
   }

   for(unsigned int j = 0 ; j < K ; j++)
   {
      STDOUT << timestamp() << "eval(" << (j + 1) 
	 << "): " << eval(j) << ", sum squared error: "
	 << err(j) << std::endl;
   }

   double mse = err.array().sum() / (dat.N * K);
   double rmse = sqrt(mse);

   STDOUT << timestamp() << "Mean squared error: " << mse
      << ", Root mean squared error: " << rmse
      << " (n=" << dat.N << ")" << std::endl;
   
}

MatrixXd maf2meansd(MatrixXd maf)
{
   MatrixXd X_meansd(maf.rows(), 2);
   X_meansd.col(0) = maf * 2.0;
   X_meansd.col(1) = maf.array() * 2.0 * (1.0 - maf.array());
   return X_meansd;
}

// Project new samples onto existing principal components.
//
// Doesn't do a lot of sanity checking.
void RandomPCA::project(Data& dat, unsigned int block_size,
   std::string loadings_file, std::string maf_file,
   std::string meansd_file)
{
   std::cout << "[project] loadings_file " << loadings_file << std::endl;

   // Read the PCs
   // TODO: expects no rownames etc, just numeric
   // TODO: missing values?
   // TODO: check that SNP ids match
   NamedMatrixWrapper M = read_text(loadings_file.c_str(), 2, -1, 1);
   V = M.X;

   // Read the means+sds or the MAF (and convert MAF to means+sds)
   if(maf_file != "")
   {
      // TODO: missing/non-numeric values?
      verbose && STDOUT << timestamp() << "Reading MAF file "
	 << maf_file << std::endl;
      NamedMatrixWrapper M2 = read_text(maf_file.c_str(), 2, -1, 1);   
      dat.X_meansd = maf2meansd(M2.X);
      dat.use_preloaded_maf = true;
   }
   else if(meansd_file != "")
   {
      verbose && STDOUT << timestamp()
	 << " Reading mean/stdev file " << meansd_file << std::endl;
      NamedMatrixWrapper M2 = read_text(meansd_file.c_str(), 2, -1, 1);   
      dat.X_meansd = M2.X;
      dat.use_preloaded_maf = true;
   }
   else
   {
      verbose && STDOUT << timestamp()
	 << " Using MAF from the data" << std::endl;
      dat.use_preloaded_maf = false;
   }

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

