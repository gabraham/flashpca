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
      SVDWide(const MatrixXd& mat_, bool verbose_=false): mat(mat_), n(mat_.rows())
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
	 verbose && STDOUT << timestamp() << " Matrix operation "
	    << nops << std::endl;
	 y.noalias() = mat * (mat.transpose() * x);
	 nops++;
      }
};

class SVDWideOnline
{
   private:
      Data& dat;
      const unsigned int n, p;
      unsigned int nblocks;
      unsigned int *start, *stop;
      int stand_method;
      bool verbose;
      unsigned int nops;
      unsigned int block_size, actual_block_size;

   public:
      SVDWideOnline(Data& dat_, unsigned int block_size_, int stand_method_,
	 bool verbose_): dat(dat_), n(dat_.N), p(dat_.nsnps)
      {
	 verbose = verbose_;
	 block_size = block_size_;
	 nblocks = (unsigned int)ceil((double)p / block_size);
	 STDOUT << "Using blocksize " << block_size << ", " <<
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
      }

      ~SVDWideOnline()
      {
	 delete[] start;
	 delete[] stop;
      }

      inline unsigned int rows() const { return n; }
      inline unsigned int cols() const { return n; }

      void perform_op(double *x_in, double* y_out)
      {
	 Map<VectorXd> x(x_in, n);
	 Map<VectorXd> y(y_out, n);
	 unsigned int actual_block_size;

	 verbose && STDOUT << timestamp() << " Matrix operation " << nops << std::endl;
	 verbose && STDOUT << timestamp() << "   Reading block 0" << std::endl;
	 actual_block_size = stop[0] - start[0] + 1;
	 dat.read_snp_block(start[0], stop[0], false, false);
	 y.noalias() = dat.X.leftCols(actual_block_size) *
	    (dat.X.leftCols(actual_block_size).transpose() * x);

	 for(unsigned int k = 1 ; k < nblocks ; k++)
	 {
	    verbose && STDOUT << timestamp() << "   Reading block " << k << std::endl;
	    actual_block_size = stop[k] - start[k] + 1;
	    dat.read_snp_block(start[k], stop[k], false, false);
	    y.noalias() = y + dat.X.leftCols(actual_block_size) *
	       (dat.X.leftCols(actual_block_size).transpose() * x);
	 }

	 nops++;
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
      verbose && STDOUT << timestamp() << " SVD begin" << std::endl;

      JacobiSVD<MatrixXd> svd(B, ComputeThinU | ComputeThinV);
      U = svd.matrixU();
      MatrixXd V = svd.matrixV();
      d = svd.singularValues().array().pow(2);

      verbose && STDOUT << timestamp() << " SVD done" << std::endl;
   }
   else if(method == METHOD_EIGEN)
   {
      verbose && STDOUT << timestamp() << " Eigen-decomposition begin" << std::endl;

      MatrixXd BBT = B * B.transpose();

      verbose && STDOUT << timestamp() << " dim(BBT): " << dim(BBT) << std::endl;

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
      verbose && STDOUT << timestamp() << " Sampling" << std::endl;

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

   verbose && STDOUT << timestamp() << " Median Euclidean distance: "
      << med << std::endl;

   return med;
}

MatrixXd rbf_kernel(MatrixXd& X, const double sigma, bool rbf_center,
   bool verbose)
{
   unsigned int n = X.rows();
   VectorXd norms = X.array().square().rowwise().sum();
   VectorXd ones = VectorXd::Ones(n);
   MatrixXd R = norms * ones.transpose();
   MatrixXd D = R + R.transpose() - 2 * X * X.transpose();
   D = D.array() / (-1 * sigma * sigma);
   MatrixXd K = D.array().exp();

   if(rbf_center)
   {
      verbose && STDOUT << timestamp() << " Centering RBF kernel" << std::endl;

      MatrixXd M = ones * ones.transpose() / n;
      MatrixXd I = ones.asDiagonal();
      K = (I - M) * K * (I - M);
   }
   return K;
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

void RandomPCA::pca_fast(MatrixXd& X, unsigned int block_size, int method, bool transpose,
   unsigned int ndim, unsigned int nextra, unsigned int maxiter, double tol,
   long seed, int kernel, double sigma, bool rbf_center,
   unsigned int rbf_sample, bool save_kernel, bool do_orth, bool do_loadings,
   int mem)
{
   unsigned int N, p;

   if(transpose)
   {
      if(stand_method_x != STANDARDIZE_NONE)
	  X_meansd = standardize_transpose(X, stand_method_x, verbose);
      N = X.cols();
      p = X.cols();
   }
   else
   {
      if(stand_method_x != STANDARDIZE_NONE)
	 X_meansd = standardize(X, stand_method_x, verbose);
      N = X.rows();
      p = X.cols();
   }

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
      d = eigs.eigenvalues().array() / div;
   }
   else
   {
      std::cerr << "Spectra eigen-decomposition was not successful, status: " 
	 << eigs.info() << std::endl;
   }
}

void RandomPCA::pca_fast(Data& dat, unsigned int block_size, int method, bool transpose,
   unsigned int ndim, unsigned int nextra, unsigned int maxiter, double tol,
   long seed, int kernel, double sigma, bool rbf_center,
   unsigned int rbf_sample, bool save_kernel, bool do_orth, bool do_loadings,
   int mem)
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
      d = eigs.eigenvalues().array() / div;
   }
   else
   {
      std::cerr << "Spectra eigen-decomposition was not successful, status: " 
	 << eigs.info() << std::endl;
   }
}

// stub for now
void RandomPCA::pca(Data& dat, int method, bool transpose,
   unsigned int ndim, unsigned int nextra, unsigned int maxiter, double tol,
   long seed, int kernel, double sigma, bool rbf_center,
   unsigned int rbf_sample, bool save_kernel, bool do_orth, bool do_loadings,
   int mem)
{
   //Spectra::DenseGenMatProd<double> op(dat.X);
   //Spectra::GenEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseGenMatProd<double> > eigs(&op, ndim, 2 * ndim + 1);
   //eigs.init();
   //int nconv = eigs.compute();

   //if(eigs.info() == Spectra::SUCCESSFUL)
   //{
   //   U = eigs.eigenvectors();
   //   d = eigs.eigenvalues().array().sqrt();
   //}
}

void RandomPCA::pca(MatrixXd &X, int method, bool transpose,
   unsigned int ndim, unsigned int nextra, unsigned int maxiter, double tol,
   long seed, int kernel, double sigma, bool rbf_center,
   unsigned int rbf_sample, bool save_kernel, bool do_orth, bool do_loadings,
   int mem)
{
   unsigned int N, p;

   if(kernel != KERNEL_LINEAR)
   {
      transpose = false;

      verbose && STDOUT << timestamp()
	 << " Kernel not linear, can't transpose" << std::endl;
   }

   verbose && STDOUT << timestamp() << " Transpose: " 
      << (transpose ? "yes" : "no") << std::endl;

   if(transpose)
   {
      if(stand_method_x != STANDARDIZE_NONE)
	  X_meansd = standardize_transpose(X, stand_method_x, verbose);
      N = X.cols();
      p = X.rows();
   }
   else
   {
      if(stand_method_x != STANDARDIZE_NONE)
	 X_meansd = standardize(X, stand_method_x, verbose);
      N = X.rows();
      p = X.cols();
   }

   unsigned int total_dim = ndim + nextra;
   MatrixXd R = make_gaussian(X.cols(), total_dim, seed);
   MatrixXd Y = X * R;

   verbose && STDOUT << timestamp() << " dim(Y): " << dim(Y) << std::endl;

   normalize(Y);
   MatrixXd Yn;

   verbose && STDOUT << timestamp() << " dim(X): " << dim(X) << std::endl;

   double div = 1;
   if(divisor == DIVISOR_N1)
      div = N - 1;
   else if(divisor == DIVISOR_P)
      div = p;

   MatrixXd K; 
   if(mem == HIGHMEM && kernel == KERNEL_RBF)
   {
      if(sigma == 0)
      {
	 unsigned int med_samples = fminl(rbf_sample, N);
      	 double med = median_dist(X, med_samples, seed, verbose);
      	 sigma = sqrt(med);
      }

      verbose && STDOUT << timestamp() << " Using RBF kernel with sigma="
	 << sigma << std::endl;

      K.noalias() = rbf_kernel(X, sigma, rbf_center, verbose);
   }
   else if(mem == HIGHMEM)
   {
      verbose && STDOUT << timestamp() << " Using linear kernel" << std::endl;

      K.noalias() = X * X.transpose() / div;
   }

   //trace = K.diagonal().array().sum() / (N - 1);
   if(mem == LOWMEM)
      trace = X.array().square().sum();
   else
   {
      verbose && STDOUT << timestamp() << " dim(K): " << dim(K) << std::endl;
      trace = K.diagonal().array().sum();
   }
   
   trace /= div;

   verbose && STDOUT << timestamp() << " Trace(K): " << trace 
      << " (N: " << N << ")" << std::endl;

   if(mem == HIGHMEM && save_kernel)
   {
      verbose && STDOUT << timestamp() << " saving K" << std::endl;

      //save_text(K, "kernel.txt");
   }

   MatrixXd Xy(X.cols(), Y.cols());

   for(unsigned int iter = 0 ; iter < maxiter ; iter++)
   {
      verbose && STDOUT << timestamp() << " iter " << iter;

      if(mem == LOWMEM)
      {
	 Xy.noalias() = X.transpose() * Y;
	 Yn.noalias() = X * Xy;
      }
      else
	 Yn.noalias() = K * Y;
      if(do_orth)
      {
	 verbose && STDOUT << " (orthogonalising)";

	 ColPivHouseholderQR<MatrixXd> qr(Yn);
	 MatrixXd I = MatrixXd::Identity(Yn.rows(), Yn.cols());
	 Yn = qr.householderQ() * I;
	 Yn.conservativeResize(NoChange, Yn.cols());
      }
      else
	 normalize(Yn);

      double diff = (Y - Yn).array().square().sum() / Y.size(); 

      verbose && STDOUT << " " << diff << std::endl;

      Y.noalias() = Yn;
      if(diff < tol)
	 break;
   }

   verbose && STDOUT << timestamp() << " QR begin" << std::endl;

   ColPivHouseholderQR<MatrixXd> qr(Y);
   MatrixXd Q = MatrixXd::Identity(Y.rows(), Y.cols());
   Q = qr.householderQ() * Q;
   Q.conservativeResize(NoChange, Y.cols());

   verbose && STDOUT << timestamp() << " dim(Q): " << dim(Q) << std::endl;
   verbose && STDOUT << timestamp() << " QR done" << std::endl;

   MatrixXd B = Q.transpose() * X;

   verbose && STDOUT << timestamp() << " dim(B): " << dim(B) << std::endl;

   MatrixXd Et;
   pca_small(B, method, Et, d, verbose);

   verbose && STDOUT << timestamp() << " dim(Et): " << dim(Et) << std::endl;

   d = d.array() / div;

   // TODO: unlike Scholkopf, when we do kernel PCA we don't divide the
   // eigenvectors by the sqrt of each eigenvalue

   if(transpose)
   {
      V.noalias() = Q * Et;
      // We divide Px by sqrt(N - 1) since X has not been divided
      // by it (but B has)
      Px.noalias() = X.transpose() * V;
      VectorXd s;
      s = 1 / (d.array().sqrt() * sqrt(div));

      MatrixXd Dinv = s.asDiagonal();
      U = Px * Dinv;
   }
   else
   {
      // Px = U D = X V
      U.noalias() = Q * Et;
      VectorXd dtmp = d.array().sqrt();
      Px.noalias() = U * dtmp.asDiagonal();
      if(do_loadings)
      {
	 VectorXd s;
	 s = 1 / (d.array().sqrt() * sqrt(div));
	 MatrixXd Dinv = s.asDiagonal();
	 V = X.transpose() * U * Dinv;
      }
   }

   Px.conservativeResize(NoChange, ndim);
   U.conservativeResize(NoChange, ndim);
   V.conservativeResize(NoChange, ndim);
   d.conservativeResize(ndim);
   pve = d.array() / trace;
}

// ZCA of genotypes
//void RandomPCA::zca_whiten(bool transpose)
//{
//   verbose && STDOUT << timestamp() << " Whitening begin" << std::endl;
//   VectorXd s = 1 / d.array();
//   MatrixXd Dinv = s.asDiagonal();
//
//   if(transpose)
//      W.noalias() = U * Dinv * U.transpose() * X.transpose();
//   else
//      W.noalias() = U * Dinv * U.transpose() * X;
//   verbose && STDOUT << timestamp() << " Whitening done (" << dim(W) << ")" << std::endl;
//}

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
	 
	 // At least up to Eigen 3.2.1, the QR produces full Q matrices, i.e.,
	 // U.rows() x U.rows() instead of U.rows() x ndim, so we use
	 // an indirect to get back the original dimensions.
	 // We purposefully don't use column pivoting QR, as this could break
	 // ordering of columns of U wrt columns of V.
	 //HouseholderQR<MatrixXd> qrU(U.leftCols(j + 1)), qrV(V.leftCols(j + 1));
	 //MatrixXd Iu = MatrixXd::Identity(U.rows(), j + 1);
	 //MatrixXd Iv = MatrixXd::Identity(V.rows(), j + 1);
	 //U.leftCols(j + 1) = qrU.householderQ() * Iu;
	 //V.leftCols(j + 1) = qrV.householderQ() * Iv;
	 //U.conservativeResize(NoChange, Iu.cols());
	 //V.conservativeResize(NoChange, Iv.cols());
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

	 if(iter > 0 && (v_old.array() - v.array()).abs().maxCoeff() < tol
	       && (u_old.array() - u.array()).abs().maxCoeff() < tol)
	 {
	    verbose && STDOUT << timestamp() << " dim " << j << " finished in "
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

      verbose && STDOUT << timestamp() << " U_" << j 
	 << " non-zeros: " << nzu << ", V_" << j
	 << " non-zeros: " << nzv << std::endl;

      // Use X and Y, not X2 and Y2
      d[j] = (X * U.col(j)).transpose() * (Y * V.col(j)); 
      verbose && STDOUT << timestamp() << " d[" << j << "]: "
	 << d[j] << std::endl;
   }
}

void scca_highmem(MatrixXd& X, MatrixXd &Y, MatrixXd& U, MatrixXd& V,
   VectorXd& d, double lambda1, double lambda2,
   unsigned int maxiter, double tol, bool verbose)
{
   verbose && STDOUT << timestamp() << " Begin computing X^T Y" << std::endl;

   MatrixXd XY = X.transpose() * Y;

   verbose && STDOUT << timestamp() << " End computing X^T Y" << std::endl;

   MatrixXd XYj;
   VectorXd u, v, u_old, v_old;

   for(unsigned int j = 0 ; j < U.cols() ; j++)
   {
      verbose && STDOUT << timestamp() << " dim " << j << std::endl;

      if(j == 0)
	 XYj = XY;
      else
      {
	 XYj = XYj - d[j-1] * U.col(j-1) * V.col(j - 1).transpose();

	 //HouseholderQR<MatrixXd> qrU(U.leftCols(j + 1)), qrV(V.leftCols(j + 1));
	 //MatrixXd Iu = MatrixXd::Identity(U.rows(), j + 1);
	 //MatrixXd Iv = MatrixXd::Identity(V.rows(), j + 1);
	 //U.leftCols(j + 1) = qrU.householderQ() * Iu;
	 //V.leftCols(j + 1) = qrV.householderQ() * Iv;
	 //U.conservativeResize(NoChange, Iu.cols());
	 //V.conservativeResize(NoChange, Iv.cols());
      }
	 
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

	 if(iter > 0 && (v_old.array() - v.array()).abs().maxCoeff() < tol
	       && (u_old.array() - u.array()).abs().maxCoeff() < tol)
	 {
	    verbose && STDOUT << timestamp() << " dim " << j << " finished in "
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

      verbose && STDOUT << timestamp() << " U_" << j 
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
   {
      X_meansd = standardize(X, stand_method_x);
      Y_meansd = standardize(Y, stand_method_x);
   }

   verbose && STDOUT << timestamp() << " dim(X): " << dim(X) << std::endl;
   verbose && STDOUT << timestamp() << " dim(Y): " << dim(Y) << std::endl;
   verbose && STDOUT << timestamp() << " lambda1: " << lambda1 
      << " lambda2: " << lambda2 << std::endl;

   unsigned int p = X.cols();

   //V = make_gaussian(k, ndim, seed);
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
   MatrixXd covYinv = svd.matrixV() * d2.asDiagonal() * svd.matrixV().transpose();
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

   std::cout << "ucca online mode, N=" << n << " p=" << p << std::endl;

   JacobiSVD<MatrixXd> svd(data.Y, ComputeThinU | ComputeThinV);
   VectorXd d = svd.singularValues() / sqrt(n - 1);
   VectorXd d2 = d.array().pow(-2);
   MatrixXd V = svd.matrixV();
   MatrixXd covYinv = svd.matrixV() * d2.asDiagonal() * svd.matrixV().transpose();
   ArrayXd r2 = ArrayXd(p);

   for(unsigned int j = 0 ; j < p ; j++)
   {
      data.read_snp_block(j, j, false, true);
      if(stand_method_x != STANDARDIZE_NONE)
	 X_meansd = standardize(data.X, stand_method_x);
      varx = var(data.X.col(0));
      covXY = cov(data.X.col(0), data.Y);

      // take absolute value to prevent numerical issues with negative numbers
      // close to zero
      r2(j) = fabs(1.0 / varx * covXY * covYinv * covXY.transpose());
   }

   res = wilks(r2, n, data.Y.cols());
}

