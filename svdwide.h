
#pragma once

#include <Eigen/Core>
#include <Eigen/Eigen>

#include "data.h"

using namespace Eigen;

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
      };

      inline unsigned int rows() const { return n; }
      inline unsigned int cols() const { return n; }

      void perform_op(double *x_in, double* y_out);
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

      ~SVDWideOnline();

      inline unsigned int rows() const { return n; }
      inline unsigned int cols() const { return n; }

      // y = X X' * x
      void perform_op(double *x_in, double* y_out);

      // y = X X' * x
      MatrixXd perform_op_mat(const MatrixXd x);

      // Like R crossprod(): y = X' * x
      // Note: size of x must be number of samples, size y must be number of SNPs
      void crossprod(double *x_in, double *y_out);

      // Like R crossprod(): y = X' * x
      // Note: size of x must be number of samples, size y must number of SNPs
      MatrixXd crossprod2(const MatrixXd& x);

      // Like y = X %*% x
      // Note: size of x must be number of SNPs,
      // size of y must be the number of samples
      void prod(double *x_in, double *y_out);

      // return Y = X * X' * x where x is a matrix (despite x being lower case)
      MatrixXd perform_op_multi(const MatrixXd& x);

      // Like Y = x' * X, where X is genotypes, x is a matrix
      MatrixXd prod2(const MatrixXd& x);

      // Like Y = X * x, where X is genotypes, x is a matrix
      MatrixXd prod3(const MatrixXd& x);
};

