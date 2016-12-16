
#include "svdtall.h"

void SVDTall::perform_op(double *x_in, double* y_out)
{
   Map<VectorXd> x(x_in, n);
   Map<VectorXd> y(y_out, n);
   verbose && STDOUT << timestamp() << "Matrix operation "
      << nops << std::endl;
   y.noalias() = mat * (mat.transpose() * x);
   nops++;
}

SVDTallOnline::~SVDTallOnline()
{
   delete[] start;
   delete[] stop;
}

// y = X X' * x
void SVDTallOnline::perform_op(double *x_in, double* y_out)
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
MatrixXd SVDTallOnline::perform_op_mat(const MatrixXd x)
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
// Note: size of x must be number of samples, size y must be number of SNPs
void SVDTallOnline::crossprod(double *x_in, double *y_out)
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
MatrixXd SVDTallOnline::crossprod2(const MatrixXd& x)
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
// size of y must be the number of samples
void SVDTallOnline::prod(double *x_in, double *y_out)
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
MatrixXd SVDTallOnline::perform_op_multi(const MatrixXd& x)
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
MatrixXd SVDTallOnline::prod2(const MatrixXd& x)
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
MatrixXd SVDTallOnline::prod3(const MatrixXd& x)
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

