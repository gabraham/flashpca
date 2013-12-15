
#include "randompca.hpp"
#include "util.hpp"

void sampleTwoGaussian(double& f1, double& f2){
  double v1 = (double)(rand() + 1.f) / ((double)RAND_MAX+2.f);
  double v2 = (double)(rand() + 1.f) / ((double)RAND_MAX+2.f);
  double len = sqrt(-2.f * log(v1));
  f1 = len * cos(2.f * M_PI * v2);
  f2 = len * sin(2.f * M_PI * v2);
}

// redsvd, replace
void sampleGaussianMat(MatrixXd& mat){
  for (int i = 0; i < mat.rows(); ++i){
    int j = 0;
    for ( ; j+1 < mat.cols(); j += 2){
      double f1, f2;
      sampleTwoGaussian(f1, f2);
      mat(i,j  ) = f1;
      mat(i,j+1) = f2;
    }
    for (; j < mat.cols(); j ++){
      double f1, f2;
      sampleTwoGaussian(f1, f2);
      mat(i, j)  = f1;
    }
  }
} 

MatrixXd make_gaussian(unsigned int rows, unsigned int cols)
{
   //boost::mt19937 rng;    // The uniform pseudo-random algorithm
   //boost::normal_distribution<double> norm;  // The gaussian combinator
   //boost::variate_generator<boost::mt19937&,boost::normal_distribution<double> >
   //    randN; // The 0-mean unit-variance normal generator

   MatrixXd G(rows, cols);
   sampleGaussianMat(G);
   return G;
}

// normalize each column of X to unit l2 norm
void normalize(MatrixXd X)
{
   unsigned int p = X.cols();
   for(unsigned int j = 0 ; j < p ; j++)
   {
      double s = 1 / sqrt(X.col(j).array().pow(2).sum());
      X.col(j) = X.col(j).array() * s;
   }
}

MatrixXd random_pca(MatrixXd X, bool transpose,
   unsigned int ndim, unsigned int nextra, unsigned int maxiter)
{
   MatrixXd M = standardize(X);
   std::cout << "transpose: " << transpose << std::endl;
   if(transpose)
      M.transposeInPlace();

   unsigned int total_dim = ndim + nextra;
   MatrixXd R = make_gaussian(M.cols(), total_dim);
   std::cout << dim(X) << std::endl;
   std::cout << dim(M) << std::endl;
   std::cout << dim(R) << std::endl;
   MatrixXd Y = M * R;

   for(unsigned int iter = 0 ; iter < maxiter ; iter++)
   {
      std::cout << "iter " << iter << std::endl;
      Y = M * (M.transpose() * Y);
      normalize(Y);
   }

   std::cout << " M: " << dim(M) << std::endl;
   std::cout << "QR ... ";
   std::cout << " Y: " << dim(Y);

   ColPivHouseholderQR<MatrixXd> qr(Y);
   MatrixXd Q = MatrixXd::Identity(Y.rows(), Y.cols());
   Q = qr.householderQ() * Q;

   Q.conservativeResize(NoChange, Y.cols());
   std::cout << " Q: " << dim(Q);
   MatrixXd B = Q.transpose() * M;
   std::cout << " done" << std::endl;

   std::cout << "SVD ... " << "B:" << dim(B) << std::endl;
   JacobiSVD<MatrixXd> svd(B, ComputeThinU | ComputeThinV);
 
   //// TODO: untested
   if(transpose)
   {
      return (svd.matrixU().leftCols(ndim).transpose() * M).transpose();
   }
   else
   {
      return M * svd.matrixV().leftCols(ndim);
   }
}


