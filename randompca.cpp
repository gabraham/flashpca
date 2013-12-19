
#include "randompca.hpp"
#include "util.hpp"

// TODO: old redsvd code, replace
void sampleTwoGaussian(double& f1, double& f2){
  double v1 = (double)(rand() + 1.f) / ((double)RAND_MAX+2.f);
  double v2 = (double)(rand() + 1.f) / ((double)RAND_MAX+2.f);
  double len = sqrt(-2.f * log(v1));
  f1 = len * cos(2.f * M_PI * v2);
  f2 = len * sin(2.f * M_PI * v2);
}

// TODO: old redsvd code, replace
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
void normalize(MatrixXd& X)
{
   unsigned int p = X.cols();
   for(unsigned int j = 0 ; j < p ; j++)
   {
      double s = 1 / sqrt(X.col(j).array().pow(2).sum());
      X.col(j) = X.col(j).array() * s;
   }
}

MatrixXd RandomPCA::pca(MatrixXd X, bool transpose,
   unsigned int ndim, unsigned int nextra, unsigned int maxiter)
{
   M = standardize(X);
   std::cout << timestamp() << " Transpose: " 
      << (transpose ? "yes" : "no") << std::endl;
   if(transpose)
   {
      std::cout << timestamp() << " transpose start" << std::endl;
      M.transposeInPlace();
      std::cout << timestamp() << " transpose done" << std::endl;
   }

   unsigned int total_dim = ndim + nextra;
   MatrixXd R = make_gaussian(M.cols(), total_dim);
   save_text("R.txt", R);
   MatrixXd Y = M * R;
   std::cout << timestamp() << " dim(Y): " << dim(Y) << std::endl;
   normalize(Y);
   MatrixXd Yn;

   for(unsigned int iter = 0 ; iter < maxiter ; iter++)
   {
      std::cout << timestamp() << " iter " << iter << " ";
      Yn = M * (M.transpose() * Y);
      normalize(Yn);
      double diff =  (Y -  Yn).array().square().sum() / Y.size(); 
      std::cout << diff << std::endl;
      Y = Yn;
      if(diff < tol)
	 break;
   }

   std::cout << timestamp() << " QR begin" << std::endl;
   ColPivHouseholderQR<MatrixXd> qr(Y);
   MatrixXd Q = MatrixXd::Identity(Y.rows(), Y.cols());
   Q = qr.householderQ() * Q;
   Q.conservativeResize(NoChange, Y.cols());
   MatrixXd B = Q.transpose() * M;
   std::cout << timestamp() << " QR done" << std::endl;

   std::cout << timestamp() << " SVD begin" << std::endl;
   std::cout << timestamp() << " dim(B): " << dim(B) << std::endl;
   JacobiSVD<MatrixXd> svd(B, ComputeThinU | ComputeThinV);
   std::cout << timestamp() << " SVD done" << std::endl;

   U = Q * svd.matrixU();
   V = svd.matrixV();
   d = svd.singularValues();

   std::cout << timestamp() << " dim(U): " << dim(U) << std::endl;
   std::cout << timestamp() << " dim(V): " << dim(V) << std::endl;

   //// TODO: untested
   if(transpose)
      P = (U.leftCols(ndim).transpose() * M).transpose();
   else
      P = M * V.leftCols(ndim);

   return P;
}

// ZCA
MatrixXd RandomPCA::zca_whiten()
{
   std::cout << timestamp() << " Whitening begin" << std::endl;
   VectorXd s = 1 / d.array();
   MatrixXd D = s.asDiagonal();
   W = U * D * U.transpose() * M;
   std::cout << timestamp() << " Whitening done (" << dim(W) << ")" << std::endl;
   return W;
}

