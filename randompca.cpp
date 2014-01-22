
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
inline void normalize(MatrixXd& X)
{
   unsigned int p = X.cols();
   for(unsigned int j = 0 ; j < p ; j++)
   {
      double s = 1 / sqrt(X.col(j).array().pow(2).sum());
      X.col(j) = X.col(j).array() * s;
   }
}

void RandomPCA::pca(MatrixXd &X, int method, bool transpose,
   unsigned int ndim, unsigned int nextra, unsigned int maxiter, double tol)
{
   M = standardize(X, true, stand_method);
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
   MatrixXd Y = M * R;
   std::cout << timestamp() << " dim(Y): " << dim(Y) << std::endl;
   normalize(Y);
   MatrixXd Yn;

   std::cout << timestamp() << " dim(M): " << dim(M) << std::endl;
   MatrixXd MMT = M * M.transpose();
   std::cout << timestamp() << " dim(MMT): " << dim(MMT) << std::endl;

   for(unsigned int iter = 0 ; iter < maxiter ; iter++)
   {
      std::cout << timestamp() << " iter " << iter << " ";
      Yn.noalias() = MMT * Y;
      normalize(Yn);
      double diff =  (Y -  Yn).array().square().sum() / Y.size(); 
      std::cout << diff << std::endl;
      Y.noalias() = Yn;
      if(diff < tol)
	 break;
   }

   std::cout << timestamp() << " QR begin" << std::endl;
   ColPivHouseholderQR<MatrixXd> qr(Y);
   MatrixXd Q = MatrixXd::Identity(Y.rows(), Y.cols());
   Q = qr.householderQ() * Q;
   Q.conservativeResize(NoChange, Y.cols());
   std::cout << timestamp() << " QR done" << std::endl;

   MatrixXd B = Q.transpose() * M;
   std::cout << timestamp() << " dim(B): " << dim(B) << std::endl;

   if(method == METHOD_SVD)
   {
      std::cout << timestamp() << " SVD begin" << std::endl;
      JacobiSVD<MatrixXd> svd(B, ComputeThinU | ComputeThinV);
      U = Q * svd.matrixU();
      V = svd.matrixV();
      d = svd.singularValues() / sqrt(X.rows() - 1);
      std::cout << timestamp() << " SVD done" << std::endl;
   }
   else if(method == METHOD_EIGEN)
   {
      std::cout << timestamp() << " Eigen-decomposition begin" << std::endl;
      MatrixXd BBT = B * B.transpose();
      std::cout << timestamp() << " dim(BBT): " << dim(BBT) << std::endl;
      SelfAdjointEigenSolver<MatrixXd> eig(BBT);

      // The eigenvalues are sorted in *increasing* order, we need decreasing
      // order
      VectorXd eval = eig.eigenvalues();
      MatrixXd evec = eig.eigenvectors();
      d.resize(eval.size());
      U.resize(BBT.rows(), BBT.rows());

      unsigned int k = 0;
      for(unsigned int i = d.size() - 1 ; i != 0 ; --i)
      {
	 // we get eigenvalues, which are the squared singular values
	 d(k) = sqrt(eval(i)) / sqrt(X.rows() - 1);
	 U.col(k) = evec.col(i);
	 k++;
      }

      U = Q * U;

      std::cout << timestamp() << " Eigen-decomposition done" << std::endl;
   }

   //// TODO: untested
   //if(transpose)
    //  P.noalias() = (U.leftCols(ndim).transpose() * M).transpose();
   //else
      //P.noalias() = M * V.leftCols(ndim);

   P.noalias() = U * d.asDiagonal();
}

// ZCA
void RandomPCA::zca_whiten()
{
   std::cout << timestamp() << " Whitening begin" << std::endl;
   VectorXd s = 1 / d.array();
   MatrixXd Dinv = s.asDiagonal();
   W.noalias() = U * Dinv * U.transpose() * M;
   std::cout << timestamp() << " Whitening done (" << dim(W) << ")" << std::endl;
}

