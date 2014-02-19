
#include "randompca.hpp"
#include "util.hpp"

MatrixXd make_gaussian(unsigned int rows, unsigned int cols, long seed)
{
   boost::mt19937 rng;
   rng.seed(seed);
   boost::normal_distribution<double> nrm;
   boost::variate_generator<boost::mt19937&,
      boost::normal_distribution<double> > rand(rng, nrm);

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

void pca_small(MatrixXd &B, int method, MatrixXd& U, VectorXd &d)
{
   if(method == METHOD_SVD)
   {
      std::cout << timestamp() << " SVD begin" << std::endl;
      JacobiSVD<MatrixXd> svd(B, ComputeThinU | ComputeThinV);
      U = svd.matrixU();
      MatrixXd V = svd.matrixV();
      d = svd.singularValues();
      std::cout << timestamp() << " SVD done" << std::endl;
   }
   else if(method == METHOD_EIGEN)
   {
      std::cout << timestamp() << " Eigen-decomposition begin" << std::endl;
      MatrixXd BBT = B * B.transpose();
      std::cout << timestamp() << " dim(BBT): " << dim(BBT) << std::endl;
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
	 d(k) = sqrt(eval(i));
	 U.col(k) = evec.col(i);
	 k++;
      }
   }
}

// Compute median of pairwise distances on sample of size n from the matrix X
double median_dist(MatrixXd& X, unsigned int n, long seed)
{
   boost::mt19937 rng;
   rng.seed(seed);
   boost::uniform_int_distribution<> randi(1, X.rows()); 

   return 0;
}

MatrixXd rbf_kernel(MatrixXd& X, const double sigma)
{
   MatrixXd K = MatrixXd::Zero(X.rows(), X.rows());
   unsigned int n = X.rows();
   int nt = nbThreads();

   setNbThreads(0);

   #pragma omp parallel for shared(K, X)
   for(unsigned int i = 0 ; i < n ; i++)
   {
      K(i, i) = 1; // exp(0)
      for(unsigned int j = 1 ; j < i ; j++)
      {
	 double z = (X.row(i).array() - X.row(j).array()).square().sum();
	 K(i, j) = K(j, i) = exp(-z / sigma);
      }
   }

   setNbThreads(nt);

   VectorXd m = VectorXd::Ones(n);
   MatrixXd M = m * m.transpose() / n;
   MatrixXd I = m.asDiagonal();
   K = (I - M) * K * (I - M);
   return K;
}

//MatrixXd rbf_kernel(MatrixXd& X)
//{
//   MatrixXd K = MatrixXd::Zero(X.rows(), X.rows());
//   unsigned int n = X.rows();
//   for(unsigned int i = 0 ; i < n ; i++)
//   {
//      K(i, i) = 1; // exp(0)
//      for(unsigned int j = 1 ; j < i ; j++)
//      {
//	 double z = (X.row(i).array() - X.row(j).array()).square().sum();
//	 K(i, j) = K(j, i) = z;
//      }
//   }
//
//   K.triangularView(Upper).data();
//
//   VectorXd m = VectorXd::Ones(n);
//   MatrixXd M = m * m.transpose() / n;
//   MatrixXd I = m.asDiagonal();
//   K = (I - M) * K * (I - M);
//   return K;
//}

void RandomPCA::pca(MatrixXd &X, int method, bool transpose,
   unsigned int ndim, unsigned int nextra, unsigned int maxiter, double tol,
   long seed, int kernel, double sigma)
{
   unsigned int N;

   if(kernel != KERNEL_LINEAR)
   {
      transpose = false;
      std::cout << timestamp()
	 << " Kernel not linear, can't transpose" << std::endl;
   }

   std::cout << timestamp() << " Transpose: " 
      << (transpose ? "yes" : "no") << std::endl;

   if(transpose)
   {
      M = standardize_transpose(X, true, stand_method);
      N = X.cols();
   }
   else
   {
      M = standardize(X, true, stand_method);
      N = X.rows();
   }

   unsigned int total_dim = ndim + nextra;
   MatrixXd R = make_gaussian(M.cols(), total_dim, seed);
   MatrixXd Y = M * R;
   std::cout << timestamp() << " dim(Y): " << dim(Y) << std::endl;
   normalize(Y);
   MatrixXd Yn;

   std::cout << timestamp() << " dim(M): " << dim(M) << std::endl;
   MatrixXd K; 
   if(kernel == KERNEL_RBF)
   {
      std::cout << timestamp() << " Using RBF kernel with sigma="
	 << sigma << std::endl;
      K.noalias() = rbf_kernel(M, sigma);
   }
   else
   {
      std::cout << timestamp() << " Using linear kernel" << std::endl;
      K.noalias() = M * M.transpose();
   }

   //save_text("K.txt", K);
   std::cout << timestamp() << " dim(K): " << dim(K) << std::endl;

   for(unsigned int iter = 0 ; iter < maxiter ; iter++)
   {
      std::cout << timestamp() << " iter " << iter << " ";
      Yn.noalias() = K * Y;
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
   std::cout << timestamp() << " dim(Q): " << dim(Q) << std::endl;
   std::cout << timestamp() << " QR done" << std::endl;

   MatrixXd B = Q.transpose() * M;
   B = B.array() / sqrt(N - 1);
   std::cout << timestamp() << " dim(B): " << dim(B) << std::endl;

   MatrixXd Ut;
   pca_small(B, method, Ut, d);
   std::cout << timestamp() << " dim(Ut): " << dim(Ut) << std::endl;
   U.noalias() = Q * Ut;
   std::cout << timestamp() << " dim(U): " << dim(U) << std::endl;

   if(transpose)
   {
      // P = X V
      // since U is actually V (eigenvectors of X^T X), since data is transposed.
      // We divide by sqrt(N - 1) since X has not been divided by it (but B
      // has)
      P.noalias() = M.transpose() * U;
      double z = 1.0 / sqrt(N - 1);
      P = P.array() * z;
   }
   else
   {
      // P = U D
      P.noalias() = U * d.asDiagonal();
   }

   P.conservativeResize(NoChange, ndim);
   d.conservativeResize(ndim);
}

// ZCA of genotypes
void RandomPCA::zca_whiten(bool transpose)
{
   std::cout << timestamp() << " Whitening begin" << std::endl;
   VectorXd s = 1 / d.array();
   MatrixXd Dinv = s.asDiagonal();

   if(transpose)
      W.noalias() = U * Dinv * U.transpose() * M.transpose();
   else
      W.noalias() = U * Dinv * U.transpose() * M;
   std::cout << timestamp() << " Whitening done (" << dim(W) << ")" << std::endl;
}

