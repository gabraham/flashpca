
#include "randompca.hpp"
#include "util.hpp"

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

void pca_small(MatrixXd &B, int method, MatrixXd& U, MatrixXd& V, VectorXd &d)
{
   if(method == METHOD_SVD)
   {
      std::cout << timestamp() << " SVD begin" << std::endl;
      JacobiSVD<MatrixXd> svd(B, ComputeThinU | ComputeThinV);
      U = svd.matrixU();
      MatrixXd V = svd.matrixV();
      d = svd.singularValues().array().pow(2);
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

      unsigned int k = 0, s = d.size();
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
// 
// Todo: use Boost accummulators for quantiles
// http://boost-sandbox.sourceforge.net/libs/accumulators/doc/html/accumulators/user_s_guide/the_statistical_accumulators_library.html#accumulators.user_s_guide.the_statistical_accumulators_library.p_square_quantile
double median_dist(MatrixXd& X, unsigned int n, long seed)
{
   boost::random::mt19937 rng;
   rng.seed(seed);
   boost::random::uniform_real_distribution<> dist(0, 1);
   double prop = (double)X.rows() / n;

   std::cout << timestamp() << 
      " Computing median Euclidean distance (" << n << " samples)" <<
      std::endl;

   MatrixXd X2(n, X.cols());
   if(n < X.rows()) 
   {
      std::cout << timestamp() << " Sampling" << std::endl;
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
   MatrixXd M = norms * ones.transpose();
   MatrixXd D = M + M.transpose() - 2 * X2 * X2.transpose();

   unsigned int m = D.size();
   double *d = D.data();
   std::sort(d, d + m); 
   double med;

   if(m % 2 == 0)
      med = (d[m / 2 - 1] + d[m / 2]) / 2;
   else
      med = d[m / 2];

   std::cout << timestamp() << " Median Euclidean distance: "
      << med << std::endl;

   return med;
}

MatrixXd rbf_kernel(MatrixXd& X, const double sigma, bool rbf_center)
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
      std::cout << timestamp() << " Centering RBF kernel" << std::endl;
      MatrixXd M = ones * ones.transpose() / n;
      MatrixXd I = ones.asDiagonal();
      K = (I - M) * K * (I - M);
   }
   return K;
}

void RandomPCA::pca(MatrixXd &X, int method, bool transpose,
   unsigned int ndim, unsigned int nextra, unsigned int maxiter, double tol,
   long seed, int kernel, double sigma, bool rbf_center,
   unsigned int rbf_sample, bool save_kernel)
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
      if(stand_method != STANDARDIZE_NONE)
	 M = standardize_transpose(X, stand_method);
      N = X.cols();
   }
   else
   {
      if(stand_method != STANDARDIZE_NONE)
	 M = standardize(X, stand_method);
      N = X.rows();
   }

   unsigned int total_dim = ndim + nextra;
   MatrixXd R = make_gaussian(M.cols(), total_dim, seed);
   MatrixXd Y = M * R; // TODO: is this right for kernel PCA?
   std::cout << timestamp() << " dim(Y): " << dim(Y) << std::endl;
   normalize(Y);
   MatrixXd Yn;

   std::cout << timestamp() << " dim(M): " << dim(M) << std::endl;
   MatrixXd K; 
   if(kernel == KERNEL_RBF)
   {
      if(sigma == 0)
      {
	 unsigned int med_samples = fminl(rbf_sample, N);
      	 double med = median_dist(M, med_samples, seed);
      	 sigma = sqrt(med);
      }
      std::cout << timestamp() << " Using RBF kernel with sigma="
	 << sigma << std::endl;
      K.noalias() = rbf_kernel(M, sigma, rbf_center);
   }
   else
   {
      std::cout << timestamp() << " Using linear kernel" << std::endl;
      K.noalias() = M * M.transpose();
   }

   trace = K.diagonal().array().sum() / (N - 1);
   std::cout << timestamp() << " Trace(K): " << trace 
      << " (N: " << N << ")" << std::endl;

   std::cout << timestamp() << " dim(K): " << dim(K) << std::endl;
   if(save_kernel)
   {
      std::cout << timestamp() << " saving K" << std::endl;
      save_text("kernel.txt", K);
   }

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
   pca_small(B, method, Ut, V, d);
   std::cout << timestamp() << " dim(Ut): " << dim(Ut) << std::endl;
   U.noalias() = Q * Ut;
   std::cout << timestamp() << " dim(U): " << dim(U) << std::endl;


   // TODO: unlike Scholkopf, when we do kernel PCA we don't divide the
   // eigenvectors by the sqrt of each eigenvalue

   if(transpose)
   {
      // P = X V
      // since U is actually V (eigenvectors of X^T X), since data is transposed.
      // We divide by sqrt(N - 1) since X has not been divided by it (but B
      // has)
      Px.noalias() = M.transpose() * U;
      double z = 1.0 / sqrt(N - 1);
      Px = Px.array() * z;
   }
   else
   {
      // P = U D
      Px.noalias() = U * d.asDiagonal();
   }

   Px.conservativeResize(NoChange, ndim);
   d.conservativeResize(ndim);
   pve = d.array() / trace;

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

void RandomPCA::cca(MatrixXd &X, MatrixXd &Y, long seed)//, int method, bool transpose,
//   unsigned int ndim, unsigned int nextra, unsigned int maxiter, double tol,
//   long seed, int kernel, double sigma, bool rbf_center,
//   unsigned int rbf_sample, bool save_kernel)
{
   std::cout << timestamp() << " Begin computing covariance matrices" << std::endl;
   MatrixXd Sx = X.transpose() * X;
   MatrixXd Sy = Y.transpose() * Y;
   MatrixXd Sxy = X.transpose() * Y;

   std::cout << timestamp() << " Covariance done" << std::endl;

   double lambda = 1e-3;

   // TODO: Add ridge regularisation to Sx and Sy
   VectorXd dx = Sx.diagonal();
   VectorXd dy = Sy.diagonal();
   Sx.diagonal() = dx.array() + lambda;
   Sy.diagonal() = dy.array() + lambda;
   
   LLT<MatrixXd> lltX(Sx);
   LLT<MatrixXd> lltY(Sy);
   MatrixXd W1 = lltX.solve(Sxy);
   MatrixXd M = lltY.solve(W1.transpose());
   std::cout << timestamp() << " dim(M): " << dim(M) << std::endl;

   unsigned int maxiter = 10;
   double tol = 1e-6;
   unsigned int ndim = 10, nextra = 10;
   unsigned int N = X.rows();

   // Necessary for getting V. Also, the correlations are the sqrt of the
   // eigenvalues of the cross-covariance matrix
   int method = METHOD_SVD;

   unsigned int total_dim = ndim + nextra;
   unsigned int m = fminl(X.cols(), Y.cols());
   if(total_dim > m)
      total_dim = m;
   MatrixXd R = make_gaussian(M.cols(), total_dim, seed);
   MatrixXd Z = M * R;
   std::cout << timestamp() << " dim(Z): " << dim(Z) << std::endl;
   normalize(Z);
   MatrixXd Zn;
   MatrixXd K = M * M.transpose();

   for(unsigned int iter = 0 ; iter < maxiter ; iter++)
   {
      std::cout << timestamp() << " iter " << iter << " ";
      Zn.noalias() = K * Z;
      normalize(Zn);
      double diff =  (Z -  Zn).array().square().sum() / Z.size(); 
      std::cout << diff << std::endl;
      Z.noalias() = Zn;
      if(diff < tol)
	 break;
   }

   std::cout << timestamp() << " QR begin" << std::endl;
   ColPivHouseholderQR<MatrixXd> qr(Z);
   MatrixXd Q = MatrixXd::Identity(Z.rows(), Z.cols());
   Q = qr.householderQ() * Q;
   Q.conservativeResize(NoChange, Z.cols());
   std::cout << timestamp() << " dim(Q): " << dim(Q) << std::endl;
   std::cout << timestamp() << " QR done" << std::endl;

   MatrixXd B = Q.transpose() * M;
   B = B.array() / sqrt(N - 1);
   std::cout << timestamp() << " dim(B): " << dim(B) << std::endl;

   MatrixXd Ut;
   pca_small(B, method, Ut, V, d);
   std::cout << timestamp() << " dim(Ut): " << dim(Ut) << std::endl;
   U.noalias() = Q * Ut;
   std::cout << timestamp() << " dim(U): " << dim(U) << std::endl;
   std::cout << timestamp() << " dim(V): " << dim(V) << std::endl;

   // U and V are switched relative to standard formulation of CCA
   // so V is now for X and U is for Y

   std::cout << timestamp() << " Begin computing Xcoef/Ycoef" << std::endl;
   MatrixXd Xcoef = lltX.matrixL().transpose().solve(V);
   MatrixXd Ycoef = lltY.matrixL().transpose().solve(U);
   std::cout << timestamp() << " End computing Xcoef/Ycoef" << std::endl;

   std::cout << timestamp() << " Begin computing Xproj/Yproj" << std::endl;
   Px = X * Xcoef;
   Py = Y * Ycoef;
   std::cout << timestamp() << " dim(Px): " << dim(Px) << " dim(Py): " <<
      dim(Py) << std::endl;
   std::cout << timestamp() << " End computing Xproj/Yproj" << std::endl;

}

