
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <Eigen/SVD>
#include <Eigen/QR>

using namespace Eigen;

MatrixXd random_pca(MatrixXd X, bool transpose=false,
   unsigned int ndim=10, unsigned int nextra=10, unsigned int maxiter=10);

