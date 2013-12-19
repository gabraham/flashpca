
#include "util.hpp"

using namespace Eigen;

// Standardize matrix column-wise to zero mean and unit variance.
// If a column is all zeros, it will remain zero.
MatrixXd standardize(const MatrixXd& X, bool scale)
{
   std::cout.setf(std::ios_base::unitbuf);

   unsigned int n = X.rows(), p = X.cols();
   MatrixXd S = MatrixXd::Zero(X.rows(), X.cols());

   if(scale)
   {
      for(unsigned int j = 0 ; j < p ; j++)
      {
         double mean = X.col(j).sum() / n;
         double sd = std::sqrt((X.col(j).array() - mean).square().sum() / (n - 1));
         if(sd > VAR_TOL)
            S.col(j) = (X.col(j).array() - mean) / sd;
      }
   }
   else
   {
      for(unsigned int j = 0 ; j < p ; j++)
      {
         double mean = X.col(j).sum() / n;
	 S.col(j) = X.col(j).array() - mean;
      }
   }
   return S;
}

void usage()
{
   std::cout << "Error: " << std::endl;
}

double myatof(char* c)
{
   if(c == 0)
      return 0;
   return atof(c);
}

// saves a compressed column sparse matrix in a *row-major* format:
//
//    <row>, <column>:<value> <column>:<value> ....
//
// dim: <nrow, ncol> of type int
// nonzero: of type int
// values: of type double
// inner indices: 
void save_sparse(const char *filename, const SparseMatrix<double>& m,
   bool beta_zerobased)
{
   std::ofstream f(filename, std::ios::out);

   // Convert to row-major
   SparseMatrix<double, RowMajor> mr(m);

   int a = beta_zerobased ? 0 : 1;

   for(unsigned int k = 0 ; k < mr.outerSize() ; k++)
   {
      SparseMatrix<double, RowMajor>::InnerIterator it(mr, k);
      if(it)
      {
	 f << it.row() + a << ",";
	 for( ; it ; ++it)
	    f << it.col() << ":" << it.value()  << " ";
	 f << std::endl;
      }
   }

   f.close();
}

std::string timestamp()
{
   time_t t = time(NULL);
   char *s = asctime(localtime(&t));
   s[strlen(s) - 1] = '\0';
   std::string str(s);
   str = std::string("[") + str + std::string("]");
   return str;
}

double sign_scalar(double x)
{
   return (0 < x) - (x < 0);
}

VectorXd zapsmall(const VectorXd& X, const double tol)
{
   VectorXd Z = (X.array().abs() > tol).select(X, 0);
   return Z; 
}

double maxabs(const SMd& B)
{
   unsigned int outerSize = B.outerSize();
   double max = 0, v;
   
   for(unsigned int k = 0 ; k < outerSize ; ++k)
   {
      for(SparseMatrix<double>::InnerIterator it(B, k) ; it ; ++it)
      {
	 v = std::abs(it.value());
	 if(max < v)
	    max = v;
      }
   }

   return max;
}

