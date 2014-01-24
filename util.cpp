
#include "util.hpp"

using namespace Eigen;

// Standardize matrix column-wise to zero mean and unit variance.
// If a column is all zeros, it will remain zero.
MatrixXd standardize(const MatrixXd& X, bool scale, int method)
{
   std::cout.setf(std::ios_base::unitbuf);

   unsigned int n = X.rows(), p = X.cols();
   MatrixXd S = MatrixXd::Zero(X.rows(), X.cols());

   if(scale)
   {
      if(method == STANDARDIZE_SD)
      {
	 std::cout << timestamp() << " standardizing matrix (SD)" 
	    << " p: " << p << std::endl;
	 double mean, sd;
	 for(unsigned int j = 0 ; j < p ; j++)
      	 {
      	    mean = X.col(j).sum() / n;
      	    sd = std::sqrt((X.col(j).array() - mean).square().sum() / (n - 1));
      	    if(sd > VAR_TOL)
      	       S.col(j) = (X.col(j).array() - mean) / sd;
      	 }
      }
      // Same as Price 2006 eqn 3
      else if(method == STANDARDIZE_BINOM)
      {
	 std::cout << timestamp() << " standardizing matrix (BINOM)" 
	    << " p: " << p << std::endl;
	 double mean, r, s;
	 for(unsigned int j = 0 ; j < p ; j++)
      	 {
      	    mean = X.col(j).sum() / n;
	    r = mean / 2.0;
	    s = sqrt(r * (1 - r));
      	    if(s > VAR_TOL)
      	       S.col(j) = (X.col(j).array() - mean) / s;
      	 }
      }
      else
	 throw std::string("unknown standardization method");
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

