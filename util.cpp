
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

MatrixXd standardize_transpose(const MatrixXd& X, bool scale, int method)
{
   std::cout.setf(std::ios_base::unitbuf);

   unsigned int n = X.cols(), p = X.rows();
   MatrixXd S = MatrixXd::Zero(X.rows(), X.cols());

   if(scale)
   {
      if(method == STANDARDIZE_SD)
      {
	 std::cout << timestamp() 
	    << " standardizing transposed matrix (SD)" 
	    << " p: " << p << std::endl;
	 double mean, sd;
	 for(unsigned int j = 0 ; j < p ; j++)
      	 {
      	    mean = X.row(j).sum() / n;
      	    sd = std::sqrt((X.row(j).array() - mean).square().sum() / (n - 1));
      	    if(sd > VAR_TOL)
      	       S.row(j) = (X.row(j).array() - mean) / sd;
      	 }
      }
      // Same as Price 2006 eqn 3
      else if(method == STANDARDIZE_BINOM)
      {
	 std::cout << timestamp() 
	    << " standardizing transposed matrix (BINOM)" 
	    << " p: " << p << std::endl;
	 double mean, r, s;
	 for(unsigned int j = 0 ; j < p ; j++)
      	 {
      	    mean = X.row(j).sum() / n;
	    r = mean / 2.0;
	    s = sqrt(r * (1 - r));
      	    if(s > VAR_TOL)
      	       S.row(j) = (X.row(j).array() - mean) / s;
      	 }
      }
      else
	 throw std::string("unknown standardization method");
   }
   else
   {
      for(unsigned int j = 0 ; j < p ; j++)
      {
         double mean = X.row(j).sum() / n;
	 S.row(j) = X.row(j).array() - mean;
      }
   }
   return S;
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

/*
 * Based on http://ndevilla.free.fr/median/median/src/torben.c
 * Algorithm by Torben Mogensen, implementation by N. Devillard.
 * This code in public domain.
 */
double median_torben(double m[], int n)
{
    int         i, less, greater, equal;
    double  min, max, guess, maxltguess, mingtguess;

    min = max = m[0] ;
    for (i=1 ; i<n ; i++) {
        if (m[i]<min) min=m[i];
        if (m[i]>max) max=m[i];
    }

    while (1) {
        guess = (min+max)/2;
        less = 0; greater = 0; equal = 0;
        maxltguess = min ;
        mingtguess = max ;
        for (i=0; i<n; i++) {
            if (m[i]<guess) {
                less++;
                if (m[i]>maxltguess) maxltguess = m[i] ;
            } else if (m[i]>guess) {
                greater++;
                if (m[i]<mingtguess) mingtguess = m[i] ;
            } else equal++;
        }
        if (less <= (n+1)/2 && greater <= (n+1)/2) break ; 
        else if (less>greater) max = maxltguess ;
        else min = mingtguess;
    }
    if (less >= (n+1)/2) return maxltguess;
    else if (less+equal >= (n+1)/2) return guess;
    else return mingtguess;
}

