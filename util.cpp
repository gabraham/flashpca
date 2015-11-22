/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2014-2015 Gad Abraham
 * All rights reserved.
 */

#include "util.hpp"

using namespace Eigen;

// Standardize matrix column-wise to zero mean and unit variance.
// *Standardizes in place*
// If a column is all zeros, it will remain zero.
// Returns p by 2 matrix [mean, sd]
MatrixXd standardize(MatrixXd& X, int method, bool verbose)
{
#ifndef RENV
   std::cout.setf(std::ios_base::unitbuf);
#endif

   unsigned int n = X.rows(), p = X.cols();
   VectorXd mean = MatrixXd::Zero(X.cols(), 1);
   VectorXd sd = MatrixXd::Ones(X.cols(), 1);

   if(method == STANDARDIZE_SD)
   {
      verbose && STDOUT << timestamp() << " standardizing matrix (SD)" 
	 << " p: " << p << std::endl;

      for(unsigned int j = 0 ; j < p ; j++)
      {
	 mean(j) = X.col(j).sum() / n;
	 sd(j) = std::sqrt((X.col(j).array() - mean(j)).square().sum() / (n - 1));
	 if(sd(j) > VAR_TOL)
	    X.col(j) = (X.col(j).array() - mean(j)) / sd(j);
      }
   }
   // Same as Price 2006 eqn 3
   else if(method == STANDARDIZE_BINOM)
   {
      verbose && STDOUT << timestamp() << " standardizing matrix (BINOM)" 
	 << " p: " << p << std::endl;

      double r;
      for(unsigned int j = 0 ; j < p ; j++)
      {
	 mean(j) = X.col(j).sum() / n;
	 r = mean(j) / 2.0;
	 sd(j) = sqrt(r * (1 - r));
	 if(sd(j) > VAR_TOL)
	    X.col(j) = (X.col(j).array() - mean(j)) / sd(j);
      }
   }
   else if(method == STANDARDIZE_CENTER)
   {
      for(unsigned int j = 0 ; j < p ; j++)
      {
         mean(j) = X.col(j).sum() / n;
	 X.col(j) = X.col(j).array() - mean(j);
      }
   }
   else
      throw std::string("unknown standardization method");

   MatrixXd P = MatrixXd::Zero(X.cols(), 2); // [mean, sd]
   P.col(0) = mean;
   P.col(1) = sd;

   return P;
}

// Expects a p times N matrix X, standardized in-place
MatrixXd standardize_transpose(MatrixXd& X, int method, bool verbose)
{
#ifndef RENV
   std::cout.setf(std::ios_base::unitbuf);
#endif

   unsigned int n = X.cols(), p = X.rows();
   VectorXd mean = MatrixXd::Zero(p, 1);
   VectorXd sd = MatrixXd::Ones(p, 1);

   if(method == STANDARDIZE_SD)
   {
      verbose && STDOUT << timestamp() 
	 << " standardizing transposed matrix (SD)" 
	 << " p: " << p << std::endl;

      for(unsigned int j = 0 ; j < p ; j++)
      {
	 mean(j) = X.row(j).sum() / n;
	 sd(j) = std::sqrt((X.row(j).array() - mean(j)).square().sum() / (n - 1));
	 if(sd(j) > VAR_TOL)
	    X.row(j) = (X.row(j).array() - mean(j)) / sd(j);
      }
   }
   // Same as Price 2006 eqn 3
   else if(method == STANDARDIZE_BINOM)
   {
      verbose && STDOUT << timestamp() 
	 << " standardizing transposed matrix (BINOM)" 
	 << " p: " << p << std::endl;

      double r;
      for(unsigned int j = 0 ; j < p ; j++)
      {
	 mean(j) = X.row(j).sum() / n;
	 r = mean(j) / 2.0;
	 sd(j) = sqrt(r * (1 - r));
	 if(sd(j) > VAR_TOL)
	    X.row(j) = (X.row(j).array() - mean(j)) / sd(j);
      }
   }
   else if(method == STANDARDIZE_CENTER)
   {
      for(unsigned int j = 0 ; j < p ; j++)
      {
	 mean(j) = X.row(j).sum() / n;
	 X.row(j) = X.row(j).array() - mean(j);
      }
   }
   else
      throw std::string("unknown standardization method");

   MatrixXd P = MatrixXd::Zero(p, 2); // [mean, sd]
   P.col(0) = mean;
   P.col(1) = sd;

   return P;
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

