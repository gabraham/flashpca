/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2014-2016 Gad Abraham
 * All rights reserved.
 */

#include "util.h"

using namespace Eigen;

bool show_timestamp;

// Standardise matrix column-wise to zero mean and unit variance.
// *Standardises in place*
// If a column is all zeros, it will remain zero.
// Returns p by 2 matrix [mean, sd]
//
// Imputes missing values (nan) to the mean, where the mean was computed over
// all non-missing values
MatrixXd standardise(MatrixXd& X, int method, bool verbose)
{
#ifndef RENV
   std::cout.setf(std::ios_base::unitbuf);
#endif

   unsigned int n = X.rows(), p = X.cols();
   VectorXd mean = MatrixXd::Zero(X.cols(), 1);
   VectorXd sd = MatrixXd::Ones(X.cols(), 1);

   // Just check for missing values and impute to mean
   if(method == STANDARDISE_NONE || method == STANDARDISE_CENTER)
   {
      for(unsigned int j = 0 ; j < p ; j++)
      {
	 mean(j) = 0;
	 unsigned int nj = 0;
	 for(unsigned int i = 0 ; i < n ; i++)
	 {
	    double xij = X(i, j);
	    if(!std::isnan(xij))
	    {
	       mean(j) += xij;
	       nj++;
	    }
	 }
	 mean(j) /= (double)nj;

	 if(method == STANDARDISE_NONE)
	 {
	    for(unsigned int i = 0 ; i < n ; i++)
	       if(std::isnan(X(i, j)))
		  X(i, j) = mean(j);
	 }
	 else
	 {
	    for(unsigned int i = 0 ; i < n ; i++)
	    {
	       if(std::isnan(X(i, j)))
		  X(i, j) = 0;
	       else
		  X(i, j) = X(i, j) - mean(j);
	    }
	 }
      }
   }
   else if(method == STANDARDISE_SD)
   {
      verbose && STDOUT << timestamp() << "standardising matrix (SD)" 
	 << " p: " << p << std::endl;

      for(unsigned int j = 0 ; j < p ; j++)
      {
	 double sum = 0;
	 double sum_sqr = 0;
	 unsigned int nj = 0;
	 // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance,
	 // shifted_data_variance algorithm
	 double K = 1; // arbitrary
	 for(unsigned int i = 0 ; i < n ; i++)
	 {
	    double xij = X(i, j);
	    if(!std::isnan(xij))
	    {
	       sum += xij - K;
	       sum_sqr += (xij - K) * (xij - K);
	       nj++;
	    }
	 }
	 double varj = (sum_sqr - (sum * sum) / nj) / (nj - 1);
	 mean(j) = (sum + K * nj) / nj;
	 sd(j) = std::sqrt(varj);

	 // Note: using the stdev estimated on the subset of non-missing
	 // observations will make the final stdev of X(_, j) not b
	 // exactly 1.
	 for(unsigned int i = 0 ; i < n ; i++)
	 {
	    if(std::isnan(X(i, j)))
	       X(i, j) = 0;
	    else
	       X(i, j) = (X(i, j) - mean(j)) / sd(j);
	 }
      }
   }
   // Same as Price 2006 eqn 3
   else if(method == STANDARDISE_BINOM || method == STANDARDISE_BINOM2)
   {
      double mult = method == STANDARDISE_BINOM ? 1 : 2;

      verbose && STDOUT << timestamp() << "standardising matrix (BINOM/BINOM2)" 
	 << " p: " << p << std::endl;

      for(unsigned int j = 0 ; j < p ; j++)
      {
	 double sum = 0;
	 unsigned int nj = 0;
	 for(unsigned int i = 0 ; i < n ; i++)
	 {
	    double xij = X(i, j);
	    if(!std::isnan(xij))
	    {
	       sum += xij;
	       nj++;
	    }
	 }
	 mean(j) = sum / nj;
	 double r = mean(j) / 2.0;
	 sd(j) = std::sqrt(mult * r * (1.0 - r));

	 // Note: using the stdev estimated on the subset of non-missing
	 // observations will make the final stdev of X(_, j) not b
	 // exactly 1.
	 for(unsigned int i = 0 ; i < n ; i++)
	 {
	    if(std::isnan(X(i, j)))
	       X(i, j) = 0;
	    else
	       X(i, j) = (X(i, j) - mean(j)) / sd(j);
	 }
      }
   }
   //else if(method == STANDARDISE_CENTER)
   //{
   //   for(unsigned int j = 0 ; j < p ; j++)
   //   {
   //      double sum = 0;
   //      double sum_sqr = 0;
   //      unsigned int nj = 0;
   //      double K = 0;
   //      for(unsigned int i = 0 ; i < n ; i++)
   //      {
   //         double xij = X(i, j);
   //         if(!std::isnan(xij))
   //         {
   //            sum += xij - K;
   //            sum_sqr += (xij - K) * (xij - K);
   //            nj++;
   //         }
   //      }
   //      double varj = (sum_sqr - (sum * sum) / nj) / (nj - 1);
   //      mean(j) = (sum + K * nj) / nj;
   //      sd(j) = std::sqrt(varj);

   //      // Note: using the stdev estimated on the subset of non-missing
   //      // observations will make the final stdev of X(_, j) not b
   //      // exactly 1.
   //      for(unsigned int i = 0 ; i < n ; i++)
   //      {
   //         if(std::isnan(X(i, j)))
   //            X(i, j) = 0;
   //         else
   //            X(i, j) = (X(i, j) - mean(j)) / sd(j);
   //      }
   //   }
   //}
   else
      throw std::runtime_error(std::string("unknown standardization method"));

   MatrixXd P = MatrixXd::Zero(X.cols(), 2); // [mean, sd]
   P.col(0) = mean;
   P.col(1) = sd;

   return P;
}

// Expects a p times N matrix X, standardised in-place
MatrixXd standardise_transpose(MatrixXd& X, int method, bool verbose)
{
#ifndef RENV
   std::cout.setf(std::ios_base::unitbuf);
#endif

   unsigned int n = X.cols(), p = X.rows();
   VectorXd mean = MatrixXd::Zero(p, 1);
   VectorXd sd = MatrixXd::Ones(p, 1);

   if(method == STANDARDISE_SD)
   {
      verbose && STDOUT << timestamp() 
	 << " standardising transposed matrix (SD)" 
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
   else if(method == STANDARDISE_BINOM)
   {
      verbose && STDOUT << timestamp() 
	 << " standardising transposed matrix (BINOM)" 
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
   else if(method == STANDARDISE_BINOM2)
   {
      verbose && STDOUT << timestamp() 
	 << " standardising transposed matrix (BINOM2)"
	 << " p: " << p << std::endl;

      double r;
      for(unsigned int j = 0 ; j < p ; j++)
      {
	 mean(j) = X.row(j).sum() / n;
	 r = mean(j) / 2.0;
	 sd(j) = sqrt(2 * r * (1 - r)); // Note the factor of 2
	 if(sd(j) > VAR_TOL)
	    X.row(j) = (X.row(j).array() - mean(j)) / sd(j);
      }
   }
   else if(method == STANDARDISE_CENTER)
   {
      for(unsigned int j = 0 ; j < p ; j++)
      {
	 mean(j) = X.row(j).sum() / n;
	 X.row(j) = X.row(j).array() - mean(j);
      }
   }
   else
      throw std::runtime_error(std::string("unknown standardization method"));

   MatrixXd P = MatrixXd::Zero(p, 2); // [mean, sd]
   P.col(0) = mean;
   P.col(1) = sd;

   return P;
}

std::string timestamp()
{
   if(show_timestamp)
   {
      time_t t = time(NULL);
      char *s = asctime(localtime(&t));
      s[strlen(s) - 1] = '\0';
      std::string str(s);
      str = std::string("[") + str + std::string("] ");
      return str;
   }
   else
      return std::string("");
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

