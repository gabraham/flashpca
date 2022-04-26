/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2014 Gad Abraham
 * All rights reserved.
 */

#pragma once

#ifdef RENV
#include <Rcpp.h>
#define STDOUT Rcpp::Rcout
#else
#define STDOUT std::cout
#endif

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigen>

#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <iomanip>

#include <time.h>
#include <sys/time.h>

#define VAR_TOL 1e-9
#define STANDARDISE_NONE 0
#define STANDARDISE_SD 1
#define STANDARDISE_BINOM 2
#define STANDARDISE_BINOM2 3
#define STANDARDISE_CENTER 4

#define TXT_SEP "\t"

using namespace Eigen;

typedef Array<bool, Dynamic, Dynamic> ArrayXXb;
typedef Array<bool, Dynamic, 1> ArrayXb;

template <typename Derived>
void save(const char *filename, const Eigen::MatrixBase<Derived>& m)
{
   unsigned int rows = m.rows(), cols = m.cols();
   std::ofstream f(filename, std::ios::out | std::ios::binary);
   f.write((char *)&rows, sizeof(rows));
   f.write((char *)&cols, sizeof(cols));
   f.write((char *)m.derived().data(), sizeof(typename Derived::Scalar) * rows * cols);
   f.close();
}

template <typename Derived>
void load(const char *filename, Eigen::MatrixBase<Derived>& m)
{
   int nrows, ncols;
   std::ifstream f(filename, std::ios::binary);
   f.read(reinterpret_cast<char*>(&nrows), sizeof(int));
   f.read(reinterpret_cast<char*>(&ncols), sizeof(int));
   m.derived().resize(nrows, ncols);
   f.read((char*)m.derived().data(), sizeof(typename Derived::Scalar) * nrows * ncols);
}

template <typename Derived, typename A, typename B>
bool save_text(Eigen::MatrixBase<Derived>& m,
   const std::vector<A>& colnames,
   const std::vector<B>& rownames,
   const char *filename,
   const unsigned int precision=7)
{
   std::ofstream out(filename, std::ofstream::out);
   out << std::setprecision(precision);
   if(!out)
   {
#ifndef RENV
      std::cerr << "Error while saving to file " << filename 
	 << ":" << strerror(errno) << std::endl;
#endif
      return false;
   }

   // Write column names
   for(unsigned int i = 0 ; i < colnames.size() ; i++)
   {
      out << colnames[i];
      if(i == colnames.size() - 1) 
	 out << std::endl;
      else
         out << TXT_SEP;
   }
   
   const IOFormat fmt(precision , DontAlignCols, TXT_SEP, "\n", "", "", "", "");

   for(unsigned int j = 0 ; j < m.rows() ; j++)
   {
      if(rownames.size() > 0)
	 out << rownames[j] << TXT_SEP;
      out << m.row(j).format(fmt) << std::endl;
   }

   out.close();
   return true;
}

template <typename Derived>
std::string dim(Eigen::MatrixBase<Derived>& m)
{
   std::stringstream ss;
   ss << m.rows() << " x " << m.cols();
   return ss.str();
}

std::string timestamp();

template <typename T>
int sign(T x)
{
   return (T(0) < x) - (x < T(0));
}

Eigen::MatrixXd read_bed(const char *filename, const unsigned int nrows);
Eigen::MatrixXd read_pheno(const char *filename, unsigned int firstcol);
Eigen::MatrixXd standardise(Eigen::MatrixXd &X, int method, bool verbose=false);
Eigen::MatrixXd standardise_transpose(Eigen::MatrixXd &X, int method,
   bool verbose=false);

