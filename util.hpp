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
#define STANDARDIZE_NONE 0
#define STANDARDIZE_SD 1
#define STANDARDIZE_BINOM 2
#define STANDARDIZE_BINOM2 3
#define STANDARDIZE_CENTER 4

using namespace Eigen;

typedef Array<bool, Dynamic, Dynamic> ArrayXXb;
typedef Array<bool, Dynamic, 1> ArrayXb;


template <typename Derived>
void save(const char *filename, const MatrixBase<Derived>& m)
{
   unsigned int rows = m.rows(), cols = m.cols();
   std::ofstream f(filename, std::ios::out | std::ios::binary);
   f.write((char *)&rows, sizeof(rows));
   f.write((char *)&cols, sizeof(cols));
   f.write((char *)m.derived().data(), sizeof(typename Derived::Scalar) * rows * cols);
   f.close();
}

template <typename Derived>
void load(const char *filename, MatrixBase<Derived>& m)
{
   int nrows, ncols;
   std::ifstream f(filename, std::ios::binary);
   f.read(reinterpret_cast<char*>(&nrows), sizeof(int));
   f.read(reinterpret_cast<char*>(&ncols), sizeof(int));
   m.derived().resize(nrows, ncols);
   f.read((char*)m.derived().data(), sizeof(typename Derived::Scalar) * nrows * ncols);
}

template <typename Derived, typename T>
bool save_text(MatrixBase<Derived>& m,
   const std::vector<T>& colnames,
   const char *filename,
   const unsigned int precision=6)
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

   for(unsigned int i = 0 ; i < colnames.size() ; i++)
   {
      out << colnames[i];
      if(i == colnames.size() - 1) 
	 out << std::endl;
      else
         out << "\t";
   }
   
   const IOFormat fmt(6, DontAlignCols, "\t", "\n", "", "", "", "");
   out << m.format(fmt) << std::endl;
   out.close();
   return true;
}

template <typename Derived>
std::string dim(MatrixBase<Derived>& m)
{
   std::stringstream ss;
   ss << m.rows() << " x " << m.cols();
   return ss.str();
}

double myatof(char* c);
std::string timestamp();

template <typename T>
int sign(T x)
{
   return (T(0) < x) - (x < T(0));
}

MatrixXd read_bed(const char *filename, const unsigned int nrows);
MatrixXd read_pheno(const char *filename, unsigned int firstcol);
MatrixXd standardize(MatrixXd &X, int method, bool verbose=false);
MatrixXd standardize_transpose(MatrixXd &X, int method,
   bool verbose=false);

