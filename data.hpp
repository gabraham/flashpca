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

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigen>

#include <boost/filesystem.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#include "util.hpp"

#define PACK_DENSITY 4
#define PLINK_NA 3

#define PHENO_BINARY_12 0
#define PHENO_CONTINUOUS 1

#define PLINK_PHENO_MISSING -9

// The BED file magic numbers
#define PLINK_OFFSET 3

#define COVAR_ACTION_TRAIN_TEST 0
#define COVAR_ACTION_TRAIN_ONLY 1

#define COVAR_ACTION_TRAIN_TEST_STR "traintest"
#define COVAR_ACTION_TRAIN_ONLY_STR "train"

/* 3 is 11 in binary, we need a 2 bit mask for each of the 4 positions */
#define MASK0 3	  /* 3 << 2 * 0 */
#define MASK1 12  /* 3 << 2 * 1 */
#define MASK2 48  /* 3 << 2 * 2 */
#define MASK3 192 /* 3 << 2 * 3 */

#define BUFSIZE 100
#define DATA_MODE_TRAIN 1
#define DATA_MODE_TEST 2

using namespace Eigen;

class Data {
   public:
      
      MatrixXd X, Y;
      unsigned int N, p, K;
      unsigned long long len, np;
      unsigned int nsnps;
      const char *geno_filename;
      boost::iostreams::mapped_file_source geno_fin;
      bool verbose;
      long seed;
      std::vector<std::string> snp_ids;
      
      Data(long seed);
      ~Data();
      void prepare();
      void read_bed(bool transpose);
      void read_snp_block(unsigned int start_idx, unsigned int stop_idx,
	 bool transpose);
      void get_size();
      void read_pheno(const char *filename, unsigned int firstcol);
      void read_plink_bim(const char *filename);
      MatrixXd read_plink_pheno(const char *filename, unsigned int firstcol);

      std::string tolower(const std::string& v);

   private:
      unsigned char *tmp, *tmp2;
      std::ifstream in;
      double* avg;
      VectorXd tmpx;
      bool* visited;

      // the standardized values for the 3 genotypes + NA, for each SNP
      ArrayXXd scaled_geno_lookup;
};


