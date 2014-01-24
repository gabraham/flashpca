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
      unsigned long long len;
      unsigned int np, nsnps;
      const char *geno_filename;
      boost::iostreams::mapped_file_source geno_fin;
      bool verbose;
      long seed;
      
      Data(long seed);
      ~Data();
      void read_bed(const char *filename);
      void read_pheno(const char *filename, unsigned int firstcol, int pheno);
      MatrixXd read_plink_pheno(const char *filename, unsigned int firstcol, int pheno);

      std::string tolower(const std::string& v);
};


