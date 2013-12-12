#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <Eigen/SVD>
#include <Eigen/Sparse>

#include <boost/filesystem.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#include "util.hpp"
#include "cache.hpp"

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
      
      MatrixXd X, X2, Y;
      MatrixXd Xtrain, Xtest; // only used if data is small enough
      MatrixXd X2train, X2test;
      MatrixXd Ytrain, Ytest;
      unsigned int N, p, K;
      unsigned long long len;
      unsigned long long cachemem;
      unsigned int np, nsnps, ncovar;
      unsigned int Ntrain, Ntest, Ncurr;
      ArrayXb mask_train, mask_test, mask_curr;;
      char *geno_filename;
      boost::iostreams::mapped_file_source geno_fin;
      //std::vector<unsigned int> covar_ignore_pred_idx;
      std::vector<unsigned int> covar_actions;
      VectorXi folds;
      unsigned int nfolds;
      Cache *cache;
      unsigned int mode;
      VectorXd ones, zeros;
      VectorXd geno;
      VectorXd *geno_ptr;
      
      Data();
      ~Data();
      void read_bed(const char *filename);
      void read_pheno(const char *filename, unsigned int firstcol, int pheno);
      void read_covar(const char *filename, unsigned int firstcol);
      MatrixXd read_plink_pheno(const char *filename, unsigned int firstcol, int pheno);
      void set_mode(unsigned int mode);
      void load_snp_double(unsigned int j, double *geno);
      VectorXd load_snp(unsigned int j);

      std::string tolower(const std::string& v);
      void read_covar_actions(const char* filename);
      void mmap_bed(char *filename);
      VectorXd get_coordinate(unsigned int j);
      VectorXd get_snp(unsigned int j);
      void split_data(unsigned int fold);
      void make_folds(unsigned int rep);
};


