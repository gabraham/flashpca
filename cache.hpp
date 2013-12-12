#pragma once

#include <iostream>
#include <vector>
#include <algorithm>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <Eigen/SVD>
#include <Eigen/Sparse>

#define CACHE_NOT_EXISTS -1

using namespace Eigen;

class Cache {

   public:

      int nbins, n;
      int *mapping; // must be signed
      int *revmapping; // must be signed
      int lastfree;
      unsigned int *counter;
      unsigned long long hits, misses;

      //double *x;
      //double *tmp;
      VectorXd *x;
      VectorXd *tmp;
      
      Cache();
      Cache(unsigned int n, unsigned int p, unsigned long long maxmem=12884901888);
      ~Cache();

      //bool get(unsigned int j, double **x);
      bool get(unsigned int j, VectorXd& x);
      void put(unsigned int j, VectorXd& x);
};

