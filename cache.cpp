
#include <stdexcept>
#include "cache.hpp"

////////////////////////////////////////////////////////////////////////////////
// A ring buffer for SNP data, with priority insertion based on how frequently
// data has been requested

Cache::Cache()
{
   hits = 0;
   misses = 0;
}

Cache::Cache(unsigned int n, unsigned int p, unsigned long long maxmem)
{
   nbins = (unsigned int)(maxmem / sizeof(double) / n);
   this->n = n;
   lastfree = 0;

   if(nbins < 1)
      throw std::runtime_error("Not enough cache memory to load data");
   if(nbins > p || nbins == -1)
      nbins = p;

   std::cout << ">>> Initialising cache, maxmem: " << maxmem
      << " n: " << n << ", p: " << p 
      << ", nbins: " << nbins << std::endl;

   mapping = new int[p]; // map SNP to a bin
   revmapping = new int[nbins]; // map bin to a SNP
   counter = new unsigned int[p]();
   x = new VectorXd[(unsigned int)nbins * n];
   tmp = new VectorXd[n];

   for(unsigned int j = 0 ; j < p ; j++)
      mapping[j] = CACHE_NOT_EXISTS;

   for(unsigned int j = 0 ; j < nbins ; j++)
      revmapping[j] = CACHE_NOT_EXISTS;
   hits = 0;
   misses = 0;
}

Cache::~Cache()
{
   std::cout << ">>> Destroying cache" << std::endl;
   delete[] mapping;
   delete[] revmapping;
   delete[] counter;
   delete[] x;
   delete[] tmp;
}

bool Cache::get(unsigned int j, VectorXd& x)
{
   int m = mapping[j];
   counter[j]++;

   if(m == CACHE_NOT_EXISTS)
   {
      misses++;
      return false;
   }

   x = this->x[m];
   hits++;
   return true;
}

// Tries to insert a vector at the next free spot. If that spot is taken, see
// whether the new vector has been requested more times than the one already
// there.
void Cache::put(unsigned int j, VectorXd& x)
{
   int m = mapping[j];
   int r = revmapping[lastfree];
   int k;

   // jth vector not already there, so insert it
   if(r == CACHE_NOT_EXISTS)
      k = j;
   // jth vector already there, check the usage counter
   // to determine whether to keep old one
   else if(counter[j] > counter[r])
   {
      k = j;
      mapping[r] = CACHE_NOT_EXISTS;
   }
   else // do nothing
      return;

   // insert the chosen vector
   mapping[k] = lastfree;
   this->x[mapping[k]] = x;
   revmapping[lastfree] = k;
   lastfree++;
   
   if(lastfree >= nbins || lastfree < 0)
      lastfree = 0;
}

