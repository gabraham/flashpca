
//#include <Eigen/Core>
//#include <Eigen/Dense>
//#include <Eigen/Eigen>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <string>
#include <fstream>
#include <sstream>

#include "data.hpp"
#include "randompca.hpp"

using namespace Eigen;

int main(int argc, char * argv[])
{
   std::cout << timestamp() << " Start" << std::endl;
   long seed = 1L;
   Data data(seed);
   data.verbose = false;
   std::cout << timestamp() << " seed: " << data.seed << std::endl;
   data.included_snps_filename = "flashpca_include_snps.txt";

   std::string pheno_file = std::string("data.fam");
   std::string geno_file = std::string("data.bed");
   std::string bim_file = std::string("data.bim");

   //data.regions.reserve(4);
   //region r;
   //r.chr = 5;
   //r.begin_bp = 44000000;
   //r.end_bp = 51500000;
   //data.regions.push_back(r);
   //
   //r.chr = 6;
   //r.begin_bp = 25000000;
   //r.end_bp = 33500000;
   //data.regions.push_back(r);

   //r.chr = 8;
   //r.begin_bp = 8000000;
   //r.end_bp = 12000000;
   //data.regions.push_back(r);

   //r.chr = 11;
   //r.begin_bp = 45000000;
   //r.end_bp = 57000000;
   //data.regions.push_back(r);

   data.bim_filename = bim_file.c_str();
   data.read_plink_bim();
   //data.nsnps_sampling = fminl(2e4, data.snps.size());
   data.nsnps_sampling = data.snps.size();
   std::cout << timestamp() << " Sampling " << data.nsnps_sampling << " SNPs" << std::endl;
   data.map_regions();
   if(data.nsnps_post_removal == 0)
   {
      std::cerr << "Error: no SNPs left for analysis" << std::endl;
      return EXIT_FAILURE;
   }
      
   data.read_pheno(pheno_file.c_str(), 6, PHENO_BINARY_12);
   data.read_bed(geno_file.c_str());

   std::cout << timestamp() << " Begin PCA" << std::endl;

   //bool transpose = data.X.rows() < data.X.cols();
   RandomPCA rpca;
   bool transpose = false;
   rpca.stand_method = STANDARDIZE_BINOMIAL;
   unsigned int max_dim = fminl(data.X.rows(), data.X.cols());
   unsigned int n_dim = fminl(max_dim, 100);
   unsigned int n_extra = fminl(max_dim - n_dim, 50);
   int method = METHOD_EIGEN;
   rpca.pca(data.X, method, transpose, n_dim, n_extra);
   //std::cout << timestamp() << " Writing matrix V" << std::endl;
  // save_text("V.txt", rpca.V);
   std::cout << timestamp() << " Writing eigenvectors" << std::endl;
   save_text("eigenvectors.txt", rpca.U);
   std::cout << timestamp() << " Writing PCs" << std::endl;
   save_text("pcs.txt", rpca.P);
   save_text("eigenvalues.txt", rpca.d);

   bool whiten = false;
   if(whiten)
   {
      std::cout << timestamp() << " Loading all SNPs" << std::endl;
      data.reset_regions();
      data.read_bed(geno_file.c_str());
      rpca.zca_whiten();
      std::cout << timestamp() << " Loading SNPs done" << std::endl;
      std::cout << timestamp() << " Writing whitened data" << std::endl;
      save("whitened.bin", rpca.W);
   }

   std::cout << timestamp() << " Goodbye!" << std::endl;

   return EXIT_SUCCESS;
}

