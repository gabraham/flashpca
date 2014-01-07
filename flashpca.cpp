
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <Eigen/SVD>

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
   Data data;

   data.verbose = false;
   //data.seed = time(NULL);
   data.seed = 1;
   std::cout << timestamp() << " seed: " << data.seed << std::endl;
   data.included_snps_filename = "flashpca_include_snps.txt";

   std::string pheno_file = std::string("data.fam");
   std::string geno_file = std::string("data.bed");
   std::string bim_file = std::string("data.bim");

   data.regions.reserve(4);
   region r;
   r.chr = 5;
   r.begin_bp = 44000000;
   r.end_bp = 51500000;
   data.regions.push_back(r);
   
   r.chr = 6;
   r.begin_bp = 25000000;
   r.end_bp = 33500000;
   data.regions.push_back(r);

   r.chr = 8;
   r.begin_bp = 8000000;
   r.end_bp = 12000000;
   data.regions.push_back(r);

   r.chr = 11;
   r.begin_bp = 45000000;
   r.end_bp = 57000000;
   data.regions.push_back(r);

   data.bim_filename = bim_file.c_str();
   data.read_plink_bim();
   data.nsnps_sampling = fminl(1e5, data.snps.size());
   std::cout << timestamp() << " Sampling " << data.nsnps_sampling << " SNPs" << std::endl;
   data.map_regions();
   if(data.nsnps_post_removal == 0)
   {
      std::cerr << "Error: no SNPs left for analysis" << std::endl;
      return EXIT_FAILURE;
   }
      
   data.read_pheno(pheno_file.c_str(), 6, PHENO_BINARY_12);
   data.read_bed(geno_file.c_str());

   std::cout << timestamp() << " Begin SVD" << std::endl;

   //bool transpose = data.X.rows() < data.X.cols();
   RandomPCA rpca;
   bool transpose = false;
   rpca.pca(data.X, transpose, 50);
   std::cout << timestamp() << " Writing PCs" << std::endl;
   save_text("pcs.txt", rpca.P);

   bool whiten = false;
   if(whiten)
   {
      rpca.zca_whiten();
      // TODO: we want to load ALL SNPs and whiten them, not just the sampled
      // ones
      std::cout << timestamp() << " Loading all SNPs" << std::endl;
      data.reset_regions();
      data.read_bed(geno_file.c_str());
      std::cout << timestamp() << " Loading SNPs done" << std::endl;
      std::cout << timestamp() << " Writing whitened data" << std::endl;
      save("whitened.bin", rpca.W);
   }

   std::cout << timestamp() << " Goodbye!" << std::endl;

   return EXIT_SUCCESS;
}

