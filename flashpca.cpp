
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
   Data data;

   std::string pheno_file = std::string("data.fam");
   std::string geno_file = std::string("data.bed");
   
   data.read_pheno(pheno_file.c_str(), 6, PHENO_BINARY_12);
   data.read_bed(geno_file.c_str());

   std::cout << ">>> Begin SVD" << std::endl;
   // TODO: sample SNPs
   // TODO: filter by LD
   
   MatrixXd P = random_pca(data.X, false);

   save_text("pcs.txt", P);
   
   std::cout << ">>> Goodbye!" << std::endl;

   return EXIT_SUCCESS;
}

