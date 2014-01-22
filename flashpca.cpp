
//#include <Eigen/Core>
//#include <Eigen/Dense>
//#include <Eigen/Eigen>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <boost/program_options.hpp>

#include <string>
#include <fstream>
#include <sstream>

#include "data.hpp"
#include "randompca.hpp"

using namespace Eigen;
namespace po = boost::program_options;

int main(int argc, char * argv[])
{

   ////////////////////////////////////////////////////////////////////////////////
   // Parse commandline arguments
   po::options_description desc("Options");
   desc.add_options()
      ("help", "produce help message")
      ("numthreads", po::value<int>(), "set number of OpenMP threads")
      ("seed", po::value<long>(), "set random seed")
      ("bed", po::value<std::string>(), "PLINK bed file")
      ("bim", po::value<std::string>(), "PLINK bim file")
      ("fam", po::value<std::string>(), "PLINK fam file")
      ("bfile", po::value<std::string>(), "PLINK root name")
      ("ndim", po::value<unsigned int>(), "number of PCs to output")
      ("nextra", po::value<unsigned int>(),
	 "number of extra dimensions to use in randomized PCA")
      ("stand", po::value<std::string>(), "standardization method")
   ;

   po::variables_map vm;
   po::store(po::parse_command_line(argc, argv, desc), vm);
   po::notify(vm);

   if(vm.count("help"))
   {
      std::cerr << desc << std::endl;
      return EXIT_FAILURE;
   }

   int num_threads = 1;
   if(vm.count("numthreads"))
      num_threads = vm["numthreads"].as<int>();


   long seed = 1L;
   if(vm.count("seed"))
      seed = vm["seed"].as<long>();

   std::string pheno_file, geno_file, bim_file;

   if(vm.count("bfile"))
   {
      geno_file = vm["bfile"].as<std::string>() + std::string(".bed");
      bim_file = vm["bfile"].as<std::string>() + std::string(".bim");
      pheno_file = vm["bfile"].as<std::string>() + std::string(".fam");
   }
   else
   {
      bool good = true;
      if(vm.count("bed"))
	 geno_file = vm["bed"].as<std::string>();
      else
	 good = false;

      if(good && vm.count("bim"))
	 bim_file = vm["bim"].as<std::string>();
      else
	 good = false;

      if(good && vm.count("fam"))
	 pheno_file = vm["fam"].as<std::string>();
      else
	 good = false;
      
      if(!good)
      {
	 std::cerr << "Error: you must specify either --bfile "
	    << "or --bed / --fam / --bim" << std::endl;
	 return EXIT_FAILURE;
      }
   }

   unsigned int n_dim = 0;
   if(vm.count("ndim"))
      n_dim = vm["ndim"].as<unsigned int>();

   unsigned int n_extra = 0;
   if(vm.count("nextra"))
      n_extra = vm["nextra"].as<unsigned int>();
   
   int stand_method;
   if(vm.count("stand"))
   {
      std::string m = vm["stand"].as<std::string>();
      if(m == "price")
	 stand_method = STANDARDIZE_BINOMIAL;
      else if(m == "sd")
	 stand_method = STANDARDIZE_SD;
      else
      {
	 std::cerr << "Error: unknown standardization method (-stand): "
	    << m << std::endl;
	 return EXIT_FAILURE;
      }
   }


   ////////////////////////////////////////////////////////////////////////////////
   // End command line parsing
      
   std::cout << timestamp() << " Start" << std::endl;
   setNbThreads(num_threads);
   std::cout << timestamp() << " Using " << num_threads 
      << " OpenMP threads" << std::endl;

   Data data(seed);
   data.verbose = false;
   std::cout << timestamp() << " seed: " << data.seed << std::endl;
   data.included_snps_filename = "flashpca_include_snps.txt";

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
   data.nsnps_sampling = data.snps.size();
   //std::cout << timestamp() << " Sampling " << data.nsnps_sampling << " SNPs" << std::endl;
   //data.map_regions();
   //if(data.nsnps_post_removal == 0)
   //{
   //   std::cerr << "Error: no SNPs left for analysis" << std::endl;
   //   return EXIT_FAILURE;
   //}
      
   data.read_pheno(pheno_file.c_str(), 6, PHENO_BINARY_12);
   data.read_bed(geno_file.c_str());

   std::cout << timestamp() << " Begin PCA" << std::endl;

   //bool transpose = data.X.rows() < data.X.cols();
   RandomPCA rpca;
   bool transpose = false;
   rpca.stand_method = stand_method;
   unsigned int max_dim = fminl(data.X.rows(), data.X.cols());
   
   if(n_dim == 0)
      n_dim = fminl(max_dim, 100);

   if(n_extra == 0)
      n_extra = fminl(max_dim - n_dim, 100);
   int method = METHOD_EIGEN;
   rpca.pca(data.X, method, transpose, n_dim, n_extra);
   //std::cout << timestamp() << " Writing matrix V" << std::endl;
  // save_text("V.txt", rpca.V);
   std::cout << timestamp() << " Writing " << n_dim << " eigenvectors" << std::endl;
   save_text("eigenvectors.txt", rpca.U);
   std::cout << timestamp() << " Writing " << n_dim << " PCs" << std::endl;
   save_text("pcs.txt", rpca.P);
   save_text("eigenvalues.txt", rpca.d);

   bool whiten = false;
   if(whiten)
   {
      std::cout << timestamp() << " Loading all SNPs" << std::endl;
      //data.reset_regions();
      data.read_bed(geno_file.c_str());
      rpca.zca_whiten();
      std::cout << timestamp() << " Loading SNPs done" << std::endl;
      std::cout << timestamp() << " Writing whitened data" << std::endl;
      save("whitened.bin", rpca.W);
   }

   std::cout << timestamp() << " Goodbye!" << std::endl;

   return EXIT_SUCCESS;
}

