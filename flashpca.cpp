
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

#ifndef NDEBUG
   std::cout << "******* Running in Debug mode *******" << std::endl;
#endif

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
      ("ndim", po::value<int>(), "number of PCs to output")
      ("nextra", po::value<int>(),
	 "number of extra dimensions to use in randomized PCA")
      ("stand", po::value<std::string>(), "standardization method")
      ("method", po::value<std::string>(), "PCA method (svd/eigen)")
      ("outpc", po::value<std::string>(), "PC output file")
      ("outvec", po::value<std::string>(), "Eigenvector output file")
      ("outval", po::value<std::string>(), "Eigenvalue output file")
      ("whiten", "whiten the data")
      ("outwhite", po::value<std::string>(), "whitened data output file")
      ("v", "verbose")
      ("maxiter", po::value<int>(), "maximum number of randomized PCA iterations")
      ("tol", po::value<double>(), "tolerance for randomized PCA iterations")
      ("transpose", "force a transpose of the data")
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

   int n_dim = 0;
   if(vm.count("ndim"))
   {
      n_dim = vm["ndim"].as<int>();
      if(n_dim <= 1)
      {
	 std::cerr << "Error: --ndim can't be less than 2" << std::endl;
	 return EXIT_FAILURE;
      }
   }

   int n_extra = 0;
   if(vm.count("nextra"))
   {
      n_extra = vm["nextra"].as<int>();
      if(n_extra < 1)
      {
	 std::cerr << "Error: --n_extra can't be less than 1" << std::endl;
	 return EXIT_FAILURE;
      }
   }
   
   int stand_method = STANDARDIZE_BINOM;
   if(vm.count("stand"))
   {
      std::string m = vm["stand"].as<std::string>();
      if(m == "binom")
	 stand_method = STANDARDIZE_BINOM;
      else if(m == "sd")
	 stand_method = STANDARDIZE_SD;
      else
      {
	 std::cerr << "Error: unknown standardization method (--stand): "
	    << m << std::endl;
	 return EXIT_FAILURE;
      }
   }

   int method = METHOD_EIGEN;
   if(vm.count("method"))
   {
      std::string m = vm["method"].as<std::string>();
      if(m == "svd")
	 method = METHOD_SVD;
      else if(m == "eigen")
	 method = METHOD_EIGEN;
      else
      {
	 std::cerr << "Error: unknown PCA method (--method): "
	    << m << std::endl;
	 return EXIT_FAILURE;
      }
   }

   std::string pcfile = "pcs.txt";
   if(vm.count("outpc"))
      pcfile = vm["outpc"].as<std::string>();

   std::string eigvecfile = "eigenvectors.txt";
   if(vm.count("outvec"))
      eigvecfile = vm["outvec"].as<std::string>();

   std::string eigvalfile = "eigenvalues.txt";
   if(vm.count("outval"))
      eigvalfile = vm["outval"].as<std::string>();

   std::string whitefile = "whitened.bin";
   if(vm.count("outwhite"))
      whitefile = vm["outwhite"].as<std::string>();

   bool whiten = vm.count("whiten");
   bool verbose = vm.count("v");
   bool transpose = vm.count("transpose");
   
   int maxiter = 10;
   if(vm.count("maxiter"))
   {
      maxiter = vm["maxiter"].as<int>();
      if(maxiter <= 0)
      {
	 std::cerr << "Error: --maxiter can't be less than 1"
	    << std::endl;
	 return EXIT_FAILURE;
      }
   }

   double tol = 1e-7;
   if(vm.count("tol"))
   {
      tol = vm["tol"].as<double>();
      if(tol <= 0)
      {
	 std::cerr << "Error: --tol can't be zero or negative"
	    << std::endl;
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
   data.verbose = verbose;
   std::cout << timestamp() << " seed: " << data.seed << std::endl;

   data.read_pheno(pheno_file.c_str(), 6, PHENO_BINARY_12);
   data.geno_filename = geno_file.c_str();
   data.get_size();
   transpose = transpose || data.N > data.nsnps;
   data.read_bed(transpose);

   RandomPCA rpca;
   rpca.stand_method = stand_method;
   unsigned int max_dim = fminl(data.X.rows(), data.X.cols());
   
   if(n_dim == 0)
      n_dim = fminl(max_dim, 10);

   if(n_extra == 0)
      n_extra = fminl(max_dim - n_dim, 190);

   ////////////////////////////////////////////////////////////////////////////////
   // Do the PCA
   std::cout << timestamp() << " PCA begin" << std::endl;
   
   rpca.pca(data.X, method, transpose, n_dim, n_extra, maxiter, tol);

   std::cout << timestamp() << " PCA done" << std::endl;

   ////////////////////////////////////////////////////////////////////////////////
   // Write out results

   std::cout << timestamp() << " Writing " << n_dim << 
      " eigenvectors to file " << eigvecfile << std::endl;
   save_text(eigvecfile.c_str(), rpca.U);

   std::cout << timestamp() << " Writing " << n_dim << 
      " eigenvalues to file " << eigvalfile << std::endl;
   save_text(eigvalfile.c_str(), rpca.d);

   std::cout << timestamp() << " Writing " << n_dim <<
      " PCs to file " << pcfile << std::endl;
   save_text(pcfile.c_str(), rpca.P);

   ////////////////////////////////////////////////////////////////////////////////
   // Whiten if required

   if(whiten)
   {
      std::cout << timestamp() << " ZCA whitening data" << std::endl;
      rpca.zca_whiten();
      std::cout << timestamp() << " Writing whitened data to file "
	 << whitefile << std::endl;
      save(whitefile.c_str(), rpca.W);
   }

   std::cout << timestamp() << " Goodbye!" << std::endl;

   return EXIT_SUCCESS;
}

