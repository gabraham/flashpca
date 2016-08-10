
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2014-2016 Gad Abraham
 * All rights reserved.
 */

//#include <Eigen/Core>
//#include <Eigen/Dense>
//#include <Eigen/Eigen>

#include <boost/program_options.hpp>

#include <string>
#include <fstream>
#include <sstream>

#include "data.hpp"
#include "randompca.hpp"

using namespace Eigen;
namespace po = boost::program_options;

extern bool show_timestamp;

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
      ("scca", "perform sparse canonical correlation analysis (SCCA)")
      ("ucca", "perform per-SNP canonical correlation analysis")
      ("project,p", "project new samples onto existing principal components")
      ("batch,b", "load all genotypes into RAM at once")
      ("rand,r", "use the legacy randomised algorithm")
      ("memory,m", po::value<int>(), "size of block for online algorithm, in GB")
      ("blocksize,b", po::value<int>(),
	 "size of block for online algorithm, in number of SNPs")
      ("numthreads,n", po::value<int>(), "set number of OpenMP threads")
      ("seed", po::value<long>(), "set random seed")
      ("bed", po::value<std::string>(), "PLINK bed file")
      ("bim", po::value<std::string>(), "PLINK bim file")
      ("fam", po::value<std::string>(), "PLINK fam file")
      ("pheno", po::value<std::string>(), "PLINK phenotype file")
      ("bfile", po::value<std::string>(), "PLINK root name")
      ("ndim,d", po::value<int>(), "number of PCs to output")
      ("nextra", po::value<int>(),
	 "number of extra dimensions to use in randomized PCA")
      ("standx,s", po::value<std::string>(),
	 "standardization method for genotypes [binom2 | binom | none | sd | center]")
      ("standy", po::value<std::string>(),
	 "standardization method for phenotypes [sd | binom2 | binom | none | center]")
      ("method", po::value<std::string>(), "PCA method [eigen | svd]")
      ("orth", po::value<std::string>(), "use orthornormalization [yes | no]")
      ("div", po::value<std::string>(),
	 "whether to divide the eigenvalues by p, n - 1, or don't divide [p | n1 | none]")
      ("outpc", po::value<std::string>(), "PC output file")
      ("outpcx", po::value<std::string>(), "X PC output file, for CCA")
      ("outpcy", po::value<std::string>(), "Y PC output file, for CCA")
      ("outvec", po::value<std::string>(), "eigenvector output file")
      ("outload", po::value<std::string>(), "SNP loadings")
      ("outvecx", po::value<std::string>(), "X eigenvector output file, for CCA")
      ("outvecy", po::value<std::string>(), "Y eigenvector output file, for CCA")
      ("outval", po::value<std::string>(), "Eigenvalue output file")
      ("outpve", po::value<std::string>(), "proportion of variance explained output file")
      ("outmeansd", po::value<std::string>(),
	 "mean+SD (used to standardize SNPs) output file")
      ("outproj", po::value<std::string>(), "PCA projection output file")
      ("inload", po::value<std::string>(), "SNP loadings input file")
      ("inmeansd", po::value<std::string>(),
	 "mean+SD (used to standardize SNPs) input file")
      ("inmaf", po::value<std::string>(), "MAF input file")
      ("whiten", "whiten the data")
      ("outwhite", po::value<std::string>(), "whitened data output file")
      ("verbose,v", "verbose")
      ("maxiter", po::value<int>(), "maximum number of randomized PCA iterations")
      ("tol", po::value<double>(), "tolerance for randomized PCA iterations")
      ("transpose", "force a transpose of the data, if possible")
      ("kernel", po::value<std::string>(), "kernel type [rbf | linear]")
      ("sigma", po::value<double>(), "sigma for RBF kernel")
      ("rbfcenter", po::value<std::string>(), "center the RBF kernel [yes | no]")
      ("rbfsample", po::value<int>(), "sample size for estimating RBF kernel width")
      ("savekernel", "save the kernel as text file")
      ("lambda1", po::value<double>(), "1st penalty for CCA/SCCA")
      ("lambda2", po::value<double>(), "2nd penalty for CCA/SCCA")
      ("debug", "debug, dumps all intermediate data (WARNING: slow, call only on small data)")
      ("suffix,f", po::value<std::string>(), "suffix for all output files")
      ("check,c", "check eigenvalues/eigenvectors")
      ("precision", po::value<int>(), "digits of precision for output")
      ("notime", "don't print timestamp in output")
      ("version,V", "version")
   ;

   po::variables_map vm;
   try
   {
      po::store(po::parse_command_line(argc, argv, desc), vm);
      po::notify(vm);
   }
   catch(std::exception& e)
   {
      std::cerr << e.what() << std::endl
	 << "Use --help to get more help"
	 << std::endl;
      return EXIT_SUCCESS;
   }

   show_timestamp = !vm.count("nostamp");

   std::cout << timestamp() << "arguments: flashpca ";
   for(int i = 0 ; i < argc ; i++)
      std::cout << argv[i] << " ";
   std::cout << std::endl;

   if(vm.count("version"))
   {
      std::cerr << "flashpca " << VERSION << std::endl;
      std::cerr << "Copyright (C) 2014-2016 Gad Abraham." << std::endl
	 << "This is free software; see the source for copying conditions.  There is NO"
	 << std::endl
	 << "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE."
	 << std::endl << std::endl;
      return EXIT_SUCCESS;
   }

   if(vm.count("help"))
   {
      std::cerr << "flashpca " << VERSION << std::endl;
      std::cerr << desc << std::endl;
      return EXIT_SUCCESS;
   }


   // Check for which mode we're working in, and verify no conlicting options
   int mode = MODE_PCA;

   std::vector<std::string>
      modes = {"cca", "ucca", "scca", "check", "project"};
   
   if(vm.count("cca"))
   {
      for(int i = 0 ; i < modes.size() ; i++)
      {
	 if(modes[i] != std::string("cca") && vm.count(modes[i]))
	 {
	    std::cerr << "Error: conflicting modes requested: --cca, --"
	       << modes[i] << std::endl
	       << "Use --help to get more help" << std::endl;
	    return EXIT_FAILURE;
	 }
      }
      mode = MODE_CCA;
      std::cerr << "Error: CCA is currently disabled" << std::endl;
      return EXIT_FAILURE;
   }
   else if(vm.count("scca"))
   {
      for(int i = 0 ; i < modes.size() ; i++)
      {
	 if(modes[i] != std::string("scca") && vm.count(modes[i]))
	 {
	    std::cerr << "Error: conflicting modes requested: --scca, --"
	       << modes[i] << std::endl
	       << "Use --help to get more help" << std::endl;
	    return EXIT_FAILURE;
	 }
      }
      mode = MODE_SCCA;
   }
   else if(vm.count("ucca"))
   {
      for(int i = 0 ; i < modes.size() ; i++)
      {
	 if(modes[i] != std::string("ucca") && vm.count(modes[i]))
	 {
	    std::cerr << "Error: conflicting modes requested: --ucca, --"
	       << modes[i] << std::endl
	       << "Use --help to get more help" << std::endl;
	    return EXIT_FAILURE;
	 }
      }
      mode = MODE_UCCA;
   }
   else if(vm.count("check"))
   {
      for(int i = 0 ; i < modes.size() ; i++)
      {
	 if(modes[i] != std::string("check") && vm.count(modes[i]))
	 {
	    std::cerr << "Error: conflicting modes requested: --check, --"
	       << modes[i] << std::endl
	       << "Use --help to get more help" << std::endl;
	    return EXIT_FAILURE;
	 }
      }
      mode = MODE_CHECK_PCA;
   }
   else if(vm.count("project"))
   {
      for(int i = 0 ; i < modes.size() ; i++)
      {
	 if(modes[i] != std::string("project") && vm.count(modes[i]))
	 {
	    std::cerr << "Error: conflicting modes requested: --project, --"
	       << modes[i] << std::endl
	       << "Use --help to get more help" << std::endl;
	    return EXIT_FAILURE;
	 }
      }
      mode = MODE_PREDICT_PCA;
      if(!vm.count("inload"))
      {
	 std::cerr << "Error: SNP-loadings must be specified using --inload"
	    << std::endl;
	 return EXIT_FAILURE;
      }

      if(!vm.count("inmaf") && !vm.count("inmeansd"))
      {
	 std::cerr
	    << "Error: one of MAF or mean/stdev must be specified using "
	    << " --inmaf or --inmeansd, respectively"
	    << std::endl;
	 return EXIT_FAILURE;
      }
   }

   int mem_mode = vm.count("batch") ? MEM_MODE_OFFLINE : MEM_MODE_ONLINE;
   bool fast_mode = !vm.count("rand");
   
   //if(mem_mode == MEM_MODE_ONLINE && !fast_mode)
   //{
   //   std::cerr
   //      << "Error: --online only supported in fast mode (--fast)"
   //      << std::endl;
   //   return EXIT_FAILURE;
   //}

   if(mode == MODE_CHECK_PCA || mode == MODE_PREDICT_PCA)
      mem_mode = MEM_MODE_ONLINE;

   //unsigned long long const MAX_MEM = _8GiB;
   int memory = 8;
   
   if(vm.count("memory"))
   {
      memory = vm["memory"].as<int>();
      if(memory < 1)
      {
	 std::cerr
	    << "Error: memory (GB) must be >=1"
	    << std::endl;
	 return EXIT_FAILURE;
      }
   }

   unsigned int block_size = 0;
   if(vm.count("blocksize"))
   {
      if(vm.count("memory"))
      {
	 std::cerr <<
	    "Error: cannot specify both --memory and --blocksize"
	    << " at the same time" << std::endl;
	 return EXIT_FAILURE;
      }

      block_size = vm["blocksize"].as<int>();
      if(block_size < 1)
      {
	 std::cerr
	    << "Error: blocksize must be >=1"
	    << std::endl;
	 return EXIT_FAILURE;
      }
   }

   int num_threads = 1;
   if(vm.count("numthreads"))
      num_threads = vm["numthreads"].as<int>();

   long seed = 1L;
   if(vm.count("seed"))
      seed = vm["seed"].as<long>();

   std::string fam_file, geno_file, bim_file, pheno_file;

   if(vm.count("bfile"))
   {
      geno_file = vm["bfile"].as<std::string>() + std::string(".bed");
      bim_file = vm["bfile"].as<std::string>() + std::string(".bim");
      fam_file = vm["bfile"].as<std::string>() + std::string(".fam");
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
	 fam_file = vm["fam"].as<std::string>();
      else
	 good = false;
      
      if(!good)
      {
	 std::cerr << "Error: you must specify either --bfile "
	    << "or --bed / --fam / --bim" << std::endl
	    << "Use --help to get more help"
	    << std::endl;
	 return EXIT_FAILURE;
      }
   }

   if(vm.count("pheno"))
      pheno_file = vm["pheno"].as<std::string>();
   else if(mode == MODE_CCA || mode == MODE_UCCA) 
   {
      std::cerr << "Error: you must specify a phenotype file "
	 "in CCA mode using --pheno" << std::endl;
      return EXIT_FAILURE;
   }

   int n_dim = 10;
   if(vm.count("ndim"))
   {
      n_dim = vm["ndim"].as<int>();
      if(n_dim < 1)
      {
	 std::cerr << "Error: --ndim can't be less than 1" << std::endl;
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
   
   int stand_method_x = STANDARDIZE_BINOM2;
   if(vm.count("standx"))
   {
      std::string m = vm["standx"].as<std::string>();
      if(m == "binom")
	 stand_method_x = STANDARDIZE_BINOM;
      else if(m == "binom2")
	 stand_method_x = STANDARDIZE_BINOM2;
      else if(m == "sd")
	 stand_method_x = STANDARDIZE_SD;
      else if(m == "center")
	 stand_method_x = STANDARDIZE_CENTER;
      else if(m == "none")
	 stand_method_x = STANDARDIZE_NONE;
      else
      {
	 std::cerr << "Error: unknown standardization method (--standx): "
	    << m << std::endl;
	 return EXIT_FAILURE;
      }
   }

   if(stand_method_x != STANDARDIZE_BINOM2
      && stand_method_x != STANDARDIZE_BINOM
      && mem_mode == MEM_MODE_ONLINE)
   {
      std::cerr << "Error: in online mode (--online), "
	 << "only --stand binom and binom2 are supported"
	 << std::endl;
      return EXIT_FAILURE;
   }

   int stand_method_y = STANDARDIZE_SD;
   if(vm.count("standy"))
   {
      std::string m = vm["standy"].as<std::string>();
      if(m == "binom")
	 stand_method_x = STANDARDIZE_BINOM;
      else if(m == "binom2")
	 stand_method_x = STANDARDIZE_BINOM2;
      else if(m == "sd")
	 stand_method_x = STANDARDIZE_SD;
      else if(m == "center")
	 stand_method_x = STANDARDIZE_CENTER;
      else if(m == "none")
	 stand_method_x = STANDARDIZE_NONE;
      else
      {
	 std::cerr << "Error: unknown standardization method (--standy): "
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

   std::string suffix = ".txt";
   if(vm.count("suffix"))
      suffix = vm["suffix"].as<std::string>();
   
   std::string pcfile = "pcs" + suffix;
   std::string pcxfile = "pcsX" + suffix;
   std::string pcyfile = "pcsY" + suffix;

   if(vm.count("outpc"))
      pcfile = vm["outpc"].as<std::string>();

   if(vm.count("outpcx"))
      pcxfile = vm["outpcx"].as<std::string>();

   if(vm.count("outpcy"))
      pcyfile = vm["outpcy"].as<std::string>();

   std::string eigvecfile = "eigenvectors" + suffix;
   if(vm.count("outvec"))
      eigvecfile = vm["outvec"].as<std::string>();

   std::string eigvecxfile = "eigenvectorsX" + suffix;
   if(vm.count("outvecx"))
      eigvecxfile = vm["outvecx"].as<std::string>();

   std::string eigvecyfile = "eigenvectorsY" + suffix;
   if(vm.count("outvecy"))
      eigvecyfile = vm["outvecy"].as<std::string>();

   std::string eigvalfile = "eigenvalues" + suffix;
   if(vm.count("outval"))
      eigvalfile = vm["outval"].as<std::string>();

   std::string eigpvefile = "pve" + suffix;
   if(vm.count("outpve"))
      eigpvefile = vm["outpve"].as<std::string>();

   std::string whitefile = "whitened" + suffix;
   if(vm.count("outwhite"))
   {
      std::cerr << "Note: whitening currently disabled" << std::endl;
      whitefile = vm["outwhite"].as<std::string>();
   }

   std::string meansdfile = "meansd" + suffix;
   bool save_meansd = false;
   if(vm.count("outmeansd"))
   {
      meansdfile = vm["outmeansd"].as<std::string>();
      save_meansd = true;
   }

   std::string projfile = "projection" + suffix;
   if(vm.count("outproj"))
      projfile = vm["outproj"].as<std::string>();

   bool whiten = vm.count("whiten");
   bool verbose = vm.count("verbose");
   bool transpose = vm.count("transpose");

   std::cout << "verbose: " << verbose << std::endl;

   bool do_orth = true;
   if(vm.count("orth"))
   {
      std::string s = vm["orth"].as<std::string>();
      if(s == "yes")
	 do_orth = true;
      else if(s == "no")
	 do_orth = false;
      else
      {
	 std::cerr << "Error: unknown option for orth " << s << std::endl;
	 return EXIT_FAILURE;
      }
   }
  
   int maxiter = 500;
   bool debug = vm.count("debug");
   
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

   double tol = 1e-6;
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

   int kernel = KERNEL_LINEAR;
   if(vm.count("kernel"))
   {
      std::string s = vm["kernel"].as<std::string>();
      if(s == "rbf")
	 kernel = KERNEL_RBF;
      else if(s == "linear")
	 kernel = KERNEL_LINEAR;
      else
      {
	 std::cerr << "Error: unknown kernel " << s << std::endl;
	 return EXIT_FAILURE;
      }
   }

   bool do_loadings = false;
   std::string loadingsfile = "";
   if(vm.count("outload"))
   {
      if(kernel == KERNEL_LINEAR)
      {
	 loadingsfile = vm["outload"].as<std::string>();
	 do_loadings = true;
      }
      else
      {
	 std::cerr
	    << "Error: won't save SNP loadings with non-linear kernel "
	    << std::endl;
	 return EXIT_FAILURE;
      }
   }

   double sigma = 0;
   if(vm.count("sigma"))
   {
      sigma = vm["sigma"].as<double>();
      if(sigma <= 0)
      {
	 std::cerr << "Error: --sigma can't be zero or negative"
	    << std::endl;
	 return EXIT_FAILURE;
      }
   }

   bool rbf_center = true;
   if(vm.count("rbfcenter"))
   {
      std::string rc = vm["rbfcenter"].as<std::string>();
      if(rc == "yes")
	 rbf_center = true;
      else if(rc == "no")
	 rbf_center = false;
      else
      {
	 std::cerr << "Error: --rbfcenter must be either 'yes' or 'no'" <<
	    std::endl;
	 return EXIT_FAILURE;
      }
   }

   unsigned int rbf_sample = 1000;
   if(vm.count("rbfsample"))
   {
      rbf_sample = vm["rbfsample"].as<int>();
      if(rbf_sample <= 1)
      {
	 std::cerr << "Error: --rbfsample too small, must be >1" << std::endl;
	 return EXIT_FAILURE;
      }
   }

   bool save_kernel = vm.count("savekernel");

   double lambda1 = 0;
   if(vm.count("lambda1"))
   {
      lambda1 = vm["lambda1"].as<double>();
      if(lambda1 < 0)
      {
	 std::cerr << "Error: --lambda1 can't be negative"
	    << std::endl;
	 return EXIT_FAILURE;
      }
   }

   double lambda2 = 0;
   if(vm.count("lambda2"))
   {
      lambda2 = vm["lambda2"].as<double>();
      if(lambda2 < 0)
      {
	 std::cerr << "Error: --lambda2 can't be negative"
	    << std::endl;
	 return EXIT_FAILURE;
      }
   }

   int divisor = DIVISOR_P;
   if(vm.count("div"))
   {
      std::string m = vm["div"].as<std::string>();
      if(m == "none")
	 divisor = DIVISOR_NONE;
      else if(m == "n1")
	 divisor = DIVISOR_N1;
      else if(m == "p")
	 divisor = DIVISOR_P;
      else
      {
	 std::cerr << "Error: unknown divisor (--div): "
	    << m << std::endl;
	 return EXIT_FAILURE;
      }
   }

   // Options relevant to prediction mode
   std::string in_meansd_file = "";
   std::string in_maf_file = "";
   if(vm.count("inmeansd"))
   {
      if(vm.count("inmaf"))
      {
	 std::cerr <<
	    "Error: conflicting options requested --inmeansd, --inmaf"
	    << std::endl;
	 return EXIT_FAILURE;
      }
	 
      in_meansd_file = vm["inmeansd"].as<std::string>();
      if(in_meansd_file == "")
      {
	 std::cerr << "Error: no file specified for --inmeansd"
	    << std::endl;
	 return EXIT_FAILURE;
      }
   }
   else if(vm.count("inmaf"))
   {
      if(vm.count("inmeansd"))
      {
	 std::cerr <<
	    "Error: conflicting options requested --inmeansd, --inmaf"
	    << std::endl;
	 return EXIT_FAILURE;
      }
	 
      in_maf_file = vm["inmaf"].as<std::string>();
      if(in_maf_file == "")
      {
	 std::cerr << "Error: no file specified for --inmaf"
	    << std::endl;
	 return EXIT_FAILURE;
      }
   }

   std::string in_load_file = "";
   if(vm.count("inload"))
   {
      in_load_file = vm["inload"].as<std::string>();
      if(in_load_file == "")
      {
	 std::cerr << "Error: no file specified for --inload"
	    << std::endl;
	 return EXIT_FAILURE;
      }
   }

   int precision = 7;
   if(vm.count("precision"))
   {
      precision = vm["precision"].as<int>();
      if(precision <= 1)
      {
	 std::cerr << "Error: output --precision too low"
	    << std::endl;
	 return EXIT_FAILURE;
      }
   }

   ////////////////////////////////////////////////////////////////////////////////
   // End command line parsing
      
   std::cout << timestamp() << "Start flashpca (version " << VERSION
      << ")" << std::endl;
#ifdef _OPENMP
#ifdef _EIGEN_HAS_OPENMP
   omp_set_num_threads(num_threads);
   std::cout << timestamp() << "Using " << num_threads 
      << " OpenMP threads" << std::endl;
#endif
#endif

   try
   {
      Data data(seed);
      data.verbose = verbose;
      data.stand_method_x = stand_method_x; //TODO: duplication with RandomPCA
      std::cout << timestamp() << "seed: " << data.seed << std::endl;
   
      if(mode == MODE_CCA || mode == MODE_SCCA || mode == MODE_UCCA)
         data.read_pheno(pheno_file.c_str(), 3);
      else
         data.read_pheno(fam_file.c_str(), 6);
   
      data.read_plink_bim(bim_file.c_str());
      data.read_plink_fam(fam_file.c_str());
         
      data.geno_filename = geno_file.c_str();
      data.get_size();
      //transpose = mode == MODE_PCA && (transpose || data.N > data.nsnps);
      transpose = false;
   
      std::cout << "warning: transpose disabled" << std::endl;
   
      if(mem_mode == MEM_MODE_OFFLINE)
      {
         data.prepare();
         data.read_bed(transpose);
      }
      else if(mem_mode == MEM_MODE_ONLINE)
      {
         data.prepare();
      }
      
      RandomPCA rpca;
      rpca.verbose = verbose;
      rpca.debug = debug;
      rpca.stand_method_x = stand_method_x;
      rpca.stand_method_y = stand_method_y;
      rpca.divisor = divisor;

      // Spectra recommends to run with
      //    1 <= nev < n
      //    nev < ncv < n
      //    ncv >= 2 nev
      // where nev is number of requested eigenvectors, ncv is the extra
      // dimensions required for the computation.
      // see
      // http://yixuan.cos.name/spectra/doc/classSpectra_1_1SymEigsSolver.html
      unsigned int max_dim = fminl(data.N, data.nsnps) / 3.0;
      
      if(n_dim > max_dim)
      {
	 std::cout << timestamp() << "You asked for " 
	    << n_dim << " dimensions, but only "
	    << max_dim << " allowed, using " << max_dim
	    << " instead" << std::endl;
         n_dim = max_dim;
      }
   
      // Only relevant for randomised algorithm, Spectra code uses
      // n_extra = 2 n_dim + 1
      if(n_extra == 0)
         n_extra = fminl(max_dim - n_dim, 200 - n_dim);

      double mem = (double)memory * 1073741824;

      if(block_size == 0)
	 block_size = mem / ((double)data.N * 8.0);

      block_size = fminl(block_size, data.nsnps);

      std::cout << timestamp() << "blocksize: " << block_size 
	 << " (" << (double)block_size * 8.0 * data.N << " bytes)" << std::endl;
   
      ////////////////////////////////////////////////////////////////////////////////
      // The main analysis
      if(mode == MODE_PCA)
      {
         std::cout << timestamp() << "PCA begin" << std::endl;
   
         if(mem_mode == MEM_MODE_OFFLINE)
         {
	    if(fast_mode)
   	    {
	       // New Spectra algorithm
   	       rpca.pca_fast(data.X, block_size, method, transpose, n_dim,
		  n_extra, maxiter, tol, seed, kernel, sigma,
		  rbf_center, rbf_sample, save_kernel, do_orth,
		  do_loadings, mem);
   	    }
   	    else
   	    {
	       // Old randomised algorithm
   	       rpca.pca(data.X, method, transpose, n_dim, n_extra, maxiter,
   	          tol, seed, kernel, sigma, rbf_center, rbf_sample, save_kernel,
   	          do_orth, do_loadings, HIGHMEM);
   	    }
         }
	 else
	 {
	    if(fast_mode)
	    {
	       // New Spectra algorithm
	       rpca.pca_fast(data, block_size, method, transpose,
		  n_dim, n_extra, maxiter, tol, seed, kernel, sigma,
		  rbf_center, rbf_sample, save_kernel, do_orth, do_loadings, mem);
	    }
	    else
	    {
	       // Old randomised algorithm
	       rpca.pca(data, method, transpose, n_dim, n_extra, maxiter,
		  tol, seed, kernel, sigma, rbf_center, rbf_sample, save_kernel,
		  do_orth, do_loadings, block_size);
	    }
	 }
	 std::cout << timestamp() << "PCA done" << std::endl;
      }
      //else if(mode == MODE_CCA)
      //{
      //   std::cout << timestamp() << "CCA begin" << std::endl;
      //   rpca.cca(data.X, data.Y, lambda1, lambda2, seed);
      //   std::cout << timestamp() << "CCA done" << std::endl;
      //}
      else if(mode == MODE_SCCA)
      {
         std::cout << timestamp() << "SCCA begin" << std::endl;
         rpca.scca(data.X, data.Y, lambda1, lambda2, seed, n_dim, mem,
	    maxiter, tol);
         std::cout << timestamp() << "SCCA done" << std::endl;
         save_text(rpca.V0,
	    std::vector<std::string>(),
	    std::vector<std::string>(),
	    std::string("scca_v0.txt").c_str(), precision);
      }
      else if(mode == MODE_UCCA)
      {
         std::cout << timestamp() << "UCCA begin" << std::endl;
         if(mem_mode == MEM_MODE_OFFLINE)
	    rpca.ucca(data.X, data.Y);
         else
	    rpca.ucca(data);
         std::cout << timestamp() << "UCCA done" << std::endl;
      }
      else if(mode == MODE_CHECK_PCA)
      {
	 rpca.check(data, block_size, eigvecfile, eigvalfile);
      }
      else if(mode == MODE_PREDICT_PCA)
      {
	 rpca.project(data, block_size,
	    in_load_file, in_maf_file, in_meansd_file);
      }
      else
      {
	 throw std::runtime_error("Unknown mode");
      }
   
   
      ////////////////////////////////////////////////////////////////////////////////
      // Write out results
      if(mode == MODE_PCA || mode == MODE_SCCA)
      {
         std::cout << timestamp() << "Writing " << n_dim << 
	    " eigenvalues to file " << eigvalfile << std::endl;
         save_text(rpca.d,
	    std::vector<std::string>(),
	    std::vector<std::string>(),
	    eigvalfile.c_str(), precision);
      }
   
      if(mode == MODE_PCA)
      {
         std::cout << timestamp() << "Writing " << n_dim << 
   	 " eigenvectors to file " << eigvecfile << std::endl;
         save_text(rpca.U,
	    std::vector<std::string>(),
   	    std::vector<std::string>(),
   	    eigvecfile.c_str(), precision);
   
         std::cout << timestamp() << "Writing " << n_dim <<
	    " PCs to file " << pcfile << std::endl;
         std::vector<std::string> v(n_dim);
         for(unsigned int i = 0 ; i < n_dim ; i++)
	    v[i] = "PC" + std::to_string(i + 1);
         save_text(rpca.Px, v,	 
	    std::vector<std::string>(),
	    pcfile.c_str(),
	    precision);
   
         std::cout << timestamp() << "Writing " << n_dim << 
   	 " proportion variance explained to file " << eigpvefile << std::endl;
         save_text(rpca.pve,
	    std::vector<std::string>(),
	    std::vector<std::string>(),
	    eigpvefile.c_str(),
	    precision);
   
         if(do_loadings)
         {
	    std::cout << timestamp() << "Writing" <<
	       " SNP loadings to file " << loadingsfile << std::endl;
	    std::cout << rpca.V.rows() << " x " << rpca.V.cols() << std::endl;
	    std::cout << data.snp_ids.size() << std::endl;
	    save_text(rpca.V, v, data.snp_ids, loadingsfile.c_str(), precision); 
         }
      }
      else if(mode == MODE_CCA || mode == MODE_SCCA)
      {
         std::cout << timestamp() << "Writing " << n_dim << 
   	 " X eigenvectors to file " << eigvecxfile << std::endl;
         save_text(rpca.U,
	    std::vector<std::string>(),
   	    std::vector<std::string>(),
   	    eigvecxfile.c_str(), precision);
   
         std::cout << timestamp() << "Writing " << n_dim << 
   	 " Y eigenvectors to file " << eigvecyfile << std::endl;
         save_text(rpca.V,
	    std::vector<std::string>(),
	    std::vector<std::string>(),
	    eigvecyfile.c_str(), precision);
   
         std::cout << timestamp() << "Writing " << n_dim <<
   	 " PCs to file " << pcxfile << std::endl;
         save_text(rpca.Px,
	    std::vector<std::string>(),
   	    std::vector<std::string>(),
   	    pcxfile.c_str(), precision);
   
         std::cout << timestamp() << "Writing " << n_dim <<
   	 " PCs to file " << pcyfile << std::endl;
   
         save_text(rpca.Py,
	    std::vector<std::string>(),
   	    std::vector<std::string>(),
   	    pcyfile.c_str(), precision);
      }
      else if(mode == MODE_UCCA)
      {
         MatrixXd res(rpca.res);
         std::string str[] = {"SNP", "R", "Fstat", "P"};
         std::vector<std::string> v(str, str + 4);
         save_text(res, v, data.snp_ids, std::string("ucca.txt").c_str(),
	    precision);
      }
      else if(mode == MODE_PREDICT_PCA)
      {
         std::vector<std::string> v(rpca.Px.cols());
         for(unsigned int i = 0 ; i < rpca.Px.cols() ; i++)
	    v[i] = "PC" + std::to_string(i + 1);
	 save_text(rpca.Px, v, std::vector<std::string>(),
	    projfile.c_str(), precision);
      }
   
      if(save_meansd)
      {
         std::cout << timestamp() << "Writing mean + sd file "
	    << meansdfile << std::endl;
	 std::vector<std::string> v = {"SNP", "Mean", "SD"};
         save_text(rpca.X_meansd, v, data.snp_ids, meansdfile.c_str(),
	    precision);
      }
   
      ////////////////////////////////////////////////////////////////////////////////
      // Whiten if required
   
      //if(mode == MODE_PCA && whiten)
      //{
      //   std::cout << timestamp() << "ZCA whitening data" << std::endl;
      //   rpca.zca_whiten(transpose);
      //   std::cout << timestamp() << "Writing whitened data to file "
      //      << whitefile << std::endl;
      //   save_text(whitefile.c_str(), rpca.W);
      //}
   
      std::cout << timestamp() << "Goodbye!" << std::endl;
   }
   catch(std::exception& e)
   {
      std::cerr << "Terminating" << std::endl;
      return EXIT_FAILURE;
   }
   catch(...)
   {
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}

