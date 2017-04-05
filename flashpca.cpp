
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

#include "data.h"
#include "randompca.h"

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
      ("scca", "perform sparse canonical correlation analysis (SCCA) [EXPERIMENTAL]")
      ("ucca", "perform per-SNP canonical correlation analysis [EXPERIMENTAL]")
      ("project,p", "project new samples onto existing principal components")
      ("batch", "load all genotypes into RAM at once")
      ("memory,m", po::value<int>(), "size of block, in MB")
      ("blocksize,b", po::value<int>(),
	 "size of block for, in number of SNPs")
      ("numthreads,n", po::value<int>(), "set number of OpenMP threads")
      ("seed", po::value<long>(), "set random seed")
      ("bed", po::value<std::string>(), "PLINK bed file")
      ("bim", po::value<std::string>(), "PLINK bim file")
      ("fam", po::value<std::string>(), "PLINK fam file")
      ("pheno", po::value<std::string>(), "PLINK phenotype file")
      ("bfile", po::value<std::string>(), "PLINK root name")
      ("ndim,d", po::value<int>(), "number of PCs to output")
      ("standx,s", po::value<std::string>(),
	 "standardization method for genotypes [binom2 | binom]")
      ("standy", po::value<std::string>(),
	 "standardization method for phenotypes [sd | binom2 | binom | none | center]")
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
      ("verbose,v", "verbose")
      ("tol", po::value<double>(), "tolerance for PCA iterations")
      ("lambda1", po::value<double>(), "1st penalty for CCA/SCCA")
      ("lambda2", po::value<double>(), "2nd penalty for CCA/SCCA")
      ("maxiter", po::value<int>(), "maximum number of SCCA iterations")
      ("debug", "debug, dumps all intermediate data (WARNING: slow, call only on small data)")
      ("suffix,f", po::value<std::string>(), "suffix for all output files")
      ("check,c", "check eigenvalues/eigenvectors")
      ("precision", po::value<int>(), "digits of precision for output")
      ("notime", "don't print timestamp in output")
      ("save-vinit", "saves the initial v eigenvector for SCCA")
      ("version", "version")
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

   show_timestamp = !vm.count("notime");
   bool save_vinit = vm.count("save-vinit");
   bool verbose = vm.count("verbose");

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

   if(mode == MODE_CHECK_PCA || mode == MODE_PREDICT_PCA)
      mem_mode = MEM_MODE_ONLINE;
   else if(mode == MODE_SCCA) // TODO: SCCA only runs in batch mode currently
      mem_mode = MEM_MODE_OFFLINE;

   int memory = 2048; // Megabytes

   if(vm.count("memory"))
   {
      memory = vm["memory"].as<int>();
      if(memory < 1)
      {
	 std::cerr
	    << "Error: memory (MB) must be >=1"
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
   else if(mode == MODE_CCA || mode == MODE_UCCA || mode == MODE_SCCA)
   {
      std::cerr << "Error: you must specify a phenotype file "
	 "in CCA/UCCA/SCCA mode using --pheno" << std::endl;
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

   int stand_method_x = STANDARDISE_BINOM2;
   if(vm.count("standx"))
   {
      std::string m = vm["standx"].as<std::string>();
      if(m == "binom")
	 stand_method_x = STANDARDISE_BINOM;
      else if(m == "binom2")
	 stand_method_x = STANDARDISE_BINOM2;
      else
      {
	 std::cerr << "Error: unknown standardization method (--standx): "
	    << m << std::endl;
	 return EXIT_FAILURE;
      }
   }

   int stand_method_y = STANDARDISE_SD;
   if(vm.count("standy"))
   {
      std::string m = vm["standy"].as<std::string>();
      if(m == "binom")
	 stand_method_x = STANDARDISE_BINOM;
      else if(m == "binom2")
	 stand_method_x = STANDARDISE_BINOM2;
      else if(m == "sd")
	 stand_method_x = STANDARDISE_SD;
      else if(m == "center")
	 stand_method_x = STANDARDISE_CENTER;
      else if(m == "none")
	 stand_method_x = STANDARDISE_NONE;
      else
      {
	 std::cerr << "Error: unknown standardization method (--standy): "
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

   std::string uccafile = "ucca" + suffix;


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

   bool do_loadings = false;
   std::string loadingsfile = "";
   if(vm.count("outload"))
   {
      loadingsfile = vm["outload"].as<std::string>();
      do_loadings = true;
   }

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
      Data data;
      data.verbose = verbose;
      data.stand_method_x = stand_method_x; //TODO: duplication with RandomPCA
      verbose && std::cout << timestamp() << "seed: " << seed << std::endl;

      if(mode == MODE_CCA || mode == MODE_SCCA || mode == MODE_UCCA)
         data.read_pheno(pheno_file.c_str(), 3);
      else
         data.read_pheno(fam_file.c_str(), 6);

      data.read_plink_bim(bim_file.c_str());
      data.read_plink_fam(fam_file.c_str());

      data.geno_filename = geno_file.c_str();
      data.get_size();

      if(mem_mode == MEM_MODE_OFFLINE)
      {
         data.prepare();
         data.read_bed(false);
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

      //double mem = (double)memory * 1073741824;
      long long mem = (long long)memory * 1048576;

      // Memory required:
      // 0. block of SNPs: block_size * N (N=#samples), not counted here
      // 1. Average + stdev for each SNP, twice (...)
      // 2. Scaled genotypes for each SNP
      // 3. N * ndim left eigenvectors U
      // 4. p * ndim right eigenvectors V (if computing loadings)
      // 5. N/4-sized char buffer + N sized char buffer for, round up to 2x
      // 6. Misc overheads, some auxiliary Spectra vectors of size N?
      if(block_size == 0)
      {
	 //block_size = mem / ((double)data.N * 8.0);
	 long long mem_req_bytes =
	      2 * (long long)data.nsnps * 8 * 2			 // avg+stdev
	    + 3 * (long long)data.nsnps * 8			 // genotypes
	    + (long long) data.N * n_dim * 8			 // left eigenvectors U
	    + (do_loadings ? data.nsnps * n_dim * 8 : 0) // eigenvectors V
	    + 2 * data.N				 // PLINK buffers
	    + 2 * (long long)(data.N + data.nsnps) * n_dim * 8      // Spectra overheads?
	    + 2 * 1024 * 1024 + (long long)data.N * 8;		 // extra space
	 long long mem_remain_bytes = mem - mem_req_bytes;

	 verbose && STDOUT << timestamp() << "mem: "
	    << mem << " mem_req_bytes: " << mem_req_bytes
	    << " mem_remain_bytes: " << mem_remain_bytes << std::endl;

	 if(mem_remain_bytes <= 0)
	 {
	    std::cerr <<
	       "The memory specified using --memory is not sufficient, try"
	       << " increasing it to at least " << (mem_req_bytes + data.N * 8) / 1048576
	       << " MB" << std::endl;
	    return EXIT_FAILURE;
	 }

	 verbose && STDOUT << timestamp() << "Reserving " << mem_req_bytes <<
	    " bytes, leaving " << mem_remain_bytes
	    << " bytes for SNP blocks" << std::endl;

	 block_size = (unsigned int)floor(mem_remain_bytes / ((double)data.N * 8.0));
	 if(block_size < 1)
	 {
	    std::cerr <<
	       "The memory specified using --memory is not sufficient, try"
	       << " increasing it" << std::endl;
	    return EXIT_FAILURE;
	 }
      }

      block_size = fminl(block_size, data.nsnps);

      std::cout << timestamp() << "blocksize: " << block_size
	 << " (" << (long long)block_size * 8 * data.N
	 << " bytes per block)" << std::endl;

      ////////////////////////////////////////////////////////////////////////////////
      // The main analysis
      if(mode == MODE_PCA)
      {
         std::cout << timestamp() << "PCA begin" << std::endl;

         if(mem_mode == MEM_MODE_OFFLINE)
         {
	    // New Spectra algorithm
	    rpca.pca_fast(data.X, block_size, n_dim,
	       maxiter, tol, seed, do_loadings);
         }
	 else
	 {
	    // New Spectra algorithm
	    rpca.pca_fast(data, block_size, 
	       n_dim, maxiter, tol, seed, do_loadings);
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
         //rpca.scca(data.X, data.Y, lambda1, lambda2, seed, n_dim, mem,
	 //   maxiter, tol);
         rpca.scca(data, lambda1, lambda2, seed, n_dim, mem,
	    maxiter, tol, block_size);
         std::cout << timestamp() << "SCCA done" << std::endl;
	 if(save_vinit)
	 {
	    std::cout << timestamp() << "Saving initial V0 vector" << std::endl;
	    save_text(rpca.V0,
	       std::vector<std::string>(),
	       std::vector<std::string>(),
	       std::string("scca_v0.txt").c_str(), precision);
	 }
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

	 std::vector<std::string> rownames(rpca.Px.rows());
	 for(int i = 0 ; i < rpca.Px.rows() ; i++)
	    rownames[i] = data.fam_ids[i] + TXT_SEP + data.indiv_ids[i];

         std::vector<std::string> colnames(rpca.Px.cols() + 1);
	 colnames[0] = std::string("FID") + TXT_SEP + "IID";
         for(int i = 0 ; i < rpca.Px.cols() ; i++)
	    colnames[i + 1] = "U" + std::to_string(i + 1);

         save_text(rpca.U, colnames, rownames, eigvecfile.c_str(), precision);

         std::cout << timestamp() << "Writing " << n_dim <<
	    " PCs to file " << pcfile << std::endl;
         for(int i = 0 ; i < rpca.Px.cols() ; i++)
	    colnames[i + 1] = "PC" + std::to_string(i + 1);

         save_text(rpca.Px, colnames, rownames, pcfile.c_str(), precision);

         std::cout << timestamp() << "Writing " << n_dim
	    << " proportion variance explained to file "
	    << eigpvefile << std::endl;
         save_text(rpca.pve,
	    std::vector<std::string>(),
	    std::vector<std::string>(),
	    eigpvefile.c_str(),
	    precision);

	 // Write out PCA SNP loadings, i.e., the V matrix
         if(do_loadings)
         {
	    std::cout << timestamp() << "Writing" <<
	       " SNP loadings to file " << loadingsfile << std::endl;

	    std::vector<std::string> colnames =
	       {std::string("SNP") + TXT_SEP + "RefAllele"};

	    for(int i = 0 ; i < rpca.V.cols() ; i++)
	       colnames.push_back(std::string("V") + std::to_string(i + 1));

	    std::vector<std::string> rownames(data.snp_ids.size());
	    for(int i = 0 ; i < rownames.size() ; i++)
	       rownames[i] = data.snp_ids[i] + TXT_SEP + data.ref_alleles[i];
	    save_text(rpca.V, colnames, rownames, loadingsfile.c_str(), precision);
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
         save_text(res, v, data.snp_ids, uccafile.c_str(), precision);
      }
      else if(mode == MODE_PREDICT_PCA)
      {
	 std::vector<std::string> rownames(rpca.Px.rows());
	 for(int i = 0 ; i < rpca.Px.rows() ; i++)
	    rownames[i] = data.fam_ids[i] + TXT_SEP + data.indiv_ids[i];

         std::vector<std::string> colnames(rpca.Px.cols() + 1);
	 colnames[0] = std::string("FID") + TXT_SEP + "IID";
         for(int i = 0 ; i < rpca.Px.cols() ; i++)
	    colnames[i + 1] = "PC" + std::to_string(i + 1);

	 save_text(rpca.Px, colnames, rownames, projfile.c_str(), precision);
      }

      if(save_meansd)
      {
         std::cout << timestamp() << "Writing mean + sd file "
	    << meansdfile << std::endl;
	 std::vector<std::string> v =
	    {std::string("SNP") + TXT_SEP + "RefAllele", "Mean", "SD"};
	 std::vector<std::string> rownames(data.snp_ids.size());
	 for(int i = 0 ; i < rownames.size() ; i++)
	    rownames[i] = data.snp_ids[i] + TXT_SEP + data.ref_alleles[i];
         save_text(rpca.X_meansd, v, rownames, meansdfile.c_str(),
	    precision);
      }

      std::cout << timestamp() << "Goodbye!" << std::endl;
   }
   catch(std::exception& e)
   {
      std::cerr << timestamp() << "Exception: " << e.what() << std::endl;
      std::cerr << timestamp() << "Terminating" << std::endl;
      return EXIT_FAILURE;
   }
   catch(...)
   {
      std::cerr << timestamp() << "Caught unknown exception, terminating " << std::endl;
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
