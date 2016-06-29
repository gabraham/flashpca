
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

int main(int argc, char * argv[])
{

#ifndef NDEBUG
   std::cout << "******* Running in Debug mode *******" << std::endl;
#endif

   std::cout << timestamp() << " arguments: flashpca ";
   for(int i = 0 ; i < argc ; i++)
      std::cout << argv[i] << " ";
   std::cout << std::endl;

   ////////////////////////////////////////////////////////////////////////////////
   // Parse commandline arguments

   po::options_description desc("Options");
   desc.add_options()
      ("help", "produce help message")
      ("scca", "perform sparse canonical correlation analysis (SCCA)")
      ("ucca", "perform per-SNP canonical correlation analysis")
      ("online", "don't load all genotypes into RAM at once")
      ("fast", "use the fast Spectra algorithm")
      ("blocksize", po::value<int>(), "size of block for online algorith")
      ("numthreads", po::value<int>(), "set number of OpenMP threads")
      ("seed", po::value<long>(), "set random seed")
      ("bed", po::value<std::string>(), "PLINK bed file")
      ("bim", po::value<std::string>(), "PLINK bim file")
      ("fam", po::value<std::string>(), "PLINK fam file")
      ("pheno", po::value<std::string>(), "PLINK phenotype file")
      ("bfile", po::value<std::string>(), "PLINK root name")
      ("ndim", po::value<int>(), "number of PCs to output")
      ("nextra", po::value<int>(),
	 "number of extra dimensions to use in randomized PCA")
      ("standx,stand", po::value<std::string>(),
	 "standardization method for genotypes [binom | binom2 | none | sd | center] (alias: stand)")
      ("standy", po::value<std::string>(),
	 "standardization method for phenotypes [binom | binom2 | none | sd | center]")
      ("method", po::value<std::string>(), "PCA method [eigen | svd]")
      ("orth", po::value<std::string>(), "use orthornormalization [yes | no]")
      ("mem", po::value<std::string>(), "SCCA/PCA method [low | high]")
      ("no_divide_n", "whether to divide X'X by n - 1")
      ("outpc", po::value<std::string>(), "PC output file")
      ("outpcx", po::value<std::string>(), "X PC output file, for CCA")
      ("outpcy", po::value<std::string>(), "Y PC output file, for CCA")
      ("outvec", po::value<std::string>(), "Eigenvector output file")
      ("outload", po::value<std::string>(), "SNP loadings")
      ("outvecx", po::value<std::string>(), "X eigenvector output file, for CCA")
      ("outvecy", po::value<std::string>(), "Y eigenvector output file, for CCA")
      ("outval", po::value<std::string>(), "Eigenvalue output file")
      ("outpve", po::value<std::string>(), "Proportion of variance explained output file")
      ("outmeansd", po::value<std::string>(), "Mean+SD (used to standardize SNPs) output file")
      ("whiten", "whiten the data")
      ("outwhite", po::value<std::string>(), "whitened data output file")
      ("verbose", "verbose")
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
      ("debug", "debug, dumps all intermdiate data (WARNING: slow, call only on small data)")
      ("suffix", po::value<std::string>(), "suffix for all output files")
      ("version", "version")
   ;

   po::variables_map vm;
   po::store(po::parse_command_line(argc, argv, desc), vm);
   po::notify(vm);

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

   int mode = MODE_PCA;
   if(vm.count("cca"))
   {
      mode = MODE_CCA;
      std::cerr << "Error: CCA is currently disabled" << std::endl;
      return EXIT_FAILURE;
   }
   else if(vm.count("scca"))
      mode = MODE_SCCA;
   else if(vm.count("ucca"))
      mode = MODE_UCCA;

   int mem_mode = vm.count("online") ? MEM_MODE_ONLINE : MEM_MODE_OFFLINE;
   bool fast_mode = vm.count("fast");

   int block_size = 100;
   if(vm.count("blocksize"))
   {
      block_size = vm["blocksize"].as<int>();
      if(block_size <= 1)
      {
	 std::cerr
	    << "Error: blocksize must be >=1 and <= num_snps"
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
   
   int stand_method_x = STANDARDIZE_BINOM;
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

   int mem = LOWMEM;
   if(vm.count("mem"))
   {
      std::string m = vm["mem"].as<std::string>();
      if(m == "low")
	 mem = LOWMEM;
      else if(m == "high")
	 mem = HIGHMEM;
      else
      {
	 std::cerr << "Error: unknown argument (--mem): "
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
      whitefile = vm["outwhite"].as<std::string>();

   std::string meansdfile = "";
   bool save_meansd = false;
   if(vm.count("outmeansd"))
   {
      meansdfile = vm["outmeansd"].as<std::string>();
      save_meansd = true;
   }

   bool whiten = vm.count("whiten");
   bool verbose = vm.count("verbose");
   bool transpose = vm.count("transpose");

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

   double tol = 1e-4;
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

   bool divide_n = !vm.count("no_divide_n");

   ////////////////////////////////////////////////////////////////////////////////
   // End command line parsing
      
   std::cout << timestamp() << " Start flashpca (version " << VERSION
      << ")" << std::endl;
#ifdef _OPENMP
#ifdef _EIGEN_HAS_OPENMP
   omp_set_num_threads(num_threads);
   std::cout << timestamp() << " Using " << num_threads 
      << " OpenMP threads" << std::endl;
#endif
#endif

   Data data(seed);
   data.verbose = verbose;
   std::cout << timestamp() << " seed: " << data.seed << std::endl;

   if(mode == MODE_CCA || mode == MODE_SCCA || mode == MODE_UCCA)
      data.read_pheno(pheno_file.c_str(), 3);
   else
      data.read_pheno(fam_file.c_str(), 6);

   data.read_plink_bim(bim_file.c_str());
      
   data.geno_filename = geno_file.c_str();
   data.get_size();
   transpose = mode == MODE_PCA && (transpose || data.N > data.nsnps);

   if(mem_mode == MEM_MODE_OFFLINE)
   {
      data.prepare();
      data.read_bed(transpose);
      data.finish();
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
   unsigned int max_dim = fminl(data.X.rows(), data.X.cols());
   
   if(n_dim == 0)
      n_dim = fminl(max_dim, 10);

   if(n_extra == 0)
      n_extra = fminl(max_dim - n_dim, 190);

   ////////////////////////////////////////////////////////////////////////////////
   // The main analysis
   if(mode == MODE_PCA)
   {
      std::cout << timestamp() << " PCA begin" << std::endl;

      if(mem_mode == MEM_MODE_OFFLINE)
      {
	 if(fast_mode)
	 {
	    rpca.pca_fast(data.X, block_size, method, transpose, n_dim, n_extra, maxiter,
	       tol, seed, kernel, sigma, rbf_center, rbf_sample, save_kernel,
	       do_orth, do_loadings, mem, divide_n);
	 }
	 else
	 {
	    rpca.pca(data.X, method, transpose, n_dim, n_extra, maxiter,
	       tol, seed, kernel, sigma, rbf_center, rbf_sample, save_kernel,
	       do_orth, do_loadings, mem, divide_n);
	 }
      }
      else
      {
	 if(fast_mode)
	 {
	    rpca.pca_fast(data, block_size, method, transpose, n_dim, n_extra, maxiter,
	       tol, seed, kernel, sigma, rbf_center, rbf_sample, save_kernel,
	       do_orth, do_loadings, mem, divide_n);
	 }
	 else
	 {
	    rpca.pca(data, method, transpose, n_dim, n_extra, maxiter,
	       tol, seed, kernel, sigma, rbf_center, rbf_sample, save_kernel,
	       do_orth, do_loadings, mem, divide_n);
	 }
      }
      std::cout << timestamp() << " PCA done" << std::endl;
   }
   //else if(mode == MODE_CCA)
   //{
   //   std::cout << timestamp() << " CCA begin" << std::endl;
   //   rpca.cca(data.X, data.Y, lambda1, lambda2, seed);
   //   std::cout << timestamp() << " CCA done" << std::endl;
   //}
   else if(mode == MODE_SCCA)
   {
      std::cout << timestamp() << " SCCA begin" << std::endl;
      rpca.scca(data.X, data.Y, lambda1, lambda2, seed, n_dim, mem,
	 maxiter, tol);
      std::cout << timestamp() << " SCCA done" << std::endl;
   }
   else if(mode == MODE_UCCA)
   {
      std::cout << timestamp() << " UCCA begin" << std::endl;
      if(mem_mode == MEM_MODE_OFFLINE)
	 rpca.ucca(data.X, data.Y);
      else
	 rpca.ucca(data);
      std::cout << timestamp() << " UCCA done" << std::endl;
   }


   ////////////////////////////////////////////////////////////////////////////////
   //  Close open files if in online mode
   if(mem_mode == MEM_MODE_ONLINE)
   {
      data.finish();
   }

   ////////////////////////////////////////////////////////////////////////////////
   // Write out results


   // Common to all decompositions, except UCCA
   if(mode != MODE_UCCA)
   {
      std::cout << timestamp() << " Writing " << n_dim << 
	 " eigenvalues to file " << eigvalfile << std::endl;
      save_text(rpca.d,
	 std::vector<std::string>(),
	 std::vector<std::string>(),
	 eigvalfile.c_str());
   }

   if(mode == MODE_PCA)
   {
      std::cout << timestamp() << " Writing " << n_dim << 
	 " eigenvectors to file " << eigvecfile << std::endl;
      save_text(rpca.U,
	 std::vector<std::string>(),
	 std::vector<std::string>(),
	 eigvecfile.c_str());

      std::cout << timestamp() << " Writing " << n_dim <<
	 " PCs to file " << pcfile << std::endl;
      std::vector<std::string> v(n_dim);
      for(unsigned int i = 0 ; i < n_dim ; i++)
	 v[i] = "PC" + std::to_string(i + 1);
      save_text(rpca.Px, v, std::vector<std::string>(), pcfile.c_str());

      std::cout << timestamp() << " Writing " << n_dim << 
	 " proportion variance explained to file " << eigpvefile << std::endl;
      save_text(rpca.pve,
	 std::vector<std::string>(),
	 std::vector<std::string>(),
	 eigpvefile.c_str());

      if(do_loadings)
      {
	 std::cout << timestamp() << " Writing" <<
	    " SNP loadings to file " << loadingsfile << std::endl;
	 save_text(rpca.V,
	    std::vector<std::string>(),
	    std::vector<std::string>(),
	    loadingsfile.c_str()); 
      }
   }
   else if(mode == MODE_CCA || mode == MODE_SCCA)
   {
      std::cout << timestamp() << " Writing " << n_dim << 
	 " X eigenvectors to file " << eigvecxfile << std::endl;
      save_text(rpca.U,
	 std::vector<std::string>(),
	 std::vector<std::string>(),
	 eigvecxfile.c_str());

      std::cout << timestamp() << " Writing " << n_dim << 
	 " Y eigenvectors to file " << eigvecyfile << std::endl;
      save_text(rpca.V,
	 std::vector<std::string>(),
	 std::vector<std::string>(),
	 eigvecyfile.c_str());

      std::cout << timestamp() << " Writing " << n_dim <<
	 " PCs to file " << pcxfile << std::endl;
      save_text(rpca.Px,
	 std::vector<std::string>(),
	 std::vector<std::string>(),
	 pcxfile.c_str());

      std::cout << timestamp() << " Writing " << n_dim <<
	 " PCs to file " << pcyfile << std::endl;

      save_text(rpca.Py,
	 std::vector<std::string>(),
	 std::vector<std::string>(),
	 pcyfile.c_str());
   }
   else if(mode == MODE_UCCA)
   {
      MatrixXd res(rpca.res);
      std::string str[] = {"SNP", "R", "Fstat", "P"};
      std::vector<std::string> v(str, str + 4);
      save_text(res, v, data.snp_ids, std::string("ucca.txt").c_str());
   }

   if(save_meansd)
   {
      std::cout << timestamp() << " Writing mean + sd file "
	 << meansdfile << std::endl;
      save_text(rpca.X_meansd,
	 std::vector<std::string>(),
	 std::vector<std::string>(),
	 meansdfile.c_str());
   }

   ////////////////////////////////////////////////////////////////////////////////
   // Whiten if required

   //if(mode == MODE_PCA && whiten)
   //{
   //   std::cout << timestamp() << " ZCA whitening data" << std::endl;
   //   rpca.zca_whiten(transpose);
   //   std::cout << timestamp() << " Writing whitened data to file "
   //      << whitefile << std::endl;
   //   save_text(whitefile.c_str(), rpca.W);
   //}

   std::cout << timestamp() << " Goodbye!" << std::endl;

   return EXIT_SUCCESS;
}

