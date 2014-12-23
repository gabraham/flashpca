/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2014 Gad Abraham
 * All rights reserved.
 */

#include <boost/program_options.hpp>

#include <string>
#include <fstream>
#include <sstream>

#include "data.hpp"
#include "randompca.hpp"

namespace po = boost::program_options;

int main(int argc, char * argv[])
{
   po::options_description desc("Options");
   desc.add_options()
      ("help", "produce help message")
      ("bfile", po::value<std::string>(), "PLINK root name")
      ("w", po::value<std::string>(), "SNP weights filename")
      ("stand", po::value<std::string>(), "standardization method [none | binom | sd | center]")
      ("out", po::value<std::string>(), "output filename")
      ("numthreads", po::value<int>(), "set number of OpenMP threads")
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


   int stand_method = STANDARDIZE_BINOM;
   if(vm.count("stand"))
   {
      std::string m = vm["stand"].as<std::string>();
      if(m == "binom")
	 stand_method = STANDARDIZE_BINOM;
      else if(m == "sd")
	 stand_method = STANDARDIZE_SD;
      else if(m == "center")
	 stand_method = STANDARDIZE_CENTER;
      else if(m == "none")
	 stand_method = STANDARDIZE_NONE;
      else
      {
	 std::cerr << "Error: unknown standardization method (--stand): "
	    << m << std::endl;
	 return EXIT_FAILURE;
      }
   }

   std::string geno_file, fam_file;
   if(vm.count("bfile"))
   {
      geno_file = vm["bfile"].as<std::string>() + std::string(".bed");
      fam_file = vm["bfile"].as<std::string>() + std::string(".fam");
      std::cout << ">>> genotype file: " << geno_file << std::endl;
   }
   else
   {
      std::cerr << "Error: bfile not specified" << std::endl;
      return EXIT_FAILURE;
   }

   std::string w_file;
   if(vm.count("w"))
   {
      w_file = vm["w"].as<std::string>();
   }

   std::string out_file = "PredX.txt";
   if(vm.count("out"))
   {
      out_file = vm["out"].as<std::string>();
   }

   omp_set_num_threads(num_threads);
   std::cout << timestamp() << " Using " << num_threads 
      << " OpenMP threads" << std::endl;

   Data data(1);
   data.verbose = true;
   data.read_pheno(fam_file.c_str(), 6);
   data.geno_filename = geno_file.c_str();
   data.get_size();
   data.read_bed(false);

   standardize(data.X, stand_method, false);
   MatrixXd W = data.read_plink_pheno(w_file.c_str(), 1);
   MatrixXd P = data.X * W;

   save_text(out_file.c_str(), P);

   return EXIT_SUCCESS;
}

