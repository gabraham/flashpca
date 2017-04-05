/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2014-2016 Gad Abraham
 * All rights reserved.
 */


#include "data.h"

Data::Data()
{
   N = 0;
   p = 0;
   K = 0;
   nsnps = 0;
   visited = NULL;
   tmp = NULL;
   tmp2 = NULL;
   avg = NULL;
   verbose = false;
   use_preloaded_maf = false;
}

Data::~Data()
{
   if(visited)
      delete[] visited;
   if(tmp)
      delete[] tmp;
   if(tmp2)
      delete[] tmp2;
   if(avg)
      delete[] avg;
   in.close();
}

/* 
 *                   plink BED           sparsnp
 * minor homozyous:  00 => numeric 0     10 => numeric 2
 * heterozygous:     10 => numeric 2     01 => numeric 1
 * major homozygous: 11 => numeric 3     00 => numeric 0
 * missing:          01 => numeric 1     11 => numeric 3
 *
 *
 * http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml says,
 * The bytes in plink are read backwards HGFEDCBA, not GHEFCDAB, but we read
 * them forwards as a character (a proper byte)
 *
 * By default, plink usage dosage of the *major* allele, since allele A1 is
 * usually the minor allele and the code "1" refers to the second allele A2,
 * so that "11" is A2/A2 or major/major.
 *
 * We always use minor allele dosage, to be consistent with the output from
 * plink --recodeA which used minor allele dosage by default.
 *
 * out: array of genotypes
 * in: array of packed genotypes (bytes)
 * n: number of bytes in input
 * 
 */
void decode_plink(unsigned char * __restrict__ out,
   const unsigned char * __restrict__ in,
   const unsigned int n)
{
   unsigned int i, k;
   unsigned char tmp, geno1, geno2, geno3, geno4;
   unsigned int a1, a2;

   for(i = 0 ; i < n ; i++)
   {
      tmp = in[i];
      k = PACK_DENSITY * i;

      /* geno is interpreted as a char, however a1 and a2 are bits for allele 1 and
       * allele 2. The final genotype is the sum of the alleles, except for 01
       * which denotes missing.
       */

      geno1 = (tmp & MASK0);
      if(geno1 == 1)
	 out[k] = 3;
      else
      {
	 a1 = !(geno1 & 1);
	 a2 = !(geno1 >> 1);
	 out[k] = a1 + a2;
      }
      k++;

      geno2 = (tmp & MASK1) >> 2;
      if(geno2 == 1)
	 out[k] = 3;
      else
      {
	 a1 = !(geno2 & 1);
	 a2 = !(geno2 >> 1);
	 out[k] = a1 + a2;
      }
      k++;

      geno3 = (tmp & MASK2) >> 4;
      if(geno3 == 1)
	 out[k] = 3;
      else
      {
	 a1 = !(geno3 & 1);
	 a2 = !(geno3 >> 1);
	 out[k] = a1 + a2;
      }
      k++;

      geno4 = (tmp & MASK3) >> 6;
      if(geno4 == 1)
	 out[k] = 3;
      else
      {
	 a1 = !(geno4 & 1);
	 a2 = !(geno4 >> 1);
	 out[k] = a1 + a2;
      }
   }
}

void decode_plink_simple(unsigned char * __restrict__ out,
   const unsigned char * __restrict__ in,
   const unsigned int n)
{
   unsigned int i, k;

   for(i = 0 ; i < n ; i++)
   {
      k = PACK_DENSITY * i;

      /* geno is interpreted as a char, however a1 and a2 are bits for allele 1 and
       * allele 2. The final genotype is the sum of the alleles, except for 01
       * which denotes missing.
       */

      out[k] =   (in[i] & MASK0);
      out[k+1] = (in[i] & MASK1) >> 2;
      out[k+2] = (in[i] & MASK2) >> 4;
      out[k+3] = (in[i] & MASK3) >> 6;
   }
}

void Data::get_size()
{
   verbose && STDOUT << timestamp() << "Analyzing BED file '" 
      << geno_filename << "'";
   std::ifstream in(geno_filename, std::ios::in | std::ios::binary);

   if(!in)
   {
      std::string err = std::string("[Data::read_bed] Error reading file ")
	  + geno_filename + ", error " + strerror(errno);
      throw std::runtime_error(err);
   }

   in.seekg(0, std::ifstream::end);

   // file size in bytes, ignoring first 3 bytes (2byte magic number + 1byte mode)
   len = (unsigned long long)in.tellg() - 3;

   // size of packed data, in bytes, per SNP
   np = (unsigned long long)ceil((double)N / PACK_DENSITY);
   nsnps = (unsigned int)(len / np);
   in.seekg(3, std::ifstream::beg);
   in.close();

   verbose && STDOUT << ", found " << (len + 3) << " bytes, "
      << nsnps << " SNPs" << std::endl;
}

// Prepare input stream etc before reading in SNP blocks
void Data::prepare()
{
   in.open(geno_filename, std::ios::in | std::ios::binary);
   in.seekg(3, std::ifstream::beg);

   if(!in)
   {
      std::string err = std::string("[Data::read_bed] Error reading file ")
	 + geno_filename;
      throw std::runtime_error(err);
   }

   tmp = new unsigned char[np];

   // Allocate more than the sample size since data must take up whole bytes
   tmp2 = new unsigned char[np * PACK_DENSITY];

   avg = new double[nsnps](); 
   visited = new bool[nsnps]();
   X_meansd = MatrixXd::Zero(nsnps, 2); // TODO: duplication here with avg

   scaled_geno_lookup = ArrayXXd::Zero(4, nsnps);

   verbose && STDOUT << timestamp() << "Detected BED file: "
      << geno_filename << " with " << (len + 3)
      << " bytes, " << N << " samples, " << nsnps 
      << " SNPs." << std::endl;
}

// Reads a _contiguous_ block of SNPs [start, stop] at a time.
// The block will contain standardised genotypes already, no need to
// standardise them again.
//
// If resize is false, then the calling code is responsible for ensuring that
// X is handled accordingly with the block size (X may be bigger than the
// block).
void Data::read_snp_block(unsigned int start_idx, unsigned int stop_idx,
   bool transpose, bool resize)
{
   in.seekg(3 + np * start_idx);

   unsigned int actual_block_size = stop_idx - start_idx + 1;

   // Resize matrix, e.g., with final block that may be smaller than
   // $blocksize$
   if(transpose)
   {
      if(X.rows() == 0 || (resize && X.rows() != actual_block_size))
      {
         verbose && STDOUT << timestamp()
	    << "Reallocating memory: " << X.rows() << " -> " <<
	    actual_block_size << std::endl;
         if(X.rows() > actual_block_size)
	 {
            X = MatrixXd(actual_block_size, N);
	 }
      }
   }
   else if(X.cols() == 0 || (resize && X.cols() != actual_block_size))
   {
      verbose && STDOUT << timestamp()
	 << "Reallocating memory: " << X.cols() << " -> " <<
	 actual_block_size << std::endl;
      X = MatrixXd(N, actual_block_size);
   }

   for(unsigned int j = 0; j < actual_block_size; j++)
   {
      unsigned int k = start_idx + j;

      // read raw genotypes
      in.read((char*)tmp, sizeof(char) * np);

      // Compute average per SNP, excluding missing values
      double snp_avg = 0;
      unsigned int ngood = 0;

      // We've seen this SNP, don't need to compute its average again
      if(!visited[k])
      {
	 // decode the genotypes and convert to 0/1/2/NA
	 decode_plink(tmp2, tmp, np);

	 double P, sd;

	 if(!use_preloaded_maf)
	 {
	    for(unsigned int i = 0 ; i < N ; i++)
      	    {
      	       double s = (double)tmp2[i];
      	       if(tmp2[i] != PLINK_NA)
      	       {
      	          snp_avg += s;
      	          ngood++;
      	       }
      	    }
      	    snp_avg /= ngood;

	    // Store the 4 possible standardised genotypes for each SNP
	    P = snp_avg / 2.0;
	    if(stand_method_x == STANDARDISE_BINOM)
	       sd = sqrt(P * (1 - P));
	    else if(stand_method_x == STANDARDISE_BINOM2)
	       sd = sqrt(2.0 * P * (1 - P));
	    else
	    {
	       std::string err = std::string("unknown standardisation method: ")
		  + std::to_string(stand_method_x);
	       throw std::runtime_error(err);
	    }

	    X_meansd(k, 0) = snp_avg;
	    X_meansd(k, 1) = sd;
	 }
	 else
	 {
	    snp_avg = X_meansd(k, 0);
	    sd = X_meansd(k, 1);
	 }

	 // scaled genotyped initialised to zero
	 if(sd > VAR_TOL)
	 {
	    // Note thet scaled values for the genotypes are stored based on
	    // the PLINK indexing rather than the actual dosage indexing,
	    // which lets us later just read the PLINK data and not have to
	    // convert the dosages.
	    //
	    // plink '3' -> actual dosage '0'
	    // plink '2' -> actual dosage '1'
	    // plink '0' -> actual dosage '2'
	    // plink '1' -> actual dosage '3' (NA)
	    //*                   plink BED           sparsnp
 	    //* minor homozyous:  00 => numeric 0     10 => numeric 2
 	    //* heterozygous:     10 => numeric 2     01 => numeric 1
 	    //* major homozygous: 11 => numeric 3     00 => numeric 0
 	    //* missing:          01 => numeric 1     11 => numeric 3
	    scaled_geno_lookup(3, k) = (0 - snp_avg) / sd;
	    scaled_geno_lookup(2, k) = (1 - snp_avg) / sd;
	    scaled_geno_lookup(0, k) = (2 - snp_avg) / sd;
	    scaled_geno_lookup(1, k) = 0; // impute to average
	 }
	 visited[k] = true;
      }

      // Unpack the genotypes, but don't convert to 0/1/2/NA, keep in
      // original format (see comments for decode_plink).
      // There is a bit of waste here in the first time the SNP is visited, as
      // we unpack the data twice, once with decoding and once without.
      decode_plink_simple(tmp2, tmp, np);

      for(unsigned int i = 0 ; i < N ; i++)
      {
	 X(i, j) = scaled_geno_lookup(tmp2[i], k);
      }
   }
}

// Reads an entire bed file into memory
// Expects PLINK bed in SNP-major format
void Data::read_bed(bool transpose)
{
   if(transpose)
      X = MatrixXd(nsnps, N);
   else
      X = MatrixXd(N, nsnps);

   unsigned int md = nsnps / 50;

   // iterate over all SNPs
   for(unsigned int j = 0 ; j < nsnps; j++)
   {
      // read raw genotypes
      in.read((char*)tmp, sizeof(char) * np);

      // decode the genotypes
      decode_plink(tmp2, tmp, np);

      // Compute average per SNP, excluding missing values
      avg[j] = 0;
      unsigned int ngood = 0;
      for(unsigned int i = 0 ; i < N ; i++)
      {
	 double s = (double)tmp2[i];
	 if(tmp2[i] != PLINK_NA)
	 {
	    avg[j] += s;
	    ngood++;
	 }
      }
      avg[j] /= ngood;

      // Impute using average per SNP
      for(unsigned int i = 0 ; i < N ; i++)
      {
	 double s = (double)tmp2[i];
	 {
	    if(transpose)
	    {
	       if(s != PLINK_NA)
		  X(j, i) = s;
	       else
		  X(j, i) = avg[j];
	    }
	    else
	    {
	       if(s != PLINK_NA)
		  X(i, j) = s;
	       else
		  X(i, j) = avg[j];
	    }
	 }
      }

      if(verbose && j % md == md - 1)
	 STDOUT << timestamp() << "Reading genotypes, "
	    << roundl(((double)j / nsnps) * 100) << "% done" 
	    << std::endl;
   }

   if(transpose)
      p = X.rows();
   else
      p = X.cols();

   verbose && STDOUT << timestamp() << "Loaded genotypes: "
      << N << " samples, " << p << " SNPs" << std::endl;
}

void Data::read_pheno(const char *filename, unsigned int firstcol)
{
    NamedMatrixWrapper M = read_text(filename, firstcol);
    Y = M.X;
    N = M.X.rows();
}

// Reads PLINK phenotype files:
// FID IID pheno1 pheno2 ...
// Need to be able to read continuous phenotypes
// 
// firstcol is _one-based_, 3 for pheno file, 6 for FAM file (ignoring sex),
// 5 for FAM file (with gender)
NamedMatrixWrapper read_text(const char *filename,
   unsigned int firstcol, unsigned int nrows, unsigned int skip,
   bool verbose)
{
   NamedMatrixWrapper M;

   unsigned int line_num = 0;

   std::ifstream in(filename, std::ios::in);

   if(!in)
   {
      std::string err = std::string("Error reading file '")
	 + filename + "': " + strerror(errno);
      throw std::runtime_error(err);
   }
   std::vector<std::string> lines;
   
   while(in)
   {
      std::string line;
      std::getline(in, line);
      if(!in.eof() && (nrows == -1 || line_num < nrows))
      {
	 if(line_num >= skip)
	    lines.push_back(line);
	 line_num++;
      }
   }

   verbose && STDOUT << timestamp() << "Detected text file " <<
      filename << ", " << lines.size() << " rows" << std::endl;

   in.close();

   unsigned int numtok = 0, numfields, numfields_1st = 0;

   M.X = MatrixXd(0, 0);

   for(unsigned int i = 0 ; i < lines.size() ; i++)
   {
      std::stringstream ss(lines[i]);
      std::string s;
      std::vector<std::string> tokens;

      while(ss >> s)
	 tokens.push_back(s);

      numtok = tokens.size();
      numfields = numtok - firstcol + 1;
      if(i == 0)
      {
	 M.X.resize(lines.size(), numfields);
	 numfields_1st = numfields;
      }
      else if(numfields_1st != numfields)
      {
	 std::string err = std::string("Error reading file '")
	    + filename + "': inconsistent number of columns";
	 throw std::runtime_error(err);
      }

      VectorXd y(numfields);
      char* err;
      errno = 0;
      for(unsigned int j = 0 ; j < numfields ; j++)
      {
	 //y(j) = std::atof(tokens[j + firstcol - 1].c_str());
	 double m = std::strtod(tokens[j + firstcol - 1].c_str(), &err);
	 if(*err != '\0' || errno != 0)
	 {
	    std::string err = std::string("Error reading file '")
	       + filename + "', line " + std::to_string(i + 1)
	       + ": '" + tokens[j + firstcol - 1] + "'"
	       + " cannot be parsed as a number";
	    throw std::runtime_error(err);
	 }
	 y(j) = m;
      }
      M.X.row(i) = y;
   }

   return M;
}

void Data::read_plink_bim(const char *filename)
{
   std::ifstream in(filename, std::ios::in);

   if(!in)
   {
      std::string err = std::string("Error reading file ")
	 + filename;
      throw std::runtime_error(err);
   }
   std::vector<std::string> lines;
   
   while(in)
   {
      std::string line;
      std::getline(in, line);
      if(!in.eof())
	 lines.push_back(line);
   }

   verbose && STDOUT << timestamp() << "Detected bim file " <<
      filename << ", " << lines.size() << " SNPs" << std::endl;
   in.close();

   for(unsigned int i = 0 ; i < lines.size() ; i++)
   {
      std::stringstream ss(lines[i]);
      std::string s;
      std::vector<std::string> tokens;

      while(ss >> s)
	 tokens.push_back(s);
      snp_ids.push_back(tokens[1]);
      ref_alleles.push_back(tokens[4]);
      alt_alleles.push_back(tokens[5]);

      char* err;
      errno = 0;
      unsigned long long m = std::strtol(tokens[3].c_str(), &err, 10);
      if(*err != '\0' || errno != 0)
      {
	 std::string err = std::string("Error reading file '")
	    + filename + "', line " + std::to_string(i + 1)
	    + ": '" + tokens[3] + "' cannot be parsed as a number";
	 throw std::runtime_error(err);
      }
      bp.push_back(m);
   }
}

void Data::read_plink_fam(const char *filename)
{
   std::ifstream in(filename, std::ios::in);

   if(!in)
   {
      std::string err = std::string(
	 "[Data::read_plink_fam] Error reading file ") + filename;
      throw std::runtime_error(err);
   }
   std::vector<std::string> lines;
   
   while(in)
   {
      std::string line;
      std::getline(in, line);
      if(!in.eof())
	 lines.push_back(line);
   }

   in.close();

   for(unsigned int i = 0 ; i < lines.size() ; i++)
   {
      std::stringstream ss(lines[i]);
      std::string s;
      std::vector<std::string> tokens;

      while(ss >> s)
	 tokens.push_back(s);
      fam_ids.push_back(tokens[0]);
      indiv_ids.push_back(tokens[1]);
   }
}

std::string Data::tolower(const std::string& v)
{
   std::string r = v;
   std::transform(r.begin(), r.end(), r.begin(), ::tolower);
   return r;
}

