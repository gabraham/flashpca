
#include "data.hpp"

Data::Data(long seed)
{
   N = 0;
   p = 0;
   K = 0;
   nsnps = 0;
   //cache = NULL;
   this->seed = seed;
   srand48(seed);
}

Data::~Data()
{
   geno_fin.close();
 //  delete cache;
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
void decode_plink(unsigned char *out,
      const unsigned char *in, const unsigned int n)
{
   unsigned int i, k;
   unsigned char tmp, geno;
   unsigned int a1, a2;

   for(i = 0 ; i < n ; ++i)
   {
      tmp = in[i];
      k = PACK_DENSITY * i;
      
      /* geno is interpreted as a char, however a1 and a2 are bits for allele 1 and
       * allele 2. The final genotype is the sum of the alleles, except for 01
       * which denotes missing.
       */
      geno = (tmp & MASK0);
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : a1 + a2;
      k++;

      geno = (tmp & MASK1) >> 2; 
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : a1 + a2;
      k++;

      geno = (tmp & MASK2) >> 4; 
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : a1 + a2;
      k++;

      geno = (tmp & MASK3) >> 6; 
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : a1 + a2;
   }
}

// Expects PLINK BED in SNP-major format
void Data::read_bed(const char *filename)
{
   std::cout << timestamp() << " Reading BED file '" << filename << "'" << std::endl;
   std::ifstream in(filename, std::ios::in | std::ios::binary);

   if(!in)
   {
      std::cerr << "[Data::read_bed] Error reading file " << filename << std::endl;
      throw std::runtime_error("io error");
   }

   in.seekg(0, std::ifstream::end);

   // file size in bytes, ignoring first 3 bytes (2byte magic number + 1byte mode)
   len = (unsigned int)in.tellg() - 3;

   // size of packed data, in bytes, per SNP
   np = (unsigned int)ceil((double)N / PACK_DENSITY);
   nsnps = len / np;
   in.seekg(3, std::ifstream::beg);

   unsigned char* tmp = new unsigned char[np];

   // Allocate more than the sample size since data must take up whole bytes
   unsigned char* tmp2 = new unsigned char[np * PACK_DENSITY];

   X = MatrixXd(N, nsnps);

   std::cout << timestamp() << " Detected BED file: " << filename <<
      " with " << len << " bytes, " << N << " samples, " << nsnps 
      << " SNPs." << std::endl;
   VectorXd tmp3(N);

   double* avg = new double[nsnps]; 

   unsigned int idx = 0;

   // iterate over all SNPs, only decode those that are included in analysis
   for(unsigned int i = 0 ; i < nsnps; i++)
   {
      // read raw genotypes
      in.read((char*)tmp, sizeof(char) * np);

      if(!snps[i].included)
	 continue;

      // decode the genotypes
      decode_plink(tmp2, tmp, np);

      // Compute average per SNP, excluding missing values
      avg[idx] = 0;
      unsigned int ngood = 0;
      for(unsigned int j = 0 ; j < N ; j++)
      {
	 double s = (double)tmp2[j];
	 if(s != PLINK_NA)
	 {
	    avg[idx] += s;
	    ngood++;
	 }
      }
      avg[idx] /= ngood;

      // Impute using average per SNP
      for(unsigned int j = 0 ; j < N ; j++)
      {
	 double s = (double)tmp2[j];
	 if(s != PLINK_NA)
	    tmp3(j) = s;
	 else
	    tmp3(j) = avg[idx];
      }

      X.col(idx) = tmp3;
      idx++;
   }

   p = X.cols();

   std::cout << timestamp() << " Loaded genotypes: " << X.rows() << " samples, " <<
      X.cols() << " SNPs post-region-removal" << std::endl;

   delete[] tmp;
   delete[] tmp2;
   delete[] avg;

   in.close();
}

void Data::read_pheno(const char *filename, unsigned int firstcol,
   int pheno)
{
   Y = read_plink_pheno(filename, firstcol, pheno);
}

// Reads PLINK phenotype files:
// FID IID pheno1 pheno2 ...
// Need to be able to read continuous phenotypes
// 
// firstcol is one-based, 3 for pheno file, 6 for FAM file (ignoring gender),
// 5 for FAM file (with gender)
MatrixXd Data::read_plink_pheno(const char *filename, unsigned int firstcol,
   int pheno)
{
   std::ifstream in(filename, std::ios::in);

   if(!in)
   {
      std::cerr << "[Data::read_plink_pheno] Error reading file " 
	 << filename << std::endl;
      throw std::string("io error");
   }
   std::vector<std::string> lines;
   
   while(in)
   {
      std::string line;
      std::getline(in, line);
      if(!in.eof())
	 lines.push_back(line);
   }

   std::cout << timestamp() << " Detected pheno file " <<
      filename << ", " << lines.size() << " samples";

   in.close();

   unsigned int numtok = 0, numfields;

   MatrixXd Z(0, 0);

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
	 Z.resize(lines.size(), numfields);

      VectorXd y(numfields);
      for(unsigned int j = 0 ; j < numfields ; j++)
	 y(j) = std::atof(tokens[j + firstcol - 1].c_str());
      Z.row(i) = y;
   }

   std::cout << ", " << Z.cols() << " columns (ex. FAM+INDIV IDs)" << std::endl;
   
   N = Z.rows();

   return Z;
}

std::string Data::tolower(const std::string& v)
{
   std::string r = v;
   std::transform(r.begin(), r.end(), r.begin(), ::tolower);
   return r;
}

//void Data::mmap_bed(char *filename)
//{
//   geno_filename = filename;
//   boost::iostreams::mapped_file_params params;
//   params.path = std::string(filename, strlen(filename));
//   params.offset = 0;
//   params.length = -1;
//   geno_fin = boost::iostreams::mapped_file_source(params);
//
//   if(N == 0)
//      throw std::runtime_error(
//	 "haven't read a FAM/PHENO file so don't know what sample size is");
//
//   len = boost::filesystem::file_size(geno_filename) - PLINK_OFFSET;
//   std::cout << geno_filename << " len: " << len << " bytes" << std::endl;
//   np = (unsigned int)ceil((double)N / PACK_DENSITY);
//   nsnps = (unsigned int)(len / np);
//
//   std::cout << geno_filename << " np: " << np << std::endl;
//   std::cout << geno_filename << " nsnps: " << nsnps << std::endl;
//}

//void Data::set_mode(unsigned int mode)
//{
//   this->mode = mode;
//   Ncurr = (mode == DATA_MODE_TRAIN) ? Ntrain : Ntest;
//   mask_curr = (mode == DATA_MODE_TRAIN) ? mask_train : mask_test;
//   ones = VectorXd::Ones(Ncurr);
//   zeros = VectorXd::Zero(Ncurr);
//}

//VectorXd Data::get_coordinate(unsigned int j)
//{
//   throw std::string("get_coordinate: unused code");
//   
//   // intercept
//   if(j == 0)
//      return ones;
//
//   // SNPs
//   if(j <= nsnps)
//      return get_snp(j - 1);
//
//   // Covariables
//   unsigned int cidx = j - nsnps - 1;
//
//   if(mode == DATA_MODE_TEST && covar_actions[cidx] == COVAR_ACTION_TRAIN_ONLY)
//   {
//      std::cout << timestamp() << " Ignoring covariable " << cidx 
//	 << " (variable " << j << ") in prediction" << std::endl;
//      return zeros;
//   }
//   
//   VectorXd x(Ncurr);
//   for(unsigned int i = 0, k = 0; i < N ; i++)
//      if(mask_curr(i))
//	 x(k++) = X2(i, cidx);
//   return x;
//}

//// Assumes data is SNP-major ordered
//VectorXd Data::get_snp(unsigned int j)
//{
//   if(mode == DATA_MODE_TRAIN)
//   {
//      // return SNP from cache if it exists
//      if(!cache->get(j, geno))
//      {
//	 geno = load_snp(j);
//	 cache->put(j, geno);
//      }
//      return geno;
//   }
//
//   // In test mode, don't try to get from cache, just load from disk. This way
//   // we don't invalidate the cache for the rest of the fitting process.
//   geno = load_snp(j);
//   return geno;
//}
//
//VectorXd Data::load_snp(unsigned int j)
//{
//   double geno_dat[Ncurr];
//   load_snp_double(j, geno_dat);
//   VectorXd geno = Map<VectorXd>(geno_dat, Ncurr, 1);
//   return geno;
//}
//
//// Read a SNP from disk, for the correct set (training or test),
//// impute missing and standardise
//void Data::load_snp_double(unsigned int j, double *geno)
//{
//   unsigned int i;
//
//   // Not found in cache, load from disk
//   unsigned char* tmp2 = new unsigned char[np * PACK_DENSITY];
//   unsigned char* data = (unsigned char*)geno_fin.data();
//   decode_plink(tmp2, data + PLINK_OFFSET + j * np, np);
//
//   // Compute average over non-missing genotypes
//   unsigned int k = 0, ngood = 0;
//   double sum = 0;
//   for(i = 0 ; i < N ; i++)
//   {
//      if(mask_curr(i))
//      {
//	 double v = (double)tmp2[i];
//         geno[k] = v;
//	 if(v != PLINK_NA)
//	 {
//	    ngood++;
//	    sum += v;
//	 }
//	 k++;
//      }
//   }
//
//   delete[] tmp2;
//
//   // Compute standard dev over non-missing genotypes
//   double mean = sum / k;
//   double sum2 = 0, v;
//   for(i = Ncurr - 1 ; i != -1 ; --i)
//   {
//      v = geno[i];
//      if(v != PLINK_NA)
//      {
//	 v -= mean;
//	 sum2 += v * v;
//      }
//   }
//   double sd = sqrt(sum2 / (ngood - 1));
//   double mean_sd = mean / sd;
//
//   if(ngood == Ncurr)
//   {
//      for(i = Ncurr - 1 ; i != -1 ; --i)
//	 geno[i] = (geno[i] - mean) / sd;
//   }
//   else
//   {
//      // Impute missing genotypes and standardise to zero-mean unit-variance
//      for(i = Ncurr - 1 ; i != -1 ; --i)
//      {
//         v = geno[i];
//         if(v == PLINK_NA)
//            geno[i] = mean_sd;
//         else
//            geno[i] = (v - mean) / sd;
//      }
//   }
//}

//// Split and standardise
//void Data::split_data(unsigned int fold)
//{
//   mask_test = (folds.array() == fold);
//   mask_train = (folds.array() != fold);
//   Ntest = mask_test.count();
//   Ntrain = N - Ntest;
//
//   Ytrain = MatrixXd::Zero(Ntrain, Y.cols());
//   Ytest = MatrixXd::Zero(Ntest, Y.cols());
//
//   unsigned int itrain = 0, itest = 0;
//   for(unsigned int r = 0 ; r < N ; r++)
//   {
//      if(mask_train(r))
//      {
//	 Ytrain.row(itrain++) = Y.row(r);
//      }
//      else
//      {
//	 Ytest.row(itest++) = Y.row(r);
//      }
//   }
//
//   std::cout << timestamp() << " Data::split_data(): Ntrain: " << Ntrain << " Ntest: " << Ntest <<
//      std::endl;
//
//   delete cache;
//
//   cache = new Cache(Ntrain, nsnps, cachemem);
//}

//void Data::make_folds(unsigned int rep)
//{
//   char buf[BUFSIZE];
//   VectorXd f = VectorXd::Random(N);
//   folds = ((f.array() + 1.0) / 2.0 * nfolds).cast<int>();
//   VectorXd foldsd = folds.cast<double>();
//   sprintf(buf, "folds_%d.txt", rep);
//   save_text(buf, foldsd);
//}

// Map the required regions in terms of genomic basepairs to indices in the
// data matrix, and mark each SNPs as included/excluded.
// In addition, each putatively included SNP is only loaded based on a random
// coin flip, to sample SNPs and reduce computational burden of computing
// correlation later
//void Data::map_regions()
//{
//   std::cout << timestamp() << " Begin mapping regions (" << regions.size()
//      << " regions, " << snps.size() << " SNPs)" << std::endl;
//
//   nsnps_post_removal = 0;
//   for(unsigned int i = 0 ; i < snps.size(); i++)
//   {
//      snps[i].included = true;
//      for(unsigned j = 0 ; j < regions.size(); j++)
//      {
//	 if(snps[i].chr == regions[j].chr 
//	    && snps[i].bp >= regions[j].begin_bp
//	    && snps[i].bp <= regions[j].end_bp)
//	 {
//	    snps[i].included = false;
//	    break;
//	 }
//      }
//      nsnps_post_removal += snps[i].included;
//      if(verbose)
//	 std::cout << "[" << i << "] "
//	    << snps[i].chr << " " << snps[i].rsid
//	    << " included: " << snps[i].included << std::endl;
//   }
//
//   std::cout << timestamp() << " Total SNPs remaining after region removal: "
//      << nsnps_post_removal << " (out of " << snps.size() << " SNPs)" <<  std::endl;
//
//   std::ofstream out(included_snps_filename, std::ios::out);
//   double prob = (double)nsnps_sampling / nsnps_post_removal;
//   if(prob < 1)
//   {
//      std::cout << timestamp() << " Sampling SNPs (prob=" << prob << ")" << std::endl;
//      nsnps_post_removal = 0;
//      for(unsigned int i = 0 ; i < snps.size(); i++)
//      {
//         if(snps[i].included)
//         {
//            if(drand48() > prob)
//            {
//               snps[i].included = false;
//            }
//            else
//            {
//               nsnps_post_removal++;
//               out << snps[i].rsid << std::endl;
//            }
//         }
//      }
//   }
//
//   out.close();
//
//   std::cout << timestamp() << " Total SNPs remaining after region removal and sampling: "
//      << nsnps_post_removal << " (out of " << snps.size() << " SNPs)" <<  std::endl;
//
//   std::cout << timestamp() << " end mapping regions" << std::endl;
//}

//void Data::reset_regions()
//{
//   nsnps_post_removal = snps.size();
//   for(unsigned int i = 0 ; i < snps.size(); i++)
//      snps[i].included = true;
//}

// BIM format:
// chr     rsID            cM	   BP      Minor   Major
// 10      rs12255619      0       88481   C       A
void Data::read_plink_bim()
{
   std::ifstream in(bim_filename, std::ios::in);

   if(!in)
   {
      std::cerr << "[Data::read_plink_bim] Error reading file "
	 << bim_filename << std::endl;
      throw std::string("io error");
   }
   std::vector<std::string> lines;
   
   while(in)
   {
      std::string line;
      std::getline(in, line);
      if(!in.eof())
	 lines.push_back(line);
   }

   std::cout << timestamp() << " Detected bim file " <<
      bim_filename << ", " << lines.size() << " SNPs" << std::endl;

   in.close();

   snps.reserve(lines.size());

   for(unsigned int i = 0 ; i < lines.size() ; i++)
   {
      std::stringstream ss(lines[i]);
      std::string s;
      std::vector<std::string> tokens;

      while(ss >> s)
	 tokens.push_back(s);

      snp sn;
      sn.chr = std::atoi(tokens[0].c_str()); 
      sn.rsid = tokens[1];
      sn.bp = std::atoi(tokens[3].c_str());
      snps.push_back(sn);
   }

   std::cout << timestamp() << " Read " << snps.size() << " SNPs" << std::endl;
}

