// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// flashpca_internal
List flashpca_internal(const Eigen::Map<Eigen::MatrixXd> X, const int stand, const unsigned int ndim, const unsigned int divisor, const unsigned int maxiter, const double tol, const long seed, const bool verbose, const bool do_loadings, const bool return_scale);
RcppExport SEXP _flashpcaR_flashpca_internal(SEXP XSEXP, SEXP standSEXP, SEXP ndimSEXP, SEXP divisorSEXP, SEXP maxiterSEXP, SEXP tolSEXP, SEXP seedSEXP, SEXP verboseSEXP, SEXP do_loadingsSEXP, SEXP return_scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int >::type stand(standSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type ndim(ndimSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type divisor(divisorSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const long >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const bool >::type do_loadings(do_loadingsSEXP);
    Rcpp::traits::input_parameter< const bool >::type return_scale(return_scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(flashpca_internal(X, stand, ndim, divisor, maxiter, tol, seed, verbose, do_loadings, return_scale));
    return rcpp_result_gen;
END_RCPP
}
// flashpca_plink_internal
List flashpca_plink_internal(const std::string fn, const int stand, const unsigned int ndim, const unsigned int divisor, const unsigned int maxiter, const unsigned int block_size, const double tol, const long seed, const bool verbose, const bool do_loadings, const bool return_scale);
RcppExport SEXP _flashpcaR_flashpca_plink_internal(SEXP fnSEXP, SEXP standSEXP, SEXP ndimSEXP, SEXP divisorSEXP, SEXP maxiterSEXP, SEXP block_sizeSEXP, SEXP tolSEXP, SEXP seedSEXP, SEXP verboseSEXP, SEXP do_loadingsSEXP, SEXP return_scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type fn(fnSEXP);
    Rcpp::traits::input_parameter< const int >::type stand(standSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type ndim(ndimSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type divisor(divisorSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type block_size(block_sizeSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const long >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const bool >::type do_loadings(do_loadingsSEXP);
    Rcpp::traits::input_parameter< const bool >::type return_scale(return_scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(flashpca_plink_internal(fn, stand, ndim, divisor, maxiter, block_size, tol, seed, verbose, do_loadings, return_scale));
    return rcpp_result_gen;
END_RCPP
}
// scca_internal
List scca_internal(const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::MatrixXd> Y, const double lambda1, const double lambda2, const unsigned int ndim, const int stand_x, const int stand_y, const int divisor, const long seed, const int maxiter, const double tol, const bool verbose, const unsigned int num_threads, const bool useV, const Eigen::Map<Eigen::MatrixXd> Vinit);
RcppExport SEXP _flashpcaR_scca_internal(SEXP XSEXP, SEXP YSEXP, SEXP lambda1SEXP, SEXP lambda2SEXP, SEXP ndimSEXP, SEXP stand_xSEXP, SEXP stand_ySEXP, SEXP divisorSEXP, SEXP seedSEXP, SEXP maxiterSEXP, SEXP tolSEXP, SEXP verboseSEXP, SEXP num_threadsSEXP, SEXP useVSEXP, SEXP VinitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda1(lambda1SEXP);
    Rcpp::traits::input_parameter< const double >::type lambda2(lambda2SEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type ndim(ndimSEXP);
    Rcpp::traits::input_parameter< const int >::type stand_x(stand_xSEXP);
    Rcpp::traits::input_parameter< const int >::type stand_y(stand_ySEXP);
    Rcpp::traits::input_parameter< const int >::type divisor(divisorSEXP);
    Rcpp::traits::input_parameter< const long >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type num_threads(num_threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type useV(useVSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Vinit(VinitSEXP);
    rcpp_result_gen = Rcpp::wrap(scca_internal(X, Y, lambda1, lambda2, ndim, stand_x, stand_y, divisor, seed, maxiter, tol, verbose, num_threads, useV, Vinit));
    return rcpp_result_gen;
END_RCPP
}
// ucca_plink_internal
List ucca_plink_internal(const std::string fn, const Eigen::Map<Eigen::MatrixXd> Y, const int stand_x, const int stand_y, const bool verbose);
RcppExport SEXP _flashpcaR_ucca_plink_internal(SEXP fnSEXP, SEXP YSEXP, SEXP stand_xSEXP, SEXP stand_ySEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type fn(fnSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const int >::type stand_x(stand_xSEXP);
    Rcpp::traits::input_parameter< const int >::type stand_y(stand_ySEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(ucca_plink_internal(fn, Y, stand_x, stand_y, verbose));
    return rcpp_result_gen;
END_RCPP
}
// scca_plink_internal
List scca_plink_internal(const std::string fn, const Eigen::Map<Eigen::MatrixXd> Y, const double lambda1, const double lambda2, const unsigned int ndim, const int stand_x, const int stand_y, const int divisor, const long seed, const int maxiter, const double tol, const bool verbose, const unsigned int num_threads, const unsigned int block_size, const bool useV, const Eigen::Map<Eigen::MatrixXd> Vinit);
RcppExport SEXP _flashpcaR_scca_plink_internal(SEXP fnSEXP, SEXP YSEXP, SEXP lambda1SEXP, SEXP lambda2SEXP, SEXP ndimSEXP, SEXP stand_xSEXP, SEXP stand_ySEXP, SEXP divisorSEXP, SEXP seedSEXP, SEXP maxiterSEXP, SEXP tolSEXP, SEXP verboseSEXP, SEXP num_threadsSEXP, SEXP block_sizeSEXP, SEXP useVSEXP, SEXP VinitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type fn(fnSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda1(lambda1SEXP);
    Rcpp::traits::input_parameter< const double >::type lambda2(lambda2SEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type ndim(ndimSEXP);
    Rcpp::traits::input_parameter< const int >::type stand_x(stand_xSEXP);
    Rcpp::traits::input_parameter< const int >::type stand_y(stand_ySEXP);
    Rcpp::traits::input_parameter< const int >::type divisor(divisorSEXP);
    Rcpp::traits::input_parameter< const long >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type num_threads(num_threadsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type block_size(block_sizeSEXP);
    Rcpp::traits::input_parameter< const bool >::type useV(useVSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Vinit(VinitSEXP);
    rcpp_result_gen = Rcpp::wrap(scca_plink_internal(fn, Y, lambda1, lambda2, ndim, stand_x, stand_y, divisor, seed, maxiter, tol, verbose, num_threads, block_size, useV, Vinit));
    return rcpp_result_gen;
END_RCPP
}
// ucca_internal
List ucca_internal(const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::MatrixXd> Y, const int stand_x, const int stand_y, const bool verbose);
RcppExport SEXP _flashpcaR_ucca_internal(SEXP XSEXP, SEXP YSEXP, SEXP stand_xSEXP, SEXP stand_ySEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const int >::type stand_x(stand_xSEXP);
    Rcpp::traits::input_parameter< const int >::type stand_y(stand_ySEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(ucca_internal(X, Y, stand_x, stand_y, verbose));
    return rcpp_result_gen;
END_RCPP
}
// standardise_impute
NumericMatrix standardise_impute(const Eigen::Map<Eigen::MatrixXd> XX, const int method);
RcppExport SEXP _flashpcaR_standardise_impute(SEXP XXSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< const int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(standardise_impute(XX, method));
    return rcpp_result_gen;
END_RCPP
}
// check_internal
List check_internal(const Eigen::MatrixXd& X, const int stand, const Eigen::MatrixXd& evec, const Eigen::VectorXd& eval, const unsigned int divisor, const bool verbose);
RcppExport SEXP _flashpcaR_check_internal(SEXP XSEXP, SEXP standSEXP, SEXP evecSEXP, SEXP evalSEXP, SEXP divisorSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int >::type stand(standSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type evec(evecSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type eval(evalSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type divisor(divisorSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(check_internal(X, stand, evec, eval, divisor, verbose));
    return rcpp_result_gen;
END_RCPP
}
// check_plink_internal
List check_plink_internal(const std::string fn, const int stand, const Eigen::MatrixXd evec, const Eigen::VectorXd eval, const unsigned int block_size, const unsigned int divisor, const bool verbose);
RcppExport SEXP _flashpcaR_check_plink_internal(SEXP fnSEXP, SEXP standSEXP, SEXP evecSEXP, SEXP evalSEXP, SEXP block_sizeSEXP, SEXP divisorSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type fn(fnSEXP);
    Rcpp::traits::input_parameter< const int >::type stand(standSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type evec(evecSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type eval(evalSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type block_size(block_sizeSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type divisor(divisorSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(check_plink_internal(fn, stand, evec, eval, block_size, divisor, verbose));
    return rcpp_result_gen;
END_RCPP
}
// project_plink_internal
List project_plink_internal(const std::string fn, const Eigen::MatrixXd& loadings, const std::vector<std::string> ref_alleles, const Eigen::VectorXd& orig_mean, const Eigen::VectorXd& orig_sd, const unsigned int block_size, const unsigned int divisor, const bool verbose);
RcppExport SEXP _flashpcaR_project_plink_internal(SEXP fnSEXP, SEXP loadingsSEXP, SEXP ref_allelesSEXP, SEXP orig_meanSEXP, SEXP orig_sdSEXP, SEXP block_sizeSEXP, SEXP divisorSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type fn(fnSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type loadings(loadingsSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string> >::type ref_alleles(ref_allelesSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type orig_mean(orig_meanSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type orig_sd(orig_sdSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type block_size(block_sizeSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type divisor(divisorSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(project_plink_internal(fn, loadings, ref_alleles, orig_mean, orig_sd, block_size, divisor, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_flashpcaR_flashpca_internal", (DL_FUNC) &_flashpcaR_flashpca_internal, 10},
    {"_flashpcaR_flashpca_plink_internal", (DL_FUNC) &_flashpcaR_flashpca_plink_internal, 11},
    {"_flashpcaR_scca_internal", (DL_FUNC) &_flashpcaR_scca_internal, 15},
    {"_flashpcaR_ucca_plink_internal", (DL_FUNC) &_flashpcaR_ucca_plink_internal, 5},
    {"_flashpcaR_scca_plink_internal", (DL_FUNC) &_flashpcaR_scca_plink_internal, 16},
    {"_flashpcaR_ucca_internal", (DL_FUNC) &_flashpcaR_ucca_internal, 5},
    {"_flashpcaR_standardise_impute", (DL_FUNC) &_flashpcaR_standardise_impute, 2},
    {"_flashpcaR_check_internal", (DL_FUNC) &_flashpcaR_check_internal, 6},
    {"_flashpcaR_check_plink_internal", (DL_FUNC) &_flashpcaR_check_plink_internal, 7},
    {"_flashpcaR_project_plink_internal", (DL_FUNC) &_flashpcaR_project_plink_internal, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_flashpcaR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
