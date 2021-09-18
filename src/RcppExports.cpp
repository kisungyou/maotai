// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// compute_SSR
double compute_SSR(arma::mat& D, arma::mat& Delta);
RcppExport SEXP _maotai_compute_SSR(SEXP DSEXP, SEXP DeltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Delta(DeltaSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_SSR(D, Delta));
    return rcpp_result_gen;
END_RCPP
}
// compute_stress
double compute_stress(arma::mat& D, arma::mat& Dhat);
RcppExport SEXP _maotai_compute_stress(SEXP DSEXP, SEXP DhatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Dhat(DhatSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_stress(D, Dhat));
    return rcpp_result_gen;
END_RCPP
}
// main_bmds
Rcpp::List main_bmds(arma::mat D, arma::mat X0, double sigg0, double a, double alpha, int maxiter, double constant, bool verbose, arma::vec betas);
RcppExport SEXP _maotai_main_bmds(SEXP DSEXP, SEXP X0SEXP, SEXP sigg0SEXP, SEXP aSEXP, SEXP alphaSEXP, SEXP maxiterSEXP, SEXP constantSEXP, SEXP verboseSEXP, SEXP betasSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X0(X0SEXP);
    Rcpp::traits::input_parameter< double >::type sigg0(sigg0SEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type constant(constantSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type betas(betasSEXP);
    rcpp_result_gen = Rcpp::wrap(main_bmds(D, X0, sigg0, a, alpha, maxiter, constant, verbose, betas));
    return rcpp_result_gen;
END_RCPP
}
// aux_shortestpath
Rcpp::NumericMatrix aux_shortestpath(NumericMatrix& wmat);
RcppExport SEXP _maotai_aux_shortestpath(SEXP wmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type wmat(wmatSEXP);
    rcpp_result_gen = Rcpp::wrap(aux_shortestpath(wmat));
    return rcpp_result_gen;
END_RCPP
}
// cppsub_2007Wang
arma::mat cppsub_2007Wang(arma::mat V0, int mm, int d, arma::mat Spu, arma::mat Stu, int maxiter, double eps);
RcppExport SEXP _maotai_cppsub_2007Wang(SEXP V0SEXP, SEXP mmSEXP, SEXP dSEXP, SEXP SpuSEXP, SEXP StuSEXP, SEXP maxiterSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type V0(V0SEXP);
    Rcpp::traits::input_parameter< int >::type mm(mmSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Spu(SpuSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Stu(StuSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(cppsub_2007Wang(V0, mm, d, Spu, Stu, maxiter, eps));
    return rcpp_result_gen;
END_RCPP
}
// gradF
arma::mat gradF(Function func, arma::mat xnow, double h);
RcppExport SEXP _maotai_gradF(SEXP funcSEXP, SEXP xnowSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Function >::type func(funcSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xnow(xnowSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(gradF(func, xnow, h));
    return rcpp_result_gen;
END_RCPP
}
// dat2centers
arma::vec dat2centers(arma::rowvec data, arma::mat& centers);
RcppExport SEXP _maotai_dat2centers(SEXP dataSEXP, SEXP centersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type centers(centersSEXP);
    rcpp_result_gen = Rcpp::wrap(dat2centers(data, centers));
    return rcpp_result_gen;
END_RCPP
}
// cpp_sylvester
arma::mat cpp_sylvester(arma::mat A, arma::mat B, arma::mat C);
RcppExport SEXP _maotai_cpp_sylvester(SEXP ASEXP, SEXP BSEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_sylvester(A, B, C));
    return rcpp_result_gen;
END_RCPP
}
// solve_lyapunov
arma::mat solve_lyapunov(arma::mat A, arma::mat B, arma::mat C);
RcppExport SEXP _maotai_solve_lyapunov(SEXP ASEXP, SEXP BSEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_lyapunov(A, B, C));
    return rcpp_result_gen;
END_RCPP
}
// cpp_weiszfeld
arma::rowvec cpp_weiszfeld(arma::mat X, double abstol, int maxiter, arma::rowvec xinit, arma::vec weights, double epsnum);
RcppExport SEXP _maotai_cpp_weiszfeld(SEXP XSEXP, SEXP abstolSEXP, SEXP maxiterSEXP, SEXP xinitSEXP, SEXP weightsSEXP, SEXP epsnumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type xinit(xinitSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< double >::type epsnum(epsnumSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_weiszfeld(X, abstol, maxiter, xinit, weights, epsnum));
    return rcpp_result_gen;
END_RCPP
}
// cpp_kmeans
Rcpp::List cpp_kmeans(arma::mat data, int k);
RcppExport SEXP _maotai_cpp_kmeans(SEXP dataSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_kmeans(data, k));
    return rcpp_result_gen;
END_RCPP
}
// emds_gamma0
double emds_gamma0(arma::mat dmat);
RcppExport SEXP _maotai_emds_gamma0(SEXP dmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type dmat(dmatSEXP);
    rcpp_result_gen = Rcpp::wrap(emds_gamma0(dmat));
    return rcpp_result_gen;
END_RCPP
}
// cpp_pairwise_L2
Rcpp::List cpp_pairwise_L2(arma::mat muA, arma::mat muB, arma::cube covA, arma::cube covB);
RcppExport SEXP _maotai_cpp_pairwise_L2(SEXP muASEXP, SEXP muBSEXP, SEXP covASEXP, SEXP covBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type muA(muASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type muB(muBSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type covA(covASEXP);
    Rcpp::traits::input_parameter< arma::cube >::type covB(covBSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_pairwise_L2(muA, muB, covA, covB));
    return rcpp_result_gen;
END_RCPP
}
// integrate_1d
double integrate_1d(arma::vec& tseq, arma::vec& fval);
RcppExport SEXP _maotai_integrate_1d(SEXP tseqSEXP, SEXP fvalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type tseq(tseqSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type fval(fvalSEXP);
    rcpp_result_gen = Rcpp::wrap(integrate_1d(tseq, fval));
    return rcpp_result_gen;
END_RCPP
}
// cpp_pdist
arma::mat cpp_pdist(arma::mat X);
RcppExport SEXP _maotai_cpp_pdist(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_pdist(X));
    return rcpp_result_gen;
END_RCPP
}
// cpp_geigen
Rcpp::List cpp_geigen(arma::mat& A, arma::mat& B);
RcppExport SEXP _maotai_cpp_geigen(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_geigen(A, B));
    return rcpp_result_gen;
END_RCPP
}
// cpp_triangle
bool cpp_triangle(arma::mat& D);
RcppExport SEXP _maotai_cpp_triangle(SEXP DSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type D(DSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_triangle(D));
    return rcpp_result_gen;
END_RCPP
}
// cpp_mmds
arma::mat cpp_mmds(arma::mat& D, int ndim, int maxiter, double abstol);
RcppExport SEXP _maotai_cpp_mmds(SEXP DSEXP, SEXP ndimSEXP, SEXP maxiterSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type ndim(ndimSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_mmds(D, ndim, maxiter, abstol));
    return rcpp_result_gen;
END_RCPP
}
// eval_gaussian
double eval_gaussian(arma::vec x, arma::vec mu, arma::mat cov);
RcppExport SEXP _maotai_eval_gaussian(SEXP xSEXP, SEXP muSEXP, SEXP covSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cov(covSEXP);
    rcpp_result_gen = Rcpp::wrap(eval_gaussian(x, mu, cov));
    return rcpp_result_gen;
END_RCPP
}
// eval_gaussian_data
arma::vec eval_gaussian_data(arma::mat X, arma::vec mu, arma::mat cov);
RcppExport SEXP _maotai_eval_gaussian_data(SEXP XSEXP, SEXP muSEXP, SEXP covSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cov(covSEXP);
    rcpp_result_gen = Rcpp::wrap(eval_gaussian_data(X, mu, cov));
    return rcpp_result_gen;
END_RCPP
}
// eval_gmm_data
arma::vec eval_gmm_data(arma::mat X, arma::mat mus, arma::cube covs, arma::vec weight);
RcppExport SEXP _maotai_eval_gmm_data(SEXP XSEXP, SEXP musSEXP, SEXP covsSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mus(musSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type covs(covsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(eval_gmm_data(X, mus, covs, weight));
    return rcpp_result_gen;
END_RCPP
}
// eval_gmm
double eval_gmm(arma::vec x, arma::mat mus, arma::cube covs, arma::vec weight);
RcppExport SEXP _maotai_eval_gmm(SEXP xSEXP, SEXP musSEXP, SEXP covsSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mus(musSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type covs(covsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(eval_gmm(x, mus, covs, weight));
    return rcpp_result_gen;
END_RCPP
}
// src_construct_by_knn
arma::sp_umat src_construct_by_knn(arma::umat& nn_idx, bool intersection);
RcppExport SEXP _maotai_src_construct_by_knn(SEXP nn_idxSEXP, SEXP intersectionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::umat& >::type nn_idx(nn_idxSEXP);
    Rcpp::traits::input_parameter< bool >::type intersection(intersectionSEXP);
    rcpp_result_gen = Rcpp::wrap(src_construct_by_knn(nn_idx, intersection));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_maotai_compute_SSR", (DL_FUNC) &_maotai_compute_SSR, 2},
    {"_maotai_compute_stress", (DL_FUNC) &_maotai_compute_stress, 2},
    {"_maotai_main_bmds", (DL_FUNC) &_maotai_main_bmds, 9},
    {"_maotai_aux_shortestpath", (DL_FUNC) &_maotai_aux_shortestpath, 1},
    {"_maotai_cppsub_2007Wang", (DL_FUNC) &_maotai_cppsub_2007Wang, 7},
    {"_maotai_gradF", (DL_FUNC) &_maotai_gradF, 3},
    {"_maotai_dat2centers", (DL_FUNC) &_maotai_dat2centers, 2},
    {"_maotai_cpp_sylvester", (DL_FUNC) &_maotai_cpp_sylvester, 3},
    {"_maotai_solve_lyapunov", (DL_FUNC) &_maotai_solve_lyapunov, 3},
    {"_maotai_cpp_weiszfeld", (DL_FUNC) &_maotai_cpp_weiszfeld, 6},
    {"_maotai_cpp_kmeans", (DL_FUNC) &_maotai_cpp_kmeans, 2},
    {"_maotai_emds_gamma0", (DL_FUNC) &_maotai_emds_gamma0, 1},
    {"_maotai_cpp_pairwise_L2", (DL_FUNC) &_maotai_cpp_pairwise_L2, 4},
    {"_maotai_integrate_1d", (DL_FUNC) &_maotai_integrate_1d, 2},
    {"_maotai_cpp_pdist", (DL_FUNC) &_maotai_cpp_pdist, 1},
    {"_maotai_cpp_geigen", (DL_FUNC) &_maotai_cpp_geigen, 2},
    {"_maotai_cpp_triangle", (DL_FUNC) &_maotai_cpp_triangle, 1},
    {"_maotai_cpp_mmds", (DL_FUNC) &_maotai_cpp_mmds, 4},
    {"_maotai_eval_gaussian", (DL_FUNC) &_maotai_eval_gaussian, 3},
    {"_maotai_eval_gaussian_data", (DL_FUNC) &_maotai_eval_gaussian_data, 3},
    {"_maotai_eval_gmm_data", (DL_FUNC) &_maotai_eval_gmm_data, 4},
    {"_maotai_eval_gmm", (DL_FUNC) &_maotai_eval_gmm, 4},
    {"_maotai_src_construct_by_knn", (DL_FUNC) &_maotai_src_construct_by_knn, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_maotai(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
