// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// idw
arma::mat idw(const arma::mat& coords_pred, const arma::mat& coords_vals, const arma::mat& vals, const int& p);
RcppExport SEXP _SpatialBSS_idw(SEXP coords_predSEXP, SEXP coords_valsSEXP, SEXP valsSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type coords_pred(coords_predSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type coords_vals(coords_valsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type vals(valsSEXP);
    Rcpp::traits::input_parameter< const int& >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(idw(coords_pred, coords_vals, vals, p));
    return rcpp_result_gen;
END_RCPP
}
// k_mat_ball
arma::mat k_mat_ball(const arma::mat& coords, const double& h);
RcppExport SEXP _SpatialBSS_k_mat_ball(SEXP coordsSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< const double& >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(k_mat_ball(coords, h));
    return rcpp_result_gen;
END_RCPP
}
// k_mat_ring
arma::mat k_mat_ring(const arma::mat& coords, const double& h1, const double& h2);
RcppExport SEXP _SpatialBSS_k_mat_ring(SEXP coordsSEXP, SEXP h1SEXP, SEXP h2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< const double& >::type h1(h1SEXP);
    Rcpp::traits::input_parameter< const double& >::type h2(h2SEXP);
    rcpp_result_gen = Rcpp::wrap(k_mat_ring(coords, h1, h2));
    return rcpp_result_gen;
END_RCPP
}
// k_mat_exp
arma::mat k_mat_exp(const arma::mat& coords, const double& h);
RcppExport SEXP _SpatialBSS_k_mat_exp(SEXP coordsSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< const double& >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(k_mat_exp(coords, h));
    return rcpp_result_gen;
END_RCPP
}
// k_mat_ball_angle
arma::mat k_mat_ball_angle(const arma::mat& coords, const double& h, const double& am, const double& tol);
RcppExport SEXP _SpatialBSS_k_mat_ball_angle(SEXP coordsSEXP, SEXP hSEXP, SEXP amSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< const double& >::type h(hSEXP);
    Rcpp::traits::input_parameter< const double& >::type am(amSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(k_mat_ball_angle(coords, h, am, tol));
    return rcpp_result_gen;
END_RCPP
}
// k_mat_ring_angle
arma::mat k_mat_ring_angle(const arma::mat& coords, const double& h1, const double& h2, const double& am, const double& tol);
RcppExport SEXP _SpatialBSS_k_mat_ring_angle(SEXP coordsSEXP, SEXP h1SEXP, SEXP h2SEXP, SEXP amSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< const double& >::type h1(h1SEXP);
    Rcpp::traits::input_parameter< const double& >::type h2(h2SEXP);
    Rcpp::traits::input_parameter< const double& >::type am(amSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(k_mat_ring_angle(coords, h1, h2, am, tol));
    return rcpp_result_gen;
END_RCPP
}
// k_mat_exp_angle
arma::mat k_mat_exp_angle(const arma::mat& coords, const double& h, const double& am, const double& tol);
RcppExport SEXP _SpatialBSS_k_mat_exp_angle(SEXP coordsSEXP, SEXP hSEXP, SEXP amSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< const double& >::type h(hSEXP);
    Rcpp::traits::input_parameter< const double& >::type am(amSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(k_mat_exp_angle(coords, h, am, tol));
    return rcpp_result_gen;
END_RCPP
}
// idx_per_block
Rcpp::List idx_per_block(const arma::mat& coords, const arma::mat& coords_block, const int d);
RcppExport SEXP _SpatialBSS_idx_per_block(SEXP coordsSEXP, SEXP coords_blockSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type coords_block(coords_blockSEXP);
    Rcpp::traits::input_parameter< const int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(idx_per_block(coords, coords_block, d));
    return rcpp_result_gen;
END_RCPP
}
// sp_lcov_sparse
arma::mat sp_lcov_sparse(const arma::mat& x, const arma::mat& k);
RcppExport SEXP _SpatialBSS_sp_lcov_sparse(SEXP xSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(sp_lcov_sparse(x, k));
    return rcpp_result_gen;
END_RCPP
}
// sp_ldiff_sparse
arma::mat sp_ldiff_sparse(const arma::mat& x, const arma::mat& k);
RcppExport SEXP _SpatialBSS_sp_ldiff_sparse(SEXP xSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(sp_ldiff_sparse(x, k));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SpatialBSS_idw", (DL_FUNC) &_SpatialBSS_idw, 4},
    {"_SpatialBSS_k_mat_ball", (DL_FUNC) &_SpatialBSS_k_mat_ball, 2},
    {"_SpatialBSS_k_mat_ring", (DL_FUNC) &_SpatialBSS_k_mat_ring, 3},
    {"_SpatialBSS_k_mat_exp", (DL_FUNC) &_SpatialBSS_k_mat_exp, 2},
    {"_SpatialBSS_k_mat_ball_angle", (DL_FUNC) &_SpatialBSS_k_mat_ball_angle, 4},
    {"_SpatialBSS_k_mat_ring_angle", (DL_FUNC) &_SpatialBSS_k_mat_ring_angle, 5},
    {"_SpatialBSS_k_mat_exp_angle", (DL_FUNC) &_SpatialBSS_k_mat_exp_angle, 4},
    {"_SpatialBSS_idx_per_block", (DL_FUNC) &_SpatialBSS_idx_per_block, 3},
    {"_SpatialBSS_sp_lcov_sparse", (DL_FUNC) &_SpatialBSS_sp_lcov_sparse, 2},
    {"_SpatialBSS_sp_ldiff_sparse", (DL_FUNC) &_SpatialBSS_sp_ldiff_sparse, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_SpatialBSS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
