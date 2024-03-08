// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "atlasqtl_types.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// coreDualLoop
void coreDualLoop(const MapMat X, const MapMat Y, MapArr2D gam_vb, const MapArr2D log_Phi_theta_plus_zeta, const MapArr2D log_1_min_Phi_theta_plus_zeta, const double log_sig2_inv_vb, const MapArr1D log_tau_vb, MapMat m1_beta, MapMat X_beta_vb, MapArr2D mu_beta_vb, const MapArr1D sig2_beta_vb, const MapArr1D tau_vb, const Eigen::VectorXi shuffled_ind, const Eigen::VectorXi sample_q, const double c);
RcppExport SEXP _atlasqtl_coreDualLoop(SEXP XSEXP, SEXP YSEXP, SEXP gam_vbSEXP, SEXP log_Phi_theta_plus_zetaSEXP, SEXP log_1_min_Phi_theta_plus_zetaSEXP, SEXP log_sig2_inv_vbSEXP, SEXP log_tau_vbSEXP, SEXP m1_betaSEXP, SEXP X_beta_vbSEXP, SEXP mu_beta_vbSEXP, SEXP sig2_beta_vbSEXP, SEXP tau_vbSEXP, SEXP shuffled_indSEXP, SEXP sample_qSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MapMat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const MapMat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< MapArr2D >::type gam_vb(gam_vbSEXP);
    Rcpp::traits::input_parameter< const MapArr2D >::type log_Phi_theta_plus_zeta(log_Phi_theta_plus_zetaSEXP);
    Rcpp::traits::input_parameter< const MapArr2D >::type log_1_min_Phi_theta_plus_zeta(log_1_min_Phi_theta_plus_zetaSEXP);
    Rcpp::traits::input_parameter< const double >::type log_sig2_inv_vb(log_sig2_inv_vbSEXP);
    Rcpp::traits::input_parameter< const MapArr1D >::type log_tau_vb(log_tau_vbSEXP);
    Rcpp::traits::input_parameter< MapMat >::type m1_beta(m1_betaSEXP);
    Rcpp::traits::input_parameter< MapMat >::type X_beta_vb(X_beta_vbSEXP);
    Rcpp::traits::input_parameter< MapArr2D >::type mu_beta_vb(mu_beta_vbSEXP);
    Rcpp::traits::input_parameter< const MapArr1D >::type sig2_beta_vb(sig2_beta_vbSEXP);
    Rcpp::traits::input_parameter< const MapArr1D >::type tau_vb(tau_vbSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi >::type shuffled_ind(shuffled_indSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi >::type sample_q(sample_qSEXP);
    Rcpp::traits::input_parameter< const double >::type c(cSEXP);
    coreDualLoop(X, Y, gam_vb, log_Phi_theta_plus_zeta, log_1_min_Phi_theta_plus_zeta, log_sig2_inv_vb, log_tau_vb, m1_beta, X_beta_vb, mu_beta_vb, sig2_beta_vb, tau_vb, shuffled_ind, sample_q, c);
    return R_NilValue;
END_RCPP
}
// coreDualMisLoop
void coreDualMisLoop(const MapMat X, const MapMat Y, MapArr2D gam_vb, const MapArr2D log_Phi_theta_plus_zeta, const MapArr2D log_1_min_Phi_theta_plus_zeta, const double log_sig2_inv_vb, const MapArr1D log_tau_vb, MapMat m1_beta, MapArr2D X_beta_vb, MapArr2D mu_beta_vb, const MapArr2D sig2_beta_vb, const MapArr1D tau_vb, const Eigen::VectorXi shuffled_ind, const MapArr2D mis_pat, const Eigen::VectorXi sample_q, const double c);
RcppExport SEXP _atlasqtl_coreDualMisLoop(SEXP XSEXP, SEXP YSEXP, SEXP gam_vbSEXP, SEXP log_Phi_theta_plus_zetaSEXP, SEXP log_1_min_Phi_theta_plus_zetaSEXP, SEXP log_sig2_inv_vbSEXP, SEXP log_tau_vbSEXP, SEXP m1_betaSEXP, SEXP X_beta_vbSEXP, SEXP mu_beta_vbSEXP, SEXP sig2_beta_vbSEXP, SEXP tau_vbSEXP, SEXP shuffled_indSEXP, SEXP mis_patSEXP, SEXP sample_qSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MapMat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const MapMat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< MapArr2D >::type gam_vb(gam_vbSEXP);
    Rcpp::traits::input_parameter< const MapArr2D >::type log_Phi_theta_plus_zeta(log_Phi_theta_plus_zetaSEXP);
    Rcpp::traits::input_parameter< const MapArr2D >::type log_1_min_Phi_theta_plus_zeta(log_1_min_Phi_theta_plus_zetaSEXP);
    Rcpp::traits::input_parameter< const double >::type log_sig2_inv_vb(log_sig2_inv_vbSEXP);
    Rcpp::traits::input_parameter< const MapArr1D >::type log_tau_vb(log_tau_vbSEXP);
    Rcpp::traits::input_parameter< MapMat >::type m1_beta(m1_betaSEXP);
    Rcpp::traits::input_parameter< MapArr2D >::type X_beta_vb(X_beta_vbSEXP);
    Rcpp::traits::input_parameter< MapArr2D >::type mu_beta_vb(mu_beta_vbSEXP);
    Rcpp::traits::input_parameter< const MapArr2D >::type sig2_beta_vb(sig2_beta_vbSEXP);
    Rcpp::traits::input_parameter< const MapArr1D >::type tau_vb(tau_vbSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi >::type shuffled_ind(shuffled_indSEXP);
    Rcpp::traits::input_parameter< const MapArr2D >::type mis_pat(mis_patSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi >::type sample_q(sample_qSEXP);
    Rcpp::traits::input_parameter< const double >::type c(cSEXP);
    coreDualMisLoop(X, Y, gam_vb, log_Phi_theta_plus_zeta, log_1_min_Phi_theta_plus_zeta, log_sig2_inv_vb, log_tau_vb, m1_beta, X_beta_vb, mu_beta_vb, sig2_beta_vb, tau_vb, shuffled_ind, mis_pat, sample_q, c);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_atlasqtl_coreDualLoop", (DL_FUNC) &_atlasqtl_coreDualLoop, 15},
    {"_atlasqtl_coreDualMisLoop", (DL_FUNC) &_atlasqtl_coreDualMisLoop, 16},
    {NULL, NULL, 0}
};

RcppExport void R_init_atlasqtl(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
