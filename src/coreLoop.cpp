/*
 *
 * This file is part of the `atlasqtl` R package:
 *     https://github.com/hruffieux/atlasqtl
 *
 * Functions for computationally expensive updates in algorithms without
 * external information.
 *
 * These functions use Eigen::Map to pass large matrices by reference from R.
 * Given dimensionalities involved in some applications, copying such matrices
 * would imply a prohibitive RAM overconsumption.
 *
 */

#include "utils.h"

double logOnePlusExp(double x) {
  double m = x;
  if (x < 0)
    m = 0;
  return log(exp(x-m) + exp(-m)) + m;
}


// for atlasqtl_core function
// [[Rcpp::export]]
void coreDualLoop(const MapMat X,
                  const MapMat Y,
                  MapArr2D gam_vb,
                  const MapArr2D log_Phi_theta_plus_zeta,
                  const MapArr2D log_1_min_Phi_theta_plus_zeta,
                  const double log_sig2_inv_vb,
                  const MapArr1D log_tau_vb,
                  MapMat m1_beta,
                  MapMat X_beta_vb,
                  MapArr2D mu_beta_vb,
                  const MapArr1D sig2_beta_vb,
                  const MapArr1D tau_vb,
                  const Eigen::VectorXi shuffled_ind,
                  const Eigen::VectorXi sample_q,
                  const double c = 1) {
  
  const Arr1D cst = -(log_tau_vb + log_sig2_inv_vb + log(sig2_beta_vb) )/ 2;
  
  for (int a = 0; a < sample_q.size(); a++) {
    int k = sample_q[a];
    
    for (int b = 0; b < shuffled_ind.size(); b++) {
      int j = shuffled_ind[b];
      
      //X_beta_vb.col(k) -= X.col(j) * m1_beta(j, k); 
      double m1_beta_jk = m1_beta(j, k);
      
      
      mu_beta_vb(j, k) = c * sig2_beta_vb[k] * tau_vb[k] *
        (Y.col(k) - (X_beta_vb.col(k) - X.col(j) * m1_beta_jk)).dot(X.col(j));

      
      gam_vb(j, k) = exp(-logOnePlusExp(c * (log_1_min_Phi_theta_plus_zeta(j, k) - log_Phi_theta_plus_zeta(j, k)
                                               - mu_beta_vb(j, k)*mu_beta_vb(j, k) / (2 * sig2_beta_vb[k])
                                               + cst[k])));
                                               
                                               
                                                m1_beta(j, k) = mu_beta_vb(j, k) * gam_vb(j, k);
                                               
                                               X_beta_vb.col(k) += X.col(j) * (m1_beta(j, k) - m1_beta_jk); 
                                               
    }
  }
  
}


// for atlasqtl_core function handling missing values in Y
// [[Rcpp::export]]
void coreDualMisLoop(const MapMat X,
                     const MapMat Y,
                     MapArr2D gam_vb,
                     const MapArr2D log_Phi_theta_plus_zeta,
                     const MapArr2D log_1_min_Phi_theta_plus_zeta,
                     const double log_sig2_inv_vb,
                     const MapArr1D log_tau_vb,
                     MapMat m1_beta,
                     MapArr2D X_beta_vb,
                     MapArr2D mu_beta_vb,
                     const MapArr2D sig2_beta_vb,
                     const MapArr1D tau_vb,
                     const Eigen::VectorXi shuffled_ind,
                     const MapArr2D mis_pat,
                     const Eigen::VectorXi sample_q,
                     const double c = 1) {
  
  const Arr1D cst = -(log_tau_vb + log_sig2_inv_vb)/ 2;
  
  X_beta_vb *= mis_pat;
  
  for (int a = 0; a < sample_q.size(); a++) {
    int k = sample_q[a];
    
    for (int b = 0; b < shuffled_ind.size(); b++) {
      int j = shuffled_ind[b];
      

      // X_beta_vb.col(k) -= (m1_beta(j, k)*X.col(j)).array() * mis_pat.col(k);
      // mu_beta_vb(j, k) = c * sig2_beta_vb(j, k) * tau_vb[k] *
      //   ((Y.col(k) - X_beta_vb.matrix().col(k)).transpose() * X.col(j)).sum();
      
      double m1_beta_jk = m1_beta(j, k);
      mu_beta_vb(j, k) = c * sig2_beta_vb(j, k) * tau_vb[k] *
        ((Y.col(k) - (X_beta_vb.col(k) - (m1_beta_jk * X.col(j)).array() * mis_pat.col(k)).matrix()).transpose() * X.col(j)).sum();
      
      

      
      //std::cout << "The value of mu_beta_vb(j, k) is: " << mu_beta_vb(j, k) << std::endl;
      
      gam_vb(j, k) = exp(-logOnePlusExp(c * (log_1_min_Phi_theta_plus_zeta(j, k) -
        log_Phi_theta_plus_zeta(j, k) - mu_beta_vb(j, k)*mu_beta_vb(j, k) / (2 * sig2_beta_vb(j, k)) -
        log(sig2_beta_vb(j, k)) / 2 + cst[k])));
      
      m1_beta(j, k) = mu_beta_vb(j, k) * gam_vb(j, k);
      
      X_beta_vb.col(k) += X.col(j).array() * (m1_beta(j, k)-m1_beta_jk) * mis_pat.col(k);
      // X_beta_vb.col(k) += (m1_beta(j, k)*X.col(j)).array() * mis_pat.col(k);
        
    }
  }
  
}
