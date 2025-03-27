// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// computes the kernel matrix for an isotropic ball kernel
//
// [[Rcpp::export]]
arma::mat k_mat_ball(const arma::mat & coords,
                     const double & h) {
  const int n = coords.n_rows;
  double distance;
  arma::mat k(n, n);
  k.ones();
  
  for(int i=0; i < n; ++i) {
    for(int j=i+1; j < n; ++j) {
      distance = norm(coords.row(i) - coords.row(j));
      if(distance > h){
        k(i,j) = k(j,i) = 0;
      }
    }
  }
  
  return k;
}

// computes the kernel matrix for an isotropic ring kernel
//
// [[Rcpp::export]]
arma::mat k_mat_ring(const arma::mat & coords,
                     const double & h1,
                     const double & h2) {
  const int n = coords.n_rows;
  double distance;
  arma::mat k(n, n);
  k.zeros();
  
  for(int i=0; i < n; ++i) {
    for(int j=i+1; j < n; ++j) {
      distance = norm(coords.row(i) - coords.row(j));
      if(h1 < distance && distance <= h2){
        k(i,j) = k(j,i) = 1;
      }
    }
  }
  
  return k;
}

// computes the kernel matrix for an isotropic gauss kernel
//
// [[Rcpp::export]]
arma::mat k_mat_exp(const arma::mat & coords,
                    const double & h) {
  const int n = coords.n_rows;
  double distance;
  arma::mat k(n, n);
  k.ones();
  
  for(int i=0; i < n; ++i) {
    for(int j=i+1; j < n; ++j) {
      distance = norm(coords.row(i) - coords.row(j));
      k(i,j) = k(j,i) = exp(- 0.5 * pow(distance, 2) / h);
    }
  }
  
  return k;
}

// computes the kernel matrix for a ball kernel with angle
//
// [[Rcpp::export]]
arma::mat k_mat_ball_angle(const arma::mat & coords,
                     const double & h,
                     const double & am,
                     const double & tol) {
  const int n = coords.n_rows;
  double distance, aij, aji;
  arma::rowvec dif_vec(2);
  arma::mat k(n, n);
  k.zeros();
  
  for(int i=0; i < n; ++i) {
    for(int j=i; j < n; ++j) {
      dif_vec = coords.row(i) - coords.row(j);
      distance = norm(dif_vec);
      if(distance <= h){
        if(distance == 0) {
          k(i,j) = k(j,i) = 1;          
        } else {
          aij = acos((dif_vec(0) * cos(am) + dif_vec(1) * sin(am)) / distance);
          aji = acos((- dif_vec(0) * cos(am) - dif_vec(1) * sin(am)) / distance);
          if (aij <= tol || aji <= tol) {
            k(i,j) = k(j,i) = 1;
          }   
        }
      }
    }
  }
  
  return k;
}

// computes the kernel matrix for a ring kernel with angle
//
// [[Rcpp::export]]
arma::mat k_mat_ring_angle(const arma::mat & coords,
                     const double & h1,
                     const double & h2,
                     const double & am,
                     const double & tol) {
  const int n = coords.n_rows;
  double distance, aij, aji;
  arma::rowvec dif_vec(2);
  arma::mat k(n, n);
  k.zeros();
  
  for(int i=0; i < n; ++i) {
    for(int j=i+1; j < n; ++j) {
      dif_vec = coords.row(i) - coords.row(j);
      distance = norm(dif_vec);
      if(h1 < distance && distance <= h2){
        aij = acos((dif_vec(0) * cos(am) + dif_vec(1) * sin(am)) / distance);
        aji = acos((- dif_vec(0) * cos(am) - dif_vec(1) * sin(am)) / distance);
        if (aij <= tol || aji <= tol) {
          k(i,j) = k(j,i) = 1;
        }
      }
    }
  }
  
  return k;
}

// computes the kernel matrix for a gauss kernel  with angle
//
// [[Rcpp::export]]
arma::mat k_mat_exp_angle(const arma::mat & coords,
                    const double & h,
                    const double & am,
                    const double & tol) {
  const int n = coords.n_rows;
  double distance, aij, aji;
  arma::rowvec dif_vec(2);
  arma::mat k(n, n);
  k.zeros();
  
  for(int i=0; i < n; ++i) {
    for(int j=i; j < n; ++j) {
      dif_vec = coords.row(i) - coords.row(j);
      distance = norm(dif_vec);
      if(distance == 0) {
        k(i,j) = k(j,i) = exp(- 0.5 * pow(distance, 2) / h);       
      } else {
        aij = acos((dif_vec(0) * cos(am) + dif_vec(1) * sin(am)) / distance);
        aji = acos((- dif_vec(0) * cos(am) - dif_vec(1) * sin(am)) / distance);
        if (aij <= tol || aji <= tol) {
          k(i,j) = k(j,i) = exp(- 0.5 * pow(distance, 2) / h);
        }   
      }
    }
  }
  
  return k;
}


