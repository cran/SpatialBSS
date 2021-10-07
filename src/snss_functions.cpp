// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// functions that gets center coordinates for each block and derives
// for each coordinate the corresponding block index
// the numbers of observations per block
// [[Rcpp::export]]
Rcpp::List idx_per_block(const arma::mat& coords,
                         const arma::mat& coords_block,
                         const int d) {
  
  const int n = coords.n_rows;
  const int m = coords_block.n_rows;
  int idx_min = 0;
  
  arma::vec distances = arma::zeros<arma::vec>(m);
  arma::vec block_idx = arma::zeros<arma::vec>(n);
  arma::vec n_per_block = arma::zeros<arma::vec>(m);
  
  for(int i=0; i < n; ++i) {
    for(int j=0; j < m; ++j) {
      distances(j) = arma::norm(coords.row(i) - coords_block.row(j), d);
    }
    idx_min = distances.index_min();
    block_idx(i) = idx_min;
    n_per_block(idx_min) += 1;  
  }
  
  return Rcpp::List::create(Rcpp::Named("block_idx") = block_idx,
                            Rcpp::Named("n_per_block") = n_per_block,
                            Rcpp::Named("n_blocks") = m);
}


