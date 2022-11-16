#------------------------------------------------#
#  local covariance matrix computation
#------------------------------------------------#
local_covariance_matrix <- function(x, kernel_list, lcov = c('lcov', 'ldiff', 'lcov_norm'), 
                                    center = TRUE) {
  lcov <- match.arg(lcov)
  
  if (center) {
    x <- white_data(x)$x_0
  }
  
  cov_sp_list <- switch(lcov,
                        'lcov'      = lapply(kernel_list, function(k_mat) sp_lcov_sparse(x = x, k = k_mat)),
                        'ldiff'     = lapply(kernel_list, function(k_mat) sp_ldiff_sparse(x = x, k = k_mat)),
                        'lcov_norm' = lapply(kernel_list, function(k_mat) sp_lcov_sparse(x = x, k = k_mat) * sqrt(nrow(k_mat) / sum(k_mat ^ 2)))
  )
  
  cov_sp_list <- lapply(cov_sp_list, function(x) (x + t(x)) / 2)
  
  attr(cov_sp_list, 'lcov') <- lcov
  
  return(cov_sp_list)
}

#------------------------------------------------#
#  robust local covariance
#------------------------------------------------#
local_gss_covariance_matrix <- function(x, kernel_list, 
                                       lcov = c('norm', 'winsor', 'qwinsor'), 
                                       center = TRUE) {
  lcov <- match.arg(lcov)
  
  if (center) {
    x <- white_data(x, whitening = "hr")$x_0
  }
  
  x <- switch(lcov,
              'norm' = norm_transform(x),
              'winsor' = winsor_transform(x),
              'qwinsor' = qwinsor_transform(x))
  
  cov_sp_list <- lapply(kernel_list, function(k_mat) sp_lcov_sparse(x = x$x, k = k_mat) * sqrt(nrow(k_mat) / sum(k_mat ^ 2)))
  cov_sp_list <- lapply(cov_sp_list, function(x) (x + t(x)) / 2)
  
  attr(cov_sp_list, 'lcov') <- lcov
  
  return(list(cov_sp_list = cov_sp_list, weights = x$weights))
}

norm_transform <- function(X) {
  n <- dim(X)[1] 
  p <- dim(X)[2]
  dists <- sqrt(rowSums(X^2))
  X_norm <- X / dists
  return(list(x = X_norm, weights = rep(1, length(dists))))
}

winsor_transform <- function(X) {
  n <- dim(X)[1] 
  p <- dim(X)[2]
  dists <- sqrt(rowSums(X^2))
  d.hmed  <- sort(dists)[floor((length(dists) + p + 1) / 2)]
  idx     <- which(dists > d.hmed)
  xi      <- rep(1, length(dists))
  xi[idx] <-  1 / dists[idx] * d.hmed
  X_winsor <- X * xi
  return(list(x = X_winsor, weights = xi))
}

qwinsor_transform <- function(X) {
  n <- dim(X)[1] 
  p <- dim(X)[2]
  dists <- sqrt(rowSums(X^2))
  d.hmed  <- sort(dists)[floor((length(dists) + p + 1) / 2)]
  idx     <- which(dists > d.hmed)
  xi      <- rep(1, length(dists))
  xi[idx] <- 1 / dists[idx]^2 * d.hmed^2
  X_qwinsor <- X * xi
  return(list(x = X_qwinsor, weights = xi))
}