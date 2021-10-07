#------------------------------------------------#
#   data whitening
#------------------------------------------------#
white_data <- function(x, rob_whitening = FALSE, lcov = c('lcov', 'ldiff', 'lcov_norm'), kernel_mat = numeric(0)) {
  lcov <- match.arg(lcov)
  
  mu <- colMeans(x)
  x_0 <- sweep(x, MARGIN = 2, STATS = mu, FUN = '-')
  
  if (rob_whitening) {
    s <- local_covariance_matrix(x, list(kernel_mat), lcov = lcov, center = TRUE)[[1]]
  } else {
    s <- crossprod(x_0) / (nrow(x) - 1)
  }

  s_evd <- eigen(s, symmetric = TRUE)
  s_inv_sqrt <- s_evd$vectors %*% tcrossprod(diag(1 / sqrt(s_evd$values)), s_evd$vectors)
  s_sqrt <- s_evd$vectors %*% tcrossprod(diag(sqrt(s_evd$values)), s_evd$vectors)
  
  x_w <- tcrossprod(x_0, s_inv_sqrt)
  colnames(x_w) <- colnames(x_0)
  
  return(list(mu = mu, x_0 = x_0, x_w = x_w, s = s, s_inv_sqrt = s_inv_sqrt, s_sqrt = s_sqrt))
}

#------------------------------------------------#
#   spatial kernel function computation
#------------------------------------------------#
spatial_kernel_matrix <- function(coords, kernel_type = c('ring', 'ball', 'gauss'), kernel_parameters) {
  kernel_type <- match.arg(kernel_type)
  
  if (kernel_type == 'ball') {
    kernel_list <- lapply(kernel_parameters, 
                          function(h) if (h >= 0) k_mat_ball(coords = coords, h = h)
                                      else stop('Radius must be zero or positive.'))    
  } else if (kernel_type == 'ring' && (length(kernel_parameters) %% 2 == 0)) {
    kernel_list <- lapply(seq(from = 1, to = length(kernel_parameters), by = 2), 
                          function(idx) if (kernel_parameters[idx] >= kernel_parameters[idx + 1]) stop('Inner radius must be smaller than outer radius.') 
                                        else k_mat_ring(coords = coords, h1 = kernel_parameters[idx], h2 = kernel_parameters[idx + 1]))  
  } else if (kernel_type == 'gauss') {
    kernel_list <- lapply(kernel_parameters, 
                          function(h) if (h >= 0) k_mat_exp(coords = coords, h = h)
                                      else stop('Parameter must be zero or positive.'))
  } else {
    stop('Invalid input. Note that the length(kernel_parameters) must be an even number for the ring kernel.')
  }
  
  return(kernel_list)
}

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

predict_idw <- function(vals, coords, p, n_grid) {
  coords_pred <- as.matrix(expand.grid(seq(from = floor(min(coords[, 1])), to = ceiling(max(coords[, 1])), length.out = n_grid), 
                                       seq(from = floor(min(coords[, 2])), to = ceiling(max(coords[, 2])), length.out = n_grid)))
  colnames(coords_pred) <- colnames(coords)
  vals_pred <- idw(coords_pred = coords_pred, coords_vals = coords, vals = vals, p = p)
  colnames(vals_pred) <- paste0(colnames(vals), '.pred')
  return(list(vals_pred_idw = vals_pred, coords_pred_idw = coords_pred))
}

#------------------------------------------------#
#  function for scatter diagonalization
#------------------------------------------------#
diag_scatters <- function(cov_list, ordered, ...) {
  decr <- if (attr(cov_list, 'lcov') == 'ldiff' && !is.null(attr(cov_list, 'lcov'))) FALSE else TRUE
  k <- length(cov_list)
  if (k == 1) {
    cov_evd <- eigen(cov_list[[1]], symmetric = TRUE)
    u <- cov_evd$vectors
    d <- diag(cov_evd$values)
  } else {
    jade <- JADE::frjd(do.call(rbind, cov_list), ...)
    u <- jade$V
    d <- jade$D
  }
  
  p <- ncol(d)
  if (ordered) {
    diags_mat <- matrix(0, nrow = k, ncol = p)
    for (idx in 1:k) {
      diags_mat[idx, ] <- diag(d[(1:p) + (idx - 1) * p, ])
    }
    diag_order <- order(colSums(diags_mat ^ 2), decreasing = decr)
    u <- u[, diag_order]
    for (idx in 1:k) {
      d[(1:p) + (idx - 1) * p, ] <- d[(1:p) + (idx - 1) * p, ][diag_order, diag_order]
    }
  } 
  
  return(list(u = u, d = d))
}

#------------------------------------------------#
#   test statistic for noise dimension tests
#------------------------------------------------#
test_stat <- function(d, q, n) {
  p <- ncol(d)
  k <- nrow(d) / p
  t <- 0
  
  for (idx in 1:k) {
    t <- t + sum(d[((q + 1):p) + (idx - 1) * p, (q + 1):p] ^ 2)
  }
  
  t <- t * n / 2
  return(t)
}

#------------------------------------------------#
#   determines the block center coordinates
#------------------------------------------------#
block_center_coords <- function(coords, n_block) {
  x_min <- min(coords[, 1])
  x_max <- max(coords[, 1])
  y_min <- min(coords[, 2])
  y_max <- max(coords[, 2])
  
  if (n_block == 'x') {
    dx <- (x_max - x_min) / 2
    block_center <- as.matrix(expand.grid((x_min + dx / 2) + (0:1) * dx, 
                                          (y_min + y_max) / 2))
  } else if (n_block == 'y') {
    dy <- (y_max - y_min) / 2
    block_center <- as.matrix(expand.grid((x_min + x_max) / 2,
                                          (y_min + dy / 2) + (0:1) * dy))    
  } else {
    dx <- (x_max - x_min) / n_block
    dy <- (y_max - y_min) / n_block
    block_center <- as.matrix(expand.grid((x_min + dx / 2) + (0:(n_block - 1)) * dx, 
                                          (y_min + dy / 2) + (0:(n_block - 1)) * dy))
  }
  
  colnames(block_center) <- c('x', 'y')
  return(block_center)
}

#------------------------------------------------#
#   prepares for each coordinate the block index
#------------------------------------------------#
make_blocks <- function(x, coords, n_block) {
  block_center <- block_center_coords(coords, n_block)
  block_idx <- idx_per_block(coords, block_center, 2L)$block_idx + 1
  
  unique_block_idx <- unique(block_idx)
  
  x_list <- vector(mode = 'list', length = length(unique_block_idx))
  coords_list <- vector(mode = 'list', length = length(unique_block_idx))
  for (idx in 1:length(unique_block_idx)) {
    x_list[[idx]] <- x[block_idx == unique_block_idx[idx], , drop = FALSE]
    coords_list[[idx]] <- coords[block_idx == unique_block_idx[idx], , drop = FALSE]
  }
  
  return(list(x_list = x_list, coords_list = coords_list))
}
