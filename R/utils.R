white_data <- function(x) {
  n <- nrow(x)
  
  mu <- colMeans(x)
  x_0 <- sweep(x, MARGIN = 2, STATS = mu, FUN = '-')
  
  s <- crossprod(x_0) / (n - 1)
  
  s_evd <- eigen(s, symmetric = TRUE)
  s_inv_sqrt <- s_evd$vectors %*% tcrossprod(diag(1 / sqrt(s_evd$values)), s_evd$vectors)
  s_sqrt <- s_evd$vectors %*% tcrossprod(diag(sqrt(s_evd$values)), s_evd$vectors)
  
  x_w <- tcrossprod(x_0, s_inv_sqrt)
  colnames(x_w) <- colnames(x_0)
  
  return(list(mu = mu, x_0 = x_0, x_w = x_w, s_inv_sqrt = s_inv_sqrt, s_sqrt = s_sqrt))
}

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

local_covariance_matrix <- function(x, kernel_list, whitening = TRUE) {
  if (whitening) {
    x_w <- white_data(x)
    cov_sp_list <- lapply(kernel_list, sp_cov_mat_sparse, x = x_w$x_w)
  } else {
    cov_sp_list <- lapply(kernel_list, sp_cov_mat_sparse, x = x)
  }
  cov_sp_list <- lapply(cov_sp_list, function(x) (x + t(x)) / 2)
  
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
