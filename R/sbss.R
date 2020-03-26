#-------------------------------------------#
# sbss generic and methods
#-------------------------------------------#
sbss <- function(x, ...) UseMethod("sbss")

sbss.default <- function(x, coords, kernel_type = c('ring', 'ball', 'gauss'), kernel_parameters, ordered = TRUE, kernel_list = NULL, ...) {
  # kernel matrix
  kernel_type <- match.arg(kernel_type)
  
  if (!missing(coords) && !missing(kernel_parameters) && is.vector(kernel_parameters)) {
    kernel_list <- spatial_kernel_matrix(coords, kernel_type = kernel_type, kernel_parameters = kernel_parameters)
  } else if (!is.null(kernel_list)) {
    coords <- NULL
  } else {
    stop('Invalid input for kernels. Either coords (or a spatial object) kernel_type and kernel_parameters (as vector) or kernel_list needs to be given.')
  }
  
  k <- length(kernel_list)
  
  # standardize data
  x_w <- white_data(x)
  
  # spatial covariance matrices
  cov_sp_list <- local_covariance_matrix(x = x_w$x_w, kernel_list = kernel_list, whitening = FALSE)
  
  # diagonalization
  if (k == 1) {
    cov_sp_evd <- eigen(cov_sp_list[[1]], symmetric = TRUE)
    u <- cov_sp_evd$vectors
    d <- diag(cov_sp_evd$values)
    
  } else {
    jade <- JADE::frjd(do.call(rbind, cov_sp_list), ...)
    u <- jade$V
    d <- jade$D
  }
  
  # unmixing matrix
  w <- crossprod(u, x_w$s_inv_sqrt)
  w_inv <- crossprod(x_w$s_sqrt, u)
  s <- tcrossprod(x_w$x_0, w)
  
  # ordering by squared (pseudo) eigenvalues
  p <- ncol(x)
  if (ordered) {
    diags_mat <- matrix(0, nrow = k, ncol = p)
    for (idx in 1:k) {
      diags_mat[idx, ] <- diag(d[(1:p) + (idx - 1) * p, ])
    }
    diag_order <- order(colSums(diags_mat ^ 2), decreasing = TRUE)
    
    u <- u[, diag_order]
    w <- w[diag_order, ]
    w_inv <- w_inv[, diag_order]
    s <- s[, diag_order]
    
    for (idx in 1:k) {
      d[(1:p) + (idx - 1) * p, ] <- d[(1:p) + (idx - 1) * p, ][diag_order, diag_order]
    }
  }
  colnames(s) <- paste0('IC.', 1:p)
  
  # results
  return(structure(list(s = s, coords = coords, w = w, w_inv = w_inv, d = d, x_mu = x_w$mu, cov_inv_sqrt = x_w$s_inv_sqrt), class = "sbss"))
}

sbss.SpatialPointsDataFrame <- function(x, ...) {
  result <- sbss.default(x = as.matrix(x@data), coords = x@coords, ...)
  x@data <- data.frame(result$s)
  result$s <- x
  result$coords <- NULL
  return(result)
}

sbss.sf <- function(x, ...) {
  if (!requireNamespace('sf', quietly = TRUE)) {
    stop('Please install the package sf to use this function.')
  } else {
    result <- sbss.default(x = as.matrix(sf::st_drop_geometry(x)), coords = sf::st_coordinates(x), ...)
    result$s <- sf::st_set_geometry(x = data.frame(result$s), value = sf::st_geometry(x))
    result$coords <- NULL
    return(result)    
  }
}



