#################################################################
##             Covariance matrix utility functions             ##
#################################################################

.to_diag_block <- function(x, n) {
  block_mat <- as(matrix(0, nrow = nrow(x) * n, ncol = ncol(x) * n), "sparseMatrix")
  for (i in 1:n) {
    ind <- (i - 1) * nrow(x) + 1
    block_mat[ind:(ind + nrow(x) - 1), ind:(ind + ncol(x) - 1)] <- x
  }
  block_mat
}

.get_Zt <- function(object) {
  inter_obs <- interaction(object$data[[object$time]], object$data[[object$subject]])
  missed <- which(!sapply(levels(inter_obs), function(x) x %in% inter_obs))

  dims = nlme::Dim(object$modelStruct$corStruct, nlme::getGroups(object))
  # dims <- attr(object$modelStruct$corStruct, "Dim")
  grps <- nlme::getGroups(
    object$data,
    eval(substitute(~ 1 | SUBJ, list(SUBJ = object$subjects)))
  )
  lvls <- get_levels(object)

  rn <- 1:dims$N
  for (i in missed) rn[rn >= i] <- rn[rn >= i] + 1

  Zt <- Matrix::sparseMatrix(rn, order(grps), x = 1, dims = c(dims$M * lvls, dims$N))
  return(Zt)
}

.get_Sigma <- function(object) {
  cov_chol <- as.matrix(chol(Matrix::forceSymmetric(covMat(object))))
  Lamt <- .to_diag_block(cov_chol, n = nlme::Dim(object$modelStruct$corStruct, nlme::getGroups(object))$M)

  Zt <- .get_Zt(object)
  Sigma <- t(Zt) %*% t(Lamt) %*% Lamt %*% Zt
  return(Sigma)
}

.get_dSigma <- function(object) {
  cor_struct <- class(object$modelStruct$corStruct)[1]
  dsigma_fn <- try(get(paste0(".dSigma_", cor_struct)), silent = T)
  if (!inherits(dsigma_fn, "function")) {
    stop(
      "Derivatives for correlation structure ",
      cor_struct,
      " is not supported!"
    )
  }

  dSig <- dsigma_fn(object)
  dSig_block <- lapply(dSig, .to_diag_block, n = nlme::Dim(object$modelStruct$corStruct, nlme::getGroups(object))$M)

  Zt <- .get_Zt(object)
  dSigma <- lapply(dSig_block, function(x) t(Zt) %*% x %*% Zt)
  return(dSigma)
}

.dSigma_list <- function(par_mat) {
  res_list <- list()
  m <- matrix(0, nrow = nrow(par_mat), ncol = ncol(par_mat))
  for (i in 1:max(par_mat)) {
    this_m <- m
    this_m[par_mat == i] <- 1
    res_list[[i]] <- this_m
  }
  res_list
}

.dSigma_corSymm <- function(object) {
  lvls <- get_levels(object)
  par_mat <- matrix(NA, nrow = lvls, ncol = lvls)

  n_var_pars <- length(.vars(object))
  diag(par_mat) <- 1:n_var_pars
  par_mat[lower.tri(par_mat)] <- n_var_pars + 1:sum(lower.tri(par_mat))
  par_mat[upper.tri(par_mat)] <- t(par_mat)[upper.tri(par_mat)]
  .dSigma_list(par_mat)
}

.dSigma_corToep <- function(object) {
  lvls <- get_levels(object)
  var_pars <- .vars(object)
  cor_pars <- coef(object$modelStruct$corStruct, uncons = F)

  par_ind_mat <- sapply(1:lvls, function(x) abs(1:lvls - x))
  v_mat <- c(1, cor_pars)[par_ind_mat + 1]
  dim(v_mat) <- dim(par_ind_mat)

  res_list <- list()
  if (object$heterogenous) {
    for (i in 1:length(var_pars)) {
      d_mat <- matrix(0, nrow = lvls, ncol = lvls)
      d_mat[, i] <- d_mat[i, ] <- sqrt(var_pars) * 1 / (2 * sqrt(var_pars[i]))
      d_mat[i, i] <- 1
      res_list <- c(res_list, list(v_mat * d_mat))
    }
    rho_mat <- tcrossprod(sqrt(var_pars), sqrt(var_pars))
  } else {
    res_list <- c(res_list, list(v_mat))
    rho_mat <- matrix(var_pars, nrow = lvls, ncol = lvls)
  }

  for (i in 1:length(cor_pars)) {
    res_list <- c(res_list, list(as.integer(par_ind_mat == i) * rho_mat))
  }

  return(res_list)
}

.dSigma_corAR1 <- .dSigma_corARMA <- function(object) {
  lvls <- get_levels(object)
  var_pars <- .vars(object)
  cor_par <- coef(object$modelStruct$corStruct, uncons = F)

  rho_exp_mat <- sapply(1:lvls, function(x) abs(1:lvls - x))
  rho_mat <- rho_exp_mat * cor_par^(abs(rho_exp_mat - 1))
  v_mat <- cor_par^rho_exp_mat
  if (object$heterogenous) {
    res_list <- list()
    for (i in 1:length(var_pars)) {
      d_mat <- matrix(0, nrow = lvls, ncol = lvls)
      d_mat[, i] <- d_mat[i, ] <- sqrt(var_pars) * 1 / (2 * sqrt(var_pars[i]))
      d_mat[i, i] <- 1
      res_list <- c(res_list, list(v_mat * d_mat))
    }
    var_prod <- tcrossprod(sqrt(var_pars), sqrt(var_pars))
    res_list <- c(res_list, list(var_prod * rho_mat))
  } else {
    res_list <- list(v_mat, (var_pars * rho_mat))
  }

  return(res_list)
}

.dSigma_corCompSymm <- function(object) {
  lvls <- get_levels(object)
  n_var_pars <- length(.vars(object))
  par_mat <- matrix(n_var_pars + 1, nrow = lvls, ncol = lvls)
  diag(par_mat) <- 1:n_var_pars
  .dSigma_list(par_mat)
}

.dSigma_corIdent <- function(object) {
  n_tp <- length(unique(object$data[[object$time]]))
  if (object$heterogenous) {
    par_mat <- diag(1:n_tp)
  } else {
    par_mat <- diag(1, n_tp, n_tp)
  }
  .dSigma_list(par_mat)
}
