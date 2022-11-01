#################################################################
##    Utils for Kenward-Roger Degrees of Freedom Adjustment    ##
##                    Using pbkrtest package                   ##
#################################################################

## ---------------------------------------------------------------
##  Adapted from pbkrtest and merDeriv packages
##  https://github.com/hojsgaard/pbkrtest
##  https://github.com/nctingwang/merDeriv
## ---------------------------------------------------------------

tr <- function(x) sum(Matrix::diag(x))

#' @importFrom pbkrtest vcovAdj
#' @import Matrix
#' @export
vcovAdj.mmrm <- function(object, information = "expected") {
  Phi <- vcov(object)
  Sigma <- .get_Sigma(object)
  dSigma <- .get_dSigma(object)
  n.ggamma <- length(dSigma)
  X <- model.matrix(object, object$data)
  SigmaInv <- chol2inv(chol(forceSymmetric(Sigma)))

  dSigmaInv <- lapply(dSigma, function(x) -SigmaInv %*% x %*% SigmaInv)
  P <- lapply(dSigmaInv, function(x) t(X) %*% x %*% X)
  Q <- list()
  for (i in 1:length(dSigmaInv)) {
    dSigmaInv_i <- dSigmaInv[[i]]
    Q[[i]] <- lapply(dSigmaInv, function(x) t(X) %*% dSigmaInv[[i]] %*% Sigma %*% x %*% X)
  }

  info_exp2 <- info_avg2 <- matrix(0, nrow = length(dSigmaInv), ncol = length(dSigmaInv))
  if (information %in% c("expected", "observed")) {
    for (i in 1:nrow(info_exp2)) {
      for (j in 1:ncol(info_exp2)) {
        info_exp2[i, j] <- tr(dSigmaInv[[i]] %*% Sigma %*% dSigmaInv[[j]] %*% Sigma) -
          tr(2 * Phi %*% Q[[i]][[j]] - Phi %*% P[[i]] %*% Phi %*% P[[j]])
      }
    }
  }

  if (information %in% c("average", "observed")) {
    r_hat <- object$residuals
    otherP <- SigmaInv -
      SigmaInv %*% X %*% solve(t(X) %*% SigmaInv %*% X) %*% t(X) %*% SigmaInv
    for (i in 1:nrow(info_avg2)) {
      for (j in 1:ncol(info_avg2)) {
        info_avg2[i, j] <- 2 * as.numeric(
          r_hat %*% otherP %*% dSigma[[i]] %*% otherP %*% dSigma[[j]] %*% otherP %*% r_hat
        )
      }
    }
  }

  if (information == "observed") {
    info2 <- -info_exp2 + info_avg2
  } else {
    info2 <- info_exp2 + info_avg2
  }

  eig_info2 <- eigen(info2, only.values = TRUE)$values
  condi <- min(abs(eig_info2))
  W <- if (condi > 1e-10) forceSymmetric(2 * solve(info2)) else forceSymmetric(2 * MASS::ginv(info2))

  U <- matrix(0, nrow = ncol(X), ncol = ncol(X))
  for (i in 1:ncol(W)) {
    for (j in 1:ncol(W)) {
      U <- U + W[i, j] * (Q[[i]][[j]] - P[[i]] %*% Phi %*% P[[j]])
    }
  }

  Lambda <- Phi %*% U %*% Phi
  PhiA <- Phi + 2 * Lambda
  attr(PhiA, "P") <- P
  attr(PhiA, "W") <- W
  attr(PhiA, "condi") <- condi
  PhiA
}
