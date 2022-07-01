#' Get Variance-Covariance matrix of MMRM model
#'
#' Copied from \link[SPR]{VarCor_gls}
#'
#' @param model a glsObject (output of \link[MMRM]{mmrm})
#'
#' @returns matrix, the variance-covoariance matrix
#'
#' @export
varcov <- function (model) {
  vars <- coef(model$modelStruct$varStruct, uncons = FALSE, allCoef = TRUE)^2 *
    model$sigma^2
  r = coef(model$modelStruct$corStruct, uncons = FALSE, allCoef = TRUE)
  cors = matrix(NA, ncol = length(vars), nrow = length(vars))
  cors[lower.tri(cors)] = r
  cors[upper.tri(cors)] = t(cors)[upper.tri(t(cors))]
  diag(cors) = rep(1, length(vars))
  covs <- diag(sqrt(vars)) %*% cors %*% diag(sqrt(vars))
  rownames(covs) <- names(vars)
  colnames(covs) <- names(vars)
  return(covs)
}
