#' @export
covMat <- function(object) {
  vars <- .vars(object)
  cors <- .cormat(object)
  vdiag <- diag(x = sqrt(vars), nrow(cors), ncol(cors))
  covs <- vdiag %*% cors %*% vdiag
  rownames(covs) <- names(vars)
  colnames(covs) <- names(vars)
  return(covs)
}

.vars <- function(object) {
  if (object$heterogenous) {
    var_order <- order(attr(object$modelStruct$varStruct, "groupNames"))
    vars <- (coef(object$modelStruct$varStruct, uncons = FALSE, allCoef = TRUE)^2 *
      object$sigma^2)[var_order]
    vars
  } else {
    sigma(object)^2
  }
}

.cormat <- function(object) {
  covs <- sort(unique(unlist(attr(object$modelStruct$corStruct, "covariate"))))
  if (is.null(covs)) covs <- sort(unique(object$data[[object$time]]))
  nlme::corMatrix(object$modelStruct$corStruct, covariate = covs)
}
