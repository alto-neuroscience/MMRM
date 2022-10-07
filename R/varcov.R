#' Get Variance-Covariance matrix of MMRM model
#'
#' @description
#' Compute the variance-covariance matrix of an mmrmObject.
#' Copy of [VarCov_gls from the SPR package](https://cjangelo.github.io/SPR/reference/VarCov_gls.html)
#'
#' @param model mmrmObject (output of [mmrm])
#'
#' @returns matrix, the variance-covoariance matrix
#'
#' @seealso [SPR::VarCov_gls]
#'
#' @export
varcov <- function(model) {
  vars = .vars(model)
  cors = .varcor(model)
  r <- coef(model$modelStruct$corStruct, uncons = FALSE, allCoef = TRUE)
  covs <- diag(sqrt(vars)) %*% cors %*% diag(sqrt(vars))
  rownames(covs) <- names(vars)
  colnames(covs) <- names(vars)
  return(covs)
}

.varcor <- function(model) {
  r <- coef(model$modelStruct$corStruct, uncons = FALSE, allCoef = TRUE)
  N <- .get_levels(model)
  cors <- diag(1, N)
  cors[lower.tri(cors)] <- r
  cors[upper.tri(cors)] <- t(cors)[upper.tri(t(cors))]
  nms = names(.vars(model))
  rownames(cors) = nms
  colnames(cors) = nms
  cors
}

.vars <- function(model) {
  var_order = order(attr(model$modelStruct$varStruct, "groupNames"))
  vars <- (coef(model$modelStruct$varStruct, uncons = FALSE, allCoef = TRUE)^2 *
             model$sigma^2)[var_order]
  vars
}

.get_levels <- function(model) {
  if (is(model, "mmrm")) {
    model = model$modelStruct
  }
  length(attr(model$varStruct, "groupNames"))
}
