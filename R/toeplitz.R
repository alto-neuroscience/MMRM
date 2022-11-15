#' Toeplitz Correlation Structure
#'
#' @description
#' This function is a constructor for the `corToep` class,
#' representing a toeplitz correlation structure.
#' Objects created using this constructor must later be initialized
#' using the appropriate `Initialize` method.
#'
#' @param value an optional vector of parameter values.
#'                Must have length number of levels - 1.
#' @param form a one sided formula of the form ~ t, or ~ t | g,
#'               specifying a time covariate and, optionally,
#'               a grouping factor g.
#' @param fixed an optional logical value indicating whether the coefficients
#'                should be allowed to vary in the optimization, or kept fixed
#'                at their initial value. Defaults to FALSE, in which case the
#'                coefficients are allowed to vary.
#'
#' @export
corToep <- function(value = numeric(0), form = ~1, fixed = FALSE) {
  attr(value, "formula") <- form
  attr(value, "fixed") <- fixed
  class(value) <- c("corToep", "corStruct")
  value
}

#' @importFrom nlme corMatrix
#' @export
corMatrix.corToep <- nlme:::corMatrix.corARMA

#' @export
coef.corToep <- function (object, unconstrained = TRUE, ...) {

  if (unconstrained) {
    return (nlme:::coef.corARMA(object, unconstrained, ...))
  }

  covs = sort(unique(unlist(attr(object, "covariate"))))
  cor_mat = corMatrix(object, covariate=covs)
  cor_pars = cor_mat[2:nrow(cor_mat), 1]
  names(cor_pars) = paste0("Phi", 1:length(cor_pars))
  cor_pars
}

#' @importFrom nlme corFactor
#' @export
corFactor.corToep <- nlme:::corFactor.corARMA

#' @importFrom nlme Initialize
#' @export
Initialize.corToep <- function(object, data, ...) {
  object <- NextMethod()
  covar <- attr(object, "covariate")
  if (!is.list(covar)) {
    covar <- list(covar)
  }
  if (any(unlist(lapply(covar, duplicated)))) {
    stop("covariate must have unique values within groups for \"corToep\" objects")
  }
  covar <- unlist(covar) - 1
  if (nlme::Dim(object)[["M"]] > 1) {
    attr(object, "covariate") <- split(covar, nlme::getGroups(object))
  } else {
    attr(object, "covariate") <- covar
  }
  attr(object, "maxLag") <- length(unique(covar)) - 1
  attr(object, "p") <- length(unique(covar)) - 1
  attr(object, "q") <- 0

  params <- as.vector(object)
  if (length(params) == 0) {
    oldAttr <- attributes(object)
    object <- double(attr(object, "maxLag"))
    attributes(object) <- oldAttr
    attr(object, "factor") <- corFactor(object)
    attr(object, "logDet") <- -attr(
      attr(object, "factor"),
      "logDet"
    )
    attr(object, "class") <- class(object)
  }
  object
}
