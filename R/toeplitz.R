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
#'               specifying a time covariate t and, optionally,
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
corMatrix.corToep <- function(object, covariate = nlme::getCovariate(object), corr = TRUE, ...) {
  corD <- nlme::Dim(object, if (is.list(covariate)) {
    if (is.null(names(covariate))) {
      names(covariate) <- seq_along(covariate)
    }
    rep(names(covariate), lengths(covariate))
  } else {
    rep(1, length(covariate))
  })

  raw_params <- as.vector(object)
  params <- sign(raw_params) * (exp(raw_params) - 1) / (1 + exp(raw_params))

  fac_inv_mat = matrix(0, nrow=length(params)+1, ncol=length(params)+1)
  fac_inv_mat[, 1] = c(1, params)
  for (rr in 2:nrow(fac_inv_mat)) {
    for (cc in 2:rr) {
      if (rr == cc) {
        fac_inv_mat[rr, cc] = sqrt(1 - sum(fac_inv_mat[rr, 1:(rr-1)]^2))
      } else {
        fac_inv_mat[rr, cc] = (
          fac_inv_mat[rr-cc+1, 1] -
            sum(fac_inv_mat[rr, 1:(cc-1)] * fac_inv_mat[cc, 1:(cc-1)])
        ) / fac_inv_mat[cc, cc]
      }
    }
  }

  toep_mat = fac_inv_mat %*% t(fac_inv_mat)
  fac_mat = solve(fac_inv_mat)
  all_covs <- sort(unique(unlist(covariate)))
  if (is.list(covariate)) {
    toep_list <- list()
    for (i in 1:length(covariate)) {
      this_covs <- sapply(all_covs, function(x) x %in% covariate[[i]])
      if (corr) {
        toep_list[[i]] <- toep_mat[this_covs, this_covs]
      } else {
        toep_list[[i]] = fac_mat[this_covs, this_covs]
      }
    }
    lD <- if (corr) NULL else sum(log(sapply(toep_list, function(x) {
      if (inherits(x, "matrix")) det(x) else 1
    })))
    val <- toep_list
  } else {
    this_covs <- sapply(all_covs, function(x) x %in% covariate)
    if (corr) {
      val <- toep_mat[this_covs, this_covs]
    } else {
      val <- fac_mat[this_covs, this_covs]
    }

    lD <- if (corr) NULL else try(log(det(val)))
  }

  if (corD[["M"]] > 1) {
    names(val) <- names(corD[["len"]])
    val <- as.list(val)
  }
  attr(val, "logDet") <- lD
  val
}

#' @export
coef.corToep <- function(object, unconstrained = TRUE, ...) {
  if (unconstrained) {
    if (attr(object, "fixed")) {
      return(numeric(0))
    } else {
      return(as.vector(object))
    }
  }
  aux <- exp(as.vector(object))
  aux <- c((aux - 1) / (aux + 1))
  names(aux) <- paste0("Rho", 1:length(aux))
  aux
}

#' @importFrom nlme corFactor
#' @export
corFactor.corToep <- function(object, ...) {
  corD <- nlme::Dim(object)
  if (corD[["sumLenSq"]] > .Machine$integer.max) {
    stop(gettextf(
      "'sumLenSq' = %g is too large (larger than maximal integer)",
      corD[["sumLenSq"]]
    ), domain = NA)
  }
  cmats <- nlme::corMatrix(object, corr = F)
  val <- unlist(cmats)
  attr(val, "logDet") <- attr(cmats, "logDet")
  val
}

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
  attr(object, "maxCov") <- length(unique(covar))

  params <- as.vector(object)
  if (length(params) == 0) {
    oldAttr <- attributes(object)
    object <- double(attr(object, "maxCov") - 1)
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
