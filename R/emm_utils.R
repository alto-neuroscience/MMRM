##################################################################
##                  emmeans emm_basis functions                 ##
##           to support kenward-rogers with glsObject           ##
##################################################################

.emm_basis_mmrm_kr <- function(object, trms, xlev, grid, mode = "kenward-rogers",
                               extra.iter = 0, options, misc,
                               pbkrtest.limit = emmeans::get_emm_option("pbkrtest.limit"),
                               ...) {
  contrasts <- object$contrasts
  m <- stats::model.frame(trms, grid, na.action = stats::na.pass, xlev = xlev)
  X <- stats::model.matrix(trms, m, contrasts.arg = contrasts)

  mf <- stats::model.frame(trms, attr(object, "data"), na.action = na.pass, xlev = xlev)
  mm <- stats::model.matrix(trms, mf, contrasts.arg = contrasts)
  mm <- emmeans::.cmpMM(mm, assign = attr(mm, "assign"))

  tmp <- stats::coef(object)
  bhat <- rep(NA, ncol(X))
  bhat[match(names(tmp), colnames(X), nomatch = 0)] <- tmp
  V <- emmeans::.my.vcov(object, ...)
  if (any(is.na(bhat))) {
    nbasis <- estimability::nonest.basis(mm)
  } else {
    nbasis <- estimability::all.estble
  }

  objN <- object$dims$N

  if (!requireNamespace("pbkrtest", quietly = TRUE)) {
    stop(
      "pbkrtest package is not installed! ",
      "Please use a different mode (e.g., mode = 'satterthwaite') ",
      "or install the package."
    )
  } else if (object$dims$N > pbkrtest.limit) {
    stop(
      "The number of observations exceeds ", pbkrtest.limit, ".\n",
      "To proceed anyways, add the argument pbkrtest.limit = ", object$dims$N, " (or larger)\n",
      "[or, globally, 'set emm_options(pbkrtest.limit = ", object$dims$N, ")' or larger];\n",
      "but be warned that this may result in large computation time and memory use."
    )
  }

  dfargs <- list(
    unadjV = V,
    adjV = pbkrtest::vcovAdj.gls(object, 0)
  )
  V <- as.matrix(dfargs$adjV)
  tst <- try(pbkrtest::Lb_ddf)
  if (!inherits(tst, "try-error")) {
    dffun <- function(k, dfargs) {
      pbkrtest::Lb_ddf(k, dfargs$unadjV, dfargs$adjV)
    }
  } else {
    stop("Failed to load pbkrtest routines!")
  }
  attr(dffun, "mesg") <- mode

  list(
    X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, dfargs = dfargs, misc = misc,
    model.matrix = mm
  )
}

emm_basis.mmrm <- function(object, trms, xlev, grid, mode="kenward",
                           extra.iter = 0, options, misc, ...) {
  if (grepl("kenward", tolower(mode)) | tolower(mode) == "KR") {
    mode <- "kenward-rogers"
    .emm_basis_mmrm_kr(object, trms, xlev, grid, mode,
                       extra.iter = 0, options, misc, ...
    )
  } else {
    emmeans:::emm_basis.gls(object, trms, xlev, grid, mode,
                            extra.iter = 0, options, misc, ...
    )
  }
}

recover_data.mmrm <- function(object, data, ...)
{
  fcall = object$call
  if (!is.null(wts <- fcall$weights)) {
    wts = nlme::varWeights(object$modelStruct)
    fcall$weights = NULL
  }
  trms = delete.response(terms(nlme::getCovariateFormula(object)))
  result = emmeans:::recover_data.call(fcall, trms, object$na.action,
                              data = object$data, ...)
  if (!is.null(wts))
    result[["(weights)"]] = wts
  if (!missing(data))
    attr(result, "misc") = list(data = data)
  attr(result, "pass.it.on") = TRUE
  result
}
