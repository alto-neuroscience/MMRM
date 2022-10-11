##################################################################
##                  emmeans emm_basis functions                 ##
##          to support kenward-rogers with mmrmObject           ##
##################################################################


## -----------------------------------------------
##  Adapted from emmeans:::emm_basis.gls
## -----------------------------------------------

#' @importFrom emmeans emm_basis
#' @export
emm_basis.mmrm <- function(object, trms, xlev, grid, mode = c(
                             "auto", "kenward-rogers", "df.error",
                             "satterthwaite", "appx-satterthwaite", "boot-satterthwaite",
                             "asymptotic"
                           ),
                           extra.iter = 0, options, misc, information = "expected",
                           pbkrtest.limit = emmeans::get_emm_option("pbkrtest.limit"),
                           force_mode = FALSE,
                           ...) {
  misc$data <- object$data
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

  # determine mode, run checks

  mode <- match.arg(mode)
  if (mode == "boot-satterthwaite") mode <- "appx-satterthwaite"
  if (!is.null(options$df)) mode <- "df.error"
  if (mode == "auto") mode <- "kenward"

  if ((!is.matrix(object$apVar)) & (!force_mode)) {
    warning(
      m$apVar,
      ". Using mode = ",
      "df.error! If you would still like to try mode = ",
      mode,
      ", please set the argument 'force_mode=TRUE'"
    )
    mode <- "df.error"
  }

  # apply mode

  if (mode == "kenward-rogers") {
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
      adjV = vcovAdj.mmrm(object, information = information)
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
  } else if (mode %in% c("satterthwaite", "appx-satterthwaite")) {
    data <- misc$data
    misc <- list()
    chk <- attr(object$apVar, "Pars")
    if ((max(abs(coef(object$modelStruct) - chk[-length(chk)])) >
      0.001) & mode == "satterthwaite") {
      message("Analytical Satterthwaite method not available; using appx-satterthwaite")
      mode <- "appx-satterthwaite"
    }
    if (mode == "appx-satterthwaite") {
      G <- try(emmeans:::gradV.kludge(object, "varBeta",
        call = object$call$model,
        data = data, extra.iter = extra.iter
      ), silent = TRUE)
    } else {
      G <- try(emmeans:::gls_grad(object, object$call, data, V))
    }
    if (inherits(G, "try-error")) {
      sugg <- ifelse(mode == "satterthwaite", "appx-satterthwaite",
        "df.error"
      )
      stop("Can't estimate Satterthwaite parameters.\n",
        "  Try adding the argument 'mode = \"", sugg,
        "\"'",
        call. = FALSE
      )
    }
    dfargs <- list(V = V, A = object$apVar, G = G)
    dffun <- function(k, dfargs) {
      est <- tcrossprod(crossprod(k, dfargs$V), k)
      g <- sapply(dfargs$G, function(M) {
        tcrossprod(crossprod(
          k,
          M
        ), k)
      })
      varest <- tcrossprod(crossprod(g, dfargs$A), g)
      2 * est^2 / varest
    }
  } else if (mode %in% c("df.error", "asymptotic")) {
    df <- ifelse(mode == "asymptotic", Inf, object$dims$N -
      object$dims$p - length(unlist(object$modelStruct)))
    dfargs <- list(df = df)
    dffun <- function(k, dfargs) dfargs$df
  }

  # return result

  attr(dffun, "mesg") <- mode

  list(
    X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, dfargs = dfargs, misc = misc,
    model.matrix = mm
  )
}

recover_data.mmrm <- function(object, data, ...) {
  fcall <- object$call
  if (!is.null(wts <- fcall$weights)) {
    wts <- nlme::varWeights(object$modelStruct)
    fcall$weights <- NULL
  }
  trms <- delete.response(terms(nlme::getCovariateFormula(object)))
  result <- emmeans:::recover_data.call(fcall, trms, object$na.action,
    data = object$data, ...
  )
  if (!is.null(wts)) {
    result[["(weights)"]] <- wts
  }
  if (!missing(data)) {
    attr(result, "misc") <- list(data = data)
  }
  attr(result, "pass.it.on") <- TRUE
  result
}
