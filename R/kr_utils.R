#################################################################
##    Utils for Kenward-Roger Degrees of Freedom Adjustment    ##
##                    Using pbkrtest package                   ##
#################################################################

#' @importFrom nlme VarCorr
#' @export
VarCorr.mmrm <- function(model){
  vars = .vars(model)
  cors = .varcor(model)

  vc = list()
  vc[[1]] = varcov(model)
  dimnames(vc[[1]]) = dimnames(cors)
  attr(vc[[1]], "stddev") = sqrt(vars)
  attr(vc[[1]], "correlation") = cors
  vc
}


##---------------------------------------------------------------
##  Adapted from https://github.com/hojsgaard/pbkrtest/pull/2   -
##---------------------------------------------------------------

#' @importFrom pbkrtest vcovAdj
#' @export
vcovAdj.mmrm <- pbkrtest:::vcovAdj.lmerMod

#' @importFrom lme4 getME
#' @export
getME.mmrm <- function(object, name, ...){
  groups <- object$groups
  glsSt <- object$modelStruct$corStruct
  model <- object$model
  mfArgs <- list(formula = nlme::asOneFormula(formula(glsSt), model, groups),
                 data = nlme::getData(object))
  mfArgs$drop.unused.levels <- TRUE
  dataMod <- do.call(model.frame, mfArgs)
  origOrder <- row.names(dataMod)	# preserve the original order
  if (!is.null(groups)) {
    grps <- groups
    ## ordering data by groups
    ord <- order(grps)
    grps <- grps[ord]
    dataMod <- dataMod[ord, ,drop = FALSE]
    revOrder <- match(origOrder, row.names(dataMod)) # putting in orig. order
    ugroups <- unique(grps)
  } else grps <- NULL

  X_raw <- model.matrix(formula(object), data=nlme::getData(object))
  X_sorted <- X_raw[ord,]

  if(name=='X'){
    return(X_raw)
  }
  if(name=='X_sorted'){
    return(X_sorted)
  }
  if(name=='is_REML'){
    return(object$method=='REML')
  }
  if(name=='Zt'){
    inter_obs = interaction(dataMod)
    missed = which(!sapply(levels(inter_obs), function(x) x %in% inter_obs))

    dims = attr(glsSt, "Dim")
    cn = order(revOrder)
    rn = 1:dims$N
    for (i in missed) rn[rn>=i] = rn[rn>=i] + 1

    Zt = Matrix::sparseMatrix(rn, cn, x=1, dims=c(dims$M*dims$maxLen, dims$N))
    return(Zt)
  }
  if(name=='X_star'){
    # Pinheiro & Bates p 202
    invsqrtLambda <- lapply(ugroups, function(i) solve(pbkrtest:::.sqrtMat(nlme::getVarCov(object, individual = i)/(sigma( object )^2))))
    X_star   <- matrix(0, nrow=nrow(X_raw), ncol=ncol(X_raw))
    for(i in 1:length(ugroups)){
      X_star[groups==ugroups[i], ] <- t(invsqrtLambda[[i]]) %*% X_sorted[grps==ugroups[i],]
    }
    return(Matrix(X_star, sparse=TRUE))
  }
  if(name=='Zt_star'){
    # Pinheiro & Bates p 202
    # Zt is (n_reff x n_subjects) rows by N cols, e.g. 1940 x 4790
    invsqrtLambda <- lapply(ugroups, function(i) solve(pbkrtest:::.sqrtMat(nlme::getVarCov(object, individual = i)/(sigma( object )^2))))
    # Pinheiro & Bates p 202
    Zt_star <- matrix(0, nrow=length(ugroups), ncol=nrow(X_raw))
    for(i in 1:length(ugroups)){
      Zt_star[i, groups==ugroups[i]] <- t(as.matrix(rep(1, sum(groups==ugroups[i])))) %*% invsqrtLambda[[i]]
    }
    return(Matrix(Zt_star, sparse=TRUE))
  }
}

#' @importFrom pbkrtest get_SigmaG
#' @export
get_SigmaG.mmrm <- function(object, details=0) {

  DB <- details > 0 ## For debugging only

  # variance-covariance of random effects
  GGamma <- nlme::VarCorr(object)

  ## Put covariance parameters for the random effects into a vector:
  ## Fixme: It is a bit ugly to throw everything into one long vector here; a list would be more elegant
  Lii <- GGamma[[1]]
  ggamma = Lii[ lower.tri( Lii, diag=TRUE ) ]
  n.ggamma <- length(ggamma)

  ## Find G_r:
  G  <- NULL
  ZZ <- getME(object, "Zt")
  n.lev = length(unique(object$groups))
  Ig    <- Matrix::sparseMatrix(1:n.lev, 1:n.lev, x=1)

  n.comp.by.RT = attr(object$modelStruct$corStruct, "Dim")$maxLen
  n.parm.by.RT = (n.comp.by.RT + 1) * n.comp.by.RT / 2
  for (rr in 1:n.parm.by.RT) {
    ii.jj = pbkrtest:::.index2UpperTriEntry(rr, n.comp.by.RT)
    ii.jj = unique(ii.jj)
    if (length(ii.jj)==1){
      EE <- Matrix::sparseMatrix(ii.jj, ii.jj, x=1, dims=rep(n.comp.by.RT, 2))
    } else {
      EE <- Matrix::sparseMatrix(ii.jj, ii.jj[2:1], dims=rep(n.comp.by.RT, 2))
    }
    EE <- Ig %x% EE  ## Kronecker product
    G  <- c( G, list( Matrix::t(ZZ) %*% EE %*% ZZ ) )
  }

  Sigma <- ggamma[1] * G[[1]]
  for (ii in 2:n.ggamma) {
    Sigma <- Sigma + ggamma[ii] * G[[ii]]
  }

  SigmaG <- list(Sigma=Sigma, G=G, n.ggamma=n.ggamma)
  SigmaG
}
