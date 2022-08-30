#' Calculate effect size for contrasts using estimated marginal means
#'
#' @description
#' Calculate timepoint by timepoint effect size for MMRM models as the
#' least-squares means / sqrt(residual variance) at each timepoint.
#' Calculated by running [emmeans::eff_size] for each timepoint contrast,
#' where the timepoint by timepoint residuals come from the diagonals of [varcov].
#' Uncertainty is estimated by [emmeans::eff_size].
#'
#' @param mmrm_model the mmrm model(s) of type `mmrm`, `mmrmList`, or `mmrmCV`
#'                     used to calculate effect size
#' @param mmrm_emm estimated marginal means from and MMRM.
#'                     output of [mmrm_emmeans]
#' @param edf equivalent degrees of freedom for the residual variance.
#'                See [emmeans::eff_size] for details, by default NULL
#' @param ... additional arguments passed to [emmeans::eff_size]
#'
#' @return data.frame with effect sizes calculated by emmeans
#'
#'
#' @export
mmrm_eff_size <- function(mmrm_model, mmrm_emm, residual_error = TRUE, ...) {
  UseMethod("mmrm_eff_size")
}

#' @exportS3Method
#' @importFrom foreach %do%
mmrm_eff_size.default <- function(mmrm_model, mmrm_emm, edf = NULL, ...) {
  if (!"contrasts" %in% names(mmrm_emm)) {
    stop("MMRM estimated marginal means object must include contrasts!")
  }

  res_var <- diag(varcov(mmrm_model))
  dfs <- if (is.null(edf)) as.data.frame(mmrm_emm$contrasts)[, "df"] else rep(edf, length(res_var))

  foreach::foreach(i = 1:nrow(mmrm_emm$contrasts@grid), .combine = rbind) %do% {
    this_emm <- mmrm_emm
    this_emm$contrasts <- mmrm_emm$contrasts[i]
    this_eff <- suppressMessages(emmeans::eff_size(this_emm,
      sigma = sqrt(res_var[i]),
      edf = dfs[i]
    ))
    this_eff@levels <- this_emm$contrasts@levels
    this_eff@grid <- this_emm$contrasts@grid
    as.data.frame(this_eff)
  }
}

#' @importFrom foreach %dopar%
#' @exportS3Method
mmrm_eff_size.mmrmList <- function(mmrm_model, mmrm_emm, residual_error = TRUE, ...) {
  .check_foreach_backend()
  eff_list <- foreach::foreach(i = 1:length(mmrm_model)) %dopar% {
    mmrm_eff_size(mmrm_model[[i]], mmrm_emm[[i]], ...)
  }
  names(eff_list) <- names(mmrm_model)
  eff_list
}

#' @importFrom foreach %dopar%
#' @exportS3Method
mmrm_eff_size.mmrmCV <- function(mmrm_model, mmrm_emm, residual_error = TRUE, ...) {
  .check_foreach_backend()
  eff_list <- foreach::foreach(i = 1:length(mmrm_model)) %dopar% {
    mmrm_eff_size(mmrm_model[[i]], mmrm_emm[[i]], ...)
  }
  names(eff_list) <- names(mmrm_model)
  eff_list
}

#' @importFrom foreach %dopar%
mmrm_emmeans.mmrmCV <- function(obj, ...) {
  .check_foreach_backend()
  em_list <- foreach::foreach(i = obj) %dopar% mmrm_emmeans(i, ...)
  names(em_list) <- names(obj)
  em_list
}