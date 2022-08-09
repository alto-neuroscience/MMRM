#' Calculate effect size for contrasts using estimated marginal means
#'
#' @description
#' Wrapper for [emmeans::eff_size] to run for MMRM models.
#'
#' @param mmrm_model the mmrm model(s) of type `mmrm`, `mmrmList`, or `mmrmCV`
#'                     used to calculate effect size
#' @param mmrm_emm estimated marginal means from and MMRM.
#'                     output of [mmrm_emmeans]
#' @param ... additional arguments passed to [emmeans::eff_size]
#'
#' @details
#' The `sigma` argument (the population standard deviation) for
#' [emmeans::eff_size] is calculated using [stats::sigma]. The `edf` argument
#' (the equivalent degrees of freedom) is estimated as the mean of the
#' degrees of freedom  across all contrasts calculated in [mmrm_emmeans]
#' (typically using Kenward-Rogers approximation).
#'
#' @export
mmrm_eff_size <- function(mmrm_model, mmrm_emm, ...) {
  UseMethod("mmrm_eff_size")
}

#' @exportS3Method
mmrm_eff_size.default <- function(mmrm_model, mmrm_emm, ...) {
  if (!"contrasts" %in% names(mmrm_emm)) {
    stop("MMRM estimated marginal means object must include contrasts!")
  }

  df <- mean(as.data.frame(mmrm_emm$contrasts)[, "df"])
  suppressMessages(emmeans::eff_size(mmrm_emm,
    sigma = sigma(mmrm_model),
    edf = df,
    ...
  ))
}

#' @importFrom foreach %dopar%
#' @exportS3Method
mmrm_eff_size.mmrmList <- function(mmrm_model, mmrm_emm, ...) {
  .check_foreach_backend()
  eff_list <- foreach::foreach(i = 1:length(mmrm_model)) %dopar% {
    mmrm_eff_size(mmrm_model[[i]], mmrm_emm[[i]], ...)
  }
  names(eff_list) <- names(mmrm_model)
  eff_list
}

#' @importFrom foreach %dopar%
#' @exportS3Method
mmrm_eff_size.mmrmCV <- function(mmrm_model, mmrm_emm, ...) {
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
