#' Calculate effect size for contrasts using estimated marginal means
#'
#' @description
#' Wrapper for [emmeans::eff_size] to run for MMRM models.
#'
#' @param mmrm_model the mmrm model used to calculate effect size
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
  if (!"contrasts" %in% names(mmrm_emm)) {
    stop("MMRM estimated marginal means object must include contrasts!")
  }

  df <- mean(as.data.frame(mmrm_emm$contrasts)[, "df"])
  emmeans::eff_size(mmrm_emm,
    sigma = sigma(mmrm_model),
    edf = df,
    ...
  )
}
