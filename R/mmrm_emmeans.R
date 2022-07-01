#' Estimated Marignal Means of MMRM models
#'
#' Alias for \link[emmeans]{emmeans} that provides support
#' for using Kenward-Rogers degrees of freedom calculation.
#' Please see \link[emmeans]{emmeans} for further documentation.
#'
#' @examples
#'
#' mmrm_emmeans(object = mmrm_object,
#'              specs = pairwise ~ time | treatment,
#'              mode = "kenward",
#'              pbkrtest.limit = emmeans::get_emm_option("pbkrtest.limit"))
#'
#' @param ... see /link[emmeans]{emmeans}
#'
#' @details Support for Kenward-Rogers degrees of freedom calculation
#'           provided with i) a custom emm_basis function that overwrites the
#'           emm_basis.glm provided by the emmeans package and
#'           ii) addition of support for glsObject in the pbkrtest package
#'           (currently in a development version
#'           \link[gkane26/pbkrtest@nlme]{https://github.com/gkane26/pbkrtest/tree/nlme}).
#'
#'           Kenward-Rogers d.f. calculation is memory and computation time intensive.
#'           By default, running this method on a large dataset will throw an error to warn users.
#'           To bypass this error, please set the pbkrtest.limit argument to a number larger than
#'           the number of observations in the model.
#'
#' @returns see \link[emmeans]{emmeans}
#'
#' @export
mmrm_emmeans <- function(...) {
  emmeans::emmeans(...)
}
