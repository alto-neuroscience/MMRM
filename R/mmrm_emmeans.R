#' Estimated Marignal Means of MMRM models
#'
#' Alias for [emmeans::emmeans] that provides support
#' for using Kenward-Rogers degrees of freedom calculation.
#' Please see [emmeans::emmeans] for further documentation.
#'
#' @details Support for Kenward-Rogers degrees of freedom calculation
#'           is made possible with i) a custom emm_basis function and
#'           ii) addition of support for glsObject in the pbkrtest package
#'           (currently in a development version
#'           [gkane26/pbkrtest@nlme](https://github.com/gkane26/pbkrtest/tree/nlme)).
#'           \cr\cr
#'           Kenward-Rogers d.f. calculation is memory and computation time intensive.
#'           By default, running this method on a large dataset will throw an error to warn users.
#'           To bypass this error, please set the pbkrtest.limit argument to a number larger than
#'           the number of observations in the model.
#'
#' @returns see [emmeans::emmeans]
#'
#' @examples
#' \dontrun{
#' mmrm_emmeans(object = mmrm_object,
#'              specs = pairwise ~ time | treatment,
#'              mode = "kenward",
#'              pbkrtest.limit = emmeans::get_emm_option("pbkrtest.limit"))
#' }
#'
#' @export
mmrm_emmeans <- emmeans::emmeans
