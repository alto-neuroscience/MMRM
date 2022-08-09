#' Estimated Marignal Means of MMRM models
#'
#' Generic method to run [emmeans::emmeans] on all [MMRM] package objects using
#' Kenward-Rogers degrees of freedom calculation.
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
#' mmrm_emmeans(
#'   object = mmrm_object,
#'   specs = pairwise ~ time | treatment,
#'   mode = "kenward",
#'   pbkrtest.limit = emmeans::get_emm_option("pbkrtest.limit")
#' )
#' }
#'
#' @export
mmrm_emmeans <- function(obj, ...) {
  UseMethod("mmrm_emmeans")
}

#' @exportS3Method
mmrm_emmeans.default <- emmeans::emmeans

#' @importFrom foreach %dopar%
#' @exportS3Method
mmrm_emmeans.mmrmList <- function(obj, ...) {
  .check_foreach_backend()
  em_list <- foreach::foreach(i = obj) %dopar% mmrm_emmeans(obj, ...)
  names(em_list) <- names(obj)
  em_list
}

#' @importFrom foreach %dopar%
#' @exportS3Method
mmrm_emmeans.mmrmCV <- function(obj, ...) {
  .check_foreach_backend()
  em_list <- foreach::foreach(i = obj) %dopar% mmrm_emmeans(i, ...)
  names(em_list) <- names(obj)
  em_list
}
