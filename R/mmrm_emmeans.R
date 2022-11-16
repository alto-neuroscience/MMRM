#' Estimated Marignal Means of MMRM models
#'
#' Generic method to run [emmeans::emmeans] on all [MMRM] package objects using
#' Kenward-Rogers degrees of freedom calculation.
#' Please see [emmeans::emmeans] for further documentation.
#'
#' @param obj a fitted MMRM model
#' @param specs A character vector specifying the names of the predictors over
#'                which EMMs are desired. See [emmeans::emmeans] for more details.
#'                If omitted, group must be specified, and emmeans and
#'                contrasts will be calculated across groups at each time point.
#' @param group the group over which to calculate emmeans and contrasts.
#'                If specs is not specified, group must be specified.
#' @param cov.reduce logical indicating whether to average over a continuous
#'                     variable (usually time points) when calculating
#'                     emmeans and contrasts. By default, FALSE
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
mmrm_emmeans <- function(obj,
                         specs=NULL,
                         group=NULL,
                         cov.reduce=FALSE,
                         nesting=NULL,
                         ...) {
  UseMethod("mmrm_emmeans")
}

#' @exportS3Method
mmrm_emmeans.default <- function(obj,
                                 specs=NULL,
                                 group=NULL,
                                 cov.reduce=FALSE,
                                 ...) {
  if (is.null(specs)) {
    specs = "pairwise ~ "
    if (!is.null(group)) {
      specs = paste0(specs, group, " | ", obj$time)
    } else {
      stop("must either specify emmeans spec or ",
           "provide the group over which to calculate contrasts")
    }
  }

  emmeans::emmeans(obj,
                   as.formula(specs),
                   cov.reduce=cov.reduce,
                   nesting=NULL,
                   ...)
}

#' @importFrom foreach %dopar%
#' @exportS3Method
mmrm_emmeans.mmrmList <- function(obj,
                                  specs=NULL,
                                  group=NULL,
                                  cov.reduce=FALSE,
                                  ...) {
  .check_foreach_backend()
  em_list <- foreach::foreach(i = obj) %dopar% mmrm_emmeans(i,
                                                            specs=specs,
                                                            group=group,
                                                            cov.reduce=cov.reduce,
                                                            ...)
  names(em_list) <- names(obj)
  em_list
}

#' @importFrom foreach %dopar%
#' @exportS3Method
mmrm_emmeans.mmrmCV <- function(obj,
                                specs=NULL,
                                group=NULL,
                                cov.reduce=FALSE,
                                ...) {
  .check_foreach_backend()
  em_list <- foreach::foreach(i = obj) %dopar% mmrm_emmeans(i,
                                                            specs=specs,
                                                            group=group,
                                                            cov.reduce=cov.reduce,
                                                            ...)
  names(em_list) <- names(obj)
  em_list
}
