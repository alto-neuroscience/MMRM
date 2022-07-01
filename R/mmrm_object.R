#' Fitted mmrm Object
#'
#' @description
#' An object returned by the [mmrm] function,
#' representing a fitted Mixed Model Repeated Measures model.
#' mmrm objects inherit from class "mmrm" and "gls", abd are a simple
#' extension of [nlme::glsObject].
#' Please see [nlme::glsObject] for full details.
#'
#' @return
#' "mmrm" objects have all components of \link[nlme]{glsObject},
#' plus the following:
#' \item{data}{the data structure used to fit the model}
#' \item{warnings}{warning messages from the model fitting procedure}
#'
#' @name mmrmObject
#'
#' @seealso [nlme::glsObject] [mmrm()]
NULL
