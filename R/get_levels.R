#' @export
get_levels <- function(object) {
  UseMethod("get_levels")
}

#' @export
get_levels.mmrm <- function(object) {
  length(unique(object$data[[object$time]]))
}

#' @export
get_levels.glsStruct <- function(object) {
  object <- object$corStruct
  get_levels(object)
}

#' @export
get_levels.corStruct <- function(object) {
  length(unique(unlist(attr(object, "covariate"))))
}
