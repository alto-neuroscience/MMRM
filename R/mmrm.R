#' Fits MMRM model using nlme::gls
#'
#' @description
#' Fits a Mixed Model Repeated Measures model (see Details).
#' The fixed effects structure is flexible -- user defined using a formula.
#' The random effects structure is hard-coded:
#' independent variance at each time point,
#' but residuals of individual subjects are correlated across timepoints.
#'
#' @param formula formula for the fixed effects structure of the model
#' @param time the time variable of the model
#' @param subjects the variable indicating unique subjects
#' @param data the data structure
#' @param cov_struct a sequence of covariance matrix specifications;
#'                       will fit a model starting with the first covariance
#'                       and check for convergence. If there are convergence
#'                       issues, it will iterate through the rest of the list
#'                       until the model converges.
#' @param control a glsControl object, see \link[nlme]{glsControl}
#' @param stop_on_convergence if TRUE, stop fitting models and
# "                              return results once one model converges.
#'                              if FALSE, fits models with all specified
#'                              covariance structures regardless of
#'                              convergence of earlier model
#' @param return_all logical; if TRUE, will return all models that were fit,
#'                    if FALSE, only returns the first model that converged
#'                    or last model attempted (if all failed).
#' @param ... additional arguments passed to \link[nlme]{gls}
#'             (e.g., na.action)
#'
#' @details TODO
#'
#' @returns glsObject or list of glsObjects
#'            if return_all = FALSE, returns only the glsObject
#'              of the first model to converge or the last attempted model.
#'            if return_all = TRUE, returns list of all attempted models
#'
#' @examples
#' \donttest{
#' mmrm(outcome ~ baseline + group + time + baseline:time + group:time,
#'      time = "time",
#'      subjects = "subjectid",
#'      data = data)
#'}
#'
#' @export
mmrm <- function(formula,
                 time,
                 subjects,
                 data,
                 cov_struct = c(
                   "unstructured",
                   "autoregressive",
                   "compound-symmetry"
                 ),
                 control = nlme::glsControl(),
                 stop_on_convergence = TRUE,
                 return_all = FALSE,
                 ...) {

  # check that time variable is a factor
  if (class(data$time) != "factor") {
    data$time = factor(data$time)
    warning("time variable is not a factor, coercing it to factor")
  }

  cov_list <- .get_cov_list(cov_struct)
  correlation_formula = stats::as.formula(paste0("~ as.numeric(", time, ") | ", subjects))
  weights_formula = stats::as.formula(paste0("~ 1 | ", time))

  res_list <- list()
  for (i in 1:length(cov_list)) {

    # It looks like information regarding model fit/convergence
    # is not returned with the glsObject.
    # Thus, the only way I can think of to evaluate convergence
    # is to catch warnings and errors when model is fit
    res <- tryCatch(
      eval(
        bquote(
          nlme::gls(
            .(formula),
            correlation = nlme::corSymm(form = .(correlation_formula)),
            weights = nlme::varIdent(form = .(weights_formula)),
            data = data,
            ...
          )
        )
      ),
        error = function(e) e,
        warning = function(w) w
    )
    res$data = data
    class(res) = c("mmrm", class(res))

    res_list[[names(cov_list)[i]]] <- res

    if (inherits(res, "error")) {
      warning("Error fitting MMRM with covariance structure = ", names(cov_list)[i], ": ", res$message)
    } else if (inherits(res, "warning")) {
      warning("Warning fitting MMRM with covariance structure = ", names(cov_list)[i], ": ", res$message)
    } else if (stop_on_convergence) {
      break
    }
      }

    if (return_all) {
      res_list
    } else {
      res
    }
  }

  .get_cov_list <- function(covariance = c(
    "unstructured",
    "autoregressive",
    "compound-symmetry"
  )) {
    cor_list <- list()
    for (c in covariance) {
      ctype <- match.arg(tolower(c), covariance)
      cor_list[[c]] <- .cov_map(ctype)
    }

    return(cor_list)
  }

  .cov_map <- function(cov_type) {
    switch(cov_type,
           "unstructured" = nlme::corSymm,
           "autoregressive" = nlme::corAR1,
           "compound-symmetry" = nlme::corCompSymm,
    )
  }
