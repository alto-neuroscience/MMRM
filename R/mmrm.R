#' Fits MMRM model using nlme::gls
#'
#' @description
#' Fits a Mixed Model Repeated Measures model (see Details).
#' In this implementation, the fixed effects structure is flexible --
#' it is user defined using a formula. This implementation does not support
#' the inculsion of random effects, and it specifies that the residual variance
#' at each time point is independent and residuals of individual subjects
#' are correlated across timepoints.
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
#' @param method character, either "REML" or "ML". If "REML" the model is fit
#'                   by maximizing the restricted log-likelihood. If "ML" the
#'                   log-likelihood is maximized. Defaults to "REML".
#' @param na.action a function that indicates what should happen when the
#'                      data contain NAs. Defaults to na.omit
#' @param control a list of control values fo the estimation algorithm to
#'                    to replace the default values returned by the function
#'                    [nlme::glsControl]
#' @param verbose logical value indicating whether to print the
#'                    evolution of the iterative algorithm. Default is FALSE
#' @param return_all logical; if TRUE, will return all models that were fit,
#'                    if FALSE, only returns the first model that converged
#'                    or last model attempted (if all failed).
#' @param stop_on_convergence ignored if return_all = FALSE.
#'                                if return_all = TRUE, stop fitting models and
#'                              return results once one model converges if
#'                              stop_on_convergence = TRUE. if return_all = TRUE
#'                              and stop_on_convergence = FALSE, fits and
#'                              returns all specified models
#'
#' @details
#' The MMRM implemented here is defined as:
#' \deqn{Y_i = X_i\beta + \epsilon_i}
#' \deqn{\epsilon_i ~ \mathcal{N}(0, \Sigma_i)}
#' * \eqn{Y_i} as the vector of outcomes with length
#' \eqn{n_{subjects}*n_{timepoints}}
#' * \eqn{X_i} as the a matrix of predictors with
#' \eqn{n_{subjects}*n_{timepoints}} rows and \eqn{n_{predictors}} columns
#' * \eqn{\beta} as the vector of coefficients with length \eqn{n_{predictors}}
#' * \eqn{\epsilon_i} as the vector of residuals
#' * \eqn{\Sigma_i} as the (\eqn{n_{subjects}*n_{timepoints} x n_{subjects}*(n_{timepoints}}) covariance matrix
#'
#' This implementation of the MMRM supports three different
#' covariance structures: unstructured, AR(1), and compound symmetry.
#' Timepoints are always treated as categorical,
#' with independent residual variance at each timepoint and correlated residuals
#' among individual subjects across timepoints.
#'
#' @returns
#' `mmrmObject` or `mmrmList`, a list of `mmrmObjects`
#'
#' If return_all = FALSE, returns only the mmrmObject of the first model to
#' converge or the last attempted model.\cr
#' If return_all = TRUE, returns list of all attempted models.
#'
#' @examples
#' \dontrun{
#' mmrm(
#'   outcome ~ baseline + group + time +
#'     baseline:time + group:time,
#'   time = "time",
#'   subjects = "subjectid",
#'   data = my_data
#' )
#' }
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
                 method = "REML",
                 na.action = na.exclude,
                 control = list(),
                 verbose = FALSE,
                 return_all = FALSE,
                 stop_on_convergence = TRUE,
                 ...) {

  # check arguments
  formula <- .check_mmrm_args(formula, time, subjects, data)

  if (!return_all && !stop_on_convergence) {
    stop_on_convergence <- TRUE
    warning(
      "return_all is FALSE, so setting stop_on_convergence = TRUE ",
      "to return the first model that converges."
    )
  }


  # ensure that time variable is a factor
  if (class(data[[time]]) != "factor") {
    data[[time]] <- factor(data[[time]])
    warning("time variable is not a factor, coercing it to factor")
  }

  cov_list <- .get_cov_list(cov_struct)
  correlation_formula <- stats::as.formula(paste0("~ as.numeric(", time, ") | ", subjects))
  weights_formula <- stats::as.formula(paste0("~ 1 | ", time))

  res_list <- list()
  for (i in 1:length(cov_list)) {

    # It looks like information regarding model fit/convergence
    # is not returned with the glsObject.
    # To evaluate convergence, catch warnings/errors during model fit
    res <- tryCatch(
      withCallingHandlers(
        {
          eval(
            bquote(
              nlme::gls(
                .(formula),
                data = data,
                correlation = .(cov_list[[i]])(form = .(correlation_formula)),
                weights = nlme::varIdent(form = .(weights_formula)),
                method = .(method),
                na.action = .(na.action),
                control = .(control),
                verbose = .(verbose)
              )
            )
          )
        },
        warning = function(w) {
          wenv <<- new.env(parent = parent.env(environment()))
          wenv$warning_log <- w
        }
      ),
      error = function(e) e
    )
    if (inherits(res, "error")) {
      warning("Error fitting MMRM with covariance structure = ", names(cov_list)[i], ": ", res$message)
    } else {
      res$call$correlation[1] <- .get_cov_call(names(cov_list)[i])
      res$data <- na.action(data)
      res$time <- time
      res$subjects <- subjects
      # res$na.action = na.action
      if (exists("wenv", mode = "environment")) {
        if (exists("warning_log", envir = wenv)) {
          res$warnings <- warning_log
        }
        rm(wenv, envir = .GlobalEnv)
      }
      if (!("warnings" %in% names(res))) {
        res$warnings <- NA
      }
      class(res) <- c("mmrm", class(res))

      res_list[[names(cov_list)[i]]] <- res
      class(res_list) <- c("mmrmList", class(res))

      if (is.na(res$warnings) && stop_on_convergence) break
    }
  }

  if (inherits(res, "error")) {
    stop(res)
  } else if (return_all) {
    res_list
  } else {
    res
  }
}

.check_mmrm_args <- function(formula, time, subjects, data) {
  if (missing(formula) ||
    missing(time) ||
    missing(subjects) ||
    missing(data)) {
    stop(
      "missing argument! ",
      "The following arguments must be specified: ",
      "formula, time, subjects, data"
    )
  }
  if (class(formula) == "character") formula <- as.formula(formula)
  if (class(formula) != "formula") {
    stop(
      "formula argument must be a formula ",
      "or a string that represents a formula"
    )
  }
  if (class(time) != "character" || class(subjects) != "character") {
    stop(
      "the time and subjects arguments must be strings -- ",
      "the column names of the time and subjects variables ",
      "in the data structure"
    )
  }
  if (!(time %in% names(data))) {
    stop(
      "the time variable (", time, ") ",
      "was not found in data"
    )
  }
  if (!(subjects %in% names(data))) {
    stop(
      "the subjects variable (", subjects, ") ",
      "was not found in data"
    )
  }
  if (!inherits(data, "data.frame")) {
    stop(
      "data is not a valid data structure ",
      "(e.g., a data.frame, data.table, or tibble)"
    )
  }

  return(formula)
}

COV_TYPES <- c(
  "unstructured",
  "autoregressive",
  "compound-symmetry"
)
.get_cov_list <- function(cov_struct = COV_TYPES) {
  cov_list <- list()
  for (c in cov_struct) {
    ctype <- match.arg(tolower(c), COV_TYPES)
    cov_list[[ctype]] <- .cov_map(ctype)
  }

  return(cov_list)
}

.cov_map <- function(cov_type) {
  switch(cov_type,
    "unstructured" = nlme::corSymm,
    "autoregressive" = nlme::corAR1,
    "compound-symmetry" = nlme::corCompSymm,
  )
}

.get_cov_call <- function(cov_struct) {
  switch(cov_struct,
    "unstructured" = quote(nlme::corSymm()),
    "autoregressive" = quote(nlme::corAR1()),
    "compound-symmetry" = quote(nlme::corCompSymm())
  )
}
