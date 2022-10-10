#' Simulate MMRM data
#'
#' @description
#' Simulate a "clinical" dataset according to assumptions of the MMRM.
#' Datasets should have "observations" at a selected number of discrete timepoints.
#' Observations are predicted by:
#' * the mean observation at each timepoint
#' * the effect of baseline at each timepoint
#' * the effect of arm at each timepoint
#' * residual covariance matrix
#'
#' @param n_per_arm the number of subjects per arm
#' @param n_timepoints the number of timepoints.
#'                     If effects are specified as vectors, this is ignored.
#' @param effect_tp mean observation at each timepoint.
#'                  Either a vector with length equal to the
#'                  number of timepoints or a scalar.
#' @param effect_baseline effect of baseline value at each timepoint.
#'                        Either a vector with length equal to the
#'                        number of timepoints or a scalar.
#' @param effect_arm effect of arm at each timepoint.
#'                   Either a vector with length equal to the
#'                   number of timepoints or a scalar.
#' @param covariance full residual covariance matrix.
#'                   Either covariance of both of variance and
#'                   correlation should be specified.
#' @param variance variance at each timepoint.
#'                   Either a vector with length equal to the
#'                   number of timepoints or a scalar.
#'                   Either both of variance and correlation or
#'                   covariance should be specified
#' @param correlation correlation of residuals across timepoints.
#'                    Should be a scalar or vector with length:
#'                    timepoints * (timepoints - 1) / 2.
#'                    Either both of variance and correlation or
#'                    covariance should be specified.
#' @param baseline mean and variance of baseline scores
#' @param p_missing probability of missing data at each timepoint.
#'                  Either a vector with length equal to the
#'                  number of timepoints or a scalar.
#'
#' @returns data.frame with simulated data in the long format
#'          with the following fields:
#' * n_subs: the nubmer of subjects
#' * subject: subject id
#' * arm: the study arm or group
#' * base: baseline score
#' * time: time point of observation
#' * chg: the change score at the given time point (i.e., the outcome)
#'
#' @export
mmrm_simulate <- function(n_per_arm=50,
                          n_timepoints=5,
                          effect_tp=8,
                          effect_baseline=-0.5,
                          effect_arm=2,
                          covariance=NULL,
                          variance=10,
                          correlation=0,
                          baseline = c(25, 5),
                          p_missing=0.1) {

  len_tp = length(effect_tp)
  len_baseline = length(effect_baseline)
  len_arm = length(effect_arm)
  max_len = max(len_tp, len_baseline, len_arm)
  if (max_len > 1) {
    lens = c(len_tp, len_baseline, len_arm)
    if (any((lens != max_len) & (lens > 1))) {
      stop("effect_tp, effect_baseline, and effect_arm have ",
           "different lengths. They must be scalars or vectors with length equal ",
           "to the number of timepoints.")
    }
    n_timepoints = max_len
  }

  if (is.null(covariance)) {
    cor_len = (n_timepoints * (n_timepoints - 1)) / 2
    if (
      ((length(variance) > 1) & (length(variance) != n_timepoints)) |
      ((length(correlation) > 1) & (length(correlation) != cor_len))
    ) {
      stop("variance and correlation must be either scalars or length ",
           "n_timepoints and n_timepoints*(n_timepoints-1)/2, respectively")
    }

    if (length(variance) < n_timepoints) variance = rep(variance, n_timepoints)
    if (length(correlation) < cor_len) correlation = rep(correlation, cor_len)
    lower_chol = diag(1, n_timepoints)
    lower_chol[lower.tri(lower_chol)] = correlation
    lower_chol[upper.tri(lower_chol)] = t(lower_chol)[upper.tri(lower_chol)]
    lower_chol = diag(sqrt(variance)) %*% lower_chol
    covariance = lower_chol %*% t(lower_chol)
  }

  if (len_tp < n_timepoints) effect_tp = rep(effect_tp, n_timepoints)
  if (len_baseline < n_timepoints) effect_baseline = rep(effect_baseline, n_timepoints)
  if (len_arm < n_timepoints) effect_arm = rep(effect_arm, n_timepoints)

  base_vals = rnorm(n_per_arm*2, baseline[1], baseline[2])
  arm_vals = sample(c(rep(0, n_per_arm), rep(1, n_per_arm)))

  mus = matrix(rep(effect_tp, n_per_arm*2), nrow=n_per_arm*2)
  mus = mus + base_vals %*% t(effect_baseline)
  mus = mus + arm_vals %*% t(effect_arm)

  noise = MASS::mvrnorm(n_per_arm*2, rep(0, n_timepoints), covariance)
  outcomes = mus + noise

  if (length(p_missing) < n_timepoints) p_missing = rep(p_missing, n_timepoints)
  for (c in 1:ncol(outcomes)) {
    mask = sample(c(NA, 1), n_per_arm*2, replace=T, prob=c(p_missing[c], 1-p_missing[c]))
    outcomes[, c] = outcomes[, c] * mask
  }

  data.frame(n_subs = n_per_arm,
             subject = 1:(n_per_arm * 2),
             arm = arm_vals,
             base = base_vals,
             time = rep(1:n_timepoints, each=n_per_arm*2),
             chg = as.numeric(outcomes))
}

#'
#'
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
