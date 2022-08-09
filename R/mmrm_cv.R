#' K-fold cross validation of MMRM models
#'
#' @description
#' Performs k-fold cross validation with the Mixed Model Repeated Measures
#' model (see [MMRM]).
#'
#' @param formula formula for the fixed effects structure of the model
#' @param time the time variable of the model
#' @param subjects the variable indicating unique subjects
#' @param data the data structure
#' @param k the number of splits for k-fold cross-validation, default = 10
#' @param seed random seed used to generating splits,
#'              If NULL does not set the seed
#' @param in_loop function to process each split of training set data. See details
#' @param ... arguments passed to \code{in_loop} and [MMRM]
#'
#' @details
#' Please see [MMRM] for details on the model.
#'
#' [mmrm_cv] supports fitting each of the k-folds in parallel using
#' [foreach::foreach] loops. To use multiple cores, please register a parallel
#' backend prior to calling [mmrm_cv]. Here is an example:
#' ```
#' cl = parallel::makeCluster(n_cores)
#' doParallel::registerDoParallel(cl)
#' ```
#'
#' The \code{in_loop} argument allows users to provide a function that
#' transforms training set data. This function must:
#' - accept the training set data as its first argument
#' - accept pass through arguments (...)
#' - return the processed data structure with all variables
#' specified in the supplied formula
#' Here is an example \code{in_loop} that performs no transformations:
#' ```
#' in_loop = function(data, ...) {
#'   # do stuff here
#'   return(data)
#' }
#' ```
#'
#'
#' @returns
#' `mmrmCV` object: a list of outputs with length `k`,
#' each element is the output of a call to [mmrm].
#' See [mmrm] and [mmrmObject] for details
#'
#' @examples
#' \dontrun{
#' mmrm_cv(
#'   outcome ~ baseline + group + time +
#'     baseline:time + group:time,
#'   time = "time",
#'   subjects = "subjectid",
#'   data = my_data,
#'   k = 10
#' )
#' }
#'
#' @importFrom foreach %dopar%
#'
#' @export
mmrm_cv <- function(formula,
                    time,
                    subjects,
                    data,
                    k = 10,
                    seed = NULL,
                    in_loop = NULL,
                    ...) {

  # check arguments
  formula <- .check_mmrm_args(formula, time, subjects, data)

  # get train and test split
  data <- .mmrm_train_test_split(data, subjects, k, seed)

  # loop through groups and fit models
  .check_foreach_backend()

  mmrm_list <- foreach::foreach(i = 1:k) %dopar% {
    split_data <- data[data$k != i, ]
    if (!is.null(in_loop)) split_data <- in_loop(split_data, ...)

    model <- mmrm(
      formula,
      time,
      subjects,
      split_data,
      ...
    )
    model$test <- data[data$k == i, ]
    model
  }

  names(mmrm_list) <- 1:k
  class(mmrm_list) <- c("mmrmCV")

  return(mmrm_list)
}

#' @importFrom magrittr %>%
#' @importFrom foreach %do%
.mmrm_train_test_split <- function(data, subjects, k = 10, seed = NULL) {
  if (!is.null(seed)) set_seed(seed)
  subs <- unique(data[[subjects]])
  n_subs <- length(subs)
  sub_group <- sample(0:(n_subs - 1) / (n_subs - 1)) %>%
    cut(seq(0, 1, 1 / k), labels = FALSE, include.lowest = TRUE)

  dsplit <- split(data, f = data[[subjects]])
  foreach::foreach(s = dsplit, .combine = rbind) %do% {
    s[["k"]] <- sub_group[which(subs == s[[subjects]][1])]
    s
  }
}

.check_foreach_backend <- function() {
  if (!foreach::getDoParRegistered()) foreach::registerDoSEQ()
}
