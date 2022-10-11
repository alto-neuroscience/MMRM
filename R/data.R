#' Simulated MMRM data sets
#'
#' @description
#' A handful of datasets simulated from an MMRM mechanism.
#' These data sets are used for unit testing and for comparing results against
#' other software, such as Stata and SAS, as well as models from other packages.
#' All datasets have two arms and 4 timepoints, with some outcome values
#' missing at random.
#'
#' * sim_mmrm_missing: A dataset with 20 subjects per arm.
#'                         All subjects have at least one missing value.
#' * sim_mmrm_20: 20 subjects per arm with values missing at random.
#' * sim_mmrm_40: 40 subjects per arm with values missing at random.
#' * sim_mmrm_60: 60 subjects per arm with values missing at random.
#' * sim_mmrm_80: 80 subjects per arm with values missing at random.
#' * sim_mmrm_100: 100 subjects per arm with values missing at random.
#'
#' @format ## `sim_mmrm`
#' \describe{
#'  \item{n_subs}{number of subjects per arm}
#'  \item{subject}{subject ID}
#'  \item{arm}{arm code (0 or 1)}
#'  \item{base}{baseline score}
#'  \item{time}{time point of measurement (1-4)}
#'  \item{mus}{mean score at timepoint}
#'  \item{chg}{change from baseline at timepoint (chg = mus + residual error)}
#' }
#' @source internal simulation
#' @name sim_mmrm
NULL

#' @rdname sim_mmrm
"sim_mmrm_missing"

#' @rdname sim_mmrm
"sim_mmrm_20"

#' @rdname sim_mmrm
"sim_mmrm_40"

#' @rdname sim_mmrm
"sim_mmrm_60"

#' @rdname sim_mmrm
"sim_mmrm_80"

#' @rdname sim_mmrm
"sim_mmrm_100"
