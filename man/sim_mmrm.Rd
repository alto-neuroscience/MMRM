% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{sim_mmrm}
\alias{sim_mmrm}
\alias{sim_mmrm_missing}
\alias{sim_mmrm_20}
\alias{sim_mmrm_40}
\alias{sim_mmrm_60}
\alias{sim_mmrm_80}
\alias{sim_mmrm_100}
\title{Simulated MMRM data sets}
\format{
\subsection{\code{sim_mmrm}}{

\describe{
\item{n_subs}{number of subjects per arm}
\item{subject}{subject ID}
\item{arm}{arm code (0 or 1)}
\item{base}{baseline score}
\item{time}{time point of measurement (1-4)}
\item{mus}{mean score at timepoint}
\item{chg}{change from baseline at timepoint (chg = mus + residual error)}
}
}

An object of class \code{spec_tbl_df} (inherits from \code{tbl_df}, \code{tbl}, \code{data.frame}) with 160 rows and 7 columns.

An object of class \code{spec_tbl_df} (inherits from \code{tbl_df}, \code{tbl}, \code{data.frame}) with 160 rows and 7 columns.

An object of class \code{spec_tbl_df} (inherits from \code{tbl_df}, \code{tbl}, \code{data.frame}) with 320 rows and 7 columns.

An object of class \code{spec_tbl_df} (inherits from \code{tbl_df}, \code{tbl}, \code{data.frame}) with 480 rows and 7 columns.

An object of class \code{spec_tbl_df} (inherits from \code{tbl_df}, \code{tbl}, \code{data.frame}) with 640 rows and 7 columns.

An object of class \code{spec_tbl_df} (inherits from \code{tbl_df}, \code{tbl}, \code{data.frame}) with 800 rows and 7 columns.
}
\source{
internal simulation
}
\usage{
sim_mmrm_missing

sim_mmrm_20

sim_mmrm_40

sim_mmrm_60

sim_mmrm_80

sim_mmrm_100
}
\description{
A handful of datasets simulated from an MMRM mechanism.
These data sets are used for unit testing and for comparing results against
other software, such as Stata and SAS, as well as models from other packages.
All datasets have two arms and 4 timepoints, with some outcome values
missing at random.
\itemize{
\item sim_mmrm_missing: A dataset with 20 subjects per arm.
All subjects have at least one missing value.
\item sim_mmrm_20: 20 subjects per arm with values missing at random.
\item sim_mmrm_40: 40 subjects per arm with values missing at random.
\item sim_mmrm_60: 60 subjects per arm with values missing at random.
\item sim_mmrm_80: 80 subjects per arm with values missing at random.
\item sim_mmrm_100: 100 subjects per arm with values missing at random.
}
}
\keyword{datasets}
