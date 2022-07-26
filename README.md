
# MMRM

R package to fit Mixed Model for Repeated Measures as is commonly used
to analyze clinical trial data. This package uses `nlme::gls` to fit the model,
and provides support for Kenward-Rogers degrees of freedom calculation.

This package is currently in beta version -- more testing and examples to come!

## Installation

##### Note: This package requires a [fork of the `pbkrtest` package](https://github.com/gkane26/pbkrtest/tree/nlme)

``` r
if(!requireNamespace("remotes")) install.packages("remotes")
try(remotes::install_github("gkane26/pbkrtest@nlme"))
remotes::install_github("alto-neuroscience/MMRM")
```

## Example

``` r
library(MMRM)


# fit an MMRM
my_mmrm = MMRM::mmrm(outcome ~ baseline + group + time + baseline:time + group:time,
                     time = "time",
                     subjects = "subjects",
                     data = my_data)
                     
# fit MMRM using k-fold cross validation (steps below work the same)
my_mmrm = MMRM::mmrm_cv(outcome ~ baseline + group + time + baseline:time + group:time,
                     time = "time",
                     subjects = "subjects",
                     data = my_data,
                     k = 10)

# get estimated marginal means
mmrm_emm = MMRM::mmrm_emmeans(my_mmrm,
                              pairwise ~ group | time,
                              mode = "kenward")

# calculate effect size
mmrm_eff = MMRM::mmrm_eff_size(my_mmrm, mmrm_emm)
```

