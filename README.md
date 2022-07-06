
# MMRM

R package to fit Mixed Model for Repeated Measures as is commonly used
to analyze clinical trial data. This package uses `nlme::gls` to fit the model,
and provides support for Kenward-Rogers degrees of freedom calculation.

TODO: more testing and examples!

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

my_mmrm = MMRM::mmrm(outcome ~ baseline + group + time + baseline:time + group:time,
                     time = "time",
                     subjects = "subjects",
                     data = my_data)

mmrm_emm = MMRM::mmrm_emmeans(my_mmrm,
                              pairwise ~ group | time,
                              mode = "kenward")

mmrm_eff = MMRM::mmrm_eff_size(my_mmrm, mmrm_emm)
```

