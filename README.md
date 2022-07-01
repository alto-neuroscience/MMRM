
# MMRM

R package to fit Mixed Model for Repeated Measures as is commonly used
to analyze clinical trial data. This package uses `nlme::gls` to fit the model,
and provides support for Kenward-Rogers degrees of freedom calculation.

TODO: more testing, documentation, and examples!

## Installation

``` r
if(!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("alto-neuroscience/MMRM")
```

## Example

``` r
library(MMRM)

my_mmrm = MMRM::mmrm(outcome ~ baseline + group + time + baseline:time + group:time,
                     time = "time",
                     subjects = "subjects",
                     data = my_data)

mmrm_emmeans = MMRM::mmrm_emmeans(my_mmrm,
                                  pairwise ~ time | group,
                                  mode = "kenward")
```

