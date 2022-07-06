
# MMRM

R package to fit Mixed Model for Repeated Measures as is commonly used
to analyze clinical trial data. This package uses `nlme::gls` to fit the model,
and provides support for Kenward-Rogers degrees of freedom calculation.

TODO: more testing and examples!

## Installation

The easiest way to install the package from a private repo is to clone and install the package locally:
```bash
mkdir MMRM_tmp && cd MMRM_tmp
git clone git@github.com:alto-neuroscience/MMRM.git
R CMD BUILD MMRM
R CMD INSTALL MMRM_*.tar.gz

# delete temporary files created by MMRM installation
cd ..
rm -rf MMRM_tmp
```
[Typical installation from public repo]
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

mmrm_emm = MMRM::mmrm_emmeans(my_mmrm,
                              pairwise ~ group | time,
                              mode = "kenward")

mmrm_eff = MMRM::mmrm_eff_size(my_mmrm, mmrm_emm)
```

