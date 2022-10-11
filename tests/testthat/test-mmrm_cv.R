library(MMRM)

###

test_data <- sim_mmrm_100
test_data$time <- factor(test_data$time)
test_data <- as.data.frame(test_data)

###

test_that(
  "Sequential Cross Validation procedure does not throw an error",
  {
    foreach::registerDoSEQ()
    expect_error(
      {
        models <- mmrm_cv(
          chg ~ base + arm + time + base:time + arm:time,
          time = "time",
          subjects = "subject",
          data = test_data,
          k = 3
        )
        emms <- mmrm_emmeans(models, specs = pairwise ~ arm | time)
        effs <- mmrm_eff_size(models, emms)
      },
      NA
    )
  }
)
