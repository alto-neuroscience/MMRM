library(MMRM)

###

test_data <- as.data.frame(nlme::Dialyzer)
names(test_data) <- c("subject", "group", "baseline", "outcome", "time")
test_data$time <- factor(test_data$time)

###

test_that(
  "Sequential Cross Validation procedure does not throw an error",
  {
    foreach::registerDoSEQ()
    expect_error(
      mcv <- mmrm_cv(
        outcome ~ baseline + group + time + baseline:time + group:time,
        time = "time",
        subjects = "subject",
        data = test_data,
        k = 3
      ),
      NA
    )
  }
)

test_that(
  "Parallel Cross Validation procedure does not throw an error",
  {
    cl <- parallel::makeCluster(3, outfile = "")
    doParallel::registerDoParallel(cl)
    expect_error(
      mcv <- mmrm_cv(
        outcome ~ baseline + group + time + baseline:time + group:time,
        time = "time",
        subjects = "subject",
        data = test_data,
        k = 3
      ),
      NA
    )
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
  }
)
