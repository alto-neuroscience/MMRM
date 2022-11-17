library(MMRM)
library(testthat)

###

test_mmrm <- function(d) {
  expect_no_error(
    m <- mmrm(chg ~ base + arm + time + base:time + arm:time,
              time = "time",
              subjects = "subject",
              categorical_time=FALSE,
              heterogenous=FALSE,
              data = d
    )
  )
  expect_no_error(
    m_emm <- mmrm_emmeans(m, group="arm", mode = "kenward", information = "expected")
  )
  expect_no_error(
    m_eff <- mmrm_eff_size(m, m_emm)
  )
}



test_that(
  "Test continuous time/homogenous variance model doesn't throw error",
  {
    set.seed(1)
    for (i in 1:10) {
      n_subs <- sample(10:50, 1)
      n_timepoints <- sample(2:6, 1)
      p_missing <- runif(n_timepoints) / 2
      sim_data <- mmrm_simulate(n_subs, n_timepoints, p_missing = p_missing)
      test_mmrm(sim_data)
    }
  }
)

