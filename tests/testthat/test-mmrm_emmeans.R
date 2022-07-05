library(MMRM)

if (!requireNamespace("lme4")) {
  stop("lme4 package not installed and is needed to run tests!")
} else {
  library(lme4)
}

###

test_data <- as.data.frame(nlme::Dialyzer)
names(test_data) <- c("subject", "group", "baseline", "outcome", "time")
test_data$time <- factor(test_data$time)

###

test_that(
  "Fixed effects of MMRM results against equivalent lmer model",
  {
    mmrm.1 <- mmrm(
      outcome ~ baseline + group + time + baseline:time + group:time,
      time = "time",
      subjects = "subject",
      data = test_data
    )
    mmrm.emm <- mmrm_emmeans(mmrm.1, pairwise ~ time | group)
    mmrm.eff <- mmrm_eff_size(mmrm.1, mmrm.emm)

    gls.1 <- nlme::gls(
      outcome ~ baseline + group + time + baseline:time + group:time,
      correlation = nlme::corSymm(form = ~ as.numeric(time) | subject),
      weights = nlme::varIdent(form = ~ 1 | time),
      data = test_data
    )

    lmer.1 <- suppressWarnings(
      lmer(
        outcome ~ baseline + group + time + baseline:time + group:time + (0 + time | subject),
        data = test_data,
        control = lmerControl(check.nobs.vs.nRE = "ignore")
      )
    )

    testthat::expect_equal(mmrm.1$coefficients, fixef(lmer.1), tolerance = 1e-3)
    testthat::expect_equal(mmrm.1$coefficients, gls.1$coefficients, tolerance = 1e-5)
    testthat::expect_equal(varcov(mmrm.1), varcov(gls.1), tolerance = 1e-5)
  }
)
