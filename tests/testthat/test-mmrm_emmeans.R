library(MMRM)
library(testthat)

if (!requireNamespace("lme4")) {
  stop("lme4 package not installed and is needed to run tests!")
} else {
  library(lme4)
}

###

compare_mmrm_lmer <- function(d) {
  d$time <- factor(d$time)

  m <- mmrm(chg ~ base + arm + time + base:time + arm:time,
    time = "time",
    subjects = "subject",
    data = d
  )

  m_emm <- mmrm_emmeans(m, group="arm", mode = "kenward", information = "expected")
  m_eff <- mmrm_eff_size(m, m_emm)

  l <- suppressMessages(suppressWarnings(
    lmer(chg ~ base + arm + time + base:time + arm:time + (0 + time | subject),
      data = d,
      control = lmerControl(check.nobs.vs.nRE = "ignore")
    )
  ))
  l_emm <- emmeans::emmeans(l, pairwise ~ arm | time)

  g <- try(
    nlme::gls(chg ~ base + arm + time + base:time + arm:time,
      correlation = nlme::corSymm(form = ~ as.numeric(time) | subject),
      weights = nlme::varIdent(form = ~ 1 | time),
      data = d,
      na.action = na.exclude
    ),
    silent = TRUE
  )
  if (!inherits(g, "try-error")) {
    expect_equal(m$coefficients, g$coefficients, tolerance = 1e-5)
  }
  tol <- if (isSingular(l)) 1e-1 else 1e-3
  if (("corSymm" %in% attr(m$modelStruct$corStruct, "class")) & (is.matrix(m$apVar))) {
    expect_equal(data.frame(m_emm$contrasts), data.frame(l_emm$contrasts), tolerance = tol)
  }
}

test_that(
  "Compare MMRM against equivalent lmer model",
  {
    set.seed(42)
    for (i in 1:100) {
      n_subs <- sample(10:50, 1)
      n_timepoints <- sample(2:6, 1)
      p_missing <- runif(n_timepoints) / 2
      sim_data <- mmrm_simulate(n_subs, n_timepoints, p_missing = p_missing)
      compare_mmrm_lmer(sim_data)
    }
  }
)
