library(nlme)

J <- 20
beta <- 0.3
tau <- 0.1
omega <- 0.05

study_data <- data.frame(
  study = LETTERS[1:J],
  nj = rpois(J, 5),
  N = 10 + rpois(J, 20)
)

es_data <- Map(function(study, nj, N) {
  d <- beta + rnorm(1, sd = tau) + rnorm(nj, sd = omega) + rnorm(nj, sd = 1 / sqrt(N))
  data.frame(study = study, d = d, v = 1 / N)
}, study = study_data$study, nj = study_data$nj, N = study_data$N)

es_data <- do.call(rbind, es_data)
es_data$es <- 1:nrow(es_data)


test_that("Fisher_info() works for a multi-level meta-analysis model with fixed sigma.", {

  mod <- lme(d ~ 1, data = es_data,
             random = ~ 1 | study,
             weights = varFixed( ~ v),
             control = lmeControl(sigma = 1))

  I_E <- Fisher_info(mod)
  expect_is(I_E, "matrix")

  mod <- lme(d ~ 1, data = es_data,
             random = ~ 1 | es,
             weights = varFixed( ~ v),
             control = lmeControl(sigma = 1))

  I_E <- Fisher_info(mod)
  expect_is(I_E, "matrix")

  mod <- lme(d ~ 1, data = es_data,
             random = ~ 1 | study / es,
             weights = varFixed( ~ v),
             control = lmeControl(sigma = 1))

  I_E <- Fisher_info(mod)
  expect_is(I_E, "matrix")

})
