library(nlme)
data(Laski)

Laski_iid <- lme(fixed = outcome ~ treatment,
                 random = ~ treatment | case,
                 data = Laski)

Laski_het <- lme(fixed = outcome ~ treatment,
                 random = ~ treatment | case,
                 weights = varIdent(form = ~ 1 | treatment),
                 data = Laski)

Laski_AR1 <- lme(fixed = outcome ~ treatment,
                 random = ~ treatment | case,
                 correlation = corAR1(0, ~ time | case),
                 data = Laski)

Laski_hetAR1 <- lme(fixed = outcome ~ treatment,
                    random = ~ treatment | case,
                    correlation = corAR1(0, ~ time | case),
                    weights = varIdent(form = ~ 1 | treatment),
                    data = Laski)

Laski_CAR1 <- lme(fixed = outcome ~ treatment,
                  random = ~ treatment | case,
                  correlation = corCAR1(0.2, ~ time | case),
                  data = Laski)

Laski_MA1 <- lme(fixed = outcome ~ treatment,
                 random = ~ treatment | case,
                 correlation = corARMA(0, ~ time | case, p = 0, q = 1),
                 data = Laski)


extract_varcomp(Laski_iid)

mod <- Laski_iid
getGroups(mod)
grp <- Laski$case

expect_correct_dims <- function(mod, grp) {
  grp_size <- as.numeric(table(grp))
  dims <- sapply(targetVariance(mod), dim)
  expect_equal(grp_size, dims[1,], check.attributes = FALSE)
  expect_equal(grp_size, dims[2,], check.attributes = FALSE)
}

test_that("targetVariance() works with 2-level models.", {
  expect_correct_dims(Laski_iid, Laski$case)
  expect_correct_dims(Laski_het, Laski$case)
  expect_correct_dims(Laski_AR1, Laski$case)
  expect_correct_dims(Laski_hetAR1, Laski$case)
  expect_correct_dims(Laski_CAR1, Laski$case)
  expect_correct_dims(Laski_MA1, Laski$case)
})
