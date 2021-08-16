library(nlme, quietly=TRUE, warn.conflicts=FALSE)

data(star, package = "mlmRev")

star <- subset(star, gr == 3 & !is.na(math))
star <- droplevels(star)
star <- star[order(star$sch, star$tch, star$id),]

star$small <- ifelse(star$cltype == "small", 1L, 0L)
star$ses <- ifelse(is.na(star$ses), "M", as.character(star$ses))
star$eth <- ifelse(is.na(star$eth), "M", as.character(star$eth))

star_2L_basic <- lme(math ~ small,
                     random = ~ 1 | sch,
                     data = star)
star_2L_control <- lme(math ~ small + schtype + ses + sx + eth,
                       random = ~ 1 | sch,
                        data = star)

star_3L_basic <- lme(math ~ small,
                     random = ~ 1 | sch / tch,
                     data = star)
star_3L_control <- lme(math ~ small + schtype + ses + sx + eth,
                       random = ~ 1 | sch / tch,
                       data = star)


star_3L_RE <- lme(math ~ small,
                  random = list(~ 1 | sch, ~ 0 + small | sch, ~ 1 | tch),
                  data = star)

star_3L_RE_control <- lme(math ~ small + schtype + ses + sx + eth,
                          random = list(~ 1 | sch, ~ 0 + small | sch, ~ 1 | tch),
                          data = star)


test_that("Fisher information matrices can be computed for STAR models.", {

  check_info_dim(star_2L_basic, 2L)
  check_info_dim(star_2L_control, 2L)
  check_info_dim(star_3L_basic, 3L)
  check_info_dim(star_3L_control, 3L)
  check_info_dim(star_3L_RE, 4L)
  check_info_dim(star_3L_RE_control, 4L)

})

test_that("g_mlm works for STAR models.", {

  g2_basic <- g_mlm(star_2L_basic, p_const = c(0, 1), r_const = c(1, 1))
  g2_explicit <- g_mlm(star_2L_basic, p_const = c(0, 1),
                       mod_denom = star_2L_basic, r_const = c(1, 1))
  g2_conditional <- g_mlm(mod = star_2L_control, p_const = c(0, 1, rep(0,11)),
                          r_const = c(1, 1))
  g2_control <- g_mlm(mod = star_2L_control, p_const = c(0, 1, rep(0,11)),
                      mod_denom = star_2L_basic, r_const = c(1, 1))

  expect_identical(g2_basic, g2_explicit)
  expect_gt(abs(g2_basic$g_AB - g2_control$g_AB), 0)
  expect_gt(g2_control$g_AB, g2_conditional$g_AB)


  g3_basic <- g_mlm(star_3L_basic, p_const = c(0, 1), r_const = c(1, 1, 1))
  g3_explicit <- g_mlm(star_3L_basic, p_const = c(0, 1),
                       mod_denom = star_3L_basic, r_const = c(1, 1, 1))
  g3_conditional <- g_mlm(mod = star_3L_control, p_const = c(0, 1, rep(0,11)),
                          r_const = c(1, 1, 1))
  g3_control <- g_mlm(mod = star_3L_control, p_const = c(0, 1, rep(0,11)),
                      mod_denom = star_3L_basic, r_const = c(1, 1, 1))

  expect_identical(g3_basic, g3_explicit)
  expect_gt(abs(g3_basic$g_AB - g3_control$g_AB), 0)
  expect_gt(g3_control$g_AB, g3_conditional$g_AB)

  gRE_basic <- g_mlm(star_3L_RE, p_const = c(0, 1), r_const = c(0, 1, 1, 1))
  gRE_explicit <- g_mlm(star_3L_RE, p_const = c(0, 1),
                       mod_denom = star_3L_RE, r_const = c(0, 1, 1, 1))
  gRE_conditional <- g_mlm(mod = star_3L_RE_control, p_const = c(0, 1, rep(0,11)),
                          r_const = c(0, 1, 1, 1))
  gRE_control <- g_mlm(mod = star_3L_RE_control, p_const = c(0, 1, rep(0,11)),
                      mod_denom = star_3L_RE, r_const = c(0, 1, 1, 1))

  expect_identical(gRE_basic, gRE_explicit)
  expect_gt(abs(gRE_basic$g_AB - gRE_control$g_AB), 0)
  expect_gt(gRE_control$g_AB, gRE_conditional$g_AB)

})
