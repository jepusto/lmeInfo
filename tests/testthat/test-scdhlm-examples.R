context("Examples from scdhlm")

library(scdhlm)

data(Lambert)
Lambert_RML <- lme(fixed = outcome ~ treatment,
                    random = ~ 1 | case,
                    correlation = corAR1(0.1, ~ time | case),
                    data = Lambert)


data(Anglesea)
Anglesea_RML <- lme(fixed = outcome ~ condition,
                    random = ~ condition | case,
                    correlation = corAR1(0.2, ~ session | case),
                    data = Anglesea)


data(Saddler)
Saddler_quality_RML <- lme(fixed = outcome ~ treatment,
                           random = ~ 1 | case,
                           correlation = corAR1(0, ~ time | case),
                           data = subset(Saddler, measure=="writing quality"))


data(Laski)
# create trt-by-time interaction
Laski$trt_time <- with(Laski,
                       unlist(tapply((treatment=="treatment") * time,
                       list(treatment, case),
                       function(x) x - min(x))) + (treatment=="treatment"))
# time-point constants
A_Laski <- 4
B_Laski <- 12
# center at follow-up time (time 12)
Center_Laski <- B_Laski
Laski$time <- Laski$time - Center_Laski

# Varying intercepts, fixed treatment effect, fixed trends
Laski_RML3 <- lme(fixed = outcome ~ time + treatment + trt_time,
                  random = ~ 1 | case,
                  correlation = corAR1(0.1, ~ time | case),
                  data = Laski)

# Varying intercepts, fixed treatment effect, varying trends

Laski_RML4 <- suppressWarnings(
  update(Laski_RML3, random = ~ time | case,
         control=lmeControl(msMaxIter = 50, apVar=FALSE, returnObject=TRUE))
)


data(Schutte)
Schutte <- subset(Schutte, case != "Case 4")
Schutte$case <- factor(Schutte$case)
Schutte$trt_week <- with(Schutte,
                         unlist(tapply((treatment=="treatment") * week,
                         list(treatment, case),
                         function(x) x - min(x))) + (treatment=="treatment"))
Schutte$week <- Schutte$week - 9 # center at follow-up time (week 9)

# Varying intercepts, fixed treatment effect, fixed trends
Schutte_RML3 <- lme(fixed = fatigue ~ week + treatment + trt_week,
                    random = ~ 1 | case,
                    correlation = corAR1(0, ~ week | case),
                    data = Schutte,
                    method = "REML")

# Varying intercepts, fixed treatment effect, varying trends
Schutte_RML4 <- update(Schutte_RML3, random = ~ week | case)

# Varying intercepts, varying trends, varying treatment-by-time interactions
Schutte_RML5 <- suppressWarnings(
  update(Schutte_RML4, random = ~ week + trt_week | case,
         control=lmeControl(msMaxIter = 50, apVar=FALSE, returnObject=TRUE))
)


test_that("lmeinfo::g_REML returns the same result as scdhlm::g_REML.", {
  check_against_scdhlm(Lambert_RML,
                       p_lmeInfo = c(0,1), r_lmeInfo = c(1,0,1),
                       p_scdhlm = c(0,1), r_scdhlm = c(1,0,1))

  check_against_scdhlm(Anglesea_RML,
                       p_lmeInfo = c(0,1), r_lmeInfo = c(1,0,0,0,1),
                       p_scdhlm = c(0,1), r_scdhlm = c(1,0,1,0,0))

  check_against_scdhlm(Saddler_quality_RML,
                       p_lmeInfo = c(0,1), r_lmeInfo = c(1,0,1),
                       p_scdhlm = c(0,1), r_scdhlm = c(1,0,1))

  # Laski
  check_against_scdhlm(Laski_RML3,
                       p_lmeInfo = c(0,0,1,8), r_lmeInfo = c(1,0,1),
                       p_scdhlm = c(0,0,1,8), r_scdhlm = c(1,0,1))

  check_against_scdhlm(Laski_RML4,
                       p_lmeInfo = c(0,0,1,8), r_lmeInfo = c(1,0,0,0,1),
                       p_scdhlm = c(0,0,1,8), r_scdhlm = c(1,0,1,0,0))

  # Schutte
  check_against_scdhlm(Schutte_RML3,
                       p_lmeInfo = c(0,0,1,7), r_lmeInfo = c(1,0,1),
                       p_scdhlm = c(0,0,1,7), r_scdhlm = c(1,0,1))

  check_against_scdhlm(Schutte_RML4,
                       p_lmeInfo = c(0,0,1,7), r_lmeInfo = c(1,0,0,0,1),
                       p_scdhlm = c(0,0,1,7), r_scdhlm = c(1,0,1,0,0))

  check_against_scdhlm(Schutte_RML5,
                       p_lmeInfo = c(0,0,1,7), r_lmeInfo = c(1,0,0,0,0,0,0,1),
                       p_scdhlm = c(0,0,1,7), r_scdhlm = c(1,0,1,0,0,0,0,0))

})

