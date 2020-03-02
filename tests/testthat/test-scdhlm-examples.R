context("Examples from scdhlm")

library(scdhlm)

data(Anglesea)
Anglesea_RML <- lme(fixed = outcome ~ condition,
                    random = ~ condition | case,
                    correlation = corAR1(0.2, ~ session | case),
                    data = Anglesea)

data(Schutte)
Schutte$trt.week <- with(Schutte,
                         unlist(tapply((treatment=="treatment") * week,
                         list(treatment,case), function(x) x - min(x))) + (treatment=="treatment"))
Schutte$week <- Schutte$week - 9
Schutte_RML <- lme(fixed = fatigue ~ week + treatment + trt.week,
                   random = ~ week | case,
                   correlation = corAR1(0.2, ~ week | case),
                   data = subset(Schutte, case != 4))

test_that("lmeinfo::g_REML returns the same result as scdhlm::g_REML.", {

  test_gREML(Anglesea_RML,
             p_lmeInfo = c(0,1), r_lmeInfo = c(1,0,0,0,1),
             p_scdhlm = c(0,1), r_scdhlm = c(1,0,1,0,0))

  test_gREML(Schutte_RML,
             p_lmeInfo = c(0,0,1,7), r_lmeInfo = c(1,0,0,0,1),
             p_scdhlm = c(0,0,1,7), r_scdhlm = c(1,0,1,0,0))
})

