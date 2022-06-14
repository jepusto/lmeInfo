pkgname <- "scdhlm"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('scdhlm')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("CI_g")
### * CI_g

flush(stderr()); flush(stdout())

### Name: CI_g
### Title: Calculates a confidence interval for a standardized mean
###   difference effect size
### Aliases: CI_g

### ** Examples

data(Laski)
Laski_RML <- lme(fixed = outcome ~ treatment,
                 random = ~ 1 | case,
                 correlation = corAR1(0, ~ time | case),
                 data = Laski)
Laski_g_REML <- suppressWarnings(
  g_REML(Laski_RML, p_const = c(0,1), 
         r_const = c(1,0,1), returnModel = FALSE)
)
CI_g(Laski_g_REML, symmetric = TRUE)
CI_g(Laski_g_REML, symmetric = FALSE)

Laski_HPS <- with(Laski, effect_size_MB(outcome, treatment, case, time))
CI_g(Laski_HPS, symmetric = FALSE)

Laski_g_mlm <- g_mlm(Laski_RML, p_const = c(0,1), r_const = c(1,0,1), returnModel = TRUE)
CI_g(Laski_g_mlm, symmetric = FALSE)




cleanEx()
nameEx("Info_Expected_lmeAR1")
### * Info_Expected_lmeAR1

flush(stderr()); flush(stdout())

### Name: Info_Expected_lmeAR1
### Title: Calculate expected information matrix
### Aliases: Info_Expected_lmeAR1

### ** Examples

data(Laski)
Laski_RML <- lme(fixed = outcome ~ treatment, 
                 random = ~ 1 | case, 
                 correlation = corAR1(0, ~ time | case), 
                 data = Laski)
Info_Expected_lmeAR1(Laski_RML)



cleanEx()
nameEx("compare_RML_HPS")
### * compare_RML_HPS

flush(stderr()); flush(stdout())

### Name: compare_RML_HPS
### Title: Run simulation comparing REML and HPS estimates
### Aliases: compare_RML_HPS

### ** Examples

compare_RML_HPS(iterations=10, beta = c(0,1,0,0), rho = 0.3, 
                 phi = 0.5, design=design_matrix(m=3,n=8))



cleanEx()
nameEx("design_matrix")
### * design_matrix

flush(stderr()); flush(stdout())

### Name: design_matrix
### Title: Create a design matrix for a single-case design
### Aliases: design_matrix

### ** Examples

design_matrix(3, 16, c(5,9,13))



cleanEx()
nameEx("effect_size_ABk")
### * effect_size_ABk

flush(stderr()); flush(stdout())

### Name: effect_size_ABk
### Title: Calculates HPS effect size
### Aliases: effect_size_ABk

### ** Examples

data(Lambert)
effect_size_ABk(outcome = outcome, treatment = treatment, id = case, 
                phase = phase, time = time, data = Lambert)
   
data(Anglesea)
effect_size_ABk(outcome = outcome, treatment = condition, id = case, 
                phase = phase, time = session, data = Anglesea)




cleanEx()
nameEx("effect_size_MB")
### * effect_size_MB

flush(stderr()); flush(stdout())

### Name: effect_size_MB
### Title: Calculates HPS effect size
### Aliases: effect_size_MB

### ** Examples

data(Saddler)
effect_size_MB(outcome = outcome, treatment = treatment, id = case, 
               time = time, data = subset(Saddler, measure=="writing quality"))

data(Laski)
effect_size_MB(outcome = outcome, treatment = treatment, id = case, 
               time = time, data = Laski)




cleanEx()
nameEx("g_REML")
### * g_REML

flush(stderr()); flush(stdout())

### Name: g_REML
### Title: Calculates adjusted REML effect size
### Aliases: g_REML

### ** Examples

data(Laski)
Laski_RML <- lme(fixed = outcome ~ treatment, 
                 random = ~ 1 | case, 
                 correlation = corAR1(0, ~ time | case), 
                 data = Laski)
summary(Laski_RML)
g_REML(Laski_RML, p_const = c(0,1), r_const = c(1,0,1), returnModel=FALSE)

data(Schutte)
Schutte$trt.week <- with(Schutte, unlist(tapply((treatment=="treatment") * week, 
         list(treatment,case), function(x) x - min(x))) + (treatment=="treatment"))
Schutte$week <- Schutte$week - 9
Schutte_RML <- lme(fixed = fatigue ~ week + treatment + trt.week, 
                   random = ~ week | case, 
                   correlation = corAR1(0, ~ week | case), 
                   data = subset(Schutte, case != 4))
summary(Schutte_RML)
Schutte_g <- g_REML(Schutte_RML, p_const = c(0,0,1,7), r_const = c(1,0,1,0,0))
summary(Schutte_g)



cleanEx()
nameEx("graph_SCD")
### * graph_SCD

flush(stderr()); flush(stdout())

### Name: graph_SCD
### Title: Graph Single Case Design Data
### Aliases: graph_SCD

### ** Examples

data(Anglesea)
graph_SCD(case=case, phase=condition, 
          session=session, outcome=outcome, 
          design="TR", treatment_name = "treatment", 
          data=Anglesea)
          
data(BartonArwood)
graph_SCD(case=case, phase=condition, 
          session=session, outcome=outcome, 
          design="MB", treatment_name = "B",  
          data=BartonArwood)




cleanEx()
nameEx("phase_pairs")
### * phase_pairs

flush(stderr()); flush(stdout())

### Name: phase_pairs
### Title: Calculate phase-pairs for a unique case
### Aliases: phase_pairs

### ** Examples


phases <- rep(c("A","B","A","B"), each = 4)
sessions <- 1:length(phases)

phase_pairs(phases, sessions)

phases <- rep(c("A","B","C","A","B","C","D"), each = 4)
phase_pairs(phases)

phases <- rep(c("B","A","C","B","A","B","C","A"), each = 4)
phase_pairs(phases)




cleanEx()
nameEx("preprocess_SCD")
### * preprocess_SCD

flush(stderr()); flush(stdout())

### Name: preprocess_SCD
### Title: Clean Single Case Design Data
### Aliases: preprocess_SCD

### ** Examples

data(Laski)
preprocess_SCD(case = case, phase = treatment,
               session = time, outcome = outcome, 
               design = "MB", center = 4, data = Laski)

          



cleanEx()
nameEx("shine_scd")
### * shine_scd

flush(stderr()); flush(stdout())

### Name: shine_scd
### Title: A shiny interface for the scdhlm package
### Aliases: shine_scd

### ** Examples

## Not run: 
##D shine_scd()
##D data(Laski)
##D shine_scd(dataset = Laski)
##D shine_scd(dataset = "SCD_data.xlsx", sheet = "Laski")
##D shine_scd(dataset = "Laski.csv") 
## End(Not run)




cleanEx()
nameEx("simulate.g_REML")
### * simulate.g_REML

flush(stderr()); flush(stdout())

### Name: simulate.g_REML
### Title: Simulate data from a fitted 'g_REML' object
### Aliases: simulate.g_REML

### ** Examples

data(Laski)
Laski_RML <- lme(fixed = outcome ~ treatment,
                 random = ~ 1 | case,
                 correlation = corAR1(0, ~ time | case), 
                 data = Laski)
Laski_g <- g_REML(Laski_RML, p_const = c(0,1), r_const = c(1,0,1))
simulate(Laski_g, nsim = 20)




cleanEx()
nameEx("simulate_MB2")
### * simulate_MB2

flush(stderr()); flush(stdout())

### Name: simulate_MB2
### Title: Simulate Model MB2 from Pustejovsky, Hedges, & Shadish (2014)
### Aliases: simulate_MB2

### ** Examples


set.seed(8)
simulate_MB2(iterations = 5, beta = c(0,1,0,0), rho = 0.4, phi = 0.5, 
             tau1_ratio = 0.5, tau_corr = -0.4, design = design_matrix(m=3, n=8))
             
set.seed(8)
simulate_MB2(iterations = 5, beta = c(0,1,0,0), rho = 0.4, phi = 0.5, 
             tau1_ratio = 0.5, tau_corr = -0.4, m = 3, n = 8, MB = FALSE)
             



cleanEx()
nameEx("simulate_MB4")
### * simulate_MB4

flush(stderr()); flush(stdout())

### Name: simulate_MB4
### Title: Simulate Model MB4 from Pustejovsky, Hedges, & Shadish (2014)
### Aliases: simulate_MB4

### ** Examples


simulate_MB4(iterations = 5, beta = c(0,1,0,0), rho = 0.8, phi = 0.5, 
             tau2_ratio = 0.5, tau_corr = 0, 
             p_const = c(0,1,0,7), r_const = c(1,0,1,0,0), 
             design = design_matrix(3, 16, treat_times=c(5,9,13), center = 12))
             
simulate_MB4(iterations = 5, beta = c(0,1,0,0), rho = 0.8, phi = 0.5, 
             tau2_ratio = 0.5, tau_corr = 0, m = 6, n = 8)
             



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
