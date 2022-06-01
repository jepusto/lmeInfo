library(readxl)
# library(lmeInfo)
library(nlme)


#-------------------------------------------------------------------------------
# Mini Study 1

# Data cleaning
MLM1_To_A <- read_xlsx("auxiliary/MLM1_To_A.xlsx")
MLM1_To_A$Towre.c <- MLM1_To_A$Towre - mean(MLM1_To_A$Towre)
MLM1_To_A$Pretest.c <- MLM1_To_A$Pretest - mean(MLM1_To_A$Pretest)

# Fit basic model
MS1_basic <- lme(fixed = Posttest ~ CONDITION,
                 random = ~ 1 | SCHOOL,
                 data = MLM1_To_A)
summary(MS1_basic)

# Fit conditional model
MS1 <- lme(fixed = Posttest ~ CONDITION + Pretest.c + Towre.c,
           random = ~ 1 | SCHOOL,
           data = MLM1_To_A)
summary(MS1)

# Effect size calculations
MS1_g1 <- g_mlm(mod = MS1, p_const = c(0,1,0,0),
                mod_denom = MS1_basic, r_const = c(1,1),
                infotype= "expected")
summary(MS1_g1)
CI_g(MS1_g1, symmetric = TRUE)

#-------------------------------------------------------------------------------
# Mini Study 3

# Data cleaning

MLM3_T_A <- read_xlsx("auxiliary/MLM3_T_A.xlsx")
MLM3_T_A$PRETEST.c <- MLM3_T_A$PRETEST - mean(MLM3_T_A$PRETEST)
MLM3_T_A$TOWRE.c <- MLM3_T_A$TOWRE - mean(MLM3_T_A$TOWRE)

# Fit basic model
MS3_basic <- lme(fixed = POSTTEST ~ CONDITION,
                 random = ~ 1 | TEACHER,
                 data = MLM3_T_A)
summary(MS3_basic)

# Fit conditional model

MS3 <- lme(fixed = POSTTEST ~ CONDITION + PRETEST.c,
             random = ~ 1 | TEACHER,
             data = MLM3_T_A)
summary(MS3)

# Effect size calculations
MS3_g1 <- g_mlm(MS3, p_const = c(0,1,0),
                mod_denom = MS3_basic, r_const = c(1,1),
                infotype = "expected")
summary(MS3_g1)
CI_g(MS3_g1, symmetric = TRUE)
