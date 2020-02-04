# fit some example models for varPower, varExp and varConstPower

# dat_A is the model without specifying a var function
# dat_B is the model using fitted(.) as the covariate by default
# dat_C is the model using a covariate from the data set but not specifying the stratum
# dat_D is the model specifying a covariate and a stratum

library(nlme)

# for varPower

# Orthodont

Ortho_A <- lme(distance ~ age + Sex,
               data = Orthodont,
               random = ~ age | Subject)

Ortho_B_Power <- update(Ortho_A, weights = varPower()) # fitted(.) is used by default
summary(Ortho_B_Power)
struct <- Ortho_B_Power$modelStruct$varStruct
dsd_dvarStruct(struct)

Ortho_C_Power <- update(Ortho_A, weights = varPower(form = ~ age))
summary(Ortho_C_Power)
struct <- Ortho_C_Power$modelStruct$varStruct
dsd_dvarStruct(struct)

Ortho_D_Power <- update(Ortho_A, weights = varPower(form = ~ age | Sex))
summary(Ortho_D_Power)
struct <- Ortho_D_Power$modelStruct$varStruct
dsd_dvarStruct(struct)

# Dialyzer

Dialyzer_A <- lme(rate ~ (pressure + I(pressure^2)) * QB,
                  data = Dialyzer,
                  random = ~ pressure | Subject)

Dialyzer_B_Power <- update(Dialyzer_A, weights = varPower())
summary(Dialyzer_B_Power)
struct <- Dialyzer_B_Power$modelStruct$varStruct
dsd_dvarStruct(struct)

Dialyzer_C_Power <- update(Dialyzer_A, weights = varPower(form = ~ pressure))
summary(Dialyzer_C_Power)
struct <- Dialyzer_C_Power$modelStruct$varStruct
dsd_dvarStruct(struct)

Dialyzer_D_Power <- update(Dialyzer_A, weights = varPower(form = ~ pressure | QB))
summary(Dialyzer_D_Power)
struct <- Dialyzer_D_Power$modelStruct$varStruct
dsd_dvarStruct(struct)

# BodyWeight

BodyWeight_A <- lme(weight ~ Time * Diet,
                    data = BodyWeight,
                    random = ~ Time | Rat)

Bodyweight_B_Power <- update(BodyWeight_A, weights = varPower()) # fitted(.) is used by default
summary(Bodyweight_B_Power)
struct <- Bodyweight_B_Power$modelStruct$varStruct
dsd_dvarStruct(struct)

Bodyweight_C_Power <- update(BodyWeight_A, weights = varPower(form = ~ Time))
summary(Bodyweight_C_Power)
struct <- Bodyweight_C_Power$modelStruct$varStruct
dsd_dvarStruct(struct)

Bodyweight_D_Power <- update(BodyWeight_A, weights = varPower(form = ~ Time | Diet))
summary(Bodyweight_D_Power)
struct <- Bodyweight_D_Power$modelStruct$varStruct
dsd_dvarStruct(struct)

mod <- Bodyweight_B_Power
struct <- mod$modelStruct$varStruct
struct
str(struct)
attr(struct, "covariate") # the variance covariate
attr(struct, "weights") # the inverse of variance function

# for varExp

Ortho_B_exp <- update(Ortho_A, weights = varExp())
summary(Ortho_B_exp)
struct <- Ortho_B_exp$modelStruct$varStruct
dsd_dvarStruct(struct)

Ortho_C_exp <- update(Ortho_A, weights = varExp(form = ~ age))
summary(Ortho_C_exp)
struct <- Ortho_C_exp$modelStruct$varStruct
dsd_dvarStruct(struct)

Ortho_D_exp <- update(Ortho_A, weights = varExp(form = ~ age | Sex))
summary(Ortho_D_exp)
struct <- Ortho_D_exp$modelStruct$varStruct
dsd_dvarStruct(struct)

# for varConstPower

Ortho_B_Const <- update(Ortho_A, weights = varConstPower()) # fitted(.) is used by default
summary(Ortho_B_Const)
struct <- Ortho_B_Const$modelStruct$varStruct
dsd_dvarStruct(struct)

Ortho_D_Const <- update(Ortho_A, weights = varConstPower(form = ~ age | Sex))
summary(Ortho_D_Const)
struct <- Ortho_D_Const$modelStruct$varStruct
dsd_dvarStruct(struct)

mod <- Ortho_B_Const
struct <- mod$modelStruct$varStruct
struct
str(struct)
pars <- coef(struct, FALSE)
pars
attr(struct, "covariate") # the variance covariate
attr(struct, "weights") # the inverse of variance function
dsd_dvarStruct(struct) # the derivative results


