# fit some example models for varPower, varExp and varConstPower

# dat_A is the model without specifying a var function
# dat_B is the model using fitted(.) as the covariate by default
# dat_C is the model using a covariate from the data set but not specifying the stratum
# dat_D is the model specifying a covariate and a stratum

library(nlme)
data(Orthodont)
data(Dialyzer)
data(BodyWeight)

# for varIdent
Laski_het <- lme(fixed = outcome ~ treatment,
                 random = ~ treatment | case,
                 weights = varIdent(form = ~ 1 | treatment),
                 data = Laski)

Laski_het$modelStruct$varStruct
dsd_dvarStruct(Laski_het$modelStruct$varStruct)
dV_dvarStruct(Laski_het)

# for varPower

# Orthodont

Ortho_A <- lme(distance ~ age + Sex,
               data = Orthodont,
               random = ~ age | Subject)
dsd_dvarStruct(Ortho_A$modelStruct$varStruct)
dV_dvarStruct(Ortho_A)

Ortho_B_Power <- update(Ortho_A, weights = varPower()) # fitted(.) is used by default
summary(Ortho_B_Power)
dsd_dvarStruct(Ortho_B_Power$modelStruct$varStruct)
dV_dvarStruct(Ortho_B_Power)

Ortho_C_Power <- update(Ortho_A, weights = varPower(form = ~ age))
summary(Ortho_C_Power)
dsd_dvarStruct(Ortho_C_Power$modelStruct$varStruct)
dV_dvarStruct(Ortho_C_Power)

Ortho_D_Power <- update(Ortho_A, weights = varPower(form = ~ age | Sex))
summary(Ortho_D_Power)
dsd_dvarStruct(Ortho_D_Power$modelStruct$varStruct)
dV_dvarStruct(Ortho_D_Power)

# Dialyzer

Dialyzer_A <- lme(rate ~ (pressure + I(pressure^2)) * QB,
                  data = Dialyzer,
                  random = ~ pressure | Subject)
dV_dvarStruct(Dialyzer_A)


Dialyzer_B_Power <- update(Dialyzer_A, weights = varPower())
summary(Dialyzer_B_Power)
dsd_dvarStruct(Dialyzer_B_Power$modelStruct$varStruct)
dV_dvarStruct(Dialyzer_B_Power)

Dialyzer_C_Power <- update(Dialyzer_A, weights = varPower(form = ~ pressure))
summary(Dialyzer_C_Power)
dsd_dvarStruct(Dialyzer_C_Power$modelStruct$varStruct)
dV_dvarStruct(Dialyzer_C_Power)

Dialyzer_D_Power <- update(Dialyzer_A, weights = varPower(form = ~ pressure | QB))
summary(Dialyzer_D_Power)
dsd_dvarStruct(Dialyzer_D_Power$modelStruct$varStruct)
dV_dvarStruct(Dialyzer_D_Power)
mod <- Dialyzer_D_Power


# BodyWeight

BodyWeight_A <- lme(weight ~ Time * Diet,
                    data = BodyWeight,
                    random = ~ Time | Rat)
dV_dvarStruct(BodyWeight_A)

Bodyweight_B_Power <- update(BodyWeight_A, weights = varPower()) # fitted(.) is used by default
summary(Bodyweight_B_Power)
struct <- Bodyweight_B_Power$modelStruct$varStruct
dsd_dvarStruct(struct)
dV_dvarStruct(Bodyweight_B_Power)

Bodyweight_C_Power <- update(BodyWeight_A, weights = varPower(form = ~ Time))
summary(Bodyweight_C_Power)
struct <- Bodyweight_C_Power$modelStruct$varStruct
dsd_dvarStruct(struct)
dV_dvarStruct(Bodyweight_C_Power)

Bodyweight_D_Power <- update(BodyWeight_A, weights = varPower(form = ~ Time | Diet))
summary(Bodyweight_D_Power)
struct <- Bodyweight_D_Power$modelStruct$varStruct
dsd_dvarStruct(struct)
dV_dvarStruct(Bodyweight_D_Power)

mod <- Bodyweight_B_Power
struct <- mod$modelStruct$varStruct
struct
str(struct)
attr(struct, "covariate") # the variance covariate
attr(struct, "weights") # the inverse of variance function

# for varExp

Ortho_B_Exp <- update(Ortho_A, weights = varExp())
summary(Ortho_B_Exp)
dsd_dvarStruct(Ortho_B_Exp$modelStruct$varStruct)
dV_dvarStruct(Ortho_B_Exp)

Ortho_C_Exp <- update(Ortho_A, weights = varExp(form = ~ age))
summary(Ortho_C_Exp)
dsd_dvarStruct(Ortho_C_Exp$modelStruct$varStruct)
dV_dvarStruct(Ortho_C_Exp)

Ortho_D_Exp <- update(Ortho_A, weights = varExp(form = ~ age | Sex))
summary(Ortho_D_Exp)
dsd_dvarStruct(Ortho_D_Exp$modelStruct$varStruct)
dV_dvarStruct(Ortho_D_Exp)

# for varConstPower

Ortho_B_Const <- update(Ortho_A, weights = varConstPower()) # fitted(.) is used by default
summary(Ortho_B_Const)
struct <- Ortho_B_Const$modelStruct$varStruct
dsd_dvarStruct(struct)
dV_dvarStruct(Ortho_B_Const)

Ortho_D_Const <- update(Ortho_A, weights = varConstPower(form = ~ age | Sex))
summary(Ortho_D_Const)
struct <- Ortho_D_Const$modelStruct$varStruct
dsd_dvarStruct(struct)
dV_dvarStruct(Ortho_D_Const)

mod <- Ortho_B_Const
struct <- mod$modelStruct$varStruct
struct
str(struct)
pars <- coef(struct, FALSE)
pars
attr(struct, "covariate") # the variance covariate
attr(struct, "weights") # the inverse of variance function
dsd_dvarStruct(struct) # the derivative results
