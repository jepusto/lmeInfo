# fit some example models

# for varPower

# Orthodont

Ortho_A <- lme(distance ~ age + Sex,
               data = Orthodont,
               random = ~ age | Subject)

Ortho_B <- update(Ortho_A, weights = varPower(form = ~ age))
summary(Ortho_B)
struct <- Ortho_B$modelStruct$varStruct
dsd_dvarStruct(struct)

Ortho_C <- update(Ortho_A, weights = varPower(form = ~ age | Sex))
summary(Ortho_C)
struct <- Ortho_C$modelStruct$varStruct
dsd_dvarStruct(struct)

# Dialyzer

Dialyzer_A <- lme(rate ~ (pressure + I(pressure^2)) * QB,
                  data = Dialyzer,
                  random = ~ pressure | Subject)

Dialyzer_B <- update(Dialyzer_A, weights = varPower(form = ~ pressure))
summary(Dialyzer_B)
struct <- Dialyzer_B$modelStruct$varStruct
dsd_dvarStruct(struct)

Dialyzer_C <- update(Dialyzer_A, weights = varPower(form = ~ pressure | QB))
summary(Dialyzer_C)
struct <- Dialyzer_C$modelStruct$varStruct
dsd_dvarStruct(struct)

# BodyWeight

BodyWeight_A <- lme(weight ~ Time * Diet,
                    data = BodyWeight,
                    random = ~ Time | Rat)

Bodyweight_B <- update(BodyWeight_A, weights = varPower()) # fitted(.) is used by default
summary(Bodyweight_B)
struct <- Bodyweight_B$modelStruct$varStruct
dsd_dvarStruct(struct)

mod <- Bodyweight_B
struct <- mod$modelStruct$varStruct
struct
str(struct)
attr(struct, "covariate") # the variance covariate
attr(struct, "weights") # the inverse of variance function
dsd_dvarStruct(struct) # the derivative results

# for varExp

Ortho_A <- lme(distance ~ age + Sex,
               data = Orthodont,
               random = ~ age | Subject)

Ortho_B_exp <- update(Ortho_A, weights = varExp(form = ~ age))
summary(Ortho_B_exp)
struct <- Ortho_B_exp$modelStruct$varStruct
dsd_dvarStruct(struct)

Ortho_C_exp <- update(Ortho_A, weights = varExp(form = ~ age | Sex))
summary(Ortho_C_exp)
struct <- Ortho_C_exp$modelStruct$varStruct
dsd_dvarStruct(struct)

Ortho_D_exp <- update(Ortho_A, weights = varExp())
summary(Ortho_D_exp)
struct <- Ortho_D_exp$modelStruct$varStruct
dsd_dvarStruct(struct)

# for varConstPower

Ortho_C_Const <- update(Ortho_A, weights = varConstPower(form = ~ age | Sex))
summary(Ortho_C_Const)
struct <- Ortho_C_Const$modelStruct$varStruct
dsd_dvarStruct(struct)

Dialyzer_C_Const <- update(Dialyzer_A, weights = varConstPower(form = ~ pressure | QB))
summary(Dialyzer_C_Const)
struct <- Dialyzer_C_Const$modelStruct$varStruct
dsd_dvarStruct(struct)
