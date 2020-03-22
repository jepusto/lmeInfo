# fit some example models for varPower, varExp and varConstPower

library(nlme)
data(Orthodont)
data(Dialyzer)
data(BodyWeight)

# Orthodont

Ortho_A <- lme(distance ~ age + Sex,
               data = Orthodont,
               random = ~ age | Subject)

Ortho_B_Power <- update(Ortho_A, weights = varPower()) # fitted(.) is used by default
Ortho_C_Power <- update(Ortho_A, weights = varPower(form = ~ age))
Ortho_D_Power <- update(Ortho_A, weights = varPower(form = ~ age | Sex))
Ortho_B_Exp <- update(Ortho_A, weights = varExp())
Ortho_C_Exp <- update(Ortho_A, weights = varExp(form = ~ age))
Ortho_D_Exp <- update(Ortho_A, weights = varExp(form = ~ age | Sex))
Ortho_B_Const <- update(Ortho_A, weights = varConstPower()) # fitted(.) is used by default
Ortho_D_Const <- update(Ortho_A, weights = varConstPower(form = ~ age | Sex))
Ortho_D_Comb <- update(Ortho_A, weights = varComb(varIdent(form = ~1|Sex), varPower()))

test_that("targetVariance() works with Orthodont models.", {
  test_Sigma_mats(Ortho_A, Orthodont$Subject)
  test_Sigma_mats(Ortho_B_Power, Orthodont$Subject)
  test_Sigma_mats(Ortho_C_Power, Orthodont$Subject)
  test_Sigma_mats(Ortho_D_Power, Orthodont$Subject)
  test_Sigma_mats(Ortho_B_Exp, Orthodont$Subject)
  test_Sigma_mats(Ortho_C_Exp, Orthodont$Subject)
  test_Sigma_mats(Ortho_D_Exp, Orthodont$Subject)
  test_Sigma_mats(Ortho_B_Const, Orthodont$Subject)
  test_Sigma_mats(Ortho_D_Const, Orthodont$Subject)
  test_Sigma_mats(Ortho_D_Comb, Orthodont$Subject)
})

test_that("Derivative matrices are of correct dimension with Orthodont models.", {
  test_deriv_dims(Ortho_A)
  test_deriv_dims(Ortho_B_Power)
  test_deriv_dims(Ortho_C_Power)
  test_deriv_dims(Ortho_D_Power)
  test_deriv_dims(Ortho_B_Exp)
  test_deriv_dims(Ortho_C_Exp)
  test_deriv_dims(Ortho_D_Exp)
  test_deriv_dims(Ortho_B_Const)
  test_deriv_dims(Ortho_D_Const)
  expect_error(test_deriv_dims(Ortho_D_Comb))
})

test_that("Information matrices work with FIML with Orthodont models.", {
  test_with_FIML(Ortho_A)
  test_with_FIML(Ortho_B_Power)
  test_with_FIML(Ortho_C_Power)
  test_with_FIML(Ortho_D_Power)
  test_with_FIML(Ortho_B_Exp)
  test_with_FIML(Ortho_C_Exp)
  test_with_FIML(Ortho_D_Exp)
  test_with_FIML(Ortho_B_Const)
  expect_error(test_with_FIML(Ortho_D_Const))

})

test_that("Results do not depend on order of data.", {
  test_after_shuffling(Ortho_A, seed = 20)
  test_after_shuffling(Ortho_B_Power, seed = 21)
  test_after_shuffling(Ortho_C_Power, seed = 20)
  test_after_shuffling(Ortho_D_Power, seed = 20)
  test_after_shuffling(Ortho_B_Exp, seed = 20)
  test_after_shuffling(Ortho_C_Exp, seed = 20)
  test_after_shuffling(Ortho_D_Exp, seed = 20)
  test_after_shuffling(Ortho_B_Const, seed = 20)
  test_after_shuffling(Ortho_D_Const, seed = 21)

})


# Dialyzer

Dialyzer_A <- lme(rate ~ (pressure + I(pressure^2)) * QB,
                  data = Dialyzer,
                  random = ~ pressure | Subject)
Dialyzer_B_Power <- update(Dialyzer_A, weights = varPower())
Dialyzer_C_Power <- update(Dialyzer_A, weights = varPower(form = ~ pressure))
Dialyzer_D_Power <- update(Dialyzer_A, weights = varPower(form = ~ pressure | QB))

test_that("targetVariance() works with Dialyzer models.", {
  test_Sigma_mats(Dialyzer_A, Dialyzer$Subject)
  test_Sigma_mats(Dialyzer_B_Power, Dialyzer$Subject)
  test_Sigma_mats(Dialyzer_C_Power, Dialyzer$Subject)
  test_Sigma_mats(Dialyzer_D_Power, Dialyzer$Subject)
})

test_that("Derivative matrices are of correct dimension with Dialyzer models.", {
  test_deriv_dims(Dialyzer_A)
  test_deriv_dims(Dialyzer_B_Power)
  test_deriv_dims(Dialyzer_C_Power)
  test_deriv_dims(Dialyzer_D_Power)
})

test_that("Information matrices work with FIML with Dialyzer models.", {
  test_with_FIML(Dialyzer_A)
  test_with_FIML(Dialyzer_B_Power)
  test_with_FIML(Dialyzer_C_Power)
  test_with_FIML(Dialyzer_D_Power)
})

# BodyWeight

BodyWeight_A <- lme(weight ~ Time * Diet,
                    data = BodyWeight,
                    random = ~ Time | Rat)
Bodyweight_B_Power <- update(BodyWeight_A, weights = varPower()) # fitted(.) is used by default
Bodyweight_C_Power <- update(BodyWeight_A, weights = varPower(form = ~ Time))
Bodyweight_D_Power <- update(BodyWeight_A, weights = varPower(form = ~ Time | Diet))

test_that("targetVariance() works with BodyWeight models.", {
  test_Sigma_mats(BodyWeight_A, BodyWeight$Rat)
  test_Sigma_mats(Bodyweight_B_Power, BodyWeight$Rat)
  test_Sigma_mats(Bodyweight_C_Power, BodyWeight$Rat)
  test_Sigma_mats(Bodyweight_D_Power, BodyWeight$Rat)
})

test_that("Derivative matrices are of correct dimension with BodyWeight models.", {
  test_deriv_dims(BodyWeight_A)
  test_deriv_dims(Bodyweight_B_Power)
  test_deriv_dims(Bodyweight_C_Power)
  test_deriv_dims(Bodyweight_D_Power)
})

test_that("Information matrices work with FIML with Dialyzer models.", {
  test_with_FIML(BodyWeight_A)
  test_with_FIML(Bodyweight_B_Power)
  test_with_FIML(Bodyweight_C_Power)
  test_with_FIML(Bodyweight_D_Power)
})

