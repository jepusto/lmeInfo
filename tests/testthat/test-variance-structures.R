# fit some example models

Ortho_A <- lme(distance ~ age + Sex,
               data = Orthodont,
               random = ~ age | Subject)

Ortho_B <- update(Ortho_A, weights = varPower(form = ~ age | Sex))

summary(Ortho_B)

mod <- Ortho_B
struct <- mod$modelStruct$varStruct
struct
str(struct)
