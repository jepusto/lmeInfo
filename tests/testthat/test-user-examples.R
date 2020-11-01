library(nlme)

test_that("scdhlm Works with Narozanic and Blair.", {

  skip("Don't worry about this user example.")

  # Naro <- read.csv("auxiliary/Narozanic and Blair_datasets_long_academic.csv", stringsAsFactors = FALSE)
  Naro <- read.csv("../auxiliary/Narozanic and Blair_datasets_long_academic.csv", stringsAsFactors = FALSE)
  Naro$Session_int <- round(Naro$Session)

  Naro_dbl <- lme(Outcome ~ Phase,
                  random = ~ 1 | Participant,
                  correlation = corAR1(0.1, ~ Session | Participant),
                  data = Naro)
  Naro_int <- lme(Outcome ~ Phase,
                  random = ~ 1 | Participant,
                  correlation = corAR1(0.1, ~ Session_int | Participant),
                  data = Naro)

  # Errors because corAR1 _truncates_ Session, leading to zero time between sessions and NPD correlation matrices
  expect_error(Fisher_info(Naro_dbl))

  expect_identical(dim(Fisher_info(Naro_int)), c(3L,3L))

})
