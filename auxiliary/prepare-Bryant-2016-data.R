# Read in Bryant 2016

Bryant16 <- read.csv("auxiliary/Bryant2016.csv", stringsAsFactors = FALSE)
str(Bryant16)

Bryant16 <- within(Bryant16, {
  school <- factor(school)
  case <- factor(case)
  treatment <- ifelse(treatment == "A", "baseline", "treatment")
  treatment <- factor(treatment)
})

save(Bryant16, file = "data/Bryant2016.RData", compress = TRUE)
