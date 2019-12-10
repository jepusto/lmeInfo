# Read in Thiemann 2004

Thiemann2004 <- read.csv("auxiliary/Thiemann2004.csv", stringsAsFactors = FALSE)
str(Thiemann2004)

Thiemann2004 <- within(Thiemann2004, {
  case <- factor(case)
  series <- factor(series)
  treatment <- ifelse(treatment == "A", "baseline", "treatment")
  treatment <- factor(treatment)
})

save(Thiemann2004, file = "data/Thiemann2004.RData", compress = TRUE)
