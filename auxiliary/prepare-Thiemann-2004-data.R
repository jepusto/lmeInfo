# Read in Thiemann 2001

Thiemann2001 <- read.csv("auxiliary/Thiemann2001.csv", stringsAsFactors = FALSE)
str(Thiemann2001)

Thiemann2001 <- within(Thiemann2001, {
  case <- factor(case)
  series <- factor(series)
  treatment <- ifelse(treatment == "A", "baseline", "treatment")
  treatment <- factor(treatment)
})

save(Thiemann2001, file = "data/Thiemann2001.RData", compress = TRUE)
