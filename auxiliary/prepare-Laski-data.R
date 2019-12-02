Laski <- read.csv("auxiliary/Laski.csv", stringsAsFactors = FALSE)
str(Laski)
Laski <- within(Laski, {
  case <- factor(paste("Case",case))
  treatment <- factor(treatment, levels = 0:1, labels = c("baseline","treatment"))
})

save(Laski, file = "data/Laski.RData", compress = TRUE)
