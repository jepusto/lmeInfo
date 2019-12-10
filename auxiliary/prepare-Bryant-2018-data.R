# Read in Bryant 2018

Bryant18 <- read.csv("auxiliary/Bryant2018.csv", stringsAsFactors = FALSE)
str(Bryant18)

Bryant18 <- within(Bryant18, {
  Study_ID <- rep("Bryant2018", length(Bryant18$school))
  school <- factor(paste("school", school))
  case <- factor(studentID)
  outcome <- AC
  treatment <- ifelse(phase == "Baseline", "baseline", "treatment")
  treatment <- factor(treatment)
})

Bryant18 <- Bryant18[c("Study_ID", "school","case","treatment","session","session_trt","outcome")]

save(Bryant18, file = "data/Bryant2018.RData", compress = TRUE)
