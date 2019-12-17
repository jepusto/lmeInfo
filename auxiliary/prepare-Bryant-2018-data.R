# Read in Bryant 2018

Bryant2018 <- read.csv("auxiliary/Bryant2018.csv", stringsAsFactors = FALSE)
str(Bryant2018)

Bryant2018 <- within(Bryant2018, {
  Study_ID <- rep("Bryant2018", length(Bryant2018$school))
  school <- factor(paste("school", school))
  case <- factor(studentID)
  outcome <- AC
  treatment <- ifelse(phase == "Baseline", "baseline", "treatment")
  treatment <- factor(treatment)
})

Bryant2018 <- Bryant2018[c("Study_ID", "school","case","treatment","session","session_trt","outcome")]

# time-point constants
A <- 4
B <- 21

# center at follow-up time
Center <- B
Bryant2018$session_c <-Bryant2018$session - Center

save(Bryant2018, file = "data/Bryant2018.RData", compress = TRUE)
