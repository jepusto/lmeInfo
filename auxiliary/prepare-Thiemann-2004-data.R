# Read in Thiemann 2004

Thiemann2004 <- read.csv("auxiliary/Thiemann2004.csv", stringsAsFactors = FALSE)
str(Thiemann2004)

Thiemann2004 <- within(Thiemann2004, {
  case <- factor(case)
  series <- factor(series)
  treatment <- ifelse(treatment == "A", "baseline", "treatment")
  treatment <- factor(treatment)
})

Thiemann2004 <- Thiemann2004 %>%
  arrange(case, series, time) %>%
  group_by(case, series) %>%
  mutate(
    trt_time = pmax(0, time - max(time[treatment == "A"])))

# time-point constants
A <- 10
B <- 33

# center at follow-up time
Center <- B
Thiemann2004$time_c <- Thiemann2004$time - Center

save(Thiemann2004, file = "data/Thiemann2004.RData", compress = TRUE)
