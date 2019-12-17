# Read in Thiemann 2001

Thiemann2001 <- read.csv("auxiliary/Thiemann2001.csv", stringsAsFactors = FALSE)
str(Thiemann2001)

Thiemann2001 <- within(Thiemann2001, {
  case <- factor(case)
  series <- factor(series)
  treatment <- ifelse(treatment == "A", "baseline", "treatment")
  treatment <- factor(treatment)
})

Thiemann2001 <- Thiemann2001 %>%
  arrange(case, series, time) %>%
  group_by(case, series) %>%
  mutate(trt_time = pmax(0, time - max(time[treatment == "A"])))

# time-point constants
A <- 8
B <- 30

# center at follow-up time
Center <- B
Thiemann2001$time_c <- Thiemann2001$time - Center

save(Thiemann2001, file = "data/Thiemann2001.RData", compress = TRUE)
