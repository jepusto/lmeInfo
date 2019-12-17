# Read in Bryant 2016

Bryant2016 <- read.csv("auxiliary/Bryant2016.csv", stringsAsFactors = FALSE)
str(Bryant2016)

Bryant2016 <- within(Bryant2016, {
  school <- factor(school)
  case <- factor(case)
  treatment <- ifelse(treatment == "A", "baseline", "treatment")
  treatment <- factor(treatment)
})

Bryant2016 <- Bryant2016 %>%
  arrange(school, case, session) %>%
  group_by(school, case) %>%
  mutate(
    trt_time = pmax(0, session - max(session[treatment == "A"])))

# time-point constants
A <- 4
B <- 9

# center at follow-up time
Center <- B
Bryant2016$session_c <- Bryant2016$session - Center


save(Bryant2016, file = "data/Bryant2016.RData", compress = TRUE)
