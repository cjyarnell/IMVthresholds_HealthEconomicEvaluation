# Calculate the baseline weibull coefficients for the cohort

# First you must derive the cohort pdx_b using 
# the approach from the repository at doi:10.5281/zenodo.7314132

library(cmdstanr)
library(bayesplot)
library(tidyverse)
setwd("C:/Git/PAHRC/imvThresholdsCEA")

pdx_b <- readRDS("C:/Git/PAHRC/imvThresholdsTargetTrial/pdx_b.rds") %>%
  select(stay_id, imv_time, icu_discharge_time, time_to_discharge, time_to_death, eligibility_time,
         anchor_age, gender) %>%
  mutate(anchor_age = (anchor_age - 66.5)/16.1,
         icu_discharge_time = icu_discharge_time - eligibility_time,
         time_to_death = time_to_death - eligibility_time,
         time_to_discharge = time_to_discharge - eligibility_time) %>%
  # exclude patients who die before IMV on this ICU admission
  filter(is.na(time_to_death) | 
           (time_to_death > imv_time) |
           (time_to_death > icu_discharge_time))%>%
  select(-eligibility_time)

# 5 states
##############
# 1 - Oxygen therapy
# 2 - IMV
# 3 - ICU discharge
# 4 - Hospital discharge
# 5 - Death

# 1 -> 2, 3 (we rephrase 1->5 as 1->2 because everyone full code in the trial)
# 2 -> 3, 5
# 3 -> 4, 5
# 4 -> 5 (not from this model)

# so we estimate 7 weibull shape and scale parameters.
# 12, 13, 15, 23, 25, 34, 35

mf <- 28*24*60 # time conversion

pdx12 <- 
  pdx_b %>%
    group_by(stay_id) %>%
    mutate(time = min(c(imv_time, icu_discharge_time, time_to_death), na.rm = T),
           status = which.min(c(imv_time, icu_discharge_time, time_to_death)) %in% c(1,3)) %>%
  select(status, time, anchor_age, gender) %>%
  mutate(priortime = 0)

pdx13 <-  
  pdx_b %>%
  group_by(stay_id) %>%
  mutate(time = min(c(imv_time, icu_discharge_time, time_to_discharge), na.rm = T),
         status = as.integer(which.min(c(imv_time, icu_discharge_time, time_to_death)) == 2)) %>%
  select(status, time, anchor_age, gender) %>%
  mutate(priortime = 0) %>%
  filter(time > 0)

pdx23 <- 
  pdx_b %>%
    filter(stay_id %in% pdx12$stay_id[pdx12$status == 1]) %>%
    mutate(icu_discharge_time = icu_discharge_time - imv_time,
           time_to_death = time_to_death - imv_time) %>%
    group_by(stay_id) %>%
    mutate(time = min(c(icu_discharge_time, time_to_death), na.rm = T),
           status = which.min(c(icu_discharge_time, time_to_death)) == 1) %>%
  mutate(priortime = imv_time) %>%
  select(status, time, priortime, anchor_age, gender) %>%
  filter(time > 0)

pdx25 <- 
  pdx_b %>%
  filter(stay_id %in% pdx12$stay_id[pdx12$status == 1]) %>%
  mutate(icu_discharge_time = icu_discharge_time - imv_time,
         time_to_death = time_to_death - imv_time) %>%
  group_by(stay_id) %>%
  mutate(time = min(c(icu_discharge_time, time_to_death), na.rm = T),
         status = which.min(c(icu_discharge_time, time_to_death)) == 2) %>%
  mutate(priortime = imv_time) %>%
  select(status, time, priortime, anchor_age, gender) %>%
  filter(time > 0)

pdx34 <- 
  bind_rows(filter(pdx13,status == 1),
            filter(pdx23,status == 1),
            .id = "priorstate") %>%
    left_join(pdx_b, by = c("stay_id", "anchor_age", "gender")) %>%
    mutate(time = time + priortime) %>%
    select(stay_id, time, time_to_discharge, time_to_death, anchor_age, gender) %>%
    rename(priortime = time) %>%
    select(stay_id, priortime, time_to_discharge, time_to_death, anchor_age, gender) %>%
    mutate(time_to_discharge = time_to_discharge-priortime,
           time_to_death = time_to_death -priortime) %>%
  group_by(stay_id) %>%
  mutate(time = min(c(time_to_discharge, time_to_death), na.rm = T),
         status = is.na(time_to_death)) %>%
  select(status, time, priortime, anchor_age, gender) %>%
  filter(time > 0)

pdx35 <- 
  bind_rows(filter(pdx13,status == 1),
            filter(pdx23,status == 1),
            .id = "priorstate") %>%
  left_join(pdx_b, by = c("stay_id", "anchor_age", "gender")) %>%
  mutate(time = time + priortime) %>%
  select(stay_id, time, time_to_discharge, time_to_death, anchor_age, gender) %>%
  rename(priortime = time) %>%
  select(stay_id, priortime, time_to_discharge, time_to_death, anchor_age, gender) %>%
  mutate(time_to_discharge = time_to_discharge-priortime,
         time_to_death = time_to_death -priortime) %>%
  group_by(stay_id) %>%
  mutate(time = min(c(time_to_discharge, time_to_death), na.rm = T),
         status = !is.na(time_to_death)) %>%
  select(status, time, priortime, anchor_age, gender) %>%
  filter(time > 0)

#####################
# fit the models
wmod <- cmdstan_model("weibull.stan")

wfit <- function(pdx){
  standat <- 
    pdx %>%
    ungroup() %>%
    select(time, status, anchor_age, gender) %>%
    mutate(time = time/(7*24*60),
           gender = ifelse(gender == "F",1,0)) %>%
    filter(time> 0)
  
X_obs <- standat %>%
  filter(status == 1) %>%
  mutate(intercept = 1) %>%
  relocate(intercept, time, status, anchor_age, gender)

X_cen <- standat %>%
  filter(status == 0) %>%
  mutate(intercept = 1) %>%
  relocate(intercept, time, status, anchor_age, gender)

standat  <- list(P = ncol(X_obs)-2,
          N_obs = nrow(X_obs),
          N_cen = nrow(X_cen),
          X_obs = select(X_obs,-status,-time),
          X_cen = select(X_cen,-status,-time),
          y_obs = X_obs$time,
          y_cen = X_cen$time)

fit <- wmod$sample(data = standat,
            iter_warmup = 500,
            chains = 4, 
            iter_sampling = 250)

fitdraws <- fit$draws(variable = c("beta","logalpha"))

dimnames(fitdraws)$variable[-4] <- names(X_obs[c(-2,-3)])
fitdraws
}

process_fit <- function(wfit){
  data.frame(
    names = dimnames(wfit)$variable,
    mean = apply(wfit, 3, mean),
    sd = apply(wfit, 3, sd)
  )
}


wfit12 <- wfit(pdx12)
wfit13 <- wfit(pdx13)
wfit23 <- wfit(pdx23)
wfit25 <- wfit(pdx25)
wfit34 <- wfit(pdx34)
wfit35 <- wfit(pdx35)

saveRDS(list(wfit12, 
             wfit13,
             wfit23,
             wfit25,
             wfit34,
             wfit35),
        "weibulls.rds")

weibulls <- readRDS("weibulls.rds")

weibull_coefficients <- bind_rows(process_fit(wfit12),
          process_fit(wfit13),
          process_fit(wfit23),
          process_fit(wfit25),
          process_fit(wfit34),
          process_fit(wfit35),
          .id = "transition") %>%
  mutate(transition = factor(transition,
                             levels = 1:6,
                             labels = c("12","13",
                                        "23","25",
                                        "34","35")))

weibull_coefficients %>% 
  write_csv("weibull_coefficients.csv")


########################################################

# weibull for ventilator duration

vents_cea <- readRDS("C:/Git/PAHRC/imvThresholdsTargetTrial/vents_cea.rds")

stanvents <- 
  vents_cea %>% select(-stay_id,-oxygen_duration,-hospital_duration,-imv28) %>%
  filter(vent_duration > 0) %>%
  mutate(sex =ifelse(gender == "F",1,0)) %>%
  select(-gender) %>%
    mutate(anchor_age = (anchor_age - 65)/15) %>%
  mutate(intercept = 1) %>%
  relocate(vent_duration, intercept, death28, anchor_age, sex)

standat <- list(
  P = ncol(stanvents)-1,
  N_obs = nrow(stanvents),
  X_obs = select(stanvents, -vent_duration),
  y_obs = stanvents$vent_duration
)

wmod_nc <- cmdstan_model("weibull_no_censored.stan")

w_vents <- wmod_nc$sample(data = standat, 
                          iter_warmup = 750,
                          iter_sampling = 250,
                          adapt_delta = 0.85)

mcmc_trace(w_vents$draws())

vent_coefficients <- data.frame( names = c(names(stanvents)[-1],"logalpha"),
                                 mean = apply(w_vents$draws(variables = c("logalpha", "beta")),
                                              3, mean),
                                 sd = apply(w_vents$draws(variables = c("logalpha", "beta")),
                                              3, sd))
                                 
saveRDS(vent_coefficients, "vent_coefficients.rds")

stano2 <- 
  vents_cea %>% select(-stay_id,-vent_duration,-hospital_duration) %>%
  filter(oxygen_duration > 0) %>%
  mutate(sex =ifelse(gender == "F",1,0)) %>%
  select(-gender) %>%
  mutate(anchor_age = (anchor_age - 65)/15) %>%
  mutate(intercept = 1) %>%
  select(oxygen_duration, intercept, death28, imv28, anchor_age, sex)

standat <- list(
  P = ncol(stano2)-1,
  N_obs = nrow(stano2),
  X_obs = select(stano2, -oxygen_duration),
  y_obs = stano2$oxygen_duration
)

wmod_o2 <- cmdstan_model("weibull_no_censored.stan")

w_o2 <- wmod_nc$sample(data = standat, 
                          iter_warmup = 750,
                          iter_sampling = 250,
                          adapt_delta = 0.85)

mcmc_trace(w_o2$draws())

o2_coefficients <- data.frame( names = c(names(stano2)[-1],"logalpha"),
                                 mean = apply(w_o2$draws(variables = c("logalpha", "beta")),
                                              3, mean),
                                 sd = apply(w_o2$draws(variables = c("logalpha", "beta")),
                                            3, sd))

saveRDS(o2_coefficients, "o2_coefficients.rds")


stanhosp <- 
  vents_cea %>% select(-stay_id) %>%
  filter(hospital_duration > 0) %>%
  mutate(sex =ifelse(gender == "F",1,0)) %>%
  select(-gender) %>%
  mutate(anchor_age = (anchor_age - 65)/15) %>%
  mutate(intercept = 1) %>%
  select(hospital_duration, intercept, death28, imv28, anchor_age, sex)

standat <- list(
  P = ncol(stanhosp)-1,
  N_obs = nrow(stanhosp),
  X_obs = select(stanhosp, -hospital_duration),
  y_obs = stanhosp$hospital_duration
)

wmod_nc <- cmdstan_model("weibull_no_censored.stan")

w_hosp <- wmod_nc$sample(data = standat, 
                       iter_warmup = 750,
                       iter_sampling = 250,
                       adapt_delta = 0.85)

mcmc_trace(w_hosp$draws())

hosp_coefficients <- data.frame( names = c(names(stanhosp)[-1],"logalpha"),
                               mean = apply(w_hosp$draws(variables = c("logalpha", "beta")),
                                            3, mean),
                               sd = apply(w_hosp$draws(variables = c("logalpha", "beta")),
                                          3, sd))

saveRDS(hosp_coefficients, "hosp_coefficients.rds")

bind_rows(readRDS("vent_coefficients.rds"),
          readRDS("o2_coefficients.rds"),
          readRDS("hosp_coefficients.rds"),
          .id = "outcome") %>%
  mutate(rn = row_number(),
         outcome = factor(outcome, 
                          levels = 1:3,
                          labels = c("vent","oxygen", "hospital"))) %>%
  relocate(rn, mean, sd, names, outcome) %>%
  write.csv("weibull_coefficients.csv")
