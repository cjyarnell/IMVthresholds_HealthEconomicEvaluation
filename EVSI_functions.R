# EVSI functions

total_imv <- function(survival, imv_surv, imv_death){
  survival*imv_surv + (1-survival)*imv_death
}

odds <- function(p){p/(1-p)}

OR <- function(p,q){odds(p)/odds(q)}

simulate_survival_imv <- function(N,p){
  temp1 <- rbinom(1, size = N, prob = unlist(p[2]))
  temp2 <- rbinom(1, size = sum(temp1), prob = unlist(p[1]))
  temp3 <- rbinom(1, size = N-sum(temp1), prob = unlist(p[3]))
  c(temp1, temp2, temp3)
}

evsi_gam_rct <- function(df){
  gam(NMB~ s(x1,x2,x3,x4,x5,x6), data = df)
}

evsi_gam_rct_sens <- function(df){
  gam(NMB~ s(x1,x2), data = df)
}

evsi_rct_simple <- function(PSA, N){
  parameters_rct <- c(
    38:39,91 # usual care
    ,44:45,94 # SF 110
    #  ,68:72 # imv duration
    #  ,73:78 # oxygen duration
    #  ,79:84 # hospital duration
    #  ,88:92 # longterm disability
  )
  
  parameters_rct <- data.frame(PSA_names = names(PSA)[parameters_rct],
                               PSA_numbers = parameters_rct)
  
  
  
  ###################
  
  
  # observations from simulated trial
  simdata_obs_name <- c("logOR_survival_sf110",
                        "logOR_imv_sf110")
  
  simdata_obsdim <- length(simdata_obs_name)
  
  # make array
  simdata <- array(dim = c(nrow(PSA),simdata_obsdim))
  
  
  for (s in 1:nrow(PSA)){
    theta <- PSA[s, parameters_rct$PSA_numbers]
    simdata[s,] <- simulate_rct_sens(N, theta)
  }
  
  # INMB for each iteration 
  # WTP 100000
  
  WTP = 100000
  
  # NMB calculations
  nmb_uc <- PSA[,4]*WTP - PSA[,2]
  nmb_sf88 <- PSA[,5]*WTP - PSA[,3]
  
  df_uc <- cbind(nmb_uc, simdata)
  df_sf88 <- cbind(nmb_sf88, simdata)
  
  xs <- rep("",simdata_obsdim) 
  for (n in 1:simdata_obsdim){
    xs[n] <- paste0("x",n)
  }
  
  names(df_uc) <- c("NMB",xs)
  names(df_sf88) <- c("NMB",xs)
  
  fit_uc <- evsi_gam_rct_sens(df_uc)
  fit_sf88 <- evsi_gam_rct_sens(df_sf88)
  
  # check for normally-distributed residuals
  #plot(density(fit_uc$residuals))
  #plot(density(fit_sf88$residuals))
  #plot(density(fit_sf98$residuals))
  #plot(density(fit_sf110$residuals))
  
  g_k <- data.frame(iter = 1:length(fit_uc$fitted.values),
                    uc = fit_uc$fitted.values,
                    sf88  = fit_sf88$fitted.values) %>%
    mutate(max_of_means = max(mean(uc),
                              mean(sf88))) %>%
    group_by(iter) %>%
    mutate(max_per_iter = max(uc, sf88))
  
  EVSI <- mean(g_k$max_per_iter) - mean(g_k$max_of_means)
  
  EVSI
  
}

evsi_rct_sens <- function(PSA, N){
  parameters_rct <- c(
     36:37,89 # usual care
    ,38:39,90 # SF 88
    #  ,68:72 # imv duration
    #  ,73:78 # oxygen duration
    #  ,79:84 # hospital duration
    #  ,88:92 # longterm disability
  )
  
  parameters_rct <- data.frame(PSA_names = names(PSA)[parameters_rct],
                               PSA_numbers = parameters_rct)
  
  
  
  ###################
  
  
  # observations from simulated trial
  simdata_obs_name <- c("logOR_survival_sf88",
                        "logOR_imv_sf88")
  
  simdata_obsdim <- length(simdata_obs_name)
  
  # make array
  simdata <- array(dim = c(nrow(PSA),simdata_obsdim))
  
  
  for (s in 1:nrow(PSA)){
    theta <- PSA[s, parameters_rct$PSA_numbers]
    simdata[s,] <- simulate_rct_sens(N, theta)
  }
  
  # INMB for each iteration 
  # WTP 100000
  
  WTP = 100000
  
  # NMB calculations
  nmb_uc <- PSA[,4]*WTP - PSA[,2]
  nmb_sf88 <- PSA[,5]*WTP - PSA[,3]
  
  df_uc <- cbind(nmb_uc, simdata)
  df_sf88 <- cbind(nmb_sf88, simdata)
  
  xs <- rep("",simdata_obsdim) 
  for (n in 1:simdata_obsdim){
    xs[n] <- paste0("x",n)
  }
  
  names(df_uc) <- c("NMB",xs)
  names(df_sf88) <- c("NMB",xs)
  
  fit_uc <- evsi_gam_rct_sens(df_uc)
  fit_sf88 <- evsi_gam_rct_sens(df_sf88)

  # check for normally-distributed residuals
  #plot(density(fit_uc$residuals))
  #plot(density(fit_sf88$residuals))
  #plot(density(fit_sf98$residuals))
  #plot(density(fit_sf110$residuals))
  
  g_k <- data.frame(iter = 1:length(fit_uc$fitted.values),
                    uc = fit_uc$fitted.values,
                    sf88  = fit_sf88$fitted.values) %>%
    mutate(max_of_means = max(mean(uc),
                              mean(sf88))) %>%
    group_by(iter) %>%
    mutate(max_per_iter = max(uc, sf88))
  
  EVSI <- mean(g_k$max_per_iter) - mean(g_k$max_of_means)
  
  EVSI
  
}

simulate_rct_sens <- function(N,theta){
  vars <- c(
    "Survival"
    ,"IMV_Surv"
    ,"IMV_Death"
    #   ,"IMV_duration"
    #    ,"Oxygen_duration"
    #    ,"Hospital_duration"
    #    ,"LTD"
  )
  
  assignment <- sample(1:2, size = N, replace = T)
  N_uc <- sum(assignment == 1)
  N_sf88 <- sum(assignment == 2)

  x <- data.frame(outcome = rep(vars,2),
                  observed = c(
                    simulate_survival_imv(N_uc, p = theta[1:3])/N_uc,
                    simulate_survival_imv(N_sf88, p = theta[4:6])/N_sf88)
  )
  
  imv_obs <- c(total_imv(x$observed[1],x$observed[2], x$observed[3]),
               total_imv(x$observed[4],x$observed[5], x$observed[6]))
  
  temp <- log(c(logOR_survival_sf88 = OR(x$observed[4], x$observed[1]),
                logOR_imv_sf88 = OR(imv_obs[2], imv_obs[1])))
  
  temp[temp > 5] = 5
  temp[temp < -5] = -5
  temp
  }

pop_evsi <- function(I_t = 5000, EVSI, duration = 10){
  I_t <- 5000
  
  popEVSI <- 0
  
  r = 0.015
  
  for (i in 0:(duration-1)){
    print((I_t)/((1+r)^i))
    popEVSI <- popEVSI + EVSI*(I_t)/((1+r)^i)
  }
  popEVSI
}

