# psa functions

# colours
c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")
c_blue = "cornflower blue"

mean95cri <- function(x){paste0(signif(mean(x), 3),
                                " (",
                                signif(quantile(x,0.025), 3),
                                " to ",
                                signif(quantile(x,0.975), 3),
                                ")")
}

savi_parameters <- function(psa){
  params <- c(51:55,
              60:84,
              88:116)
  psa[,params]
}

savi_parameters_sensitivity <- function(psa){
  params <- c(27:31,
              36:60,
              64:92)
  psa[,params]
}


savi_costs <- function(psa){
  temp <- psa[,2:5]
  names(temp) <- gsub("Cost\n","", names(temp))
  temp
}

savi_effectiveness <- function(psa){
  temp <- psa[,6:9]
  names(temp) <- gsub("Effectiveness\n","", names(temp))
  temp
}

savi_costs_sensitivity <- function(psa){
  temp <- psa[,2:3]
  names(temp) <- gsub("Cost\n","", names(temp))
  temp
}

savi_effectiveness_sensitivity <- function(psa){
  temp <- psa[,4:5]
  names(temp) <- gsub("Effectiveness\n","", names(temp))
  temp
}


CEAC <- function(psa, cost_columns,
                 qaly_columns,
                 wtp_thresholds){

    temp <- psa[,c(1,cost_columns,
                 qaly_columns)] %>%
    pivot_longer(-Iteration) %>%
    mutate(
      outcome = case_when(
        grepl("Cost", name) ~ "Cost",
        grepl("Effectiveness", name) ~ "QALYs"
      ),
      strategy = factor(rep(1:4, times = 2000),
                        levels = 1:4,
                        labels = c("Usual care",
                                   "SF < 88",
                                   "SF < 98",
                                   "SF < 110"))) %>%  
    select(Iteration, outcome, strategy, value) %>%
    group_by(Iteration) %>%
    pivot_wider(names_from = outcome,
                values_from = value)
    
    ceaclist <- list()

    for(i in 1:length(wtp_thresholds)){
      ceaclist[[i]] <-
      temp %>%
        mutate(wtp = wtp_thresholds[i]) %>%
        mutate(nmb = QALYs*wtp - Cost) %>%
        group_by(Iteration) %>%
        mutate(optimal = ifelse(strategy == strategy[which.max(nmb)],1,0)) %>%
        ungroup() %>%
        group_by(strategy) %>%
        summarise(prob = mean(optimal),
                  wtp = mean(wtp))
    }
    
    ceacdf <- bind_rows(ceaclist) %>%
      mutate(strategy = factor(strategy, 
                               levels = c("Usual care",
                                          "SF < 110",
                                          "SF < 98",
                                          "SF < 88"),
                               ordered = T))
    
    ceacdf %>%
      rename(Strategy = strategy) %>%
      ggplot(aes(x = wtp/1000, y = prob, color = Strategy)) +
      geom_line(size = 1.75) +
      geom_point(data = filter(ceacdf, wtp %% 20000 == 0), 
                 color = "white", size = 2.25) +
      geom_point(data = filter(ceacdf, wtp %% 20000 == 0), 
                 color = "black", size = 1.3) +
      theme_minimal() +
      lims(y = c(0,1)) +
      labs(y = "Probability of having highest net monetary benefit",
           x = "Willingness-to-pay (1000's CAD)") +
      scale_color_manual(values = c(c_blue, c_dark, c_mid, c_light)) +
      theme(panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
            axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
}

CEAC_sens <- function(psa, cost_columns,
                 qaly_columns,
                 wtp_thresholds){
  
  temp <- psa[,c(1,cost_columns,
                 qaly_columns)] %>%
    pivot_longer(-Iteration) %>%
    mutate(
      outcome = case_when(
        grepl("Cost", name) ~ "Cost",
        grepl("Effectiveness", name) ~ "QALYs"
      ),
      strategy = factor(rep(1:2, times = 2000),
                        levels = 1:2,
                        labels = c("Usual care",
                                   "Hypothetical threshold"))) %>%  
    select(Iteration, outcome, strategy, value) %>%
    group_by(Iteration) %>%
    pivot_wider(names_from = outcome,
                values_from = value)
  
  ceaclist <- list()
  
  for(i in 1:length(wtp_thresholds)){
    ceaclist[[i]] <-
      temp %>%
      mutate(wtp = wtp_thresholds[i]) %>%
      mutate(nmb = QALYs*wtp - Cost) %>%
      group_by(Iteration) %>%
      mutate(optimal = ifelse(strategy == strategy[which.max(nmb)],1,0)) %>%
      ungroup() %>%
      group_by(strategy) %>%
      summarise(prob = mean(optimal),
                wtp = mean(wtp))
  }
  
  ceacdf <- bind_rows(ceaclist) %>%
    mutate(strategy = factor(strategy, 
                             levels = c("Usual care",
                                        "Hypothetical threshold"),
                             ordered = T))
  
  ceacdf %>%
    rename(Strategy = strategy) %>%
    ggplot(aes(x = wtp/1000, y = prob, color = Strategy)) +
    geom_line(size = 1.75) +
    geom_point(data = filter(ceacdf, wtp %% 20000 == 0), 
               color = "white", size = 2.25) +
    geom_point(data = filter(ceacdf, wtp %% 20000 == 0), 
               color = "black", size = 1.3) +
    theme_minimal() +
    lims(y = c(0,1)) +
    labs(y = "Probability of having highest net monetary benefit",
         x = "Willingness-to-pay (1000's CAD)") +
    scale_color_manual(values = c(c_blue, "orange2")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
}


jointCEAC <- function(psa, cost_columns,
                 qaly_columns,
                 wtp_thresholds){
  
  temp <- 
    psa[,c(1:2,cost_columns,
                 qaly_columns)] %>%
    pivot_longer(-(source:Iteration)) %>%
    mutate(
      outcome = case_when(
        grepl("Cost", name) ~ "Cost",
        grepl("Effectiveness", name) ~ "QALYs"
      ),
      strategy = factor(rep(1:4, times = 4000),
                        levels = 1:4,
                        labels = c("Usual care",
                                   "SF < 88",
                                   "SF < 98",
                                   "SF < 110"))) %>%  
    select(source,Iteration, outcome, strategy, value) %>%
    group_by(source,Iteration) %>%
    pivot_wider(names_from = outcome,
                values_from = value)
  
  ceaclist <- list()
  
  for(i in 1:length(wtp_thresholds)){
    ceaclist[[i]] <-
      temp %>%
      mutate(wtp = wtp_thresholds[i]) %>%
      mutate(nmb = QALYs*wtp - Cost) %>%
      group_by(source,Iteration) %>%
      mutate(optimal = ifelse(strategy == strategy[which.max(nmb)],1,0)) %>%
      ungroup() %>%
      group_by(source,strategy) %>%
      summarise(prob = mean(optimal),
                wtp = mean(wtp))
  }
  
  ceacdf <- bind_rows(ceaclist) %>%
    mutate(strategy = factor(strategy, 
                             levels = c("Usual care",
                                        "SF < 110",
                                        "SF < 98",
                                        "SF < 88"),
                             ordered = T),
           source = factor(source,
                           levels = c("MIMIC","AMDS"),
                           labels = c("Primary Analysis","Sensitivity Analysis")))
  
  ceacdf %>%
    rename(Strategy = strategy) %>%
    ggplot(aes(x = wtp/1000, y = prob, color = Strategy)) +
    geom_line(size = 1.75) +
    geom_point(data = filter(ceacdf, wtp %% 20000 == 0), 
               color = "white", size = 2.25) +
    geom_point(data = filter(ceacdf, wtp %% 20000 == 0), 
               color = "black", size = 1.3) +
    theme_minimal() +
    lims(y = c(0,1)) +
    labs(y = "Probability of having highest net monetary benefit",
         x = "Willingness-to-pay (1000's CAD)") +
    scale_color_manual(values = c(c_blue, c_light, c_mid, c_dark)) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          legend.position = "right") +
    facet_wrap(source~., nrow = 2, scales = "free")
}
