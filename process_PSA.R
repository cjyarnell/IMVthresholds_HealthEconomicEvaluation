# Process the PSA results
library(tidyverse)
library(ggridges)
setwd("C:/Git/PAHRC/imvThresholdsCEA")

WTP = 100000

source("psa_functions.r")
psa <- read_csv("PSA.csv",skip = 2)


# mortality and imv 

mortality_columns <- which(grepl("t_survive", names(psa)))
imv_columns <- 27:30
ltd_columns <- which(grepl("Longterm disability", names(psa)))
le_columns <-  which(grepl("Life expectancy", names(psa)))
imvdur_columns <-  which(grepl("t_imvduration", names(psa)))
o2dur_columns <-  which(grepl("t_oxygenduration", names(psa)))
hospdur_columns <- which(grepl("t_hospitalduration", names(psa)))
nmb_columns <-  which(grepl("NMB", names(psa)))
qaly_columns <-  which(grepl("Effectiveness", names(psa)))
cost_columns <-  2:5

psa[,c(mortality_columns, imv_columns)] %>%
  pivot_longer(cols = everything()) %>%
  mutate(outcome = ifelse(grepl("survive", name), "survival","imv"),
         strategy = factor(rep(1:4, times = 2000),
                           levels = 1:4,
                           labels = c("Usual care",
                                      "SF < 88",
                                      "SF < 98",
                                      "SF < 110"))) %>%
  select(outcome, strategy, value) %>%
  ggplot(aes(x = value, y = strategy)) +
  geom_density_ridges() +
  facet_grid(outcome~.) +
  theme_minimal()

tbl2 <-
  psa[,c(mortality_columns, imv_columns,
       ltd_columns, le_columns, 
       imvdur_columns, o2dur_columns,
       hospdur_columns,cost_columns,
       qaly_columns)] %>%
  pivot_longer(cols = everything()) %>%
  mutate(outcome = case_when(
    grepl("oxygenduration", name) ~ "ICU, non-IMV duration",
    grepl("imvduration", name) ~ "IMV duration",
    grepl("hospital", name) ~ "Ward duration",
    grepl("imv", name) ~ "Invasive ventilation",
    grepl("survive", name) ~ "Survival",
    grepl("Longterm", name) ~ "Longterm disability",
    grepl("Life", name) ~ "Life expectancy",
    grepl("Cost", name) ~ "Cost",
    grepl("Effectiveness", name) ~ "QALYs"
    ),
         strategy = factor(rep(1:4, times = 9000),
                           levels = 1:4,
                           labels = c("Usual care",
                                      "SF < 88",
                                      "SF < 98",
                                      "SF < 110"))) %>%
  mutate(outcome = factor(outcome,
                          levels = c("ICU, non-IMV duration",
                                     "IMV duration",
                                     "Ward duration",
                                     "Invasive ventilation",
                                     "Survival",
                                     "Longterm disability",
                                     "Life expectancy",
                                     "QALYs",
                                     "Cost"),
                          ordered = T)) %>%
  select(outcome, strategy, value) %>%
  group_by(outcome, strategy) %>%
  mutate(iter = 1:1000,
         value = ifelse(outcome == "Cost", value/1000, value)) %>%
  mutate(value = ifelse(outcome %in% c("Invasive ventilation",
                                       "Survival",
                                       "Longterm disability"),
                        value*100, value)) %>%
  pivot_wider(names_from = outcome,
              values_from = value) %>%
  mutate(NMB = QALYs*WTP/1000 - Cost,
         `IMV duration` = 100*`IMV duration`/`Invasive ventilation`) %>%
  group_by(iter) %>%
  mutate(Optimal = ifelse(strategy == strategy[which.max(NMB)], 1, 0))  %>%
  ungroup() %>%
  select(-iter) %>%
  pivot_longer(-strategy,
               names_to = "outcome") %>%  
  group_by(strategy, outcome) %>%
    mutate(outcome = factor(outcome,
                            levels = c("ICU, non-IMV duration",
                                       "IMV duration",
                                       "Ward duration",
                                       "Invasive ventilation",
                                       "Survival",
                                       "Longterm disability",
                                       "Life expectancy",
                                       "QALYs",
                                       "Cost",
                                       "NMB",
                                       "Optimal"),
                            ordered = T)) %>% 
  summarise(output = mean95cri(value)) %>%
    pivot_wider(names_from = strategy,
                values_from = output)
  
icur <-  
  psa[,c(cost_columns,
         qaly_columns)] %>%
    pivot_longer(cols = everything()) %>%
    mutate(outcome = case_when(
      grepl("Cost", name) ~ "Cost",
      grepl("Effectiveness", name) ~ "QALYs"
    ),
    strategy = factor(rep(1:4, times = 2000),
                      levels = 1:4,
                      labels = c("Usual care",
                                 "SF < 88",
                                 "SF < 98",
                                 "SF < 110"))) %>%
    mutate(outcome = factor(outcome,
                            levels = c("ICU, non-IMV duration",
                                       "IMV duration",
                                       "Ward duration",
                                       "Invasive ventilation",
                                       "Survival",
                                       "Longterm disability",
                                       "Life expectancy",
                                       "QALYs",
                                       "Cost"),
                            ordered = T)) %>%
    select(outcome, strategy, value) %>%    
    group_by(strategy, outcome) %>%
    summarise(mean = mean(value)) %>%
    pivot_wider(names_from = outcome,
                values_from = mean) %>%
    ungroup() %>%
    mutate(whichmincost = 1) %>%
    mutate(ICUR = ifelse(Cost == Cost[whichmincost],
                         0,
                         (Cost - Cost[whichmincost])/
                           (QALYs - QALYs[whichmincost]))) %>%
    pivot_longer(-strategy,
                 names_to = "outcome") %>%
    mutate(value = prettyNum(signif(value,3), big.mark = ",")) %>%
    pivot_wider(names_from = strategy) %>%
    filter(outcome == "ICUR")

bind_rows(tbl2, icur) %>% write_csv("table2.csv")

# CEAC

wtp_thresholds = seq(from = 0, to = 200000, by = 1000)

CEAC(psa, cost_columns, qaly_columns, wtp_thresholds) 

ggsave("figures/FigureCEAC.svg", width = 7, height = 5)

# joint CEAC

psa <- bind_rows(MIMIC = read_csv("PSA.csv",skip = 2),
                 AMDS  = read_csv("PSA_amd.csv",skip = 2),
                 .id = "source")

qaly_columns <-  which(grepl("Effectiveness", names(psa)))
cost_columns <-  3:6

jointCEAC(psa, cost_columns, qaly_columns, wtp_thresholds) 

ggsave("figures/Figure1jointCEAC.svg", height = 8, width = 7)


# EVPPI
psa[,-c(19:22)] %>%
savi_parameters() %>%
  write_csv("savi_parameters.csv")

psa[,-c(19:22)] %>%
savi_costs() %>%
  write_csv("savi_costs.csv")

psa[,-c(19:22)] %>%
savi_effectiveness() %>%
  write_csv("savi_effectiveness.csv")

# use SAVI web interface and save as EVPPI_mimic.csv

evppi_df <- read.csv("EVPPI_mimic.csv") %>%
  select(X, Per.Person.EVPPI...., Standard.Error) %>%
  rename(Study = X,
         EVPPI = Per.Person.EVPPI....,
         se = Standard.Error) %>%
  mutate(EVPPI = signif(EVPPI, 3),
         lb = EVPPI - 1.96*se,
         ub = EVPPI + 1.96*se) %>%
  bind_rows(data.frame(Study = "Total",
                       EVPPI = 67524,
                       se = 0,
                       lb = NA,
                       ub = NA))%>% 
  arrange(desc(EVPPI)) %>%
  mutate(Study = factor(Study,
                        labels = rev(c("Total EVPI",
                                   "Hospitalization costing study",
                                   "Observational study",
                                   "Randomized trial"
                                   ))),
         label = prettyNum(EVPPI, big.mark = ","))
evppi_df %>%
  ggplot(aes(x = fct_reorder(Study,EVPPI), y = EVPPI, fill = Study)) +
  geom_col(width = 0.8) +
  geom_errorbar(aes(ymin = lb, ymax = ub),
                width = 0.4,
                color = "gray50") +
  theme_minimal() +
  coord_flip() +
    labs(x = "",
       y = "Per-person estimated value of partial perfect information (EVPPI)") +
  scale_fill_manual(values = c(c_dark, c_mid, c_light, c_blue),
                    guide = "none") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        panel.grid.minor = element_blank()) +
  geom_label(aes(label = label), fill = "grey90", y = 500, hjust = 0)
  
ggsave("figures/Figure2EVPPI.svg", width = 8, height = 4)

# population EVPPI



EVSI <- 

I_t <- 5000

popEVSI <- 0

r = 0.015

for (i in 0:9){
  print((I_t)/((1+r)^i))
  popEVSI <- popEVSI + EVSI*(I_t)/((1+r)^i)
}

popEVSI

# population EVPPI, sensitivity


EVSI <- 71366.17

I_t <- 5000

popEVSI <- 0

r = 0.015

for (i in 0:9){
  print((I_t)/((1+r)^i))
  popEVSI <- popEVSI + EVSI*(I_t)/((1+r)^i)
}

popEVSI
