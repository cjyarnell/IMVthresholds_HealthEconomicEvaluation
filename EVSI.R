# Simulating datasets for EVSI calculation

# set working directory to source file location

source("C:/Git/PAHRC/imvThresholdsCEA/EVSI_functions.R")
setwd("C:/Git/PAHRC/imvThresholdsCEA")

c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")
c_blue = "cornflower blue"


library(mgcv)
library(tidyverse)

# MIMIC 
PSA <- read_csv("PSA.csv", skip = 2)[,-c(19:22)]

sample_sizes <- c(100, 200, 400, 800, 1600, 3200)

EVSI <- data.frame(N = sample_sizes,
                   source = "MIMIC-IV",
                   EVSI = NA)

for (n in 1:nrow(EVSI)){
  EVSI$EVSI[n] = evsi_rct(PSA, EVSI$N[n])
  saveRDS(EVSI, "evsi.rds")
}

# figure
readRDS("evsi.rds") %>%
  mutate(source = factor(source, levels = c("MIMIC-IV","AmsterdamUMCdb"),
                         labels = c("Primary Analysis", "Sensitivity Analysis"))) %>%
  mutate(popEVSI = pop_evsi(EVSI = EVSI)) %>%
  ggplot(aes(x = N, y = popEVSI/1000000)) +
  geom_col(fill = "gray50", alpha = 0.5) +
#  geom_col(aes(fill = source), position = "dodge", alpha = 0.9) +
#  geom_line(aes(color = source), size = 2) +
#  geom_point(color = "white", size = 2) +
#  geom_point(color = "black", size = 1.5) +
#  scale_color_manual(values = c(c_dark, c_blue)) +
 # scale_fill_manual(values = c("grey75", "black"),
#                    name = "Analysis") +
  theme_minimal() +
  scale_x_continuous(trans = "log2", breaks = sample_sizes) +
  scale_y_continuous(expand = c(0,0), breaks = c(2, 4, 6)) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
  labs(y = "EVSI per patient (billion CAD)",
       x = "Sample size of 4-arm trial")

ggsave(filename = "figures/EVSI.svg", width = 6, height = 6, units = "in")

# MIMIC 
PSA_s <- read_csv("PSA_sensitivity.csv", skip = 2)

sample_sizes <- c(100, 200, 400, 800, 1600, 3200)

EVSI_sens <- data.frame(N = sample_sizes,
                   source = "MIMIC-IV",
                   EVSI = NA)

for (n in 1:nrow(EVSI_sens)){
  EVSI_sens$EVSI[n] = evsi_rct_sens(PSA_s, EVSI_sens$N[n])
  saveRDS(EVSI_sens, "evsi_sens.rds")
}

# figure
readRDS("evsi_sens.rds") %>%
  mutate(source = factor(source, levels = c("MIMIC-IV","AmsterdamUMCdb"),
                         labels = c("Primary Analysis", "Sensitivity Analysis"))) %>%
  ggplot(aes(x = N, y = EVSI/1000)) +
  geom_col(fill = "gray50", alpha = 0.5) +
  #  geom_col(aes(fill = source), position = "dodge", alpha = 0.9) +
  #  geom_line(aes(color = source), size = 2) +
  #  geom_point(color = "white", size = 2) +
  #  geom_point(color = "black", size = 1.5) +
  #  scale_color_manual(values = c(c_dark, c_blue)) +
  # scale_fill_manual(values = c("grey75", "black"),
  #                    name = "Analysis") +
  facet_wrap(source~., nrow = 2, scales = "free") +
  theme_minimal() +
  scale_x_continuous(trans = "log2", breaks = sample_sizes) +
  scale_y_continuous(expand = c(0,0), breaks = c(0, 50, 100)) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
  labs(y = "EVSI per patient (1000's CAD)",
       x = "Sample size of 4-arm trial")

ggsave(filename = "figures/EVSI.svg", width = 6, height = 6, units = "in")


bind_rows(readRDS("evsi.rds"), 
          readRDS("evsi_sens.rds"),
          .id = "Scenario") %>%
  mutate(Scenario = factor(Scenario,
                           levels = c(1,2),
                           labels = c("Scenario 1: a 4-arm trial showing early invasive ventilation is beneficial",
                                      "Scenario 2: a 2-arm trial showing late invasive ventilation is safe")))  %>%
  mutate(popEVSI = pop_evsi(EVSI = EVSI)) %>%
  ggplot(aes(x = N, y = popEVSI/1e9)) +
  geom_col(aes(fill = Scenario)) +
  #  geom_col(aes(fill = source), position = "dodge", alpha = 0.9) +
  #  geom_line(aes(color = source), size = 2) +
  #  geom_point(color = "white", size = 2) +
  #  geom_point(color = "black", size = 1.5) +
  #  scale_color_manual(values = c(c_dark, c_blue)) +
   scale_fill_manual(values = c(c_dark, "orange2"),
                      name = "Scenario",
                     labels = 1:2,
                     guide = "none") +
  facet_wrap(Scenario~., nrow = 2, scales = "free_x") +
  theme_minimal() +
  scale_x_continuous(trans = "log2", breaks = sample_sizes) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  labs(y = "Population EVSI (CAD, billions)",
       x = "Total sample size of trial")

ggsave("figures/evsi_both.svg", units = "in", width = 6, height = 6)


