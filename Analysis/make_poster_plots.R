# Load necessary packages

# Read and process data
library(readr)
library(magrittr)
library(dplyr)
library(data.table)

# Plotting
library(ggplot2)
library(lattice)
library(gtable)
library(gridExtra)
library(gridGraphics)
library(ggdist)
library(ggsignif)
library(ggpubr)

# Linear models
library(afex)
library(car)
library(codingMatrices)
library(broom)
library(modelr)

# Mixed models
library(emmeans)
library(lme4)
library(performance)
library(ggsignif)
library(rtdists)

# Bayesian analyses
library(bayestestR)
library(BayesFactor)
library(brms)

# Save defaults
graphical_defaults <- par()
options_defaults <- options() 

#==========================================================================================#
#======================================= 0. Prepare data ===================================
#==========================================================================================#
# Functions for raincloud plots
#source("./Analysis/R_rainclouds.R")

source('./Analysis/Prepare_data_exp0_v2_6.R')

RT_by_condition <- ggplot(data = summaryData2 %>% 
                            group_by(ID, foreperiod, condition) %>% 
                            summarise(meanRT = mean(meanRT)),
                          aes(x = foreperiod,
                              y = meanRT,
                              color = condition)) +
  geom_jitter(height = 0, width = 0.15, size = 3.5, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 4.3, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 4.1, width = 0.1, geom = "errorbar") + 
  labs(title = paste("Experiment 2 (n = ", length(unique(summaryData2$ID)), ")", sep = ""),
       x = "FP (s)",
       y = "",
       color = "Condition") +
  theme(plot.title = element_text(size = rel(2.8), hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = rel(2.8)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 2.75, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(size = rel(2.5)),
        legend.title = element_text(size = rel(2.6)),
        legend.text = element_text(size = rel(2.4)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(5.5, 5.5, 5.5, 1),
        legend.direction = "horizontal",
        legend.position = "bottom") +
  scale_color_manual(values = c("orange", "blue"), labels = c("External", "Action"))

# Extract legend
cond_legend <- as_ggplot(get_legend(RT_by_condition))


ggsave("G:/My Drive/Post-doc/Eventos/TRF-3/Poster/RT_by_condition_exp2.tiff",
                RT_by_condition + theme(legend.position = "none"),
                width = 25,
                height = 16.66,
                units = "cm")

ggsave("G:/My Drive/Post-doc/Eventos/TRF-3/Poster/RT_by_condition_legend.tiff",
       cond_legend,
       width = 25,
       height = 16.66,
       units = "cm")

facet_labels <- c(expression("FP"[n-1]*"= 1.0"), expression("FP"[n-1]*"= 1.6"),
                  expression("FP"[n-1]*"= 2.2"), expression("FP"[n-1]*"= 2.8"))

seqEff_by_oneback <- ggplot(data = summaryData2 %>%
                              group_by(ID, foreperiod, condition, oneBackFP) %>%
                              summarise(meanRT = mean(meanRT)),
                            aes(x = foreperiod,
                                y = meanRT,
                                color=condition)) +
  #geom_violin(aes(fill = condition), alpha = 0.5) +
  geom_jitter(height = 0, width = 0.30, size = 3.5, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 4.3, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 4.1, width = 0.1, geom = "errorbar") + 
  labs(title = "Experiment 2",
       x = expression("FP"[n]*" (s)"),
       y = "",#"Mean RT (s)",
       color = "Condition") +
  ylim(min(summaryData2$meanRT), max(summaryData2$meanRT)) +
  theme(plot.title = element_text(size = rel(3.0), hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(2.8)),
        axis.title = element_text(size = rel(2.8)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 2.75, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(size = rel(2.5)),
        legend.title = element_text(size = rel(2.6)),
        legend.text = element_text(size = rel(2.4)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(5.5, 5.5, 5.5, 1)) +
  facet_wrap(~oneBackFP,
             labeller = as_labeller(c(`1` = "FP[n-1] == 1.0",
                                      `1.6` = "FP[n-1] == 1.6",
                                      `2.2` = "FP[n-1] == 2.2",
                                      `2.8` = "FP[n-1] == 2.8"),
                                    default = label_parsed)) +
  scale_color_manual(values = c("orange", "blue"), labels = c("External", "Action")) #+
  #scale_fill_manual(values = c("orange", "blue"), labels = c("External", "Action"))

ggsave("G:/My Drive/Post-doc/Eventos/TRF-3/Poster/seqEffs_exp2.tiff",
       seqEff_by_oneback,
       width = 35,
       height = 23.32,
       units = "cm")

#=============================================================================================#
#====================================== Plots em português ===================================#
#=============================================================================================#

RT_by_condition <- ggplot(data = summaryData2 %>% 
                            group_by(ID, foreperiod, condition) %>% 
                            summarise(meanRT = mean(meanRT)),
                          aes(x = foreperiod,
                              y = meanRT,
                              color = condition)) +
  geom_jitter(height = 0, width = 0.15, size = 3.5, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 3.6, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 3.4, width = 0.1, geom = "errorbar") + 
  labs(title = "TR",
       x = "",
       y = "TR médio (s)",
       color = "Condição") +
  theme(plot.title = element_text(size = rel(2.8), hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = rel(2.8)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 2.75, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(size = rel(2.5)),
        legend.title = element_text(size = rel(2.6)),
        legend.text = element_text(size = rel(2.4)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(5.5, 5.5, 5.5, 1)) +
  scale_color_manual(values = c("orange", "blue"), labels = c("Externa", "Ação"))
