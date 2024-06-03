

# Load packages

# Read and process data
library(tidyverse)
library(magrittr)
library(data.table)

# Plotting
library(ggplot2)
library(lattice)
library(gtable)
library(gridExtra)
library(gridGraphics)
library(ggdist)
library(ggpubr)
library(extrafont)
library(egg)

# Descriptives
library(Hmisc)

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

# Load fonts from extrafonts
loadfonts()

# Prepare data 
source('./Analysis/Prepare_data_6.R')

# Prepare theme for plots
source("./Analysis/plot_theme.R")
theme_set(mytheme)

#==========================================================================================#
#======================================= 0. Data quality ===================================
#==========================================================================================#
# 0.1.1. RTs across blocks (conditions aggregated)
ggplot(data=summaryData,
       aes(x=block,
           y=meanRT))+
  stat_summary(fun='mean',geom='point')+
  stat_summary(fun='mean',geom='line',linewidth=1,aes(group=1))+
  stat_summary(fun.data='mean_cl_boot',width=0.2,geom='errorbar')+
  labs(title='RT by block')

# 0.1.2. RTs across blocks (separated by condition)
ggplot(data=summaryData,
       aes(x=block,
           y=meanRT,
           color=condition))+
  stat_summary(fun='mean',geom='point')+
  stat_summary(fun='mean',geom='line',size=1,aes(group=condition))+
  stat_summary(fun.data='mean_cl_boot',width=0.2,geom='errorbar')+
  cond_cols +
  labs(title='RT by block and condition')

# 0.1.3. RTs across blocks by counterbalancing order
ggplot(data=summaryData,
       aes(x=block,
           y=meanRT,
           color=counterbalance))+
  stat_summary(fun='mean',geom='point')+
  stat_summary(fun='mean',geom='line',size=1,aes(group=counterbalance))+
  stat_summary(fun.data='mean_cl_boot',width=0.2,geom='errorbar')+
  cond_cols +
  labs(title='RT by block split by counterbalancing order')

# 0.1.4. By testing order (first or second condition)
ggplot(data=summaryData,
       aes(x=testPos,
           y=meanRT,
           color=condition))+
  stat_summary(fun='mean',geom='point')+
  stat_summary(fun='mean',geom='line',size=1,aes(group=condition))+
  stat_summary(fun.data='mean_cl_boot',width=0.2,geom='errorbar')+
  cond_cols +
  labs(title='RT by block split by counterbalancing order')

# 0.1.5. Distribution of data
dataHists <- ggplot(data=summaryData,
                    aes(x=meanRT,
                        color=foreperiod))+
  geom_histogram()+
  facet_grid(foreperiod~condition)

# 0.1.6. Individual histograms
indHistograms <- ggplot(data=data,
                        aes(x=RT))+
  geom_histogram()+
  facet_wrap(~ID)
indHistograms

# 0.1.7. Residuals distribution for each variable transformation
qqmath(~RT|ID, data=data)
qqmath(~invRT|ID, data=data)
qqmath(~logRT|ID, data=data)


# 0.1.8. Check for influence of external fixation duration
ggplot(data=filter(data,condition=='external'),
       aes(x=extFixationDuration,
           y=RT,
           color=foreperiod))+
  geom_point() +
  labs(x = "External fixation duration",
       y = "RT") +
  facet_wrap(~foreperiod)
ggsave("./Analysis/Plots/extfixduration.png",
       width = 13.4,
       height = 10)

# 0.1.9. Check for influence of latency of action key press on RT
ggplot(data=filter(data,condition=='action'),
       aes(x=action_trigger.rt,
           y=RT,
           color=foreperiod))+
  geom_point() +
  labs(x = "Action trigger delay",
       y = "RT") +
  facet_wrap(~foreperiod)
ggsave("./Analysis/Plots/actiontrigpress.png",
       width = 13.4,
       height = 10)

# 0.1.10. Correlation between RT and total ITI
ggplot(data = data,
       aes(x = ITItotal,
           y = RT,
           color = condition)) +
  geom_jitter(alpha = 0.3) +
  labs(x = "Total ITI",
       y = "RT") +
  facet_wrap(~foreperiod) +
  cond_cols


# 0.1.11. Histograms of delays between action and WS
ggplot(data = data2) +
  geom_histogram(aes(x = round(delay,3)))

# 0.1.12. Frequency of values in delay vector
table(round(data2$delay, 3))

# 0.1.13. Plot RT against delay excluding zero-delay trials
ggplot(data = filter(data2, delay > 0, delay < 0.034), aes(x = delay, y = RT)) +
  geom_point() +
  geom_abline()

# 0.1.14. Plot RT against delay with best-fitting regression line
ggplot(data = data2, aes(x = delay, y = RT)) +
  geom_point() +
  geom_abline() +
  ylim(c(0,1.0))

# 0.1.15. Correlation between RT and delay excluding zero-delay trials
cor(filter(goData2, delay > 0)$delay, filter(goData2, delay > 0)$RT) 

# 0.1.16. Plot data by foreperiod only to ascertain that there is an FP effect
ggplot(data = summaryData,
       aes(x = foreperiod,
           y = meanRT,
           group = 1)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar")

# 0.1.17. Plot RT by foreperiod in external condition only (baseline)
ggplot(data = filter(summaryData,
                     condition == "external"),
       aes(x = foreperiod,
           y = meanRT,
           group = 1,
           color = cond_cols[1])) +
  stat_summary(fun = 'mean', geom = 'point') +
  stat_summary(fun = 'mean', geom = 'line' ) +
  stat_summary(fun.data = 'mean_cl_boot', width = 0.2, geom = 'errorbar') +
  facet_wrap(~ID)

#==========================================================================================#
#================================= 1. Stopping-rule ======================================
#==========================================================================================#

#====================== 1.1. Individual linear models comparison =========================

# Variables used as predictors: numForeperiod and numOneBackFP
# Dependent variable: logRT
# Variables nested by condition and ID

buildmodel <- function(data) {
  lm(logRT ~ numForeperiod*numOneBackFP,
     data = data)
}

nested_data <- data2 %>%
  select(ID, condition, numForeperiod, numOneBackFP, logRT) %>%
  group_by(ID, condition) %>%
  nest()

fitted_data <- nested_data %>%
  mutate(fit = map(data, buildmodel),
         params = map(fit, tidy)) %>%
  ungroup() %>%
  unnest(c(params)) %>%
  select(ID, condition, term, estimate) %>%
  pivot_wider(names_from = term,
              values_from = estimate)
  

# Foreperiod
fp_bfs <- ttestBF(x = fitted_data$numForeperiod[fitted_data$condition=='external'],
        y = fitted_data$numForeperiod[fitted_data$condition=='action'],
        paired=TRUE)

onebackfp_bfs <- ttestBF(x = fitted_data$numOneBackFP[fitted_data$condition=='external'],
        y = fitted_data$numOneBackFP[fitted_data$condition=='action'],
        paired=TRUE)

interact_bfs <- ttestBF(x = fitted_data$`numForeperiod:numOneBackFP`[fitted_data$condition=='external'],
        y = fitted_data$`numForeperiod:numOneBackFP`[fitted_data$condition=='action'],
        paired=TRUE)


#============================ 1.2. Mixed models BF comparison ============================ 

library(buildmer)

# Condition x fp n
with_condition <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                          (1 + condition + numForeperiod | ID),
                        data = data2,
                        control = lmerControl(optimizer = c('bobyqa'), optCtrl = list(maxfun=2e5), calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method = 'S',
                        REML = TRUE,
                        return = 'merMod')

isSingular(with_condition)

no_condition <- mixed(formula = logRT ~ 1 + numForeperiod +
                        (1 + numForeperiod | ID),
                      data = data2,
                      control = lmerControl(optimizer = c('bobyqa'), optCtrl = list(maxfun=2e5), calc.derivs = FALSE),
                      progress = TRUE,
                      expand_re = TRUE,
                      method = 'S',
                      REML = TRUE,
                      return = 'merMod')

isSingular(no_condition)

bic_to_bf(c(BIC(no_condition),
            BIC(with_condition)),
          denominator = c(BIC(no_condition)))

# Condition x fp n x fp n-1

with_onebackfp <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                        oneBackFP + numForeperiod:oneBackFP + condition:oneBackFP + 
                        (1 + condition + numForeperiod | ID),
                      data = data2,
                      control = lmerControl(optimizer = c('bobyqa'), optCtrl = list(maxfun=2e5), calc.derivs = FALSE),
                      progress = TRUE,
                      expand_re = TRUE,
                      method = 'S',
                      REML = TRUE,
                      return = 'merMod')

isSingular(with_onebackfp)

no_onebackfp <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod +
                          (1 + condition + numForeperiod | ID),
                        data = data2,
                        control = lmerControl(optimizer = c('bobyqa'), optCtrl = list(maxfun=2e5), calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method = 'S',
                        REML = TRUE,
                        return = 'merMod')


isSingular(no_onebackfp)

BIC(no_onebackfp, with_onebackfp)

bic_to_bf(c(BIC(no_onebackfp),
            BIC(with_onebackfp)),
          denominator = c(BIC(no_onebackfp)))


# Condition x fp n x fp n-1 (with and without condition)
with_condition <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                          oneBackFP + numForeperiod:oneBackFP + condition:oneBackFP + condition:numForeperiod:oneBackFP +
                          (1 + condition + numForeperiod | ID),
                        data = data2,
                        control = lmerControl(optimizer = c('bobyqa'), optCtrl = list(maxfun=2e5), calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method = 'S',
                        REML = TRUE,
                        return = 'merMod')

isSingular(with_condition)

no_condition <- mixed(formula = logRT ~ 1 + numForeperiod + oneBackFP + numForeperiod:oneBackFP +
                      (1 + numForeperiod | ID),
                      data = data2,
                      control = lmerControl(optimizer = c('bobyqa'), optCtrl = list(maxfun=2e5), calc.derivs = FALSE),
                      progress = TRUE,
                      expand_re = TRUE,
                      method = 'S',
                      REML = TRUE,
                      return = 'merMod')

isSingular(no_condition)

bic_to_bf(c(BIC(no_condition),
            BIC(with_condition)),
          denominator = c(BIC(no_condition)))

#============================== 1.3. Sequential bayes factors ===========================

external_fits <- fitted_data[fitted_data$condition=='external',]
action_fits <- fitted_data[fitted_data$condition=='action',]

srange <- 10:nrow(external_fits)

fp_bfs <- sapply(srange, function(range) {
  extractBF(ttestBF(x = external_fits$numForeperiod[1:range],
                    y = action_fits$numForeperiod[1:range],
                    paired=TRUE),
            onlybf = TRUE)
})

plot(srange, fp_bfs)
lines(srange, fp_bfs)

onebackfp_bfs <- sapply(srange, function(range) {
  extractBF(ttestBF(x = external_fits$numOneBackFP[1:range],
                    y = action_fits$numOneBackFP[1:range],
                    paired=TRUE),
            onlybf = TRUE)
})

plot(srange, onebackfp_bfs)
lines(srange, onebackfp_bfs)

interact_bfs <- sapply(srange, function(range) {
  extractBF(ttestBF(x = external_fits$`numForeperiod:numOneBackFP`[1:range],
                    y = action_fits$`numForeperiod:numOneBackFP`[1:range],
                    paired=TRUE),
            onlybf = TRUE)
})

plot(srange, interact_bfs)
lines(srange, interact_bfs)

#============================== 1.4. Mixed models using brms ============================
b_one_back_fp <- brm(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                       numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                       (1 + condition + numForeperiod | ID),
                     data = data2)

b_one_back_fp_full <- brm(formula = logRT ~ numForeperiod * condition * numOneBackFP + 
                            (1+numForeperiod*condition*numOneBackFP|ID),
                     data = data2)


options(digits = 5)
summary(b_one_back_fp)
options(options_defaults)

#==========================================================================================#
#================================ 2. Descriptive analyses ==================================
#==========================================================================================#

#====================== 2.1. Descriptive statistics ===========================
#by(data = dataAcc$Acc, INDICES = dataAcc$condition, FUN = describe)

# Acc by condition
tapply(X = dataAcc$Acc, INDEX = dataAcc$condition, FUN = summary)

tapply(X = dataAcc$Acc, INDEX = dataAcc$condition, FUN = var)

# RT by condition
tapply(X = data2$RT, INDEX = data2$condition, FUN = summary)

tapply(X = data2$RT, INDEX = data2$condition, FUN = var)

# Acc overall
summary(dataAcc$Acc)

var(dataAcc$Acc)

# RT overall
summary(data2$RT)

var(data2$RT)

#============================== 2.2. Plots ====================================

#================== 2.2.1. RT by FP and condition ===================
# Boxplots
boxplots <- ggplot(data=summaryData,
                   aes(x=foreperiod,
                       y=meanRT,
                       color=condition))+
  geom_boxplot()+
  cond_cols

boxplots

# Histograms
histograms <- ggplot(data=summaryData,
                     aes(x=meanRT,
                         color=foreperiod))+
  geom_histogram()+
  facet_grid(foreperiod~condition)

# Histograms by participant
histograms <- ggplot(data=data,
                     aes(x=RT))+
  geom_histogram()+
  facet_wrap(~ID)

# Check for correlation between mean and SD of RT
RTstats <- data2 %>%
  group_by(ID, condition, foreperiod) %>%
  summarise(meanRT = mean(RT), sdRT = sd(RT)) %>%
  ungroup()

ggplot(data=RTstats,
       aes(x=meanRT,
           y=sdRT,
           color = condition)) +
  geom_point() +
  geom_smooth(method='lm') +
  cond_cols

View(RTstats[which.max(RTstats$sdRT),])

cor(RTstats_full$meanRT, RTstats_full$sdRT)

data2$scaledNumForeperiod <- scale(data2$numForeperiod, center = TRUE, scale = TRUE)
data2$scaledNumOneBackFP <- scale(data2$numOneBackFP, center = TRUE, scale = TRUE)



# FP and condition
rt_by_condition <- ggplot(data = summaryData2 %>% 
                               group_by(ID, foreperiod, condition) %>% 
                               summarise(meanRT = mean(meanRT)),
                             aes(x = foreperiod,
                                 y = meanRT,
                                 color = condition)) +
  geom_jitter(height = jthgt, width = jtwdt, size = psz, alpha = alp) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = lsz, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = erlsz, width = erlwth, geom = "errorbar") + 
  labs(title = "RT",
       x = "FP (s)", 
       y = "Mean RT (s)",
       color = "Condition") +
  cond_cols
ggplot2::ggsave("./Analysis/Plots/RT_by_condition.jpg",
                rt_by_condition,
                width = 16,
                height = 5,
                unit = "cm",
                dpi = 300)

# sd by FP and condition
sd_by_condition <- ggplot(data = data2 %>% 
                            group_by(ID, foreperiod, condition) %>% 
                            summarise(sdRT = sd(RT)),
                          aes(x = foreperiod,
                              y = sdRT,
                              color = condition)) +
  geom_jitter(height = jthgt, width = jtwdt, size = psz, alpha = alp) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = lsz, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = erlsz, width = erlwth, geom = "errorbar") + 
  labs(title = "RT",
       x = "FP (s)", 
       y = "SD RT (s)",
       color = "Condition") +
  cond_cols
ggplot2::ggsave("./Analysis/Plots/RT_by_condition.png",
                sd_by_condition,
                width = 16,
                height = 5)

#================== 2.2.2. Sequential effects ===================
# Sequential effects separated by FP n-1
seqEff_by_oneback <- ggplot(data = summaryData2 %>%
                              group_by(ID, foreperiod, condition, oneBackFP) %>%
                              summarise(meanRT = mean(meanRT)),
                            aes(x = foreperiod,
                                y = meanRT,
                                color=condition)) +
  geom_jitter(height = jthgt, width = jtwdt, size = psz, alpha = alp) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = lsz, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = erlsz, width = erlwth, geom = "errorbar") + 
  labs(x = "FP (s)",
       y = "Mean RT (s)",
       color = "Condition") +
  facet_wrap(~oneBackFP,
             labeller = as_labeller(c(`1` = "FP[n-1] == 1.0",
                                      `1.6` = "FP[n-1] == 1.6",
                                      `2.2` = "FP[n-1] == 2.2",
                                      `2.8` = "FP[n-1] == 2.8"),
                                    default = label_parsed)) +
  cond_cols
seqEff_by_oneback <- set_panel_size(seqEff_by_oneback, width = unit(4, "cm"),
                                    height = unit(2.6, "cm"))

ggsave("./Analysis/Plots/SeqEff.jpg",
       seqEff_by_oneback,
       width = 10,
       height = 9,
       units = "cm",
       dpi = 300)


#================== 2.2.3. Accuracy ===================
ggplot(data = summaryDataAcc) +
  geom_histogram(aes(x = meanAcc)) +
  facet_grid(foreperiod ~ condition)

# Error rates
error_by_condition <- ggplot(data = summaryDataAcc %>%
                                   group_by(ID, foreperiod, condition) %>%
                                   summarise(errorRate = mean(errorRate)),
                                 aes(x = foreperiod,
                                     y = errorRate,
                                     color = condition)) +
  labs(title = "Error Rate",
       x = "FP (s)",
       y = "Mean Error Rate",
       color = "Condition") +
  geom_jitter(height = jthgt, width = jtwdt, size = psz, alpha = alp) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = lsz, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = erlsz, width = erlwth, geom = "errorbar") + 
  theme(plot.margin = unit(c(5.5, 5.5, 5.5, 1), "pt")) +
  cond_cols

ggsave("./Analysis/Plots/errors_by_fp_condition.jpg",
       error_by_condition,
       width = 16,
       height = 5,
       unit = "cm",
       dpi = 300)

# RT and error rate in single panel
cond_legend <- ggpubr::get_legend(error_by_condition)

xaxis_title <- text_grob(error_by_condition$labels$x, # Extract x axis label to plot in the middle of both plots
                         just = "top",
                         size = mytheme$axis.title$size)
                         #size = mytheme$axis.title$size * 11) # size of grob; 11 is the base size in theme_classic
                          


xaxis_title_margin <- unit(2/5, "line") # margin to separate from other plot elements


# Visualize
grid.arrange(arrangeGrob(rt_by_condition + theme(legend.position = "none",
                                                 axis.title.x = element_blank()),
                         error_by_condition + theme(legend.position = "none",
                                                    axis.title.x = element_blank()),
                         nrow = 1,
                         widths = c(1/2, 1/2)),
             xaxis_title,
             cond_legend,
             heights = unit.c(unit(1, "null"),
                              grobHeight(xaxis_title) + xaxis_title_margin,
                              grobHeight(cond_legend)),
             nrow = 3)

# Save plots
RT_panel <- rt_by_condition + theme(legend.position = "none",
                                    axis.title.x = element_blank())
RT_panel <- set_panel_size(RT_panel,
                           width = unit(4, "cm"),
                           height = unit(2.6, "cm"))


error_panel <- error_by_condition + theme(legend.position = "none",
                                          axis.title.x = element_blank())
error_panel <- set_panel_size(error_panel, 
                              width = unit(4, "cm"),
                              height = unit(2.6, "cm"))


rt_error_plots <- arrangeGrob(arrangeGrob(RT_panel,
                                          error_panel,
                                          nrow = 1,
                                          widths = c(1/2, 1/2)),
                              xaxis_title,
                              cond_legend,
                              heights = unit.c(unit(1, "null"),
                                               grobHeight(xaxis_title) + xaxis_title_margin,
                                               grobHeight(cond_legend)),
                              nrow = 3)

ggsave("./Analysis/Plots/rt_error_plots.jpg",
       rt_error_plots,
       width = 10.5,
       height = 5,
       unit = "cm",
       dpi = 300)



#===================== 2.4. Plots for between-experiments comparisons ========================
extrFPData <- summaryData2 %>%
  filter(foreperiod %in% c("1000", "2800")) %>%
  mutate(foreperiod = fct_relevel(foreperiod, c("1000", "2800")))

rt_extr_fps_part <- ggplot(data = extrFPData,
                           aes(x = foreperiod,
                               y = meanRT,
                               color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", linewidth = 0.8, width = 0.2) +
  labs(title = "RT by condition",
       x = "Foreperiod",
       y = "Mean RT",
       color = "Condition") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(0.8)),
        axis.title = element_text(size = rel(1.2))) +
  scale_color_manual(values = c("orange", "blue"), labels = c("External", "Action")) +
  facet_wrap(~ ID)

ggsave("./Analysis/Plots/rt_extr_fps_part.png",
       rt_extr_fps_part,
       width = 6.7,
       height = 5)

#==========================================================================================#
#======================================= 3. ANOVAs =========================================
#==========================================================================================#

#==================== 3.1. FP x RT by condition ======================
#================ 3.1.1. RT =================

# Run repeated-measures anova
fpAnova <- aov_ez(id = "ID",
       dv = "meanRT",
       data = summaryData2,
       within = c("foreperiod", "condition"))

### Check assumptions

# Sphericity
check_sphericity(fpAnova)

fpAnova_plot <- afex_plot(fpAnova, x = 'foreperiod', trace = 'condition', error = 'within')

# Normality of residuals
is_norm <- check_normality(fpAnova)

plot(is_norm)

plot(is_norm, type = "qq")
plot(is_norm, type = "qq", detrend = TRUE)

testnormality = function(dfr) return(shapiro.test(dfr$invRT)$p.value)
p = as.vector(by(data, data$ID, testnormality))
names(p) = levels(data$ID)
names(p[p < 0.05])


# Try transformations
invfpAnova <- aov_ez(id = "ID",
                  dv = "meanInvRT",
                  data = summaryData2,
                  within = c("foreperiod", "condition"))

### Check assumptions

# Sphericity
check_sphericity(invfpAnova)

# Normality of residuals
is_norm <- check_normality(invfpAnova)

plot(is_norm)

plot(is_norm, type = 'qq')

plot(is_norm, type = 'qq', detrend = TRUE)

# using 1/RT does not solve the problem

logfpAnova <- aov_ez(id = "ID",
                     dv = "meanLogRT",
                     data = summaryData2,
                     within = c("foreperiod", "condition"))


### Check assumptions

# Sphericity
check_sphericity(logfpAnova)

# Normality of residuals
is_norm <- check_normality(logfpAnova)

plot(is_norm)

plot(is_norm, type = 'qq')

plot(is_norm, type = 'qq', detrend = TRUE)


# ============ 3.1.2. Fit models separately for external condition (safety check) ====================
logfpAnovaExt <- aov_ez(id = "ID",
                     dv = "meanLogRT",
                     data = summaryData2[summaryData2$condition=='external',],
                     within = c("foreperiod"))

# Sphericity
check_sphericity(logfpAnovaExt)

# Normality of residuals
is_norm <- check_normality(logfpAnovaExt)

plot(is_norm)

plot(is_norm, type = 'qq')

plot(is_norm, type = 'qq', detrend = TRUE)


# Effect is significant

logfpAnovaExt2 <- aov_ez(id = "ID",
                      dv = "meanLogRT",
                      data = summaryData2[summaryData2$condition=='external',],
                      within = c("logFP"))

# Effect is significant

#============================== 3.2. Sequential effects ================================================

#============= 3.2.1. RT ==================
# 2.2.1. Anova for FP n-1
seqEffAnova <- aov_ez(id = "ID",
                  dv = "meanRT",
                  data = summaryData2,
                  within = c("foreperiod", "condition", "oneBackFP"))

logSeqEffAnova <- aov_ez(id = "ID",
                      dv = "meanLogRT",
                      data = summaryData2,
                      within = c("foreperiod", "condition", "oneBackFP"))

nice(seqEffAnova,
     correction='none')
                  
nice(logSeqEffAnova,
     correction = "none")

# 2.2.3. FP n-2
twoBackAnova <- aov_ez(id = "ID",
                      dv = "meanRT",
                      data = summaryData2,
                      within = c("foreperiod", "condition", "twoBackFP"))

nice(twoBackAnova,
     correction = "none")

ntworegression <- lm(meanRT ~ foreperiod * oneBackFP * twoBackFP, 
                     data = summaryData)
summary(ntworegression)
anova(ntworegression)
Anova(ntworegression, type = "II")


#============= 3.2.2. Accuracy =================
AccAnova <- aov_ez(id = "ID",
                   dv = "meanAcc",
                   data = summaryDataAll,
                   within = c("foreperiod", "condition", "oneBackFP"))

AccAnova <- aov_ez(id = "ID",
                   dv = "meanAcc",
                   data = summaryDataAll,
                   within = c("foreperiod", "condition"))

acc_emm <- emmeans(AccAnova, c("condition", "foreperiod"))
pairs(acc_emm, by = "foreperiod")

# ============ 3.2.2. Fit models separately for external condition (safety check) ====================
seqEffAnovaExt <- aov_ez(id = "ID",
                     dv = "meanLogRT",
                     data = summaryData2[summaryData2$condition=='external',],
                     within = c("foreperiod", "oneBackFP"))

# Effect of FP is significant

seqEffAnovaExt2 <- aov_ez(id = "ID",
                     dv = "meanLogRT",
                     data = summaryData2[summaryData2$condition=='external',],
                     within = c("logFP", "oneBackFP"))

# Effect of FP is significant

# Compare full and reduced models
summary(logfpAnovaExt, correction = 'none')
summary(seqEffAnovaExt, correction = 'none')

summary(logfpAnovaExt2, correction = 'none')
summary(seqEffAnovaExt2, correction = 'none')

# ============ 3.3.2. Fit models separately for external condition (safety check) ====================
fpDiffAnovaExt <- aov_ez(id = "ID",
                     dv = "meanLogRT",
                     data = summaryData2[summaryData2$condition=='external',],
                     within = c("foreperiod", "oneBackFPDiff"))

# Error because there are empty cells

logfpAnova2 <- aov_ez(id = "ID",
                      dv = "meanLogRT",
                      data = summaryData2[summaryData2$condition=='external',],
                      within = c("logFP", "oneBackFPDiff"))

# Error because there are empty cells

#=============================== 3.4 Orientation ===================================
oriAnova <- aov_ez(id = 'ID',
                   dv = 'meanRT',
                   data = summaryData2,
                   within = c('orientation', 'condition', 'foreperiod'))

#============================= 4. Quadratic effects  =====================================
# fpEmmeans <- emmeans(fpAnova,
#                      pairwise ~ condition|foreperiod,
#                      adjust = 'none')


fpEmmeans <- emmeans(fpAnova,
                     pairwise ~ foreperiod|condition,
                     adjust = 'bonferroni')

fpEmmeansContrasts <- contrast(fpEmmeans[[1]],
                               interaction=c('poly'),
                               adjust='bonferroni')

fpEmmeansContrasts <- contrast(fpEmmeans[[1]],
                               interaction=c('consec'),
                               adjust='bonferroni')

contrasts(summaryData2$foreperiod) <- contr.poly(4)

# using emmeans
logfpAnova <- aov_ez(id = "ID",
                     dv = "meanLogRT",
                     data = summaryData2[summaryData2$condition=='external',],
                     within = c("foreperiod"))

fp_posthoc <- emmeans(logfpAnova,
                      specs = "foreperiod")

fp_contrasts <- contrast(fp_posthoc,
                         interaction = 'poly')

# Using lists and anova as aov
logfpAnova <- aov_ez(id = "ID",
                     dv = "meanLogRT",
                     data = summaryData2[summaryData2$condition=='external',],
                     within = c("foreperiod"),
                     return = 'aov')


summary(logfpAnova, split=list(foreperiod=list(linear=1,quadratic=2,cubic=3)))

# Using regression
logfpregression <- lm(meanLogRT ~ numForeperiod, data = summaryData2[summaryData2$condition=='external',])
logfpregression <- lm(meanLogRT ~ numForeperiod * numOneBackFPDiff, 
                      data = summaryData2[summaryData2$condition=='external',])
summary(logfpregression)
anova(logfpregression)

fp_posthoc <- emmeans(logfpregression,
                      specs = "foreperiod", 
                      adjust="tukey")

fp_contrasts <- contrast(fp_posthoc,
                         interaction = "poly",
                         adjust = "tukey")

#=================================== 5. Learning =============================================


#============================ 5.1. Effects across blocks ==============================
# Foreperiod, condition and block
blocklm <- lm(meanRT ~ foreperiod * counterbalance * block,
              data = summaryData)

anova(blocklm)
Anova(blocklm)


#========================== 5.2. Split anovas by counterbalancing order ===================

fpAnova_ae <- aov_ez(id = "ID",
                     dv = "meanRT",
                     data = summaryData2[summaryData2$counterbalance=='action-external',],
                     within = c("foreperiod", "condition"))

fpAnova_ea <- aov_ez(id = "ID",
                     dv = "meanRT",
                     data = summaryData2[summaryData2$counterbalance=='external-action',],
                     within = c("foreperiod", "condition"))


# Regression by block
blocklm <- lm(meanRT ~ foreperiod * condition * block,
              data = summaryData2)

anova(blocklm)
Anova(blocklm)

# 2.4.2. Split participants' trias in 3 bins and compare average RTs
dataBin <- data2 %>%
  group_by(ID, condition) %>%
  mutate(trialBin)



#============================= 6. Residual analysis for sequential effects ==========================

seqEffAnovares1 <- aov_ez(id = "ID",
                      dv = "meanLogRT",
                      data = filter(summaryData2, condition == 'external'),
                      within = c("oneBackFP"),
                      fun_aggregate = mean)

seqEffres <- seqEffAnovares1$lm$residuals

seqEffres <- seqEffAnovares1$lm$residuals %>%
  as_tibble() %>%
  pivot_longer(cols = 1:4, names_to = c("oneBackFP"), values_to = "residuals")
seqEffres$oneBackFP <- fct_recode(seqEffres$oneBackFP, '1000' = 'X1000', '1600' = 'X1600', '2200' = 'X2200', '2800' = 'X2800')

aovez_lm <- function(data) {
  model <- aov_ez(id = "ID",
                  dv = "meanLogRT",
                  data = data,
                  within = c("foreperiod", "oneBackFP"),
                  fun_aggregate = mean)
  modellm <- model$lm
}

aovez_lm2 <- function(data) {
  model <- aov_ez(id = "ID",
                  dv = "meanLogRT",
                  data = data,
                  within = c("oneBackFP"),
                  fun_aggregate = mean)
  modellm <- model$lm
}

seqEffResData <- summaryData2 %>%
  group_by(condition) %>%
  nest() %>%
  mutate(model = map(data, ~aovez_lm(data)))

seqEffResData <- summaryData2 %>%
  group_by(condition) %>%
  nest() %>%
  mutate(model = map(data, ~lm(meanLogRT ~ oneBackFP, data = .x)),
         predictions = map2(data, model, add_predictions))

# Add residuals
seqEffResData <- summaryData2 %>%
  group_by(condition) %>%
  nest() %>%
  mutate(model = map(data, ~aovez_lm2(data)),
         residuals = map2(data, model, add_residuals))

# Add predictions
seqEffResData <- summaryData2 %>%
  group_by(condition) %>%
  nest() %>%
  mutate(model = map(data, ~aovez_lm2(data)),
         predictions = map2(data, model, add_predictions))


seqEffResData <- summaryData2 %>%
  group_by(condition) %>%
  nest() %>%
  mutate(model = map(data, ~aovez_lm(data))) %>%
  unnest(c(data))


# Add predictions using broom
seqEffResData <- summaryData2 %>%
  group_by(condition) %>%
  nest() %>%
  mutate(model = map(data, ~aovez_lm(data)),
         model_data = map(model, broom::augment)) %>%
  unnest(model_data)

aovez_resid <- function(data) {
  model <- aov_ez(id = "ID",
                  dv = "meanLogRT",
                  data = data,
                  within = c("foreperiod", "oneBackFP"),
                  fun_aggregate = mean)
  model_res <- model$lm$residuals %>%
    as_tibble() %>%
    pivot_longer(cols = 1:16, names_to = c("foreperiod", "oneBackFP"), names_pattern = "(.*)_(.*)", values_to = "residuals")
  model_res$foreperiod <- fct_recode(model_res$foreperiod, '1000' = 'X1000', '1600' = 'X1600', '2200' = 'X2200', '2800' = 'X2800')
  model_res$oneBackFP <- fct_recode(model_res$oneBackFP, '1000' = 'X1000', '1600' = 'X1600', '2200' = 'X2200', '2800' = 'X2800')
  model_res
}

aovez_resid2 <- function(data) {
  model <- aov_ez(id = "ID",
                  dv = "meanLogRT",
                  data = data,
                  within = c("oneBackFP"),
                  fun_aggregate = mean)
  model_res <- model$lm$residuals %>%
    as_tibble() %>%
    pivot_longer(cols = 1:4, names_to = c("oneBackFP"), values_to = "residuals")
  model_res$oneBackFP <- fct_recode(model_res$oneBackFP, '1000' = 'X1000', '1600' = 'X1600', '2200' = 'X2200', '2800' = 'X2800')
  model_res
}

seqEffResData <- summaryData2 %>%
  group_by(condition) %>%
  nest() %>%
  mutate(model = map(data, ~aovez_lm(data)),
         residuals = map(data, ~aovez_resid(.x)))

aggSumData2 <- summaryData2 %>%
  group_by(ID, condition, foreperiod, oneBackFP) %>%
  summarise(meanLogRT = mean(meanLogRT)) %>%
  ungroup()

seqEffResInt <- aggSumData2 %>%
  group_by(condition, foreperiod) %>%
  nest() %>%
  mutate(residuals = map(data, ~aovez_resid2(.x))) %>%
  unnest(residuals) %>%
  ungroup()

seqEffResData <- aggSumData2 %>%
  mutate(residuals = seqEffResInt$residuals)

seqEffAnovares2 <- aov_ez(id = 'ID',
                          dv = 'residuals',
                          data = seqEffResData,
                          within = c('condition'))

summary(seqEffAnovares2)

condAnova <- aov_ez(id = "ID",
                    dv = "meanLogRT",
                    data = summaryData2,
                    within = c("condition"))

#==================================================================================#
# Functions to extract lm and residuals from aov_ez
aovez_lm <- function(data) {
  model <- aov_ez(id = "ID",
                  dv = "meanLogRT",
                  data = data,
                  within = c("foreperiod", "oneBackFP"),
                  fun_aggregate = mean)
  modellm <- model$lm
}

aovez_resid <- function(data) {
  model <- aov_ez(id = "ID",
                  dv = "meanLogRT",
                  data = data,
                  within = c("foreperiod", "oneBackFP"),
                  fun_aggregate = mean)
  model_res <- model$lm$residuals %>%
    as_tibble() %>%
    pivot_longer(cols = 1:16, names_to = c("foreperiod", "oneBackFP"), names_pattern = "(.*)_(.*)", values_to = "residuals")
  model_res$foreperiod <- fct_recode(model_res$foreperiod, '1000' = 'X1000', '1600' = 'X1600', '2200' = 'X2200', '2800' = 'X2800')
  model_res$oneBackFP <- fct_recode(model_res$oneBackFP, '1000' = 'X1000', '1600' = 'X1600', '2200' = 'X2200', '2800' = 'X2800')
  model_res
}

# Aggregate over other cells
aggSumData <- summaryData2 %>%
  group_by(ID, condition, foreperiod, oneBackFP) %>%
  summarise(meanLogRT = mean(meanLogRT)) %>%
  ungroup()

# Save residuals to aggregated dataset
seqEffResInt <- aggSumData %>%
  group_by(condition) %>%
  nest() %>%
  mutate(residuals = map(data, ~aovez_resid(.x))) %>%
  unnest(residuals) %>%
  ungroup()

seqEffResData <- aggSumData %>%
  mutate(residuals = seqEffResInt$residuals)

# Run anova with condition as IV
seqEffAnovares2 <- aov_ez(id = 'ID',
                          dv = 'residuals',
                          data = seqEffResData,
                          within = c('condition'))


condAnova <- aov_ez(id = "ID",
                    dv = "meanLogRT",
                    data = summaryData2,
                    within = c("condition"))


#=====================#
aovez_resid2 <- function(data) {
  model <- aov_ez(id = "ID",
                  dv = "meanLogRT",
                  data = data,
                  within = c("oneBackFP"),
                  fun_aggregate = mean)
  model_res <- model$lm$residuals %>%
    as_tibble() %>%
    pivot_longer(cols = 1:4, names_to = c("oneBackFP"), values_to = "residuals")
  model_res$oneBackFP <- fct_recode(model_res$oneBackFP, '1000' = 'X1000', '1600' = 'X1600', '2200' = 'X2200', '2800' = 'X2800')
  model_res
}

seqEffResData2 <- summaryData2 %>%
  group_by(condition) %>%
  nest() %>%
  mutate(model = map(data, ~aovez_lm(data)),
         residuals = map(data, ~aovez_resid2(.x)))

aggSumData2 <- summaryData2 %>%
  group_by(ID, condition, foreperiod, oneBackFP) %>%
  summarise(meanLogRT = mean(meanLogRT)) %>%
  ungroup()

seqEffResInt2 <- aggSumData2 %>%
  group_by(condition, foreperiod) %>%
  nest() %>%
  mutate(residuals = map(data, ~aovez_resid2(.x))) %>%
  unnest(residuals) %>%
  ungroup()

seqEffResData2 <- aggSumData2 %>%
  mutate(residuals = seqEffResInt2$residuals)

seqEffAnovares2 <- aov_ez(id = 'ID',
                          dv = 'residuals',
                          data = seqEffResData2,
                          within = c('condition'))


condAnova <- aov_ez(id = "ID",
                    dv = "meanLogRT",
                    data = summaryData2,
                    within = c("condition"))

#============================#

aovez_resid3 <- function(data) {
  model <- aov_ez(id = "ID",
                  dv = "meanLogRT",
                  data = data,
                  within = c("foreperiod"),
                  fun_aggregate = mean)
  model_res <- model$lm$residuals %>%
    as_tibble() %>%
    pivot_longer(cols = 1:4, names_to = c("foreperiod"), values_to = "residuals")
  model_res$foreperiod <- fct_recode(model_res$foreperiod, 
                                     '1000' = 'X1000', '1600' = 'X1600', '2200' = 'X2200', '2800' = 'X2800')
  model_res
}

aggSumData3 <- summaryData2 %>%
  group_by(ID, condition, foreperiod, oneBackFP) %>%
  summarise(meanLogRT = mean(meanLogRT)) %>%
  ungroup()

seqEffResInt3 <- aggSumData3 %>%
  group_by(condition, oneBackFP) %>%
  nest() %>%
  mutate(residuals = map(data, ~aovez_resid3(.x))) %>%
  unnest(residuals) %>%
  ungroup()

seqEffResData3 <- aggSumData3 %>%
  mutate(residuals = seqEffResInt3$residuals)

seqEffAnovares3 <- aov_ez(id = 'ID',
                          dv = 'residuals',
                          data = seqEffResData3,
                          within = c('condition'))

summary(seqEffAnovares4)