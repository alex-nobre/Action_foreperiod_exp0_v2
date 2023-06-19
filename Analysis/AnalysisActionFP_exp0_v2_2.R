

#================================================================================#
# Changes
# Run analyses on scaled predictors
#================================================================================#

# Load necessary packages
library(magrittr)
library(tidyverse)
library(broom)
library(lattice)
library(afex)
library(emmeans)
library(lme4)
library(car)
library(data.table)
library(codingMatrices)
library(performance)
library(modelr)
library(bayestestR)
library(BayesFactor)
library(brms)
library(rtdists)


# Save defaults
graphical_defaults <- par()
options_defaults <- options() 

#==========================================================================================#
#======================================= 0. Prepare data ===================================
#==========================================================================================#
# Functions for raincloud plots
#source("./Analysis/R_rainclouds.R")

source('./Analysis/Prepare_data_exp0_v2_6.R')

#==========================================================================================#
#======================================= 1. Data quality ===================================
#==========================================================================================#
# 1.1.1. RTs across blocks (conditions aggregated)
ggplot(data=summaryData,
       aes(x=block,
           y=meanRT))+
  stat_summary(fun='mean',geom='point')+
  stat_summary(fun='mean',geom='line',size=1,aes(group=1))+
  stat_summary(fun.data='mean_cl_boot',width=0.2,geom='errorbar')+
  theme(plot.title=element_text(size = rel(2), hjust = 0.5),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)))+
  labs(title='RT by block')

# 1.1.2. RTs across blocks (separated by condition)
ggplot(data=summaryData,
       aes(x=block,
           y=meanRT,
           color=condition))+
  stat_summary(fun='mean',geom='point')+
  stat_summary(fun='mean',geom='line',size=1,aes(group=condition))+
  stat_summary(fun.data='mean_cl_boot',width=0.2,geom='errorbar')+
  theme(plot.title=element_text(size = rel(2), hjust = 0.5),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)))+
  scale_color_manual(values=c('orange','blue')) +
  labs(title='RT by block and condition')

# 1.1.3. RTs across blocks by counterbalancing order
ggplot(data=summaryData2,
       aes(x=block,
           y=meanRT,
           color=counterbalance))+
  stat_summary(fun='mean',geom='point')+
  stat_summary(fun='mean',geom='line',size=1,aes(group=counterbalance))+
  stat_summary(fun.data='mean_cl_boot',width=0.2,geom='errorbar')+
  theme(plot.title=element_text(size = rel(2), hjust = 0.5),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)))+
  scale_color_manual(values=c('blue','orange')) +
  labs(title='RT by block split by counterbalancing order')

ggplot(data=summaryData,
       aes(x=testPos,
           y=meanRT,
           color=condition))+
  stat_summary(fun='mean',geom='point')+
  stat_summary(fun='mean',geom='line',size=1,aes(group=condition))+
  stat_summary(fun.data='mean_cl_boot',width=0.2,geom='errorbar')+
  theme(plot.title=element_text(size = rel(2), hjust = 0.5),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)))+
  scale_color_manual(values=c('orange','blue')) +
  labs(title='RT by block split by counterbalancing order')

# Distribution of data
dataHists <- ggplot(data=summaryData,
                    aes(x=meanRT,
                        color=foreperiod))+
  geom_histogram()+
  facet_grid(foreperiod~condition)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Individual histograms
indHistograms <- ggplot(data=data2,
                        aes(x=RT))+
  geom_histogram()+
  facet_wrap(~ID)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
indHistograms

qqmath(~RT|ID, data=data)
qqmath(~invRT|ID, data=data)

# 1.4.1. Plot RT by foreperiod only
ggplot(data = summaryData,
       aes(x = foreperiod,
           y = meanRT,
           group = 1)) +
  stat_summary(fun = 'mean', geom = 'point') +
  stat_summary(fun = 'mean', geom = 'line' ) +
  stat_summary(fun.data = 'mean_cl_boot', width = 0.2, geom = 'errorbar') +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +#,
        #axis.text = element_text(size = rel(1.5)),
        #axis.title = element_text(size = rel(1.5))) +
  facet_wrap(~ID)

# Check for influence of external fixation duration
ggplot(data=filter(data,condition=='external'),
       aes(x=extFixationDuration,
           y=RT,
           color=foreperiod))+
  geom_point() +
  facet_wrap(~foreperiod)

# Check for influence of latency of action key press on RT
ggplot(data=filter(data,condition=='action'),
       aes(x=action_trigger.rt,
           y=RT,
           color=foreperiod))+
  geom_point() +
  facet_wrap(~foreperiod)

# Plot data by foreperiod only to ascertain that there is an FP effect
ggplot(data = summaryData,
       aes(x = foreperiod,
           y = meanRT,
           group = 1)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)))



#==========================================================================================#
#================================= 1.2. Stopping-rule ======================================
#==========================================================================================#

#====================== 1.2.1. Individual linear models comparison =========================

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


#============================ 1.2.2. Mixed models BF comparison ============================ 
library(afex)
library(lme4)
library(buildmer)

with_onebackfp <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                        numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                        (1 + condition + numForeperiod | ID),
                      data = data2,
                      control = lmerControl(optimizer = c('bobyqa'), optCtrl = list(maxfun=2e5), calc.derivs = FALSE),
                      progress = TRUE,
                      expand_re = TRUE,
                      method = 'S',
                      REML = TRUE,
                      return = 'merMod')

isSingular(no_onebackfp)

no_onebackfp <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod +
                          (1 + condition + numForeperiod | ID),
                        data = data2,
                        control = lmerControl(optimizer = c('bobyqa'), optCtrl = list(maxfun=2e5), calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method = 'S',
                        REML = TRUE,
                        return = 'merMod')


isSingular(with_onebackfp)

BIC(no_onebackfp, with_onebackfp)

bic_to_bf(c(BIC(no_onebackfp),
            BIC(with_onebackfp)),
          denominator = c(BIC(no_onebackfp)))



with_condition <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                          numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                          (1 + condition + numForeperiod | ID),
                        data = data2,
                        control = lmerControl(optimizer = c('bobyqa'), optCtrl = list(maxfun=2e5), calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method = 'S',
                        REML = TRUE,
                        return = 'merMod')

isSingular(with_condition)

no_condition <- mixed(formula = logRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP +
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

#============================== 1.2.3. Sequential bayes factors ===========================

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

#============================== 1.2.4. Mixed models using brms ============================
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
#================================ 2. Descriptive analysis ==================================
#==========================================================================================#
# Boxplots
boxplots <- ggplot(data=summaryData,
                   aes(x=foreperiod,
                       y=meanRT,
                       color=condition))+
  geom_boxplot()+
  scale_color_manual(values=c('orange','blue'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

boxplots


logBoxplots <- ggplot(data=summaryData,
                      aes(x=foreperiod,
                          y=meanLogRT,
                          color=condition))+
  geom_boxplot()+
  scale_color_manual(values=c('orange','blue'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

logBoxplots


# Histograms
histograms <- ggplot(data=summaryData,
                     aes(x=meanRT,
                         color=foreperiod))+
  geom_histogram()+
  facet_grid(foreperiod~condition)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

histograms <- ggplot(data=data,
                     aes(x=RT))+
  geom_histogram()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~ID)

# Plot histograms by participant and condition

RTstats <- data2 %>%
  group_by(ID, condition) %>%
  summarise(meanRT = mean(RT), sdRT = sd(RT)) %>%
  ungroup()


ggplot(data=RTstats,
       aes(x=meanRT,
           y=sdRT)) +
  geom_point() +
  geom_smooth(method='lm')


# Do the same by participant, condition and foreperiod
RTstats_full <- data2 %>%
  group_by(ID, condition, foreperiod) %>%
  summarise(meanRT = mean(RT), sdRT = sd(RT)) %>%
  ungroup()

ggplot(data=RTstats_full,
       aes(x=meanRT,
           y=sdRT)) +
  geom_point() +
  geom_smooth(method='lm')

View(RTstats[which.max(RTstats$sdRT),])

cor(RTstats_full$meanRT, RTstats_full$sdRT)

ggplot(data=data2,
       aes(x=RT)) +
  geom_histogram() +
  facet_wrap(~ ID + condition)

data2$scaledNumForeperiod <- scale(data2$numForeperiod, center = TRUE, scale = TRUE)
data2$scaledNumOneBackFP <- scale(data2$numOneBackFP, center = TRUE, scale = TRUE)

# Try distributions
invgausfit <- brm(RT ~ condition + (1|ID),
                  dat = data2,
                  family = inverse.gaussian(),
                  file = 'invgaussianfit')

#==========================================================================================#
#======================================= 3. ANOVAs =========================================
#==========================================================================================#

# Set constrasts for variables used in ANOVAs
contrasts(summaryData$foreperiod) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)
contrasts(summaryData$condition) <- c(-1/2, 1/2)
contrasts(summaryData$oneBackFP) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)

contrasts(summaryData2$foreperiod) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)
contrasts(summaryData2$condition) <- c(-1/2, 1/2)
contrasts(summaryData2$oneBackFP) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)

contrasts(summaryDataAll$foreperiod) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)
contrasts(summaryDataAll$condition) <- c(-1/2, 1/2)
contrasts(summaryDataAll$oneBackFP) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)

#==================== 3.1. FP x RT by condition ======================
#================ 3.1.1. RT =================

# FP and condition by participante
fp_cond_part <- ggplot(data = summaryData2,
                             aes(x = foreperiod,
                                 y = meanRT,
                                 color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.1, geom = "errorbar") +
  labs(title = "RT by FP and condition",
       x = "Foreperiod", 
       y = "Mean RT",
       color = "Condition") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_color_manual(values = c("orange","blue")) +
  facet_wrap(~ ID)


# FP and condition
lines_by_condition <- ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.1, geom = "errorbar") +
  labs(title = "RT by FP and condition",
       x = "Foreperiod", 
       y = "Mean RT",
       color = "Condition") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_color_manual(values = c("orange","blue"))
ggplot2::ggsave("./Analysis/Plots/RT_by_condition.png",
                lines_by_condition,
                width = 8.5,
                height = 5.7)

# using LogRT
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanLogRT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", size = 1, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  scale_color_manual(values = c("orange", "blue"))


# Individual panels by participant for RT and log RT
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_color_manual(values = c("orange","blue")) +
  facet_wrap(~ID)

ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanLogRT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_color_manual(values = c("orange","blue")) +
  facet_wrap(~ID)

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

# The log-transform does not solve the problem either

fpregression <- lm(meanRT ~ condition * foreperiod, data = summaryData)
summary(fpregression)
anova(fpregression)

logfpregression <- lm(meanRT ~ condition * logFP, data = summaryData)
anova(logfpregression)

#================ 3.1.2. Accuracy =================
ggplot(data = summaryDataAll,
       aes(x = foreperiod,
           y = meanAcc,
           color = condition)) +
  labs(title = "Accuracy (proportion of correct responses)",
       x = "Foreperiod",
       y = "Mean Acc" ) +
  stat_summary(fun = "mean", geom = "point", size = 1.8) +
  stat_summary(fun = "mean", geom = "line", aes(group = condition), linewidth = 0.8) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", linewidth = 0.7, width = 0.05) +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_color_manual(values = c("orange", "blue"))
ggsave("./Analysis/Plots/acc_by_condition.png",
       width = 7.5,
       height = 5)

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
# Aggregated across conditions
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color=oneBackFP)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", size = 0.8, aes(group=oneBackFP)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1, size = 0.8) +
  labs(x = "Foreperiod",
       y = "Mean RT",
       color = "FP n -1") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)))+
  scale_color_manual(values = c('blue','orange','green', 'magenta'))
ggsave("./Analysis/Plots/seqeff.png",
       width = 8.5,
       height = 5.7)

# Separated by condition
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color = oneBackFP)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.6, aes(group = oneBackFP)) +
  stat_summary(fun.data = "mean_cl_boot", size = 0.3, width = 0.05, geom = "errorbar") +
  labs(x = "Foreperiod",
       y = "Mean RT",
       color = "FP n-1") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  facet_wrap(~condition) +
  scale_color_manual(values = c('blue','orange','green', 'magenta'))
ggsave("./Analysis/Plots/seqeff_by_condition.png",
       width = 8.5,
       height = 5.7)

# Differences by condition
ggplot(data = summaryData2,
       aes(x = oneBackFP,
           y = meanSeqEff,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.6, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", size = 0.4, width = 0.1, geom = "errorbar") +
  facet_wrap(~ foreperiod) +
  scale_color_manual(values = c("orange", "blue"))


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
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color = twoBackFP)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.6, aes(group = twoBackFP)) +
  stat_summary(fun.data = "mean_cl_boot", size = 0.3, width = 0.05, geom = "errorbar") +
  labs(x = "Foreperiod",
       y = "Mean RT",
       color = "FP n-2") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  facet_wrap(~condition) +
  scale_color_manual(values = c('blue','orange','green', 'magenta'))

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

# Effects of previous orientation
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color=prevOri)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group = prevOri)) +
  stat_summary(fun.data = "mean_cl_boot", size = 0.8, width = 0.2, geom = "errorbar") +
  labs(x = "Foreperiod",
       y = "Mean RT",
       color = "Previous Orientation") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  facet_wrap(~condition) +
  scale_color_manual(values = c('blue', 'orange'))

ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color=seqOri)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group = seqOri)) +
  stat_summary(fun.data = "mean_cl_boot", size = 0.8, width = 0.2, geom = "errorbar") +
  labs(x = "Foreperiod",
       y = "Mean RT",
       color = "Previous Orientation") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  facet_wrap(~condition) +
  scale_color_manual(values = c('blue', 'orange'))


#============= 3.2.2. Accuracy =================
ggplot(data = summaryDataAll,
       aes(x = foreperiod,
           y = meanAcc,
           color = oneBackFP)) +
  labs(title = "Accuracy (proportion of correct responses)",
       x = "Foreperiod",
       y = "Mean Acc",
       color = "FP n-1") +
  stat_summary(fun = "mean", geom = "point", size = 1.8) +
  stat_summary(fun = "mean", geom = "line", aes(group = oneBackFP), linewidth = 0.8) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", linewidth = 0.7, width = 0.05) +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_color_manual(values = c("orange", "blue", "green", "magenta")) +
  facet_wrap(~ condition)
ggsave("./Analysis/Plots/seqeff_by_condition.png",
       width = 7.5,
       height = 5)

ggplot(data = summaryDataAll) +
  geom_histogram(aes(x = meanAcc)) +
  facet_grid(foreperiod ~ condition)

AccAnova <- aov_ez(id = "ID",
                   dv = "meanAcc",
                   data = summaryDataAll,
                   within = c("foreperiod", "condition", "oneBackFP"))

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

#=================== 3.3 Sequential effects with difference between current and previous FP ============
ggplot(data = summaryData2,
       aes(x = numForeperiod,
           y = meanRT,
           color = oneBackFPDiff)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = oneBackFPDiff)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  #scale_color_manual(values = c("orange","blue")) +
  facet_wrap(~condition)

# using RT difference as DV
ggplot(data = summaryData2,
       aes(x = numOneBackFPDiff,
           y = meanSeqEff,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  scale_color_manual(values = c("orange","blue"))


# 2.2.2. Lm with difference between the durations of FPn and FPn-1 as regressor
fpDiffRegression <- lm(meanRT ~ foreperiod * condition * oneBackFPDiff,
                       data = summaryData)
summary(fpDiffRegression)
anova(fpDiffRegression)


ggplot(data = summaryData2,
       aes(x = oneBackFPDiff,
           y = meanSeqEff)) +
  stat_summary(fun = 'mean', geom = 'point')

ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanSeqEff)) +
  stat_summary(fun = 'mean', geom = 'point')

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
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color = seqOri)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = seqOri)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_color_manual(values = c("lightgoldenrod1","indianred2"))

ggplot(data = summaryData2,
       aes(x = orientation,
           y = meanRT,
           color = prevOri)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = prevOri)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_color_manual(values = c("lightgoldenrod1","indianred2"))


# Orientation
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color = orientation)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = orientation)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  scale_color_manual(values = c("deeppink3","chartreuse3"))

# No apparent difference
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


ggplot(data = data2,
       aes(x = trial,
           y = RT,
           color = condition)) +
  #geom_point() +
  geom_line(aes(group = condition)) +
  geom_smooth(method = 'lm') +
  facet_wrap(~ID) +
  scale_color_manual(values = c('orange', 'blue'))

#========================== 5.2. Split anovas by counterbalancing order ===================

fpAnova_ae <- aov_ez(id = "ID",
                     dv = "meanRT",
                     data = summaryData2[summaryData2$counterbalance=='action-external',],
                     within = c("foreperiod", "condition"))

fpAnova_ea <- aov_ez(id = "ID",
                     dv = "meanRT",
                     data = summaryData2[summaryData2$counterbalance=='external-action',],
                     within = c("foreperiod", "condition"))


#===================================================================================================#
#=================================== 6. Analysis with scaled predictors =============================
#===================================================================================================#

fpregression <- lm(meanRT ~ condition * foreperiod, data = summaryData)
summary(fpregression)
anova(fpregression)

logfpregression <- lm(meanRT ~ condition * logFP, data = summaryData)
anova(logfpregression)



# ============ 6.1.2. Fit models separately for external condition (safety check) ====================
logfpAnovaExt <- aov_ez(id = "ID",
                        dv = "meanLogRT",
                        data = summaryData2[summaryData2$condition=='external',],
                        within = c("foreperiod"))

# Effect is significant

logfpAnovaExt2 <- aov_ez(id = "ID",
                         dv = "meanLogRT",
                         data = summaryData2[summaryData2$condition=='external',],
                         within = c("logFP"))

# Effect is significant

#============================== 6.2. Sequential effects ================================================
# Sequential effects (aggregated across conditions)
ggplot(data = summaryData,
       aes(x = foreperiod,
           y = meanRT,
           color=oneBackFP)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", size = 0.8, aes(group=oneBackFP)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)))+
  scale_color_manual(values = c('blue','orange','green', 'pink'))


# Sequential effects (separated by condition)
ggplot(data = summaryData,
       aes(x = foreperiod,
           y = meanRT,
           color = oneBackFP)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", size = 0.8, aes(group = oneBackFP)) +
  stat_summary(fun.data = "mean_cl_boot", size = 0.8, width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  facet_wrap(~condition) +
  scale_color_manual(values = c('blue','orange','green', 'pink'))


# 2.2.1. Anova for FP n-1
seqEffAnova <- aov_ez(id = "ID",
                      dv = "meanRT",
                      data = summaryData2,
                      within = c("foreperiod", "condition", "oneBackFP"))

seqEffAnovaInv <- aov_ez(id = "ID",
                      dv = "meanInvRT",
                      data = summaryData,
                      within = c("foreperiod", "condition", "oneBackFP"))

nice(seqEffAnova,
     correction='none')

seqEffAnovaLog <- aov_ez(id = "ID",
                      dv = "meanLogRT",
                      data = summaryData,
                      within = c("foreperiod", "condition", "oneBackFP"))

nice(seqEffAnovaLog,
     correction='none')

seqEFffregression <- lm(meanRT ~ foreperiod * oneBackFP * condition, 
                        data = summaryData)
summary(seqEFffregression)
anova(seqEFffregression)

logseqEffregression <- lm(meanRT~logFP*logoneBackFP*condition,
                          data=summaryData) 
anova(logseqEffregression)

# 2.2.3. Anova for FP n-2
ntworegression <- lm(meanRT ~ foreperiod * oneBackFP * twoBackFP, 
                     data = summaryData)
summary(ntworegression)
anova(ntworegression)
Anova(ntworegression, type = "II")

# ============ 6.2.2. Fit models separately for external condition (safety check) ====================
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

#=================== 6.3 Sequential effects with difference between current and previous FP ============
ggplot(data = summaryData2,
       aes(x = numForeperiod,
           y = meanRT,
           color = oneBackFPDiff)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = oneBackFPDiff)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  #scale_color_manual(values = c("orange","blue")) +
  facet_wrap(~condition)

# using RT difference as DV
ggplot(data = summaryData2,
       aes(x = numOneBackFPDiff,
           y = meanSeqEff,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  scale_color_manual(values = c("orange","blue"))


# Lm with difference between the durations of FPn and FPn-1 as regressor
fpDiffRegression <- lm(meanRT ~ foreperiod * condition * oneBackFPDiff,
                       data = summaryData)
summary(fpDiffRegression)
anova(fpDiffRegression)


ggplot(data = summaryData2,
       aes(x = oneBackFPDiff,
           y = meanSeqEff)) +
  stat_summary(fun = 'mean', geom = 'point')

ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanSeqEff)) +
  stat_summary(fun = 'mean', geom = 'point')

# ============ 6.3.2. Fit models separately for external condition (safety check) ====================
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

#=============================== 6.4 Orientation ===================================
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color = seqOri)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = seqOri)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_color_manual(values = c("lightgoldenrod1","indianred2"))

ggplot(data = summaryData2,
       aes(x = orientation,
           y = meanRT,
           color = prevOri)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = prevOri)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_color_manual(values = c("lightgoldenrod1","indianred2"))


# Orientation
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color = orientation)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = orientation)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  scale_color_manual(values = c("deeppink3","chartreuse3"))

# No apparent difference
oriAnova <- aov_ez(id = 'ID',
                   dv = 'meanRT',
                   data = summaryData2,
                   within = c('orientation', 'condition', 'foreperiod'))

#============================= 6.5. Quadratic effects  =====================================
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

#============================= 7. Residual analysis for sequential effects ==========================

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