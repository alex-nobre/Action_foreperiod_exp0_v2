

#================================================================================#
# Changes
# Manually set contrasts for anova
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

setwd('G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_0_v2')

#==========================================================================================#
#======================================= 0. Prepare data ===================================
#==========================================================================================#
# Functions for raincloud plots
#source("./Analysis/R_rainclouds.R")

source('./Analysis/Prepare_data_exp0_v2_5.R')

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
ggplot(data=summaryData,
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

# Check for influence of external fixation duration
ggplot(data=filter(data,condition=='external'),
       aes(x=extFixationDuration,
           y=RT,
           color=foreperiod))+
  geom_point() +
  facet_wrap(~foreperiod)


# Check for influence of total fixation duration
data <- data %>%
  mutate(totalFix = extFixationDuration + numForeperiod)

ggplot(data = filter(data, condition == 'external'),
       aes(x = totalFix,
           y = RT)) +
  geom_point()


# Check for influence of latency of action key press on RT
ggplot(data=filter(data,condition=='action'),
       aes(x=action_trigger.rt,
           y=RT,
           color=foreperiod))+
  geom_point() +
  facet_wrap(~foreperiod)



# Plot data by foreperiod only to ascertain that there is an FP effect
ggplot(data = summaryData2,
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


# Fit by participant
ggplot(data = data2,
       aes(x = foreperiod,
           y = RT)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", size = 1, aes(group=1)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.0)),
        axis.title = element_text(size = rel(1.0))) +
  facet_wrap(~ID) +
  scale_color_manual(values = c("orange","blue"))


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

interact_bfs <- sapply(srange, function(range) {
  extractBF(ttestBF(x = external_fits$`numForeperiod:numOneBackFP`[1:range],
                    y = action_fits$`numForeperiod:numOneBackFP`[1:range],
                    paired=TRUE),
            onlybf = TRUE)
})

plot(srange, interact_bfs)

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

#==================== 3.1. FP x RT by condition ======================
# Lines by condition
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
  scale_color_manual(values = c("orange","blue"))

# using LogRT
ggplot(data = summaryData,
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


# Individual panels by participant
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



# ============ 3.1.2. Fit models separately for external condition (safety check) ====================
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

#============================== 3.2. Sequential effects ================================================
# Sequential effects (aggregated across conditions)
ggplot(data = summaryData,
       aes(x = foreperiod,
           y = meanRT,
           color=oneBackFP)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group=oneBackFP)) +
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

seqEffAnova <- aov_ez(id = "ID",
                      dv = "meanInvRT",
                      data = summaryData,
                      within = c("foreperiod", "condition", "oneBackFP"))

nice(seqEffAnova,
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







