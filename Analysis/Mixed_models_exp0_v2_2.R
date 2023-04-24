#================================================================================================================#
# Changes:

# To choose the dependent variable, we fit models using only random intercepts instead of the full
# random-effects structure, as suggested by Salet et al. (2022) in their model structure Rmd
# Removed old code in settings contrasts section
#================================================================================================================#

# Load necessary packages
library(magrittr)
library(plyr)
library(tidyverse)
library(lattice)
library(gridExtra)
library(afex)
library(emmeans)
library(lme4)
library(effects)
library(MuMIn)
library(car)
library(data.table)
library(buildmer)
library(janitor)
library(broom.mixed)
library(knitr)
library(brms)
library(bayestestR)
library(BayesFactor)

# Save defaults
graphical_defaults <- par()
options_defaults <- options() 

setwd('G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_0_v2')

source('G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_0_v2/Analysis/Prepare_data_exp0_v2_2.R')


# Functions
hist_resid <- function(M,ptitle='Residuals') {
  d <- data.frame(resid=residuals(M)) 
  d  %>% ggplot(aes(x=resid)) + 
    geom_histogram(aes(y=..density..), bins=75, color='black', fill='grey') + 
    geom_density(color='darkred') + 
    ggtitle(ptitle) -> pl
  return(pl)
}

fitstats = function(M,mname='M') {
  QQ<-qqnorm(residuals(M), plot.it=FALSE)
  R2qq <- cor(QQ$x,QQ$y)^2
  dfqq = data.frame(stat='R2qq', V1=R2qq)
  r2tab <- r.squaredGLMM(M)  %>% 
    t  %>% as.data.frame  %>% rownames_to_column(var='stat')  %>% 
    rbind(.,dfqq)
  r2tab$stat = c("$R^2_m$","$R^2_c$",'$R^2_{qq}$' )
  colnames(r2tab) <- c('stat',mname)
  return(r2tab)
}

#==========================================================================================#
#================================= 1. Explore individual data ==============================
#==========================================================================================#

# Plot RT by FP by participant and model using complete pooling 
FPfitAll=lm(meanRT ~ foreperiod,
            data=summaryData)

fit.params=tidy(FPfitAll)

summary(FPfitAll)


ggplot(data=summaryData,
       aes(x=foreperiod,
           y=meanRT)) +
  stat_summary(fun="mean", geom="point", size=1.5)+
  geom_abline(intercept=fit.params$estimate[1],
              slope=fit.params$estimate[2],
              color="blue")+
  facet_wrap(~ ID, ncol=6)


# Plot RT by FP by participant and model using individual data (no pooling)
dataGroupedByRT <- summaryData %>% 
  group_by(ID,foreperiod) %>% 
  summarise(meanRT=mean(meanRT)) %>%
  ungroup() %>%
  mutate(numForeperiod=as.numeric(as.character(foreperiod)))

data.no_pooling <- dataGroupedByRT %>%
  select(-foreperiod) %>%
  group_by(ID) %>%
  nest(data = c(numForeperiod, meanRT)) %>%
  mutate(fit = map(data, ~ lm(meanRT ~ numForeperiod, data = .)),
         params = map(fit, tidy)) %>%
  ungroup() %>%
  unnest(c(params)) %>%
  select(ID, term, estimate) %>%
  complete(ID, term, fill = list(estimate = 0)) %>%
  pivot_wider(names_from = term,
              values_from = estimate) %>% 
  clean_names()


data.no_pooling <- data.no_pooling %>%
  rename(ID=id,
         numForeperiod=num_foreperiod)


ggplot(data = dataGroupedByRT,
       aes(x = numForeperiod, y = meanRT)) + 
  geom_abline(data = data.no_pooling,
              aes(intercept = intercept,
                  slope = numForeperiod),
              color = "blue") +
  geom_point() +
  facet_wrap(~ID, ncol=6) + 
  scale_x_continuous(breaks = 0:4 * 2) +
  theme(strip.text = element_text(size = 12),
        axis.text.y = element_text(size = 12))

fp_no_pooling <- data.no_pooling$numForeperiod

# Compare results to see how much they differ
data_grouped_by_fp <- data %>%
  group_by(ID, foreperiod) %>%
  summarise(meanRT=mean(meanRT)) %>%
  ungroup()

fit_partial_pooling <- lmer(formula = RT ~ numForeperiod + 
                               (1 + numForeperiod|ID),
                             data = data)

data_partial_pooling <- fit_partial_pooling %>%
  augment() %>%
  select(ID, numFore, RT, .fitted) %>%
  rename(fitted=.fitted)

#================ 2.2. FP n-1 ==================
dataGroupedByRT <- summaryData %>% 
  group_by(ID,oneBackFP) %>% 
  summarise(meanRT=mean(meanRT)) %>%
  ungroup() %>%
  mutate(numOneBackFP = as.numeric(as.character(oneBackFP)))

data.no_pooling <- dataGroupedByRT %>%
  select(-oneBackFP) %>%
  group_by(ID) %>%
  nest(data = c(numOneBackFP, meanRT)) %>%
  mutate(fit = map(data, ~ lm(meanRT ~ numOneBackFP, data = .)),
         params = map(fit, tidy)) %>%
  ungroup() %>%
  unnest(c(params)) %>%
  select(ID, term, estimate) %>%
  complete(ID, term, fill = list(estimate = 0)) %>%
  pivot_wider(names_from = term,
              values_from = estimate) %>% 
  clean_names()


data.no_pooling <- data.no_pooling %>%
  rename(ID=id,
         numOneBackFP=num_one_back_fp)


ggplot(data = dataGroupedByRT,
       aes(x = numOneBackFP, y = meanRT)) + 
  geom_abline(data = data.no_pooling,
              aes(intercept = intercept,
                  slope = numOneBackFP),
              color = "blue") +
  geom_point() +
  facet_wrap(~ID, ncol=6) + 
  scale_x_continuous(breaks = 0:4 * 2) +
  theme(strip.text = element_text(size = 12),
        axis.text.y = element_text(size = 12))


#==========================================================================================#
#====================================== 2. Prepare model ====================================
#==========================================================================================#

# Set contrasts
contrasts(data$foreperiod) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)
contrasts(data$condition) <- c(-1/2, 1/2)
contrasts(data$prevOri) <- c(-1/2, 1/2)
contrasts(data$oneBackFP) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)

contrasts(data2$foreperiod) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)
contrasts(data2$condition) <- c(-1/2, 1/2)
contrasts(data2$prevOri) <- c(-1/2, 1/2)
contrasts(data2$oneBackFP) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)


#=========================== 2.1. Choose dependent variable =============================
# We use the strategy of keeping it maximal to find a model that converges and progressively
# remove terms, one of the strategies recommended to avoid overfitting:
# https://rdrr.io/cran/lme4/man/isSingular.html

#  Trying to fit the model using foreperiod and FPn-1 as factors results
# in R hanging during execution; for this reason, we use them as
# numerical variables

# To choose the data transformation that leads to the optimal random effects structure, we fit models including only 
# random intercepts and compare R2 and residuals

# Fit models with RT and inverse RT without trimming
fplmm1 <- mixed(formula = RT ~ numForeperiod*condition*numOneBackFP + 
                  (1|ID),
                data = data,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE,
                return = "merMod")


# Now we run the same model with inverse RT and logRT as outcomes
invfplmm1 <- mixed(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                     (1|ID),
                   data = data,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

logfplmm1 <- mixed(formula = logRT ~ numForeperiod*condition*numOneBackFP + 
                     (1|ID),
                   data = data,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

# Let's check that the current structure does not provide a singular fit:
isSingular(fplmm1)
isSingular(invfplmm1)
isSingular(logfplmm1)

# None of them return singular fits!

# Amount of variance accounted for by the model
var <- data.frame(dataset = 'no trim',
                  'RT' = cor(fitted(fplmm1), data$RT)^2,
                  'invRT' = cor(fitted(invfplmm1), data$invRT)^2,
                  'logRT' = cor(fitted(logfplmm1), data$logRT)^2)


# Check normality of residuals
par(mfrow=c(1,3))
qqnorm(resid(fplmm1),
       main="Normal q-qplot fplmm1")
qqnorm(resid(invfplmm1),
       main="Normal q-qplot invfplmm1")
qqnorm(resid(logfplmm1),
       main="Normal q-qplot logfplmm1")

par(graphical_defaults)

# All models show departures from normality, although this is slightly less for log RT

# Plot residuals
par(mfrow=c(3,1))
plot(resid(fplmm1), fitted(fplmm1),
     main="Residuals fplmm1")
plot(resid(invfplmm1), fitted(invfplmm1),
     main="Residuals invfplmm1")
plot(resid(logfplmm1), fitted(logfplmm1),
     main="Residuals logfplmm1")

par(graphical_defaults)
# There are no clear correlations

# Residual histograms
grid.arrange(hist_resid(fplmm1, 'RT'),
             hist_resid(invfplmm1, '1/RT'),
             hist_resid(logfplmm1, 'logRT'),
             ncol=1)

# All appear to be relatively normally distributed, with a slight positive skew

# Fit models with RT and inverse RT without trimming
trimfplmm1 <- mixed(formula = RT ~ numForeperiod*condition*numOneBackFP + 
                  (1|ID),
                data = data2,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE,
                return = "merMod")


# Now we run the same model with inverse RT and logRT as outcomes
triminvfplmm1 <- mixed(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                     (1|ID),
                   data = data2,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

# Now we run the same model with inverse RT and logRT as outcomes
trimlogfplmm1 <- mixed(formula = logRT ~ numForeperiod*condition*numOneBackFP + 
                     (1|ID),
                   data = data2,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

isSingular(trimfplmm1)
isSingular(triminvfplmm1)
isSingular(trimlogfplmm1)

# No singular fits here either

# Amount of variance accounted for by the model
var <- rbind(var,
             data.frame(dataset = 'trim',
             'RT' = cor(fitted(trimfplmm1), data2$RT)^2,
             'invRT' = cor(fitted(triminvfplmm1), data2$invRT)^2,
             'logRT'= cor(fitted(trimlogfplmm1), data2$logRT)^2))



# Again, the log model accounts for a larger amount of the variance, although this difference is pretty small

# check normality
par(mfrow=c(2,3))
qqnorm(resid(fplmm1),
       main="Normal q-qplot fplmm1")
qqnorm(resid(invfplmm1),
       main="Normal q-qplot invfplmm1")
qqnorm(resid(logfplmm1),
       main="Normal q-qplot logfplmm1")
qqnorm(resid(trimfplmm1),
       main="Normal q-qplot trimfplmm")
qqnorm(resid(triminvfplmm1),
       main="Normal q-qplot triminvfplmm")
qqnorm(resid(trimlogfplmm1),
       main="Normal q-qplot trimlogfplmm")
par(graphical_defaults)

# The log plots show the least deviations from normality, and trimming does seem to help a bit

# Plot residuals
par(mfrow=c(3,2))
plot(resid(fplmm1), fitted(fplmm1),
     main="Residuals RT")
plot(resid(invfplmm1), fitted(invfplmm1),
     main="Residuals 1/RT")
plot(resid(logfplmm1), fitted(logfplmm1),
     main="Residuals logRT")
plot(resid(trimfplmm1), fitted(trimfplmm1),
     main="Residuals trimmed RT")
plot(resid(triminvfplmm1), fitted(triminvfplmm1),
     main="Residuals trimmed 1/RT")
plot(resid(trimlogfplmm1), fitted(trimlogfplmm1),
     main="Residuals trimmed logRT")
par(graphical_defaults)

# LIttle difference after trimming, although outliers are rarer; logRT still appears to perform better
grid.arrange(hist_resid(fplmm1, 'RT'),
             hist_resid(invfplmm1, '1/RT'),
             hist_resid(logfplmm1, 'logRT'),
             hist_resid(trimfplmm1, 'trimmed RT'),
             hist_resid(triminvfplmm1, 'trimmed 1/RT'),
             hist_resid(trimlogfplmm1, 'trimmed logRT'),
             ncol=2,
             as.table=FALSE)

# The same pattern is apparent as in the model with no trimming

R2table <- fitstats(fplmm1, 'RT') %>%
  plyr::join(., fitstats(invfplmm1, '1/(RT)'), by='stat') %>%
  plyr::join(., fitstats(logfplmm1, 'log(RT)'), by='stat') %>%
  plyr::join(., fitstats(trimfplmm1, 'trim RT'), by='stat') %>%
  plyr::join(., fitstats(triminvfplmm1, 'trim 1/(RT)'), by='stat') %>%
  plyr::join(., fitstats(trimlogfplmm1, 'trim log(RT)'), by='stat') %>%
  kable(digits=4)

# Indeed, the model using trimmed logRT performs best according to fit. All R2 are very low 

# Use z-scores for centering with no outlier trimming
fplmm3 <- mixed(formula = RTzscore ~ numForeperiod*condition*numOneBackFP + 
                  (1|ID),
                data = data,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE,
                return = "merMod")

isSingular(fplmm3)

# Amount of variance accounted for by the model
cor(fitted(fplmm3), data$RT)^2

# Very low variance explained

# check normality
qqnorm(resid(fplmm3),
       main="Normal q-qplot fplmm3")

# Far from normality

# Plot residuals
plot(fplmm3, resid(.) ~ fitted(.),
     main="Residuals fplmm3")

# Correlation between residuals and fitted values seems to dissappear

qplot(resid(fplmm3),
      main="Residuals fplmm3")

# Residuals are quite assymnetric

# Use z-scores for centering with trimming
fplmm4 <- mixed(formula = RTzscore ~ numForeperiod*condition*numOneBackFP + 
                  (1|ID),
                data = data2,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE,
                return = "merMod")

isSingular(fplmm4)

# Singular fit

# Amount of variance accounted for by the model
cor(fitted(fplmm4), data2$RT)^2

# check normality
qqnorm(resid(fplmm4),
       main="Normal q-qplot fplmm4")

# Plot residuals
plot(fplmm4, resid(.) ~ fitted(.),
     main="Residuals fplmm4")

qplot(resid(fplmm4),
      main="Residuals fplmm4")

# The models with z-scores for RT performs horribly


# Variance explained is higher with trimming, although this is probably not significant
# Q-q plots are better with trimming
# Residuals are less correlated with fitted values with trimming
# Residuals are more normally distributed with trimming

# The model with inverse RT performs better than the one with RT

#==============================================================================================#
#============================== 3. Find random effects structure ===============================
#==============================================================================================#
options(scipen=1)

# For this, we use buildmer
trimfplmm1 <- buildmer(RT ~ numForeperiod * condition * numOneBackFP + 
                            (1+numForeperiod*condition*numOneBackFP|ID), 
                          data=data2,
                          buildmerControl = list(#direction='backward',
                                                 crit='LRT',#ddf = "Satterthwaite",
                                                 family=gaussian(link = 'identity'),
                                                 calc.anova = TRUE))

formula(trimfplmm1)

triminvfplmm1 <- buildmer(RT ~ numForeperiod * condition * numOneBackFP + 
                            (1+numForeperiod*condition*numOneBackFP|ID), 
                          data=data2,
                          buildmerControl = list(#direction='backward',
                                                 crit='LRT',
                                                 family=gaussian(link = 'inverse'),
                                                 calc.anova = TRUE))

formula(triminvfplmm1)
summary(triminvfplmm1)

trimlogfplmm1 <- buildmer(logRT ~ numForeperiod * condition * numOneBackFP + 
                            (1+numForeperiod*condition*numOneBackFP|ID), 
                          data=data2,
                          buildmerControl = list(direction='backward',
                                                 crit='LRT',
                                                 family=gaussian(link='identity'),
                                                 calc.anova = TRUE))
                                                            
formula(trimlogfplmm1)
formula(fulltrimlogfplmm1)

summary(trimlogfplmm1)
isSingular(trimlogfplmm1)

options(options_defaults)


# Visualize effects
trimlogfplmm1 <- mixed(logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                         numOneBackFP + numForeperiod:numOneBackFP + (1 + condition | ID),
                       data=data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod")

fp_effect <- effect(term='numForeperiod', mod=trimlogfplmm1)
summary(fp_effect)

fp_effect_df <- as.data.frame(fp_effect)

ggplot() +
  stat_summary(fun='mean', geom='point', data=data2, aes(x=numForeperiod, y=logRT)) +
  geom_line(data=fp_effect_df, aes(x=numForeperiod, y=fit), color='red') +
  geom_ribbon(data=fp_effect_df, aes(x=numForeperiod, ymin=lower, ymax=upper), alpha=0.3) +
  labs(x='Foreperiod (continuous)', y='logRT')

# Systematic comparisons via BIC
trimlogfplmm1 <- mixed(logRT ~ 1 + numForeperiod + condition + numForeperiod:condition + 
                         numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                         (1 + numForeperiod + condition | ID),
                       data=data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")


trimlogfplmm2 <- mixed(logRT ~ 1 + numForeperiod + condition + numForeperiod:condition + 
                         numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                         (1 + numForeperiod + condition || ID),
                       data=data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")

trimlogfplmm3 <- mixed(logRT ~ 1 + numForeperiod + condition + numForeperiod:condition + 
                         numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                         (1 + numForeperiod | ID),
                       data=data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")

trimlogfplmm4 <- mixed(logRT ~ 1 + numForeperiod + condition + numForeperiod:condition + 
                         numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                         (1 + condition | ID),
                       data=data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")

trimlogfplmm5 <- mixed(logRT ~ 1 + numForeperiod + condition + numForeperiod:condition + 
                         numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                         (1 | ID),
                       data=data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")

BIC(trimlogfplmm1, trimlogfplmm2, trimlogfplmm3, trimlogfplmm4, trimlogfplmm5) %>%
  kable()

# Compare result from buildmer and manual result
trimlogfplmm6 <- mixed(logRT ~ 1 + numForeperiod + condition + numForeperiod:condition + 
                         numOneBackFP + numForeperiod:numOneBackFP +
                         (1 + condition| ID),
                       data=data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")

BIC(trimlogfplmm5, trimlogfplmm6) %>%
  kable()

# Compare dependent variables using random-effects structure
fplmm1 <- mixed(formula = RT ~ 1 + numForeperiod + condition + numForeperiod:condition + 
                  numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                  (1 + numForeperiod | ID),
                data = data,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'KR',
                REML=TRUE,
                return = "merMod")


invfplmm1 <- mixed(formula = invRT ~ 1 + numForeperiod + condition + numForeperiod:condition + 
                     numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                     (1 + numForeperiod | ID),
                   data = data,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'KR',
                   REML=TRUE,
                   return = "merMod")

logfplmm1 <- mixed(formula = logRT ~ 1 + numForeperiod + condition + numForeperiod:condition + 
                     numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                     (1 + numForeperiod | ID),
                   data = data,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'KR',
                   REML=TRUE,
                   return = "merMod")

trimfplmm1 <- mixed(formula = RT ~ 1 + numForeperiod + condition + numForeperiod:condition + 
                  numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                  (1 + numForeperiod | ID),
                data = data2,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'KR',
                REML=TRUE,
                return = "merMod")


triminvfplmm1 <- mixed(formula = invRT ~ 1 + numForeperiod + condition + numForeperiod:condition + 
                     numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                     (1 + numForeperiod | ID),
                   data = data2,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'KR',
                   REML=TRUE,
                   return = "merMod")

trimlogfplmm1 <- mixed(formula = logRT ~ 1 + numForeperiod + condition + numForeperiod:condition + 
                     numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                     (1 + numForeperiod | ID),
                   data = data2,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'KR',
                   REML=TRUE,
                   return = "merMod")

# Compare model R2 and residuals
# Amount of variance accounted for by the model
var <- data.frame(dataset = 'no trim',
                  'RT' = cor(fitted(fplmm1), data$RT)^2,
                  'invRT' = cor(fitted(invfplmm1), data$invRT)^2,
                  'logRT' = cor(fitted(logfplmm1), data$logRT)^2)

var <- rbind(var,
             data.frame(dataset = 'trim',
                        'RT' = cor(fitted(trimfplmm1), data2$RT)^2,
                        'invRT' = cor(fitted(triminvfplmm1), data2$invRT)^2,
                        'logRT'= cor(fitted(trimlogfplmm1), data2$logRT)^2))

# Again, the log model accounts (barely) for a larger amount of the variance

# check normality
par(mfrow=c(2,3))
qqnorm(resid(fplmm1),
       main="Normal q-qplot fplmm1")
qqnorm(resid(invfplmm1),
       main="Normal q-qplot invfplmm1")
qqnorm(resid(logfplmm1),
       main="Normal q-qplot logfplmm1")
qqnorm(resid(trimfplmm1),
       main="Normal q-qplot trimfplmm")
qqnorm(resid(triminvfplmm1),
       main="Normal q-qplot triminvfplmm")
qqnorm(resid(trimlogfplmm1),
       main="Normal q-qplot trimlogfplmm")
par(graphical_defaults)

# The log plots show the least deviations from normality, and trimming does seem to help a bit

# Plot residuals
par(mfrow=c(3,2))
plot(resid(fplmm1), fitted(fplmm1),
     main="Residuals RT")
plot(resid(invfplmm1), fitted(invfplmm1),
     main="Residuals 1/RT")
plot(resid(logfplmm1), fitted(logfplmm1),
     main="Residuals logRT")
plot(resid(trimfplmm1), fitted(trimfplmm1),
     main="Residuals trimmed RT")
plot(resid(triminvfplmm1), fitted(triminvfplmm1),
     main="Residuals trimmed 1/RT")
plot(resid(trimlogfplmm1), fitted(trimlogfplmm1),
     main="Residuals trimmed logRT")
par(graphical_defaults)

grid.arrange(hist_resid(fplmm1, 'RT'),
             hist_resid(invfplmm1, '1/RT'),
             hist_resid(logfplmm1, 'logRT'),
             hist_resid(trimfplmm1, 'trimmed RT'),
             hist_resid(triminvfplmm1, 'trimmed 1/RT'),
             hist_resid(trimlogfplmm1, 'trimmed logRT'),
             ncol=2,
             as.table=FALSE)


R2table <- fitstats(fplmm1, 'RT') %>%
  plyr::join(., fitstats(invfplmm1, '1/(RT)'), by='stat') %>%
  plyr::join(., fitstats(logfplmm1, 'log(RT)'), by='stat') %>%
  plyr::join(., fitstats(trimfplmm1, 'trim RT'), by='stat') %>%
  plyr::join(., fitstats(triminvfplmm1, 'trim 1/(RT)'), by='stat') %>%
  plyr::join(., fitstats(trimlogfplmm1, 'trim log(RT)'), by='stat') %>%
  kable(digits=4)


# By most measures, logRT with trimming still yields the best fit

#====================== 3.2. Model with difference between FP and FP n-1 =====================
trimlogfpdifflmm1 <- buildmer(formula = logRT ~ numForeperiod + condition + numOneBackFPDiff + 
                                numForeperiod:condition + condition:numOneBackFPDiff +
                                (1 + numForeperiod + condition + numOneBackFPDiff +
                                   numForeperiod:condition + condition:numOneBackFPDiff | ID), 
                          data=data2,
                          buildmerControl = buildmerControl(ddf = "Satterthwaite",
                                                            calc.anova = TRUE))

formula(trimlogfpdifflmm1)

summary(trimlogfpdifflmm1)

trimlogfpdifflmm1 <- mixed(formula = logRT ~ 1 + condition + numForeperiod + numOneBackFPDiff + condition:numOneBackFPDiff + 
                             (1 + condition | ID),
                    data = data2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    return = "merMod")

trimlogfpdifflmm2 <- update(trimlogfpdifflmm1, formula = ~ . -condition:numOneBackFPDiff)
trimlogfpdifflmm2 <- update(trimlogfpdifflmm2, formula = ~ . -numOneBackFPDiff)

fp_effect <- effect(term='numForeperiod', mod=trimlogfpdifflmm2)
summary(fp_effect)

fp_effect_df <- as.data.frame(fp_effect)

ggplot() +
  stat_summary(fun='mean', geom='point', data=data2, aes(x=numForeperiod, y = logRT)) +
  geom_line(data=fp_effect_df, aes(x=numForeperiod, y=fit), color='red') +
  geom_ribbon(data=fp_effect_df, aes(x=numForeperiod, ymin=lower, ymax=upper), alpha=0.3) +
  labs(x='Foreperiod (continuous)', y='RT')

trimlogfpdiff_modcomp <- anova(trimlogfpdifflmm1, trimlogfpdifflmm2)

bf_trimlogfpdiff <- exp((BIC(trimlogfpdifflmm1)-BIC(trimlogfpdifflmm2))/2)

bf_trimlogfpdiff <- bic_to_bf(c(BIC(trimlogfpdifflmm2), BIC(trimlogfpdifflmm1)),
                              denominator = BIC(trimlogfpdifflmm2))
bf_trimlogfpdiff


#===================================================================================#
# Sanity check: model comparisons without trimming
#===================================================================================#
invfplmm2 <- mixed(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                         (1+numForeperiod+condition||ID),
                       data = data,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod")

anova(invfplmm2, invfplmm1, refit=FALSE)

invfplmm3 <- mixed(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                     (1+numForeperiod+numOneBackFP||ID),
                   data = data,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

anova(invfplmm3, invfplmm1, refit=FALSE)

invfplmm4 <- mixed(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                     (1+condition+numOneBackFP||ID),
                   data = data,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

anova(invfplmm4, invfplmm1, refit=FALSE)

# Amount of variance accounted for by the models
cor(fitted(invfplmm1), data$RT)^2
cor(fitted(invfplmm2), data$invRT)^2
cor(fitted(invfplmm3), data$invRT)^2
cor(fitted(invfplmm4), data$invRT)^2


#==============================================================================================#
#================================== 4. Choose link function ====================================
#==============================================================================================#

#======================== 4.1. Visualize distributions for each variable =======================
# Pooled
RTHistograms <- ggplot(data=data2,
                          aes(x=RT))+
  geom_histogram()
RTHistograms

invRTHistograms <- ggplot(data=data2,
                             aes(x=invRT)) +
  geom_histogram()
invRTHistograms

logRTHistograms <- ggplot(data=data2,
                             aes(x=logRT)) +
  geom_histogram()
logRTHistograms


# By participant
indRTHistograms <- ggplot(data=data2,
                        aes(x=RT))+
  geom_histogram()+
  facet_wrap(~ID)
indRTHistograms

indinvRTHistograms <- ggplot(data=data2,
                             aes(x=invRT)) +
  geom_histogram() +
  facet_wrap(~ID)
indinvRTHistograms

indlogRTHistograms <- ggplot(data=data2,
                             aes(x=logRT)) +
  geom_histogram() +
  facet_wrap(~ID)
indlogRTHistograms

# Try models
trimfpgauss <- mixed(formula = RT ~ 1 + numForeperiod + condition + numForeperiod:condition + 
                         numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                         (1 + numForeperiod | ID),
                       data = data2,
                       family=gaussian(link = "identity"),
                       #control = glmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       return = "merMod")

summary(trimfpgauss)

trimlogfpgauss <- mixed(formula = RT ~ 1 + numForeperiod + condition + numForeperiod:condition + 
                       numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                       (1 + numForeperiod | ID),
                     data = data2,
                     family=gaussian(link = "log"),
                     #control = glmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                     progress = TRUE,
                     expand_re = TRUE,
                     method =  'KR',
                     return = "merMod")

triminvfpgauss <- mixed(formula = logRT ~ 1 + numForeperiod + condition + numForeperiod:condition + 
                             numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                             (1 + numForeperiod | ID),
                           data = data2,
                           family=gaussian(link='inverse'),
                           #control = glmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                           progress = TRUE,
                           expand_re = TRUE,
                           method =  'KR',
                           return = "merMod")

#==============================================================================================#
#============================== 5. Bayesian mixed models using brms ============================
#==============================================================================================#

#======= 3.2.1. Models using effects from previous experiment as priors

# Baseado no exp 0 v1
prior1 <- c(set_prior('normal(3.07,0.06)', class = 'Intercept'),
            set_prior("normal(-0.09, 0.027)", class = "b", coef = "numForeperiod"),
            set_prior("normal(0.007, 0.08)", class = "b", coef = "condition1"),
            set_prior("normal(-0.034, 0.027)", class = "b", coef = "numOneBackFP")
)

# Gaussian distr
b_one_back_fp <- brm(formula = logRT ~ numForeperiod * condition * numOneBackFP + 
                       (1+numForeperiod*condition*numOneBackFP|ID),
                     data = data2,
                     family = gaussian(),
                     prior = prior1,
                     save_all_pars = TRUE,
                     warmup = 2000,
                     iter = 10000)

b_one_back_fp <- readRDS('./Analysis/b_one_back_fp.rds')


ggplot() +
  stat_summary(fun='mean', geom='point', data=data2, aes(x=numForeperiod, y=logRT)) +
  geom_point(data=fp_effect_df, aes(x=numForeperiod, y=fit), color='red') +
  geom_line(data=fp_effect_df, aes(x=numForeperiod, y=fit), color='red') +
  geom_ribbon(data=fp_effect_df, aes(x=numForeperiod, ymin=lower, ymax=upper), alpha=0.3) +
  labs(x='Foreperiod (continuous)', y='logRT')



# Inverse gaussian distr
b_one_back_fp_invgaus <- brm(formula = logRT ~ numForeperiod * condition * numOneBackFP + 
                       (1+numForeperiod*condition*numOneBackFP|ID),
                     data = data2,
                     family = gaussian(),
                     prior = prior1,
                     save_all_pars = TRUE,
                     warmup = 2000,
                     iter = 10000)


b_one_back_fp_null <- brm(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                            numOneBackFP + numForeperiod:numOneBackFP + 
                            (1 + condition + numForeperiod + numOneBackFP + 
                               condition:numForeperiod + numForeperiod:numOneBackFP| ID),
                     data = data2,
                     family = gaussian(),
                     prior = prior1,
                     save_all_pars = TRUE,
                     warmup = 2000,
                     iter = 10000,
                     file="b_one_back_fp_null.rds")



bf_one_back_brm <- bayes_factor(b_one_back_fp, b_one_back_fp_null)


#================================= 3.2.2 Using a vague prior ================================
b_one_back_fp_vagueprior <- brm(formula = logRT ~ numForeperiod * condition * numOneBackFP + 
                                  (1+numForeperiod*condition*numOneBackFP|ID),
                                data = data2,
                                family = gaussian(),
                                save_all_pars = TRUE,
                                warmup = 2000,
                                iter = 10000)


options(digits = 5)
summary(b_one_back_fp)
options(options_defaults)














fpDifflmm2 <- mixed(formula = invRT ~ foreperiod * condition * oneBackFPDiff + 
                      (1+foreperiod*condition*oneBackFPDiff|ID),
                    data = data2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    return = "merMod")

anova(fpDifflmm2)

fpDifflmm3 <- mixed(formula = invRT ~ foreperiod * condition * oneBackFPDiff + 
                      (1+foreperiod*condition*oneBackFPDiff|ID),
                    data = data2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    return = "merMod")

fpDifflmm4 <- mixed(formula = invRT ~ foreperiod * condition * oneBackFPDiff + 
                      (1+foreperiod*condition*oneBackFPDiff||ID),
                    data = data2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    return = "merMod")

fpDifflmm5 <- mixed(formula = invRT ~ foreperiod * condition * oneBackFPDiff + 
                      (1+foreperiod+condition+oneBackFPDiff||ID),
                    data = data2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    return = "merMod")

cor(fitted(fpDifflmm2), data2$invRT)^2
cor(fitted(fpDifflmm4), data2$invRT)^2
cor(fitted(fpDifflmm5), data2$invRT)^2


#============================================================================================================#
#=========================== Difference between the durations of FP n and FP n-1 as predictor ================
#============================================================================================================#

# Fullest model using foreperiod and one back fp as categorical variables
fpDifflmm1 <- mixed(formula = invRT ~ condition * oneBackFPDiff + 
                      (1+condition*oneBackFPDiff|ID),
                    data = data2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    check_contrasts=FALSE,
                    return = "merMod")

# Remove correlations
fpDifflmm2 <- mixed(formula = invRT ~ foreperiod * condition * oneBackFPDiff + 
                      (1+foreperiod*condition*oneBackFPDiff||ID),
                    data = data2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    check_contrasts=FALSE,
                    return = "merMod")

# Remove interactions
fpDifflmm3 <- mixed(formula = invRT ~ foreperiod * condition * oneBackFPDiff + 
                      (1+foreperiod+condition+oneBackFPDiff||ID),
                    data = data2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    check_contrasts=FALSE,
                    return = "merMod")

# Find optimal structure using buildmer
fpDifflmm <- buildmer(invRT ~ foreperiod * condition * oneBackFPDiff + 
                        (1+foreperiod*condition*oneBackFPDiff|ID), 
                      data=data2,
                      buildmerControl = buildmerControl(include = ~ foreperiod * condition * oneBackFPDiff, calc.anova = TRUE, ddf = "Satterthwaite"))


