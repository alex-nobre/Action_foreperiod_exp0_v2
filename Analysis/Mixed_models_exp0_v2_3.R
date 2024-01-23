#================================================================================================================#
# Changes:

# To choose the dependent variable, we fit models using only random intercepts instead of the full
# random-effects structure, as suggested by Salet et al. (2022) in their model structure Rmd
# Removed old code in settings contrasts section
#================================================================================================================#

# Load packages

# Data processing and plotting
library(magrittr)
library(tidyverse)
library(lattice)
library(gridExtra)
library(data.table)

# Simple models
library(car)
library(janitor)

# Mixed-effects modeling
library(afex)
library(emmeans)
library(lme4)
library(MuMIn)
library(buildmer)
library(broom.mixed)
library(marginaleffects)

# Bayesian models
library(brms)
library(bayestestR)
library(BayesFactor)

# Assess models and results
library(effects)
library(ggeffects)
library(performance)
library(knitr)
library(kableExtra)
library(sjPlot)
library(prediction)


# Save defaults
graphical_defaults <- par()
options_defaults <- options()

# emm options
emm_options(lmer.df = "satterthwaite", lmerTest.limit = 12000)


#================================== 0. Read data ================================
# Create dataset
source('./Analysis/Prepare_data_exp0_v2_6.R')

# Set contrasts
contrasts(data$foreperiod) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)
contrasts(data$logFP) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)
contrasts(data$condition) <- c(-1/2, 1/2)
contrasts(data$prevOri) <- c(-1/2, 1/2)
contrasts(data$oneBackFP) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)

contrasts(data2$foreperiod) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)
contrasts(data2$logFP) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)
contrasts(data2$condition) <- c(-1/2, 1/2)
contrasts(data2$prevOri) <- c(-1/2, 1/2)
contrasts(data2$oneBackFP) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)


#==========================================================================================#
#================================= 1. Explore individual data ==============================
#==========================================================================================#

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
options(scipen = 999)

fplmm1 <- mixed(formula = RT ~ numForeperiod*condition*numOneBackFP + 
                  (1|ID),
                data = data,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE,
                return = "merMod")

summary(fplmm1)


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


var
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



# Variance explained is higher with trimming, although this is probably not significant
# Q-q plots are better with trimming
# Residuals are less correlated with fitted values with trimming
# Residuals are more normally distributed with trimming

# The model with inverse RT performs better than the one with RT

#================================= 2.2. Find random effects structure ==========================

#=========================== 2.2.1. FP and FP n-1 as numerical ============================
trimlogfplmm1 <- buildmer(logRT ~ numForeperiod * condition * numOneBackFP + 
                            (1+numForeperiod*condition*numOneBackFP|ID), 
                          data=data2,
                          buildmerControl = list(direction='backward',
                                                 crit='LRT',#ddf = "Satterthwaite",
                                                 family=gaussian(link = 'identity'),
                                                 calc.anova = TRUE))

isSingular(trimlogfplmm1)
formula(trimlogfplmm1)
summary(trimlogfplmm1)
# Systematic comparisons betweeen lmm's via BIC

# Model obtained with buildmer
trimlogfplmm1 <- mixed(logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                         numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                         (1 + condition | ID),
                       data=data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")

# Random-intercept only model
trimlogfplmm1v2 <- mixed(logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                         numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                         (1 | ID),
                       data=data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")

# Uncorrelated intercepts and slopes
trimlogfplmm1v3 <- mixed(logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                         numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                         (1 | ID) + (0+condition | ID),
                       data=data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")

# Compare BICs and AICs
BIC(trimlogfplmm1, trimlogfplmm1v2, trimlogfplmm1v3) %>%
  kable()

AIC(trimlogfplmm1, trimlogfplmm1v2, trimlogfplmm1v3) %>%
  kable()

cor(fitted(trimlogfplmm1), data2$logRT)^2
cor(fitted(trimlogfplmm2), data2$logRT)^2
cor(fitted(trimlogfplmm3), data2$logRT)^2


#===================== 2.2.2. Using FP n-1 as categorical (for emm comparisons) ===================
trimlogfplmm2 <- buildmer(logRT ~ numForeperiod * condition * oneBackFP + 
                            (1+numForeperiod*condition*oneBackFP|ID), 
                          data=data2,
                          buildmerControl = list(direction='backward',
                                                 crit='LRT',#ddf = "Satterthwaite",
                                                 family=gaussian(link = 'identity'),
                                                 calc.anova = TRUE))

isSingular(trimlogfplmm2)
formula(trimlogfplmm2)


#==========================  2.2.3. Both FP and FP n-1 as categorical ===============================
trimlogfplmm3 <- buildmer(logRT ~ foreperiod * condition * oneBackFP + 
                            (1+foreperiod*condition*oneBackFP|ID), 
                          data=data2,
                          buildmerControl = list(crit='LRT',#ddf = "Satterthwaite",
                                                 family=gaussian(link = 'identity'),
                                                 calc.anova = TRUE))

isSingular(trimlogfplmm3)
formula(trimlogfplmm3)

#==============================================================================================#
#==================================== 3. Model assessment ======================================
#==============================================================================================#

#=============================== 3.1. FP and FP n-1 as numerical ==================================
trimlogfplmm1 <- mixed(logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                         numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + condition:numForeperiod:numOneBackFP +
                         (1 + condition | ID),
                       data=data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod",
                       check_contrasts = FALSE)

summary(trimlogfplmm1)
anova(trimlogfplmm1)

fp_effect <- effect(term='numForeperiod', mod=trimlogfplmm1)
summary(fp_effect)

fp_effect_df <- as.data.frame(fp_effect)

ggplot() +
  stat_summary(fun='mean', geom='point', data=data2, aes(x=numForeperiod, y=logRT)) +
  geom_line(data=fp_effect_df, aes(x=numForeperiod, y=fit), color='red') +
  geom_ribbon(data=fp_effect_df, aes(x=numForeperiod, ymin=lower, ymax=upper), alpha=0.3) +
  labs(x='Foreperiod (continuous)', y='logRT')

# Use RT
trimfplmm1 <- mixed(RT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                      numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                      condition:numForeperiod:numOneBackFP +
                      (1 + condition | ID),
                    data=data2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'KR',
                    REML=TRUE,
                    return = "merMod",
                    check_contrasts = FALSE)

summary(trimfplmm1)
anova(trimfplmm1)

# Pairwise comparisons
Fp_by_Previous=emtrends(trimfplmm1, "numOneBackFP", var = "numForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")

Fp_by_Previous=emtrends(trimfplmm1, c("condition", "numOneBackFP"), var = "numForeperiod")
emtrends(trimfplmm1, by = "numOneBackFP", var = "numForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "none")

slopes(trimfplmm1)

#====================== 3.2. Using FP n-1 as categorical for emm comparisons ==========================
#============= 3.2.1. Using RT ===============
trimfplmm <- mixed(RT ~ 1 + scaledNumForeperiod + condition + scaledNumForeperiod:condition + 
                     oneBackFP + scaledNumForeperiod:oneBackFP + scaledNumForeperiod:condition:oneBackFP +
                     (1 + condition | ID),
                   data=data2,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'KR',
                   REML=TRUE,
                   return = "merMod")

anova(trimfplmm)

# Pairwise comparisons
Fp_by_Previous=emtrends(trimfplmm, "oneBackFP", var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")

Fp_by_Previous=emtrends(trimfplmm, c("condition", "oneBackFP"), var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "none")

#============= 3.2.2. Using logRT ===============
trimlogfplmm <- mixed(logRT ~ 1 + scaledNumForeperiod + condition + scaledNumForeperiod:condition + 
                        oneBackFP + scaledNumForeperiod:oneBackFP + condition:oneBackFP +
                        scaledNumForeperiod:condition:oneBackFP + 
                        (1 + condition + scaledNumForeperiod + scaledNumForeperiod:condition | ID),
                      data=data2,
                      control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                      progress = TRUE,
                      expand_re = FALSE,
                      method =  'KR',
                      REML=TRUE,
                      return = "merMod")

anova(trimlogfplmm)

#Visualize random effects
dotplot(ranef(trimlogfplmm, condVar = TRUE))

# Visualize interactions
emmip(trimlogfplmm,
      oneBackFP ~ condition, style = "factor") # Using averaged FP

emmip(ref_grid(trimlogfplmm, at = list(scaledNumForeperiod = c(1.0, 1.6, 2.2, 2.8))),
      oneBackFP * scaledNumForeperiod~ condition, style = "factor") # Split by FP (not very informative)

# Single slopes tests
fp_by_condition <- slopes(trimlogfplmm, by = "condition", variables = "scaledNumForeperiod",
                          p_adjust = "holm")

test(emtrends(trimlogfplmm, ~ condition, var="scaledNumForeperiod")) # equivalent to slopes

fp_by_oneback <- slopes(trimlogfplmm, by = "oneBackFP", variables = "scaledNumForeperiod",
                        p_adjust = "holm")

test(emtrends(trimlogfplmm, ~ oneBackFP, var="scaledNumForeperiod")) # equivalent to slopes

threeway_int <- slopes(trimlogfplmm, by = c("oneBackFP", "condition"), variables = "scaledNumForeperiod",
                       p_adjust = "holm")




# Pairwise comparisons
fp_by_condition_comp <- emtrends(trimlogfplmm, "condition", var = "scaledNumForeperiod")
fp_by_condition_comp
update(pairs(fp_by_condition_comp), by = NULL, adjust = "holm")

Fp_by_Previous=emtrends(trimlogfplmm, "oneBackFP", var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")

threeway_int_comp = emtrends(trimlogfplmm, c("condition", "oneBackFP"), var = "scaledNumForeperiod")
threeway_int_comp
update(pairs(threeway_int_comp), by = NULL, adjust = "holm")
pairs(threeway_int_comp, simple = "condition")

# Marginal means for FP by FP n-1 and condition
oneback_by_cond = emmeans(trimlogfplmm, c("condition", "oneBackFP"))
oneback_by_cond = emmeans(trimlogfplmm, ~ oneBackFP * condition)
pairs(oneback_by_cond, simple = "condition")


#========================= 3.3. Both FP and FP n-1 as categorical ===================================

#=============== 3.3.1. using logRT =========================
trimlogfplmm3 <- mixed(logRT ~ 1 + foreperiod + condition + foreperiod:condition + oneBackFP + 
                         foreperiod:oneBackFP + condition:oneBackFP + foreperiod:condition:oneBackFP + 
                         (1 + condition | ID),
                       data=data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod",
                       check_contrasts = FALSE)


summary(trimlogfplmm3)

anova(trimlogfplmm3)

# 3.3.1.1. Pairwise comparisons by FP (estimates consecutive differences)
emmip(trimlogfplmm3, condition ~ oneBackFP|foreperiod, CIs = TRUE)

trimlogfplmm3emm <- emmeans(trimlogfplmm3, ~ oneBackFP * condition|foreperiod)

trimlogfplmm3emm <- emmeans(trimlogfplmm3, pairwise ~ oneBackFP * condition|foreperiod)
contrast(trimlogfplmm3emm[[1]], interaction = c("consec", "consec"), by = "foreperiod", adjust = "mvt")

trimlogfplmm3emm <- emmeans(trimlogfplmm3, pairwise ~ condition*oneBackFP|foreperiod)
contrast(trimlogfplmm3emm[[1]], interaction = c("consec"), by = c("foreperiod", "oneBackFP"), adjust = "holm")
contrast(trimlogfplmm3emm[[1]], interaction = c("consec"), by = c("foreperiod", "condition"), adjust = "holm")

# 3.3.1.2. Pairwise comparisons by by FP n-1 (estimate slopes)
emmip(trimlogfplmm3, condition ~ foreperiod|oneBackFP, CIs = TRUE)

trimlogfplmm3emm <- emmeans(trimlogfplmm3, ~ foreperiod * condition|oneBackFP)

trimlogfplmm3emm <- emmeans(trimlogfplmm3, pairwise ~ foreperiod * condition|oneBackFP)
contrast(trimlogfplmm3emm[[1]], interaction = c("poly", "consec"), by = "oneBackFP", adjust = "mvt")


#============= 3.3.2. Using RT instead of logRT ================
trimfplmm3 <- mixed(RT ~ 1 + foreperiod + condition + foreperiod:condition + oneBackFP + 
                      foreperiod:oneBackFP + condition:oneBackFP + foreperiod:condition:oneBackFP + 
                      (1 + condition | ID),
                    data=data2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'KR',
                    REML=TRUE,
                    return = "merMod",
                    check_contrasts = FALSE)

isSingular(trimfplmm3)
summary(trimfplmm3)
anova(trimfplmm3)

# 3.3.2.1. Pairwise comparisons by FP (estimates consecutive differences)
emmip(trimfplmm3, condition ~ oneBackFP|foreperiod, CIs = TRUE, style = "factor",
      xlab = "FP n-1") +
  labs(title = "RT pairwise comparisons") +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_color_manual(values = c("orange", "blue"))
ggsave("./Analysis/Plots/mixed_model_pairwise.png",
       width = 8.5,
       height = 5.7)

trimfplmm3emm <- emmeans(trimfplmm3, ~ oneBackFP * condition|foreperiod)
trimfplmm3emm <- emmeans(trimfplmm3, pairwise ~ oneBackFP * condition|foreperiod)

contrast(trimfplmm3emm[[1]], interaction = c("consec", "consec"), by = "foreperiod", adjust = "mvt")


# 3.3.2.2. Pairwise comparisons by FP n-1 (estimate slopes)
emmip(trimfplmm3, condition ~ foreperiod|oneBackFP, CIs = TRUE)

trimfplmm3emm <- emmeans(trimfplmm3, ~ foreperiod * condition|oneBackFP)
trimfplmm3emm <- emmeans(trimfplmm3, pairwise ~ foreperiod * condition|oneBackFP)

contrast(trimfplmm3emm[[1]], interaction = c("poly", "consec"), by = "oneBackFP", adjust = "mvt")

# 3.3.2.3. Separate models for each level of FP n
seq1000 <- mixed(RT ~ 1 + oneBackFP + condition + oneBackFP:condition +
                   (1 + condition | ID),
                 data=filter(data2, foreperiod == '1000'),
                 control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                 progress = TRUE,
                 expand_re = TRUE,
                 method =  'KR',
                 REML=TRUE,
                 return = "merMod",
                 check_contrasts = FALSE)

summary(seq1000)

anova(seq1000)

emmeans(seq1000, ~ condition|oneBackFP)
emmip(seq1000, condition ~ oneBackFP)


seq1600 <- mixed(RT ~ 1 + oneBackFP + condition + oneBackFP:condition +
                   (1 + condition | ID),
                 data=filter(data2, foreperiod == '1600'),
                 control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                 progress = TRUE,
                 expand_re = TRUE,
                 method =  'KR',
                 REML=TRUE,
                 return = "merMod",
                 check_contrasts = FALSE)

summary(seq1600)
anova(seq1600)

emmeans(seq1600, pairwise ~ condition|foreperiod)

seq2200 <- mixed(RT ~ 1 + oneBackFP + condition + oneBackFP:condition +
                   (1 + condition | ID),
                 data=filter(data2, foreperiod == '2200'),
                 control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                 progress = TRUE,
                 expand_re = TRUE,
                 method =  'KR',
                 REML=TRUE,
                 return = "merMod",
                 check_contrasts = FALSE)

summary(seq2200)
anova(seq2200)

emmeans(seq2200, pairwise ~ condition|foreperiod)

seq2800 <- mixed(RT ~ 1 + oneBackFP + condition + oneBackFP:condition +
                   (1 + condition | ID),
                 data=filter(data2, foreperiod == '2800'),
                 control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                 progress = TRUE,
                 expand_re = TRUE,
                 method =  'KR',
                 REML=TRUE,
                 return = "merMod",
                 check_contrasts = FALSE)

summary(seq2800)
anova(seq2800)

emmeans(seq2800, pairwise ~ condition|foreperiod)

#================== 3.4. Compare dependent variables using random-effects structure ==================
fplmm1 <- mixed(formula = RT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                  numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                  (1 + condition | ID),
                data = data,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'KR',
                REML=TRUE,
                return = "merMod")


invfplmm1 <- mixed(formula = invRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                     numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                     (1 + condition | ID),
                   data = data,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'KR',
                   REML=TRUE,
                   return = "merMod")

logfplmm1 <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                     numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                     (1 + condition | ID),
                   data = data,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'KR',
                   REML=TRUE,
                   return = "merMod")

trimfplmm1 <- mixed(formula = RT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                      numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                      (1 + condition | ID),
                    data = data2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'KR',
                    REML=TRUE,
                    return = "merMod")


triminvfplmm1 <- mixed(formula = invRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                         numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                         (1 + condition | ID),
                       data = data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")

trimlogfplmm1 <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                         numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                         (1 + condition | ID),
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

var
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

# Compare lm and lmm models
trimlogfplm <- lm(logRT ~ 1 + condition + numForeperiod + condition:numForeperiod +
                    numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP,
                  data = data2)

trimlogfplmmML <- mixed(logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                          numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                          (1 + condition | ID),
                        data=data2,
                        control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method =  'KR',
                        REML=FALSE,
                        return = "merMod",
                        check_contrasts = FALSE)

trimlogfplmmML <- lmer(logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                         numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                         (1 + condition | ID),
                       data=data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       method =  'KR',
                       REML=FALSE)

anova(trimlogfplmmML, trimlogfplm)

#============ 3.5. Hierarchical entry ===============
h_trimlogfplmm1 <- mixed(logRT ~ 1 + numForeperiod + (1 | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method = 'KR',
                         REML=TRUE,
                         return = "merMod")

h_trimlogfplmm2 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFP + (1 | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method =  'KR',
                         REML=TRUE,
                         return = "merMod",
                         check_contrasts = FALSE)

h_trimlogfplmm3 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + (1 | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method =  'KR',
                         REML=TRUE,
                         return = "merMod",
                         check_contrasts = FALSE)

h_trimlogfplmm4 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + condition + (1 | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method =  'KR',
                         REML=TRUE,
                         return = "merMod",
                         check_contrasts = FALSE)

h_trimlogfplmm5 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + condition + 
                           condition:numForeperiod + (1 | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method =  'KR',
                         REML=TRUE,
                         return = "merMod",
                         check_contrasts = FALSE)

h_trimlogfplmm6 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + condition + 
                           condition:numForeperiod + condition:numOneBackFP + (1 | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method =  'KR',
                         REML=TRUE,
                         return = "merMod",
                         check_contrasts = FALSE)

h_trimlogfplmm7 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + condition + 
                           condition:numForeperiod + condition:numOneBackFP + condition:numOneBackFP:numForeperiod + (1 | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method =  'KR',
                         REML=TRUE,
                         return = "merMod",
                         check_contrasts = FALSE)


h_trimlogfplmm8 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + condition + 
                           condition:numForeperiod + condition:numOneBackFP + condition:numOneBackFP:numForeperiod + (1 + condition | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method =  'KR',
                         REML=TRUE,
                         return = "merMod",
                         check_contrasts = FALSE)

anova(h_trimlogfplmm1, h_trimlogfplmm2, h_trimlogfplmm3, h_trimlogfplmm4, 
      h_trimlogfplmm5, h_trimlogfplmm6, h_trimlogfplmm7, h_trimlogfplmm8)


#==================== 3.3. Run model separately for action and external conditions ===================

#======== 3.3.1. External ============#
# FP as numeric
trimlogfplmmext1 <- mixed(logRT ~ 1 + numForeperiod + 
                            numOneBackFP + numForeperiod:numOneBackFP +
                            (1 | ID),
                          data=data2[data2$condition=='external',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmmext1)


# FP as categorical
trimlogfplmmext2 <- mixed(logRT ~ 1 + foreperiod + 
                            numOneBackFP + foreperiod:numOneBackFP +
                            (1 | ID),
                          data=data2[data2$condition=='external',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmmext2)

# logFP as numerical
trimlogfplmmext3 <- mixed(logRT ~ 1 + numLogFP + 
                            numOneBackFP + numLogFP:numOneBackFP +
                            (1 | ID),
                          data=data2[data2$condition=='external',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmmext3)

# logFP as categorical
trimlogfplmmext4 <- mixed(logRT ~ 1 + logFP + 
                            numOneBackFP + logFP:numOneBackFP +
                            (1 | ID),
                          data=data2[data2$condition=='external',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmmext4)

# FP as numerical Including only FP
trimlogfplmmext5 <- mixed(logRT ~ 1 + numForeperiod + 
                            (1 | ID),
                          data=data2[data2$condition=='external',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmmext5)
summary(trimlogfplmmext1)

BIC(trimlogfplmmext1, trimlogfplmmext5)

# Including FP n-1 improves BIC

# FP as categorical Including only FP
trimlogfplmmext6 <- mixed(logRT ~ 1 + foreperiod + 
                            (1 | ID),
                          data=data2[data2$condition=='external',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmmext6)
summary(trimlogfplmmext2)

var(lme4::fixef(trimlogfplmmext6))
var(lme4::fixef(trimlogfplmmext2))

r.squaredGLMM(trimlogfplmmext6)
r.squaredGLMM(trimlogfplmmext2)

BIC(trimlogfplmmext2, trimlogfplmmext6)

# logFP as numerical Including only FP
trimlogfplmmext7 <- mixed(logRT ~ 1 + numLogFP + 
                            (1 | ID),
                          data=data2[data2$condition=='external',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmmext7)
summary(trimlogfplmmext3)

BIC(trimlogfplmmext3, trimlogfplmmext7)

# Including FP n-1 improves BIC

# logFP as categorical Including only FP
trimlogfplmmext8 <- mixed(logRT ~ 1 + logFP + 
                            (1 | ID),
                          data=data2[data2$condition=='external',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmmext8)
summary(trimlogfplmmext4)

BIC(trimlogfplmmext4, trimlogfplmmext8)

# FP as numeric, FP n-1 as categorical
trimlogfplmmext9 <- mixed(logRT ~ 1 + numForeperiod + 
                            oneBackFP + numForeperiod:oneBackFP +
                            (1 | ID),
                          data=data2[data2$condition=='external',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmmext9)
anova(trimlogfplmmext9)

#======== 3.3.2. Action ============#
trimlogfplmm1act <- mixed(logRT ~ 1 + numForeperiod + 
                            numOneBackFP + numForeperiod:numOneBackFP +
                            (1 | ID),
                          data=data2[data2$condition=='action',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmm1act)


trimlogfplmm1 <- mixed(logRT ~ 1 + foreperiod +
                         (1 | ID),
                       data=data2[data2$condition=='external',],
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod",
                       check_contrasts = FALSE)
summary(trimlogfplmm1)

anova(trimlogfplmm1)

#==============================================================================================#
#================================== 4. Choose distribution ====================================
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
trimlogfpgauss <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                          numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                          (1 + condition | ID),
                        data = data2,
                        #family=gaussian(link = "identity"),
                        control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method =  'KR',
                        return = "merMod")

summary(trimlogfpgauss)

trimlogfpinvgauss <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                             numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                             (1 + condition | ID),
                           data = data2,
                           family=inverse.gaussian(link = "identity"),
                           control = lmerControl(optimizer = c("nloptwrap"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                           progress = TRUE,
                           expand_re = TRUE,
                           method =  'KR',
                           return = "merMod")

summary(trimlogfpinvgauss)

trimlogfpgamma <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                          numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                          (1 + condition | ID),
                        data = data2,
                        family=Gamma(link = "identity"),
                        control = lmerControl(optimizer = c("nloptwrap"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method =  'KR',
                        return = "merMod")

summary(trimlogfpgamma)

# Compare visualizations
ggpredict(model = trimlogfpgauss,
          terms = "numForeperiod",
          type = 'fe') %>%
  plot()

ggpredict(model = trimlogfpinvgauss,
          terms = "numForeperiod",
          type = 'fe') %>%
  plot()

ggpredict(model = trimlogfpgamma,
          terms = 'numForeperiod',
          type = 'fe') %>%
  plot()

ggpredict(model = trimlogfpgauss,
          terms = "condition",
          type = 'fe') %>%
  plot()

ggpredict(model = trimlogfpinvgauss,
          terms = "condition",
          type = 'fe') %>%
  plot()

ggpredict(model = trimlogfpgamma,
          terms = 'condition',
          type = 'fe') %>%
  plot()


# Compare performance across models
trimlogfpgauss %>%
  check_model()

trimlogfpinvgauss %>%
  check_model()

trimlogfpgamma %>%
  check_model()

# The gaussian model seems to provide a better fit than the alternatives


#==============================================================================================#
#======================= 5. Model with difference between FP and FP n-1 ========================
#==============================================================================================#

trimlogfpdifflmm1 <- buildmer(formula = logRT ~ numForeperiod * condition * numOneBackFPDiff +
                                (1 + numForeperiod * condition * numOneBackFPDiff | ID), 
                              data=data2,
                              buildmerControl = buildmerControl(ddf = "Satterthwaite",
                                                                calc.anova = TRUE))

formula(trimlogfpdifflmm1)

summary(trimlogfpdifflmm1)

trimlogfpdifflmm1 <- mixed(formula = logRT ~ 1 + condition + numForeperiod + numOneBackFPDiff + numForeperiod:numOneBackFPDiff + 
                             condition:numOneBackFPDiff + (1 + condition | ID),
                           data = data2,
                           control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                           progress = TRUE,
                           expand_re = TRUE,
                           method =  'S',
                           REML=TRUE,
                           return = "merMod",
                           check_contrasts = FALSE)

summary(trimlogfpdifflmm1)

trimlogfpdifflmm2 <- update(trimlogfpdifflmm1, formula = ~ . -numForeperiod:numOneBackFPDiff)
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

#=================== 5.2. Run model separately for action and external conditions ==============

#======== 5.2.1. External ============#
# FP as numerical
trimlogfpdifflmmext1 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFPDiff + numForeperiod:numOneBackFPDiff +
                                (1 | ID),
                              data=data2[data2$condition=='external',],
                              control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                              progress = TRUE,
                              expand_re = TRUE,
                              method =  'KR',
                              REML=TRUE,
                              return = "merMod",
                              check_contrasts = FALSE)

summary(trimlogfpdifflmmext1)

# FP as categorical
trimlogfpdifflmmext2 <- mixed(logRT ~ 1 + foreperiod + numOneBackFPDiff + foreperiod:numOneBackFPDiff +
                                (1 | ID),
                              data=data2[data2$condition=='external',],
                              control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                              progress = TRUE,
                              expand_re = TRUE,
                              method =  'KR',
                              REML=TRUE,
                              return = "merMod",
                              check_contrasts = FALSE)

summary(trimlogfpdifflmmext2)

# logFP as numerical
trimlogfpdifflmmext3 <- mixed(logRT ~ 1 + numLogFP + numOneBackFPDiff + numLogFP:numOneBackFPDiff +
                                (1 | ID),
                              data=data2[data2$condition=='external',],
                              control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                              progress = TRUE,
                              expand_re = TRUE,
                              method =  'KR',
                              REML=TRUE,
                              return = "merMod",
                              check_contrasts = FALSE)

summary(trimlogfpdifflmmext3)

# logFP as categorical
trimlogfpdifflmmext4 <- mixed(logRT ~ 1 + logFP + numOneBackFPDiff + logFP:numOneBackFPDiff +
                                (1 | ID),
                              data=data2[data2$condition=='external',],
                              control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                              progress = TRUE,
                              expand_re = TRUE,
                              method =  'KR',
                              REML=TRUE,
                              return = "merMod",
                              check_contrasts = FALSE)

summary(trimlogfpdifflmmext4)


# FP as numerical without interaction
trimlogfpdifflmmext5 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFPDiff +
                                (1 | ID),
                              data=data2[data2$condition=='external',],
                              control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                              progress = TRUE,
                              expand_re = TRUE,
                              method =  'KR',
                              REML=TRUE,
                              return = "merMod",
                              check_contrasts = FALSE)

summary(trimlogfpdifflmmext5)
summary(trimlogfpdifflmmext1)

# FP as categorical without interaction
trimlogfpdifflmmext6 <- mixed(logRT ~ 1 + foreperiod + numOneBackFPDiff +
                                (1 | ID),
                              data=data2[data2$condition=='external',],
                              control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                              progress = TRUE,
                              expand_re = TRUE,
                              method =  'KR',
                              REML=TRUE,
                              return = "merMod",
                              check_contrasts = FALSE)


summary(trimlogfpdifflmmext6)
summary(trimlogfpdifflmmext2)


# logFP as numerical without interaction
trimlogfpdifflmmext7 <- mixed(logRT ~ 1 + numLogFP + numOneBackFPDiff +
                                (1 | ID),
                              data=data2[data2$condition=='external',],
                              control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                              progress = TRUE,
                              expand_re = TRUE,
                              method =  'KR',
                              REML=TRUE,
                              return = "merMod",
                              check_contrasts = FALSE)

summary(trimlogfpdifflmmext7)
summary(trimlogfpdifflmmext3)


# logFP as categorical without interaction
trimlogfpdifflmmext8 <- mixed(logRT ~ 1 + logFP + numOneBackFPDiff +
                                (1 | ID),
                              data=data2[data2$condition=='external',],
                              control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                              progress = TRUE,
                              expand_re = TRUE,
                              method =  'KR',
                              REML=TRUE,
                              return = "merMod",
                              check_contrasts = FALSE)

summary(trimlogfpdifflmmext8)
summary(trimlogfpdifflmmext4)

#======== 5.2.2. Action ============#
trimlogfpdifflmm1act <- mixed(logRT ~ 1 + numForeperiod + numOneBackFPDiff + numForeperiod:numOneBackFPDiff +
                                (1 | ID),
                              data=data2[data2$condition=='action',],
                              control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                              progress = TRUE,
                              expand_re = TRUE,
                              method =  'KR',
                              REML=TRUE,
                              return = "merMod",
                              check_contrasts = FALSE)

summary(trimlogfpdifflmm1act)

#==============================================================================================#
#====== 6. Model with quadratic terms for FP and for the difference between FP and FP n-1 ======
#==============================================================================================#

# Here we use ML instead of REML to properly compare AICs/BICs
trimlogfpdifflmm1 <- mixed(formula = logRT ~ 1 + condition + numForeperiod + numOneBackFPDiff + numForeperiod:numOneBackFPDiff + 
                             condition:numOneBackFPDiff + (1 + condition | ID),
                           data = data2,
                           control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                           progress = TRUE,
                           expand_re = TRUE,
                           method =  'S',
                           REML=FALSE,
                           return = "merMod")

trimlogfpdifflmm2 <- mixed(formula = logRT ~ 1 + condition + (numForeperiod + squaredNumForeperiod) + 
                             numOneBackFPDiff + (numForeperiod + squaredNumForeperiod):numOneBackFPDiff + 
                             condition:numOneBackFPDiff + 
                             (1 + condition | ID),
                           data = data2,
                           control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                           progress = TRUE,
                           expand_re = TRUE,
                           method =  'S',
                           REML=FALSE,
                           return = "merMod")

summary(trimlogfpdifflmm2)

trimlogfpdifflmm3 <- mixed(formula = logRT ~ 1 + condition + numForeperiod + 
                             (numOneBackFPDiff + squaredNumOneBackFPDiff) + 
                             numForeperiod:(numOneBackFPDiff + squaredNumOneBackFPDiff) + 
                             condition:(numOneBackFPDiff + squaredNumOneBackFPDiff) + 
                             (1 + condition | ID),
                           data = data2,
                           control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                           progress = TRUE,
                           expand_re = TRUE,
                           method =  'S',
                           REML=FALSE,
                           return = "merMod")


summary(trimlogfpdifflmm3)


trimlogfpdifflmm4 <- mixed(formula = logRT ~ 1 + condition + (numForeperiod + squaredNumForeperiod) + 
                             (numOneBackFPDiff + squaredNumOneBackFPDiff) + 
                             (numForeperiod + squaredNumForeperiod):(numOneBackFPDiff + squaredNumOneBackFPDiff) + 
                             condition:(numOneBackFPDiff + squaredNumOneBackFPDiff) + 
                             (1 + condition | ID),
                           data = data2,
                           control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                           progress = TRUE,
                           expand_re = TRUE,
                           method =  'S',
                           REML=FALSE,
                           return = "merMod")


summary(trimlogfpdifflmm4)


AIC(trimlogfpdifflmm1, trimlogfpdifflmm2, trimlogfpdifflmm3, trimlogfpdifflmm4) %>%
  bind_cols(BIC(trimlogfpdifflmm1, trimlogfpdifflmm2, trimlogfpdifflmm3, trimlogfpdifflmm4)) %>%
  kable()

# the full model has the lowest adjustment

check_collinearity(trimlogfpdifflmm4) %>%
  kable(digits = 3) %>%
  row_spec(0, bold=TRUE) %>%
  kable_styling(position = "center")

trimlogfpdifflmm3 %>%
  check_model()

# Anovas

anova(trimlogfpdifflmm1, trimlogfpdifflmm2) %>% tidy() %>%
  kable(digits=3) %>% row_spec(0, bold=T)

anova(trimlogfpdifflmm1, trimlogfpdifflmm3) %>% tidy() %>%
  kable(digits=3) %>% row_spec(0, bold=T)

anova(trimlogfpdifflmm1, trimlogfpdifflmm4) %>% tidy() %>%
  kable(digits=3) %>% row_spec(0, bold=T)

# No significant improvement

ggplot(data=data2,
       aes(x=squaredNumForeperiod,
           y=RT)) +
  stat_summary(fun = 'mean', geom = 'point') + 
  stat_summary(fun = 'mean', geom = 'line', aes(group=1))

trimlogfpdifflmm4 %>%
  augment() %>%
  bind_cols(trimlogfpdifflmm1 %>% augment() %>% select(.fitted) %>% rename(null_fitted = .fitted)) %>%
  #clean_names() %>%
  ggplot(data = .,
         mapping = aes(x = numForeperiod,
                       y = logRT,
                       color = condition)) +
  stat_summary(fun = 'mean', geom = 'point') +
  stat_summary(fun = 'mean', geom = 'line', aes(group=condition)) +
  #geom_line(aes(y=.fitted), color = 'blue') 
  stat_summary(fun = 'mean', geom = 'line', linetype = 'dashed', aes(y = .fitted, color = condition)) +
  stat_summary(fun = 'mean', geom = 'line', linetype = 'dotted', aes(y = null_fitted, color = condition)) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank()) +
  scale_color_manual(values = c('orange', 'blue'))



# To inspect estimates, run the same models with RT instead of logRT
trimfpdifflmm1 <- mixed(formula = RT ~ 1 + condition + numForeperiod + numOneBackFPDiff + numForeperiod:numOneBackFPDiff + 
                          condition:numOneBackFPDiff + (1 + condition | ID),
                        data = data2,
                        control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method =  'S',
                        REML=TRUE,
                        return = "merMod",
                        check_contrasts = FALSE)

summary(trimfpdifflmm1)

trimfpdifflmm2 <- mixed(formula = RT ~ 1 + condition + (numForeperiod + squaredNumForeperiod) + 
                          numOneBackFPDiff + (numForeperiod + squaredNumForeperiod):numOneBackFPDiff + 
                          condition:numOneBackFPDiff + 
                          (1 + condition | ID),
                        data = data2,
                        control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method =  'S',
                        REML=TRUE,
                        return = "merMod",
                        check_contrasts = FALSE)

summary(trimfpdifflmm2)

trimfpdifflmm3 <- mixed(formula = RT ~ 1 + condition + (numForeperiod + squaredNumForeperiod) + 
                          (numOneBackFPDiff + squaredNumOneBackFPDiff) + 
                          (numForeperiod + squaredNumForeperiod):(numOneBackFPDiff + squaredNumOneBackFPDiff) + 
                          condition:(numOneBackFPDiff + squaredNumOneBackFPDiff) + 
                          (1 + condition | ID),
                        data = data2,
                        control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method =  'S',
                        REML=TRUE,
                        return = "merMod",
                        check_contrasts = FALSE)

summary(trimfpdifflmm3)

AIC(trimfpdifflmm1, trimfpdifflmm2, trimfpdifflmm3) %>%
  bind_cols(BIC(trimfpdifflmm1, trimfpdifflmm2, trimfpdifflmm3)) %>%
  kable()

# Forest plot
pred_names <- rev(c('condition', 'foreperiod', 'foreperiod^2',
                    'fp_diff', 'fp_diff^2', 'foreperiod x fp_diff', 'foreperiod x fp_diff^2',
                    'foreperiod^2 x fp_diff', 'foreperiod^2 x fp_diff^2',
                    'condition x fp_diff', 'condition x fp_diff^2'))

plot_model(trimfpdifflmm3, transform = NULL,
           axis.labels = pred_names, vline.color = 'black',
           axis.lim = c(-1/10, 1/10)) + theme_sjplot()




#==================================================================================#
#========================== 7. Model using prevFPLonger ============================
#==================================================================================#
trimprevfplonglmm1 <- buildmer(logRT ~ foreperiod * condition * prevFPLonger + 
                                 (1 + foreperiod * condition * prevFPLonger|ID), 
                               data=data2,
                               buildmerControl = list(#direction='backward',
                                                      crit='LRT',#ddf = "Satterthwaite",
                                                      family=gaussian(link = 'identity'),
                                                      calc.anova = TRUE))

formula(trimprevfplonglmm1)


trimprevfplonglmm1 <- mixed(formula = logRT ~ 1 + condition + foreperiod + prevFPLonger + 
                                   condition:foreperiod + condition:prevFPLonger + 
                                   (1 + condition | ID),
                                 data = data2,
                                 control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                                 progress = TRUE,
                                 expand_re = TRUE,
                                 method =  'S',
                                 REML=TRUE,
                                 return = "merMod")

isSingular(trimprevfplonglmm1)

summary(trimprevfplonglmm1)
anova(trimprevfplonglmm1)
