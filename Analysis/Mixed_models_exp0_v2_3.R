

# Load packages

# Data processing
library(magrittr)
library(tidyverse)

# Plotting
library(lattice)
library(gridExtra)
library(data.table)
library(extrafont)
library(egg)

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

# Load fonts from extrafonts
loadfonts()

# Prepare data 
source('./Analysis/Prepare_data_6.R')

# Prepare theme for plots
source("./Analysis/plot_theme.R")
theme_set(mytheme)

# emm options
emm_options(lmer.df = "satterthwaite", lmerTest.limit = 12000)


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

fplmm <- mixed(formula = logRT ~ 1 + foreperiod + condition + foreperiod:condition + oneBackFP + 
                 foreperiod:oneBackFP + condition:oneBackFP + foreperiod:condition:oneBackFP + 
                 (1 + condition | ID),
               data=data2,
               control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
               progress = TRUE,
               expand_re = FALSE,
               method =  'KR',
               REML=TRUE,
               return = "merMod")

isSingular(fplmm)

#============ Model for sensitivity analysis ===========
contrasts(data2$foreperiod, how.many = 1) <- contr.poly(4)
contrasts(data2$condition) <- contr.sum(2)
contrasts(data2$oneBackFP) <- contr.sum(4)

fplmm <- mixed(formula = logRT ~ 1 + foreperiod + condition + foreperiod:condition + oneBackFP + 
                 foreperiod:oneBackFP + condition:oneBackFP + foreperiod:condition:oneBackFP + 
                 (1 + condition | ID),
               data=data2,
               control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
               progress = TRUE,
               expand_re = FALSE,
               method =  'KR',
               REML=TRUE,
               return = "merMod",
               check_contrasts = FALSE)

summary(fplmm)

saveRDS(fplmm, "./Analysis/fplmm.rds")


#==============================================================================================#
#==================================== 3. Model assessment ======================================
#==============================================================================================#


#=========== 3.1. Omnibus anova =================
anova(fplmm)


#Visualize random effects
dotplot(ranef(fplmm, condVar = TRUE))

#========== 3.2. Two-way interactions ==========
# 3.2.1. Compare difference between conditions within each level of FPn
fp_cond_emm <- emmeans(fplmm, ~ condition|foreperiod)
contrast(fp_cond_emm, interaction = c("pairwise"), adjust = "holm")

# 3.2.2. Exame effect of FP within each condition
cond_emm <- emmeans(fplmm, ~ foreperiod|condition)
contrast(cond_emm, interaction = c("poly"), adjust = "holm", max.degree = 1)

# Compare FP effect between conditions
cond_emm <- emmeans(fplmm, ~ foreperiod*condition)
contrast(cond_emm, interaction = c("poly"), adjust = "holm", max.degree = 1)

# 3.2.3. Compare difference between consecutive levels of FPn-1 for each level of FPn
fp_onebackfp_emm <- emmeans(fplmm, ~ oneBackFP|foreperiod)
contrast(fp_onebackfp_emm, interaction = c("consec"), adjust = "holm")

fp_onebackfp_emm <- emmeans(fplmm, ~ foreperiod|oneBackFP, )
contrast(fp_onebackfp_emm, interaction = c("poly"), adjust = "holm", max.degree = 1)

#========== 3.3 Three-way interaction =========
# Multiple plots to visualize interactions
# Pairwise comparisons by FP (estimates consecutive differences)
emmip(fplmm, condition ~ foreperiod|oneBackFP, CIs = TRUE,
      xlab = "FP",
      facelab = "label_both") +
  labs(title = "RT pairwise comparisons") +
  theme(plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_color_manual(values = c("orange", "blue"))


emmip(fplmm, oneBackFP ~ foreperiod|condition, CIs = TRUE,
      xlab = "FP",
      facelab = "label_both") +
  labs(title = "RT pairwise comparisons") +
  theme(plot.title = element_text(size = 16, hjust = 0.5))


# First analysis: compare two-way interactions (FPn x FPn-1) between conditions
# Compare FP linear fits between consecutive levels of FPn-1, separately for each condition
fpemm <- emmeans(fplmm, ~ foreperiod*oneBackFP|condition)
contrast(fpemm, interaction = c("poly", "consec", "pairwise"), adjust = "holm", max.degree = 1)

# Now compare magnitude of those interactions between conditions
fpemm <- emmeans(fplmm, ~ foreperiod*oneBackFP*condition)
contrast(fpemm, interaction = c("poly", "consec", "pairwise"), adjust = "holm", max.degree = 1)

# Second analysis: compare FP effect between conditions, separately for each FPn-1
# Test significance of FP linear fits separately for each level of FPn-1 and each condition
fpemm <- emmeans(fplmm, ~ foreperiod|oneBackFP*condition)
contrast(fpemm, interaction = c("poly"), adjust = "holm", max.degree = 1)


# Now test difference of FP fits between conditions, separately for each level of FPn-1
fpemm <- emmeans(fplmm, ~ foreperiod*condition|oneBackFP)
contrast(fpemm, interaction = c("poly", "pairwise"), adjust = "holm", max.degree = 1)


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


#===============================================================================================#
#===================================== 8. Accuracy ==============================================
#===============================================================================================#

fpaccglmm <- mixed(formula = error_result ~ 1 + foreperiod:condition:oneBackFP + foreperiod + 
                     condition + oneBackFP + foreperiod:condition + foreperiod:oneBackFP + 
                     condition:oneBackFP + (1|ID),
                   data=dataAcc,
                   family = binomial(link = "logit"),
                   control = glmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5)),
                   progress = TRUE,
                   expand_re = FALSE,
                   method = "LRT")

isSingular(fpaccglmm)
anova(fpaccglmm)

#========== 8.2. Two-way interactions ==========
# Compare difference between conditions within each level of FPn
fp_cond_acc_emm <- emmeans(fpaccglmm, ~ condition|foreperiod)
contrast(fp_cond_acc_emm, interaction = c("pairwise"), adjust = "holm")

# 3.2.2. Exame effect of FP within each condition
fp_cond_acc_emm <- emmeans(fpaccglmm, ~ foreperiod|condition)
contrast(fp_cond_acc_emm, interaction = c("poly"), adjust = "holm", max.degree = 1)

# Compare FP effect between conditions
fp_cond_acc_emm <- emmeans(fpaccglmm, ~ foreperiod*condition)
contrast(fp_cond_acc_emm, interaction = c("poly"), adjust = "holm", max.degree = 1)


#========== 8.3. Three-way interactions ==========
# Test significance of FP linear fits separately for each level of FPn-1 and each condition
fp_acc_emm <- emmeans(fpaccglmm, ~ foreperiod|oneBackFP*condition)
contrast(fp_acc_emm, interaction = c("poly"), adjust = "holm", max.degree = 1)


# Now test difference of FP fits between conditions, separately for each level of FPn-1
fp_acc_emm <- emmeans(fpaccglmm, ~ foreperiod*condition|oneBackFP)
contrast(fp_acc_emm, interaction = c("poly", "pairwise"), adjust = "holm", max.degree = 1)

#===============================================================================================#
#============================= 9. Analyses using exp(-foreperiod) ==============================
#===============================================================================================#
trimexpfplmm1 <- buildmer(RT ~ expFP * condition * oneBackFP + 
                            (1+expFP*condition*oneBackFP|ID), 
                          data=data2,
                          buildmerControl = list(direction='backward',
                                                 crit='LRT',#ddf = "Satterthwaite",
                                                 family=gaussian(link = 'identity'),
                                                 calc.anova = TRUE))

#===============================================================================================#
#=================================== 10. Analyses with smaller N ====================================
#===============================================================================================#

nPartsExp3 <- 28
nReps <- 1000
pValues <- vector(mode = "numeric", length = nReps)

# Sample 28 participants nReps times and fit mixed model for each sample, extracting p-value of comparison of interest
for(iRep in 1:nReps) {
  
  # Sample participants
  sampParts <- sample(levels(data2$ID), nPartsExp3)
  sampData <- data2 %>%
    filter(ID %in% sampParts)
  
  # Recompute scaled variables (which were de-centered by trimming)
  sampData$squaredNumForeperiod <-  (sampData$numForeperiod)^2
  sampData$squaredNumOneBackFP <- (sampData$numOneBackFP)^2
  sampData$scaledNumForeperiod <-  scale(sampData$numForeperiod, scale = FALSE)[,1]
  sampData$squaredScaledNumForeperiod <- (sampData$scaledNumForeperiod)^2
  sampData$scaledNumOneBackFP <- scale(sampData$numOneBackFP, scale = FALSE)[,1]
  
  
  # Fit model using same structure as before
  sampLmm <- mixed(formula = logRT ~ 1 + scaledNumForeperiod + oneBackFP + scaledNumForeperiod:oneBackFP + 
                     condition + scaledNumForeperiod:condition + oneBackFP:condition + 
                     scaledNumForeperiod:oneBackFP:condition + 
                     (1 + condition + scaledNumForeperiod | ID),
                   data=sampData,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = FALSE,
                   method =  'KR',
                   REML=TRUE,
                   return = "merMod")
  
  
  # Run anova and extract p-value
  sampAnova <- anova(sampLmm)
  #pValues[iRep] <- sampAnova['scaledNumForeperiod:oneBackFP:condition',]$`Pr(>F)`
  pValues[iRep] <- sampAnova['scaledNumForeperiod:condition',]$`Pr(>F)`
  
}

# Get p-value of comparison with full sample
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

stdAnova <- anova(trimlogfplmm)

# Plot sampled values against original values
jpeg("./Analysis/Plots/Smaller_sample_analyses/28_parts_1000_sims.jpg", width = 20, height = 12, units = "cm", res = 300)
hist(pValues, breaks = 100, col = "lightblue",
     main = "Exp. 2: p-Values for sample size = 28 (1000 simulations)")
#abline(v = stdAnova['scaledNumForeperiod:oneBackFP:condition',]$`Pr(>F)`, lwd = 2)
abline(v = stdAnova['scaledNumForeperiod:condition',]$`Pr(>F)`, lwd = 2)
abline(v = 0.05, lwd = 2, lty = 2)
text(x = 0.025, y = 600, labels = "p-value for full sample")
text(x = 0.06, y = 600, labels = "p = .05")
dev.off()
