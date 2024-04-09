

# Load necessary packages
library(readr)
library(ggplot2)
library(magrittr)
library(dplyr)
library(lattice)
library(afex)
library(emmeans)
library(lme4)
library(car)
library(data.table)
library(buildmer)
library(tidyr)
library(purrr)
library(janitor)
library(broom.mixed)
library(brms)


source('./Analysis/Prepare_data_6.R')

#==========================================================================================#
#================================= 2. Explore individual data ==============================
#==========================================================================================#


# Plot RT by FP by participant and model using pooled data
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


# Plot RT by FP by participant and model using individual data
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


#==========================================================================================#
#====================================== 2. Mixed models ====================================
#==========================================================================================#

# Set contrasts for variables used in the models
contrasts(goData$foreperiod) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)
contrasts(goData$condition) <- c(-1/2, 1/2)
contrasts(goData$oneBacktrialType) <- c(-1/2, 1/2)
contrasts(goData$oneBackFP) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)

contrastCodes <- cbind(c(3/4, -1/4, -1/4, -1/4),
                       c(-1/4, 3/4, -1/4, -1/4),
                       c(-1/4, -1/4, 3/4, -1/4)) %>%
  set_colnames(c('nogoVsShorter', 'longerVsShorter', 'equalVsShorter'))

contrasts(goData$oneBackFPDiff) <- contrastCodes

contrasts(goData2$foreperiod) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)
contrasts(goData2$condition) <- c(-1/2, 1/2)
contrasts(goData2$oneBacktrialType) <- c(-1/2, 1/2)
contrasts(goData2$oneBackFP) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)

contrastCodes <- cbind(c(3/4, -1/4, -1/4, -1/4),
                       c(-1/4, 3/4, -1/4, -1/4),
                       c(-1/4, -1/4, 3/4, -1/4)) %>%
  set_colnames(c('nogoVsShorter', 'longerVsShorter', 'equalVsShorter'))

contrasts(goData2$oneBackFPDiff) <- contrastCodes

#=========================== 3.1. Foreperiod, condition and sequential effects =============================
# We use the strategy of keeping it maximal to find a model that converges and progressively
# remove terms, one of the strategies recommended to avoid overfitting:
# https://rdrr.io/cran/lme4/man/isSingular.html

#  Trying to fit the model using foreperiod and FPn-1 as factors results
# in R hanging during execution; for this reason, we use them as
# numerical variables

# We start by fitting a model using mixed, from afex.

fplmm <- mixed(formula = RT ~ numForeperiod*condition*numOneBackFP + 
                  (1+numForeperiod*condition*numOneBackFP|ID),
                data = goData,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                return = 'merMod',
                REML=TRUE)

summary(fplmm)

# Although this solves the fitting problem, it does not solve the singular fit issue, 
# which indicates that the model is too complex (i.e., it is overfitted:
# https://stats.stackexchange.com/questions/378939/dealing-with-singular-fit-in-mixed-models)

goData$scaledNumForeperiod <- scale(goData$numForeperiod)[,1]
goData$scaledNumOneBackFP <- scale(goData$numOneBackFP)[,1]

scalefplmm <- mixed(formula = RT ~ scaledNumForeperiod*condition*scaledNumOneBackFP + 
                  (1+scaledNumForeperiod*condition*scaledNumOneBackFP|ID),
                data = goData,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE)

# This did not help, so we begin removing components

# 3.1.1.2. Remove correlations of mixed part
fplmm1v2 <- mixed(formula = RT ~ numForeperiod*condition*numOneBackFP + 
                  (1+numForeperiod*condition*numOneBackFP||ID),
                data = goData,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE)
summary(fplmm1v2)
anova(fplmm1v2)

# Singular fit persists

# 3.1.1.3. Remove interactions of mixed part
fplmm1v3 <- mixed(formula = RT ~ numForeperiod*condition*numOneBackFP + 
                  (1+numForeperiod+condition+numOneBackFP||ID),
                data = goData,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE)
summary(fplmm1v3)

# We got rid of the singular fit. Additionally, this model incorporates reasonable assumptions
# about the random effects structure, since we have no a priori reasons to assume
# interactions or correlations between random effects in these data

# Choose dependent variable


# Fit models with RT and inverse RT without trimming
fplmm1 <- mixed(formula = RT ~ numForeperiod*condition*numOneBackFP + 
                  (1+numForeperiod+condition+numOneBackFP||ID),
                data = goData,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE,
                return = "merMod")

# If we ran this with the maximal structure, it does not converge and returns a
# singular fit. This is the case for all models below we ran with the maximal structure
# (not shown here).

# Let's check that the current structure does not provide a singular fit:
isSingular(fplmm1)

# Now we run the same model with inverse RT as outcome
invfplmm1 <- mixed(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                     (1+numForeperiod+condition+numOneBackFP||ID),
                   data = goData,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

isSingular(invfplmm1)


anova(fplmm1)
anova(invfplmm1)

# Amount of variance accounted for by the model
cor(fitted(fplmm1), goData$RT)^2
cor(fitted(invfplmm1), goData$invRT)^2

# The first model explains a larger amount of the variance 

# Check normality of residuals
qqnorm(resid(fplmm1),
       main="Normal q-qplot fplmm1")
qqnorm(resid(invfplmm1),
       main="Normal q-qplot invfplmm1")

# Both models show considerable departures from normality

# Plot residuals
plot(fplmm1, resid(.) ~ fitted(.),
     main="Residuals fplmm1")
plot(invfplmm1, resid(.) ~ fitted(.),
     main="Residuals invfplmm1")

# It appears that residuals correlate somewhat with fitted values; there are also outliers

# Residual histograms
qplot(resid(fplmm1),
      main="Residuals fplmm1")
qplot(resid(invfplmm1),
      main="Residuals invfplmm1")

# Both appear to be relatively normally distributed, although the first has 
# a larger positive skew

# Fit models with RT and inverse RT after outlier trimming
trimfplmm <- mixed(formula = RT ~ numForeperiod*condition*numOneBackFP + 
                  (1+numForeperiod+condition+numOneBackFP||ID),
                data = goData2,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE,
                return = "merMod")

triminvfplmm <- mixed(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                     (1+numForeperiod+condition+numOneBackFP||ID),
                   data = goData2,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

isSingular(trimfplmm)
isSingular(triminvfplmm)

anova(trimfplmm)
anova(triminvfplmm)

# Amount of variance accounted for by the model
cor(fitted(trimfplmm), goData2$RT)^2
cor(fitted(triminvfplmm), goData2$invRT)^2

# Again, the first plot accounts for a larget amount of the variance

# check normality
qqnorm(resid(trimfplmm),
       main="Normal q-qplot trimfplmm")
qqnorm(resid(triminvfplmm),
       main="Normal q-qplot triminvfplmm")

# The second model is much closer to normality

# Plot residuals
plot(trimfplmm, resid(.) ~ fitted(.),
     main="Residuals trimfplmm")
plot(triminvfplmm, resid(.) ~ fitted(.),
     main="Residuals triminvfplmm")

# Outliers are gone, but residuals still appear to correlate with fitted values in
# the first model

qplot(resid(trimfplmm),
      main="Residuals trimfplmm")
qplot(resid(triminvfplmm),
      main="Residuals triminvfplmm")

# Both still appear relatively normally distributed, with the second model
# performing better

# Use z-scores for centering with no outlier trimming
fplmm3 <- mixed(formula = RTzscore ~ numForeperiod*condition*numOneBackFP + 
                  (1+numForeperiod+condition+numOneBackFP||ID),
                data = goData,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE,
                return = "merMod")

isSingular(fplmm3)

anova(fplmm3)

# Amount of variance accounted for by the model
cor(fitted(fplmm3), goData$RT)^2

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
                  (1+numForeperiod+condition+numOneBackFP||ID),
                data = goData2,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE,
                return = "merMod")

isSingular(fplmm4)

# Singular fit

# Amount of variance accounted for by the model
cor(fitted(fplmm4), goData2$RT)^2

# check normality
qqnorm(resid(fplmm4),
       main="Normal q-qplot fplmm4")

# Plot residuals
plot(fplmm4, resid(.) ~ fitted(.),
     main="Residuals fplmm4")

qplot(resid(fplmm4),
      main="Residuals fplmm4")

# Only models with RT or inverse RT perform well

# Of those:
# Variance explained is higher with trimming
# Q-q plots are better with trimming
# Residuals are less correlated with fitted values with trimming
# Residuals are more normally distributed with trimming

# The model with inverse RT performs better than the one with RT

# Model comparisons
triminvfplmm2 <- mixed(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                         (1+numForeperiod+condition||ID),
                       data = goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod")

anova(triminvfplmm2, triminvfplmm, refit=FALSE)

triminvfplmm3 <- mixed(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                         (1+numForeperiod+numOneBackFP||ID),
                       data = goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod")

anova(triminvfplmm3, triminvfplmm, refit=FALSE)

triminvfplmm4 <- mixed(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                         (1+condition+numOneBackFP||ID),
                       data = goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod")

anova(triminvfplmm4, triminvfplmm, refit=FALSE)

# Amount of variance accounted for by the models
cor(fitted(triminvfplmm), goData2$RT)^2
cor(fitted(triminvfplmm2), goData2$invRT)^2
cor(fitted(triminvfplmm3), goData2$invRT)^2
cor(fitted(triminvfplmm4), goData2$invRT)^2

# Arrange by BIC
BIC(triminvfplmm, triminvfplmm2, triminvfplmm3, triminvfplmm4) %>%
  arrange(BIC)

# The fullest model explains most of the variance and BICs are very close between
# this and the model without a random-effect from Fpn-1. In line with the 
# "keep it maximal" strategy, we keep this full model


#===================================================================================#
# Sanity check: model comparisons without trimming
#===================================================================================#
invfplmm2 <- mixed(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                         (1+numForeperiod+condition||ID),
                       data = goData,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod")

anova(invfplmm2, invfplmm1, refit=FALSE)

invfplmm3 <- mixed(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                     (1+numForeperiod+numOneBackFP||ID),
                   data = goData,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

anova(invfplmm3, invfplmm1, refit=FALSE)

invfplmm4 <- mixed(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                     (1+condition+numOneBackFP||ID),
                   data = goData,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

anova(invfplmm4, invfplmm1, refit=FALSE)

# Amount of variance accounted for by the models
cor(fitted(invfplmm1), goData$RT)^2
cor(fitted(invfplmm2), goData$invRT)^2
cor(fitted(invfplmm3), goData$invRT)^2
cor(fitted(invfplmm4), goData$invRT)^2

# 3.1.3.1. Fullest model
fpDifflmm1 <- mixed(formula = invRT ~ numForeperiod * condition * numOneBackFPDiff + 
                      (1+numForeperiod*condition*numOneBackFPDiff|ID),
                    data = goData2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    return = "merMod")

fpDifflmm2 <- mixed(formula = invRT ~ foreperiod * condition * oneBackFPDiff + 
                      (1+foreperiod*condition*oneBackFPDiff|ID),
                    data = goData2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    return = "merMod")

anova(fpDifflmm2)

fpDifflmm3 <- mixed(formula = invRT ~ foreperiod * condition * oneBackFPDiff + 
                      (1+foreperiod*condition*oneBackFPDiff|ID),
                    data = goData2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    return = "merMod")

fpDifflmm4 <- mixed(formula = invRT ~ foreperiod * condition * oneBackFPDiff + 
                      (1+foreperiod*condition*oneBackFPDiff||ID),
                    data = goData2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    return = "merMod")

fpDifflmm5 <- mixed(formula = invRT ~ foreperiod * condition * oneBackFPDiff + 
                      (1+foreperiod+condition+oneBackFPDiff||ID),
                    data = goData2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    return = "merMod")

cor(fitted(fpDifflmm2), goData2$invRT)^2
cor(fitted(fpDifflmm4), goData2$invRT)^2
cor(fitted(fpDifflmm5), goData2$invRT)^2

#==============================================================================================#
#============================== 5. Bayesian mixed models using brms ============================
#==============================================================================================#

#========================== 5.1. For Log RT ==============================
#======= 5.1.1 Models using effects from previous experiment as priors

# Baseado no exp 0 v1
# prior1 <- c(set_prior('normal(2.55,0.009)', class = 'Intercept'),
#             set_prior("normal(-0.03, 0.002)", class = "b", coef = "foreperiod2"),
#             set_prior("normal(-0.02, 0.002)", class = "b", coef = "foreperiod3"),
#             set_prior("normal(0.01, 0.002)", class = "b", coef = "oneBackFP2"),
#             set_prior("normal(0.02, 0.002)", class = "b", coef = "oneBackFP2"),
#             set_prior("normal(0.01, 0.007)", class = "b", coef = "condition1"),
#             set_prior("normal(-0.02, 0.005)", class = "b", coef = "foreperiod2:oneBackFP2"),
#             set_prior("normal(-0.03, 0.005)", class = "b", coef = "foreperiod3:oneBackFP2"),
#             set_prior("normal(-0.03, 0.005)", class = "b", coef = "foreperiod2:oneBackFP3"),
#             set_prior("normal(-0.04, 0.005)", class = "b", coef = "foreperiod3:oneBackFP3"),
#             set_prior("normal(-0.01, 0.004)", class = "b", coef = "foreperiod2:condition1"),
#             set_prior("normal(-0.05, 0.004)", class = "b", coef = "foreperiod3:condition1"),
#             set_prior("normal(0.003, 0.004)", class = "b", coef = "oneBackFP2:condition1"),
#             set_prior("normal(0.002, 0.004)", class = "b", coef = "oneBackFP3:condition1"),
#             set_prior("normal(-0.02, 0.009)", class = "b", coef = "foreperiod2:oneBackFP2:condition1"),
#             set_prior("normal(-0.02, 0.009)", class = "b", coef = "foreperiod3:oneBackFP2:condition1"),
#             set_prior("normal(-0.02, 0.009)", class = "b", coef = "foreperiod2:oneBackFP3:condition1"),
#             set_prior("normal(-0.03, 0.009)", class = "b", coef = "foreperiod3:oneBackFP3:condition1")
# )
prior1 <- c(set_prior('normal(2.55,0.009)', class = 'Intercept'),
            set_prior("normal(-0.03, 0.002)", class = "b", coef = "foreperiod2"),
            set_prior("normal(-0.02, 0.002)", class = "b", coef = "foreperiod3"),
            set_prior("normal(0.01, 0.002)", class = "b", coef = "oneBackFP2"),
            set_prior("normal(0.02, 0.002)", class = "b", coef = "oneBackFP3"),
            set_prior("normal(0.01, 0.007)", class = "b", coef = "condition1")
)

# Gaussian distr
b_one_back_fp <- brm(formula = logRT ~ foreperiod * condition * oneBackFP + 
                       (1+numForeperiod*condition*numOneBackFP|ID),
                     data = data2,
                     family = gaussian(),
                     prior = prior1,
                     control = list(adapt_delta = 0.9),
                     save_all_pars = TRUE,
                     warmup = 2000,
                     iter = 12000,
                     cores = -1)

saveRDS(b_one_back_fp, "b_one_back_fp.rds")

#b_one_back_fp <- readRDS('./Analysis/b_one_back_fp.rds')


ggplot() +
  stat_summary(fun='mean', geom='point', data=data2, aes(x=numForeperiod, y=logRT)) +
  geom_point(data=fp_effect_df, aes(x=numForeperiod, y=fit), color='red') +
  geom_line(data=fp_effect_df, aes(x=numForeperiod, y=fit), color='red') +
  geom_ribbon(data=fp_effect_df, aes(x=numForeperiod, ymin=lower, ymax=upper), alpha=0.3) +
  labs(x='Foreperiod (continuous)', y='logRT')



# Inverse gaussian distr
b_one_back_fp_invgaus <- brm(formula = RT ~ numForeperiod * condition * numOneBackFP + 
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
                          control = list(adapt_delta = 0.9),
                          warmup = 2000,
                          iter = 10000,
                          file="b_one_back_fp_null.rds")



bf_one_back_brm <- bayes_factor(b_one_back_fp, b_one_back_fp_null)


#========== 5.1.2 Using a vague prior
b_one_back_fp_vagueprior <- brm(formula = logRT ~ numForeperiod * condition * numOneBackFP + 
                                  (1+numForeperiod*condition*numOneBackFP|ID),
                                data = data2,
                                family = gaussian(),
                                save_all_pars = TRUE,
                                control = list(adapt_delta = 0.9),
                                warmup = 2000,
                                iter = 12000,
                                cores = -1)


options(digits = 5)
summary(b_one_back_fp)
options(options_defaults)

#===================================== 5.2. For linear RT ======================================
#======= 5.2.1 Models using effects from previous experiment as priors

# Baseado no exp 0 v1
# prior2 <- c(set_prior('normal(368.41,8.34)', class = 'Intercept'),
#             set_prior("normal(-27.82, 1.78)", class = "b", coef = "foreperiod2"),
#             set_prior("normal(-24.96, 1.78)", class = "b", coef = "foreperiod3"),
#             set_prior("normal(8.57, 1.77)", class = "b", coef = "oneBackFP2"),
#             set_prior("normal(17.22, 1.77)", class = "b", coef = "oneBackFP2"),
#             set_prior("normal(12.6, 6.27)", class = "b", coef = "condition1"),
#             set_prior("normal(-14.39, 4.38)", class = "b", coef = "foreperiod2:oneBackFP2"),
#             set_prior("normal(-24.31, 4.37)", class = "b", coef = "foreperiod3:oneBackFP2"),
#             set_prior("normal(-26.36, 4.34)", class = "b", coef = "foreperiod2:oneBackFP3"),
#             set_prior("normal(-37.74, 4.37)", class = "b", coef = "foreperiod3:oneBackFP3"),
#             set_prior("normal(-15.27, 3.55)", class = "b", coef = "foreperiod2:condition1"),
#             set_prior("normal(-22.61, 3.56)", class = "b", coef = "foreperiod3:condition1"),
#             set_prior("normal(2.92, 3.55)", class = "b", coef = "oneBackFP2:condition1"),
#             set_prior("normal(1.45, 3.54)", class = "b", coef = "oneBackFP3:condition1"),
#             set_prior("normal(-12.49, 8.75)", class = "b", coef = "foreperiod2:oneBackFP2:condition1"),
#             set_prior("normal(-21.64, 8.74)", class = "b", coef = "foreperiod3:oneBackFP2:condition1"),
#             set_prior("normal(-20.01, 8.70)", class = "b", coef = "foreperiod2:oneBackFP3:condition1"),
#             set_prior("normal(-23.69, 8.75)", class = "b", coef = "foreperiod3:oneBackFP3:condition1")
# )
prior2 <- c(set_prior('normal(368.41,8.34)', class = 'Intercept'),
            set_prior("normal(-27.82, 1.78)", class = "b", coef = "foreperiod2"),
            set_prior("normal(-24.96, 1.78)", class = "b", coef = "foreperiod3"),
            set_prior("normal(8.57, 1.77)", class = "b", coef = "oneBackFP2"),
            set_prior("normal(17.22, 1.77)", class = "b", coef = "oneBackFP3"),
            set_prior("normal(12.6, 6.27)", class = "b", coef = "condition1")
)


# Gaussian distr
b_one_back_fp <- brm(formula = RT ~ foreperiod * condition * oneBackFP + 
                       (1+foreperiod*condition*oneBackFP|ID),
                     data = data2,
                     family = gaussian(),
                     prior = prior2,
                     control = list(adapt_delta = 0.9),
                     save_all_pars = TRUE,
                     warmup = 2000,
                     iter = 12000,
                     cores = -1)

saveRDS(b_one_back_fp, "b_one_back_fp.rds")

#b_one_back_fp <- readRDS('./Analysis/b_one_back_fp.rds')


ggplot() +
  stat_summary(fun='mean', geom='point', data=data2, aes(x=numForeperiod, y=logRT)) +
  geom_point(data=fp_effect_df, aes(x=numForeperiod, y=fit), color='red') +
  geom_line(data=fp_effect_df, aes(x=numForeperiod, y=fit), color='red') +
  geom_ribbon(data=fp_effect_df, aes(x=numForeperiod, ymin=lower, ymax=upper), alpha=0.3) +
  labs(x='Foreperiod (continuous)', y='logRT')



# Inverse gaussian distr
b_one_back_fp_invgaus <- brm(formula = RT ~ foreperiod * condition * oneBackFP + 
                               (1+foreperiod*condition*oneBackFP|ID),
                             data = data2,
                             family = gaussian(),
                             prior = prior1,
                             save_all_pars = TRUE,
                             warmup = 2000,
                             iter = 10000)


b_one_back_fp_null <- brm(formula = RT ~ 1 + condition + foreperiod + condition:foreperiod + 
                            oneBackFP + foreperiod:oneBackFP + 
                            (1 + condition + foreperiod + oneBackFP + 
                               condition:numForeperiod + numForeperiod:numOneBackFP| ID),
                          data = data2,
                          family = gaussian(),
                          prior = prior1,
                          save_all_pars = TRUE,
                          control = list(adapt_delta = 0.9),
                          warmup = 2000,
                          iter = 10000,
                          file="b_one_back_fp_null.rds")



bf_one_back_brm <- bayes_factor(b_one_back_fp, b_one_back_fp_null)


#========== 5.1.2 Using a vague prior
b_one_back_fp_vagueprior <- brm(formula = RT ~ foreperiod * condition * oneBackFP + 
                                  (1+foreperiod*condition*oneBackFP|ID),
                                data = data2,
                                family = gaussian(),
                                save_all_pars = TRUE,
                                control = list(adapt_delta = 0.9),
                                warmup = 2000,
                                iter = 12000,
                                cores = -1)


options(digits = 5)
summary(b_one_back_fp)
options(options_defaults)

#============================================================================================================#
#======================== 4. Difference between the durations of FP n and FP n-1 as predictor ================
#============================================================================================================#

# Fullest model using foreperiod and one back fp as categorical variables
fpDifflmm1 <- mixed(formula = invRT ~ condition * oneBackFPDiff + 
                      (1+condition*oneBackFPDiff|ID),
                    data = goData2,
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
                    data = goData2,
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
                    data = goData2,
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
                      data=goData2,
                      buildmerControl = buildmerControl(include = ~ foreperiod * condition * oneBackFPDiff, calc.anova = TRUE, ddf = "Satterthwaite"))


#==============================================================================================#
#============================== 5. Model with scaled numerical predictors ======================
#==============================================================================================#

# This allows us to see if a more complex random-effects structure can be fitted with 
# centered and standardized predictors, and also how they contribute to the variance of
# the model

prior2 <- c(set_prior('normal(2.56, 0.009)', class = 'Intercept'),
            set_prior("normal(-0.01, 0.001)", class = "b", coef = "scalednumForeperiod"),
            set_prior("normal(0.01, 0.006)", class = "b", coef = "condition1"),
            set_prior("normal(0.008, 0.001)", class = "b", coef = "scaledNumOneBackFP")
)

# Gaussian distr
b_one_back_fp <- brm(formula = logRT ~ scaledNumForeperiod * condition * scaledNumOneBackFP + 
                       (1+scaledNumForeperiod*condition*scaledNumOneBackFP|ID),
                     data = data2,
                     family = gaussian(),
                     prior = prior2,
                     control = list(adapt_delta = 0.9),
                     save_all_pars = TRUE,
                     warmup = 2000,
                     iter = 12000)
