



# Mixed-effects modeling
library(afex)
library(emmeans)
library(lme4)
library(buildmer)

library(mixedpower)

setwd("~/Dropbox/UFABC/Projetos_Alunos/AlexandreNobre/GitStuff/Action_foreperiod_exp0_v2-master")

source('./Analysis/Prepare_data_6.R')

# Save defaults
graphical_defaults <- par()
options_defaults <- options()

# emm options
emm_options(lmer.df = "satterthwaite", lmerTest.limit = 12000)

contrasts(data2$foreperiod) <- contr.sum(levels(data2$foreperiod))
contrasts(data2$condition) <-contr.sum(levels(data2$condition))
contrasts(data2$oneBackFP) <- contr.sum(levels(data2$oneBackFP))

data2$ID=as.numeric(data2$ID)

#========================= Find model

# I ran using buildmer just to see what happens and we end up with the exact same model as using heuristics
find_model <- buildmer(formula= logRT ~ 1 + foreperiod * condition * oneBackFP +(1+foreperiod*condition*oneBackFP|ID),
                       data=data2,
                       buildmerControl = buildmerControl(ddf = "Satterthwaite",calc.anova = TRUE,include='foreperiod * condition * oneBackFP'))

#Fit data with lmer
model_lmer=lmer(logRT ~ 1 + foreperiod + condition + foreperiod:condition + oneBackFP +  
                  foreperiod:oneBackFP + condition:oneBackFP + foreperiod:condition:oneBackFP +(1 + condition | ID), 
      data=data2,
      control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
      REML=TRUE)

#Set params for mixedpower
model <- model_lmer # which model do we want to simulate power for?
data <- data2 # data used to fit the model
fixed_effects <- c("foreperiod", "condition",'oneBackFP') # all fixed effects specified in FLPmodel
simvar <- "subject" # which random effect do we want to vary in the simulation?

# SIMULATION PARAMETERS
steps <- c(20, 30, 40, 50, 60) # which sample sizes do we want to look at?
critical_value <- 2 # which t/z value do we want to use to test for significance?
n_sim <- 1000 # how many single simulations should be used to estimate power?




# RUN SIMULATION
#I am using only 10 simulations just to look what the output looks like..we need way more than that (around 5000)
power_FLP <- mixedpower(model = model_lmer, data = data2,
                        fixed_effects = c("foreperiod", "condition",'oneBackFP'),
                        simvar = "ID", steps = c(20,30,40,50,60),
                        critical_value = 2, n_sim = 10)

##-----------Try with simr
library(simr)
model_lmer_extended <- extend(model_lmer, along="ID", n = 60) #extend data set to include 60 participants

getSimrOption("nsim")
oldopts <- simrOptions(nsim=10)
getSimrOption("nsim")
#simrOptions(oldopts)
#getSimrOption("nsim")

powerC <- powerCurve(fit = model_lmer, test = fixed("foreperiod"), along = "ID",
                     breaks = c(10, 15, 20))

# 

# -------------------------------------------- #
# EXTEND DATA SET
FLPmodel_extended <- extend(FLPmodel, along="subject", n = 60) #extend data set to include 60 participants




#================================== 0. Read data ================================
# Create dataset
source('./Analysis/Prepare_data_6.R')
#This model does not even run until the end...
FullModel <- mixed(logRT ~ 1 + foreperiod * condition * oneBackFP +(1+foreperiod*condition*oneBackFP|ID), 
                       data=data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod")

Model_1 <- mixed(logRT ~ 1 + foreperiod * condition * oneBackFP +(1+foreperiod*condition*oneBackFP||ID), 
                   data=data2,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

isSingular(Model_1)
anova(Model_1)

Model_2 <- mixed(logRT ~ 1 + foreperiod * condition * oneBackFP +(1+foreperiod+condition+oneBackFP||ID), 
                 data=data2,
                 control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                 progress = TRUE,
                 expand_re = TRUE,
                 method =  'S',
                 REML=TRUE,
                 return = "merMod")

isSingular(Model_2)
anova(Model_2)

Model_3 <- mixed(logRT ~ 1 + foreperiod * condition * oneBackFP +(1+foreperiod+condition||ID), 
                 data=data2,
                 control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                 progress = TRUE,
                 expand_re = TRUE,
                 method =  'S',
                 REML=TRUE,
                 return = "merMod")

isSingular(Model_3)
anova(Model_3)

Model_4 <- mixed(logRT ~ 1 + foreperiod * condition * oneBackFP +(1+condition||ID), 
                 data=data2,
                 control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                 progress = TRUE,
                 expand_re = TRUE,
                 method =  'S',
                 REML=TRUE,
                 return = "merMod")

isSingular(Model_4)
anova(Model_4)

emm_options(pbkrtest.limit = 12000)

custom_colors <- c("red", "blue")
afex_plot(Model_4, "foreperiod", "condition", id = "ID",
          data_geom = ggbeeswarm::geom_quasirandom, mapping = c("shape", "fill"),
          data_arg = list(
            dodge.width = 0.5,  ## needs to be same as dodge
            cex = 0.8,
            width = 0.1,
            color = "red"))

afex_plot(Model_4, x = "foreperiod", trace = "condition", id='ID',
                mapping = c("shape", "fill"),
                data_geom = ggplot2::geom_boxplot, 
                data_arg = list(width = 0.3))

afex_plot(Model_4, "foreperiod", "oneBackFP", "condition",id = "ID",
          data_geom = ggbeeswarm::geom_quasirandom, mapping = c("shape", "fill"),
          data_arg = list(
            dodge.width = 0.5,  ## needs to be same as dodge
            cex = 0.8,
            width = 0.1,
            color = "darkgrey"))

emm_options(lmerTest.limit = 12000)
emm_options(lmer.df = "satterthwaite") # also possible: 'satterthwaite', 'kenward-roger'
emm_fp_condition <- emmeans(Model_4, "condition", by = c("foreperiod"))
emm_fp_condition
update(pairs(emm_fp_condition), by = NULL, adjust = "holm")

emm_fp_oneback <- emmeans(Model_4, "oneBackFP", by = c("foreperiod"))
emm_fp_oneback
update(pairs(emm_fp_oneback), by = NULL, adjust = "holm")


joint_tests(Model_4, by = c("foreperiod", "condition"))
joint_tests(Model_4, by = c("foreperiod"))


