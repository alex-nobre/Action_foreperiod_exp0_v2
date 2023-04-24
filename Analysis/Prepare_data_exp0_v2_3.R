


#==============================================================================#
# Changes
# Replaced column for FP n-1 diff, since it will be used as a continuous variable
# now
# Use RT in ms instead of s, so that logRT does not assume negative values
# log of predictors is now log10 instead of log2
#==============================================================================#

# Load necessary packages
library(readr)
library(ggplot2)
library(magrittr)
library(dplyr)

source("G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_0_v2/Analysis/helper_functions.R")

# Read data
data <- read_csv("G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_0_v2/Analysis/dataActionFPAll.csv")

# Remove unnecessary columns
data <- data %>%
  dplyr::select(ID, Acc, condition, block, orientation,
                foreperiod, RT, counterbalance, 
                extFixationDuration, action_trigger.rt,
                oneBackFP, twoBackFP, oneBackEffect)  %>%
  mutate(foreperiod = foreperiod * 1000,
         RT = RT *1000,
         extFixationDuration = extFixationDuration * 1000,
         action_trigger.rt = action_trigger.rt * 1000,
         oneBackFP = oneBackFP * 1000,
         twoBackFP = twoBackFP * 1000,
         oneBackEffect = oneBackEffect * 1000)

# Coerce to factors
data$ID <- as.factor(data$ID)
data$condition <- data$condition %>%
  as.factor() %>%
  forcats::fct_relevel(c("external", "action"))
data$block <- as.factor(data$block)
data$orientation <- as.factor(data$orientation)
data$foreperiod <- as.factor(data$foreperiod)
data$counterbalance <- as.factor(data$counterbalance)
data$oneBackFP <- as.factor(data$oneBackFP)
data$twoBackFP <- as.factor(data$twoBackFP)

# Create numeric versions of foreperiod and FP n-1
data$numForeperiod <- as.numeric(as.character(data$foreperiod))
data$numOneBackFP <- as.numeric(as.character(data$oneBackFP))

# Remove practice trials
data <- data %>%
  filter(condition != 'practice')

# Create column for difference between current and previous FP as a numeric variable:
data <- data %>%
  group_by(ID, block) %>%
  mutate(numOneBackFPDiff = c(NA, diff(numForeperiod))) %>%
  mutate(oneBackFPDiff = as.factor(numOneBackFPDiff)) %>% # to summarise RTs according to value of FP n-1
  ungroup()

# Create column for previous orientation
data <- data %>%
  mutate(prevOri = ifelse(lag(orientation)==orientation, 'same', 'different')) %>%
  mutate(prevOri = as.factor(prevOri))
  

# Remove trials without n-1 FP values (i.e., first of each block)
data <- data %>%
  filter(!is.na(oneBackFP), !is.na(twoBackFP))
#filter(!is.na(oneBackFP))

# Keep only go trials with correct responses to analyze RT
data <- data %>%
  filter(!is.na(RT), Acc == 1)

# Coerce foreperiod and FP n-1 back to numeric
data$numForeperiod <- as.numeric(as.character(data$foreperiod))
data$numOneBackFP <- as.numeric(as.character(data$oneBackFP))

# Create log10 of continuous indenpendent variables
data$logFP <- log10(data$numForeperiod)
data$logOneBackFP <- log10(data$numOneBackFP)

# Remove extreme values
data <- data %>%
  filter(RT < 1000) %>%
  filter(RT > 150)

# Transform RT to reduce skew
data$logRT <- ifelse(!is.na(data$RT), log10(data$RT), NA) # log-transform
data$invRT <- ifelse(!is.na(data$RT), 1/data$RT, NA)

# Trimming
data2 <- data %>%
  group_by(ID) %>%
  mutate(RTzscore=ifelse(!is.na(RT), compute_zscore(RT), NA),
         logRTzscore=ifelse(!is.na(RT), compute_zscore(logRT), NA)) %>%
  filter(abs(logRTzscore) < 3) %>%
  ungroup()

# No trimming
data <- data %>%
  group_by(ID) %>%
  mutate(RTzscore=ifelse(!is.na(RT), compute_zscore(RT), NA),
         logRTzscore=ifelse(!is.na(RT), compute_zscore(logRT), NA)) %>%
  #filter(abs(logRTzscore) < 3) %>%
  ungroup()


# Average data
summaryData <- data %>%
  group_by(ID,foreperiod,condition,
           orientation,prevOri,
           oneBackFP,twoBackFP,oneBackFPDiff,
           block,counterbalance) %>%
  summarise(meanRT = mean(RT),
            meanLogRT = mean(logRT),
            meanRTzscore = mean(RTzscore),
            meanInvRT = mean(invRT),
            meanSeqEff = mean(oneBackEffect)) %>%
  ungroup() %>%
  mutate(numForeperiod=as.numeric(as.character(foreperiod)),
         numOneBackFPDiff=as.numeric(as.character(oneBackFPDiff)))

summaryData2 <- data2 %>%
  group_by(ID,foreperiod,condition,
           orientation,prevOri,
           oneBackFP,twoBackFP,oneBackFPDiff,
           block,counterbalance) %>%
  summarise(meanRT = mean(RT),
            meanLogRT = mean(logRT),
            meanRTzscore = mean(RTzscore),
            meanInvRT = mean(invRT),
            meanSeqEff = mean(oneBackEffect)) %>%
  ungroup() %>%
  mutate(numForeperiod=as.numeric(as.character(foreperiod)),
         numOneBackFPDiff=as.numeric(as.character(oneBackFPDiff)))