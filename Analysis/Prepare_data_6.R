
#==============================================================================#
# Changes
# Replaced column for FP n-1 diff, since it will be used as a continuous variable
# now

# Use RT in ms instead of s, so that logRT does not assume negative values
# log of predictors is now log10 instead of log2

# Created columns for quadratic trend for FP and difference between
# FP and FP n-1
# Created column for previous orientation and renamed column for comparison
# between current and previous orientation
# Created column for trial number

# Create separate columns for logFP as numerical and categorical

# Create column for scaled predictors
#==============================================================================#

# Load necessary packages
library(readr)
library(ggplot2)
library(magrittr)
library(dplyr)
library(forcats)

source("./Analysis/helper_functions.R")

# Read data
data <- read_csv("./Analysis/dataActionFPAll.csv")

# Remove unnecessary columns
data <- data %>%
  dplyr::select(ID, Acc, condition, block, orientation,
                foreperiod, RT, counterbalance, 
                extFixationDuration, action_trigger.rt, ITI,
                oneBackFP, twoBackFP, oneBackEffect)

# Coerce to factors
data <- data %>%
  mutate(across(c(ID, condition, block, orientation, foreperiod, counterbalance, oneBackFP, twoBackFP), as_factor))

data$condition <- data$condition %>%
  fct_relevel(c("external", "action"))

# Column for testing position based on counterbalancing
data <- data %>%
  mutate(testPos = case_when((condition == 'action' & counterbalance == 'action-external') ~ '1',
                             (condition == 'external' & counterbalance == 'action-external') ~ '2',
                             (condition == 'action' & counterbalance == 'external-action') ~ '2',
                             (condition == 'external' & counterbalance == 'external-action') ~ '1')) %>%
  mutate(testPos = as.factor(testPos))

# Create numeric versions of foreperiod and FP n-1
data$numForeperiod <- as.numeric(as.character(data$foreperiod))
data$numOneBackFP <- as.numeric(as.character(data$oneBackFP))

# Quadratic term for numForeperiod
data$squaredNumForeperiod <- data$numForeperiod^2

# Remove practice trials
data <- data %>%
  filter(condition != 'practice')

# Create column for trial number
data <- data %>%
  group_by(ID) %>%
  mutate(trial = seq(1,n())) %>%
  ungroup() %>%
  group_by(ID, block) %>%
  mutate(trial_bl = seq(1, n())) %>%
  ungroup()

# Create columns for numerical and categorical difference between current and previous FP, and for 
# whether the previous FP is longer than the current:
# If FPn-1 is shorter than the current FP, value is 0
# If FPn-1 is equal to the current FP, value is 0
# If FPn-1 is longer than the current FP, value is 1
data <- data %>%
  group_by(ID, block) %>%
  mutate(numOneBackFPDiff = c(NA, diff(numForeperiod))) %>%
  mutate(oneBackFPDiff = as.factor(numOneBackFPDiff), # to summarise RTs according to value of FP n-1
         squaredNumOneBackFPDiff = numOneBackFPDiff^2) %>% # quadratic term for difference between FP and FP n-1
  mutate(prevFPLonger = case_when(numOneBackFPDiff>0 ~ "0",
                                  numOneBackFPDiff<0 ~ "1",
                                  numOneBackFPDiff==0 ~ "0")
  ) %>%
  mutate(prevFPLonger = as.factor(prevFPLonger)) %>%
  mutate(prevFPLonger = fct_relevel(prevFPLonger, c("0", "1"))) %>%
  ungroup()


# Create column for previous orientation and for comparison of current and previous orientations
data <- data %>%
  mutate(seqOri = ifelse(lag(orientation)==orientation, 'same', 'different'),
         prevOri = lag(orientation)) %>%
  mutate(seqOri = as.factor(seqOri),
         prevOri = as.factor(prevOri))

# Create column for total ITI depending on condition
data <- data %>%
  mutate(ITItotal = ifelse(condition == 'action', ITI + action_trigger.rt, ITI + extFixationDuration))


# Remove trials without n-1 FP values (i.e., first of each block)
data <- data %>%
  filter(!is.na(oneBackFP)) #%>%
# filter(!is.na(oneBackFP), !is.na(twoBackFP))


# Save data with error trials to assess accuracy
dataAcc <- data

# Keep only trials with correct responses to analyze RT
data <- data %>%
  filter(!is.na(RT), Acc == 1)

# Create variable for error rate
dataAcc$Error <- abs(dataAcc$Acc - 1)

# Create log10 of continuous indenpendent variables
data$numLogFP <- log10(data$numForeperiod)
data$logFP <- as.factor(data$numLogFP)
data$logOneBackFP <- log10(data$numOneBackFP)

dataAcc$numLogFP <- log10(dataAcc$numForeperiod)
dataAcc$logFP <- as.factor(dataAcc$numLogFP)
dataAcc$logOneBackFP <- log10(dataAcc$numOneBackFP)

# Factor version of accuracy/error rate
dataAcc$acc_result <- as.factor(dataAcc$Acc)
dataAcc$error_result <- as.factor(dataAcc$Error)

# Remove extreme values
ntrials_before_extrem <- nrow(data)
data <- data %>%
  filter(RT < 1.0) %>%
  filter(RT > 0.15)
ntrials_after_extrem <- nrow(data)

# Transform RT to reduce skew
data$logRT <- ifelse(!is.na(data$RT), log10(data$RT), NA) # log-transform
data$invRT <- ifelse(!is.na(data$RT), 1/data$RT, NA)

# nrow(data)
# # Remove trials with action.trigger_rt too long
# data <- data %>%
#   filter(action_trigger.rt < 9000)
# nrow(data)

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

# Create scaled predictors
data$scaledNumForeperiod <- scale(data$numForeperiod, scale = FALSE)[,1]
data$squaredScaledNumForeperiod <- data$scaledNumForeperiod^2
data$scaledNumOneBackFP <- scale(data$numOneBackFP, scale = FALSE)[,1]
data <- data %>%
  group_by(ID, block) %>%
  mutate(scaledNumOneBackFPDiff = scale(numOneBackFPDiff, scale = FALSE)[,1]) %>%
  ungroup()
data$squaredScaledNumOneBackFPDiff = data$scaledNumOneBackFPDiff^2

data2$scaledNumForeperiod <- scale(data2$numForeperiod, scale = FALSE)[,1]
data2$squaredScaledNumForeperiod <- data2$scaledNumForeperiod^2
data2$scaledNumOneBackFP <- scale(data2$numOneBackFP, scale = FALSE)[,1]
data2 <- data2 %>%
  group_by(ID, block) %>%
  mutate(scaledNumOneBackFPDiff = scale(numOneBackFPDiff, scale = FALSE)[,1]) %>%
  ungroup()
data2$squaredScaledNumOneBackFPDiff = data2$scaledNumOneBackFPDiff^2

dataAcc$scaledNumForeperiod <- scale(dataAcc$numForeperiod, scale = FALSE)[,1]
dataAcc$squaredScaledNumForeperiod <- dataAcc$scaledNumForeperiod^2
dataAcc$scaledNumOneBackFP <- scale(dataAcc$numOneBackFP, scale = FALSE)[,1]
dataAcc <- dataAcc %>%
  group_by(ID, block) %>%
  mutate(scaledNumOneBackFPDiff = scale(numOneBackFPDiff, scale = FALSE)[,1]) %>%
  ungroup()
dataAcc$squaredScaledNumOneBackFPDiff = dataAcc$scaledNumOneBackFPDiff^2


###############################################################
# Add delay data
###############################################################
# delayData <- read_csv("./Analysis/delayDataAll.csv") %>%
#   mutate(across(c(ID, condition), as_factor)) %>%
#   select(-condition)

# data <- data %>%
#   filter(ID %!in% c("007"))
# data2 <- data2 %>%
#   filter(ID %!in% c("007"))
# dataAcc <- dataAcc %>%
#   filter(ID %!in% c("007"))

# data <- inner_join(data, delayData, by = c("trial", "ID"))
# data2 <- inner_join(data2, delayData, by = c("trial", "ID"))
# dataAcc <- inner_join(dataAcc, delayData, by = c("trial", "ID"))

################################################################

# Average data
summaryData <- data %>%
  group_by(ID,foreperiod,logFP,condition,
           orientation,prevOri,seqOri,
           oneBackFP,twoBackFP,oneBackFPDiff,prevFPLonger,
           block,counterbalance, testPos) %>%
  summarise(meanRT = mean(RT),
            meanLogRT = mean(logRT),
            meanRTzscore = mean(RTzscore),
            meanInvRT = mean(invRT),
            meanSeqEff = mean(oneBackEffect),
            meanITITotal = mean(ITItotal)) %>% #,
            #meanDelay = mean(delay)) %>%
  ungroup() %>%
  mutate(numForeperiod = as.numeric(as.character(foreperiod)),
         numOneBackFP = as.numeric(as.character(oneBackFP)),
         numOneBackFPDiff = as.numeric(as.character(oneBackFPDiff)),
         numLogFP = as.numeric(as.character(logFP))) %>%
  mutate(squaredNumForeperiod = numForeperiod^2,
         squaredNumOneBackFPDiff = numOneBackFPDiff^2,
         squaredNumLogFP = numLogFP^2,
         scaledNumForeperiod = scale(numForeperiod)[,1],
         squaredScaledNumForeperiod = scaledNumForeperiod^2,
         scaledNumOneBackFP = scale(numOneBackFP)[,1],
         scaledNumOneBackFPDiff = scale(numOneBackFPDiff)[,1],
         squaredScaledNumOneBackFPDiff = scaledNumOneBackFPDiff^2)

summaryData2 <- data2 %>%
  group_by(ID,foreperiod,logFP,condition,
           orientation,prevOri,seqOri,
           oneBackFP,twoBackFP,oneBackFPDiff,prevFPLonger,
           block,counterbalance, testPos) %>%
  summarise(meanRT = mean(RT),
            meanLogRT = mean(logRT),
            meanRTzscore = mean(RTzscore),
            meanInvRT = mean(invRT),
            meanSeqEff = mean(oneBackEffect),
            meanITITotal = mean(ITItotal)) %>% #,
            #meanDelay = mean(delay)) %>%
  ungroup() %>%
  mutate(numForeperiod = as.numeric(as.character(foreperiod)),
         numOneBackFP = as.numeric(as.character(oneBackFP)),
         numOneBackFPDiff = as.numeric(as.character(oneBackFPDiff)),
         numLogFP = as.numeric(as.character(logFP))) %>%
  mutate(squaredNumForeperiod = numForeperiod^2,
         squaredNumOneBackFPDiff = numOneBackFPDiff^2,
         squaredNumLogFP = numLogFP^2,
         scaledNumForeperiod = scale(numForeperiod)[,1],
         squaredScaledNumForeperiod = scaledNumForeperiod^2,
         scaledNumOneBackFP = scale(numOneBackFP)[,1],
         scaledNumOneBackFPDiff = scale(numOneBackFPDiff)[,1],
         squaredScaledNumOneBackFPDiff = scaledNumOneBackFPDiff^2)


summaryDataAcc <- dataAcc %>%
  group_by(ID,foreperiod,condition,
           oneBackFP) %>%
  summarise(meanRT = mean(RT),
            meanAcc = mean(Acc),
            errorRate = mean(Error),
            meanSeqEff = mean(oneBackEffect),
            meanITITotal = mean(ITItotal)) %>% #,
            #meanDelay = mean(delay)) %>%
  ungroup() %>%
  mutate(numForeperiod = as.numeric(as.character(foreperiod)),
         numOneBackFP = as.numeric(as.character(oneBackFP))) %>%
  mutate(squaredNumForeperiod = numForeperiod^2,
         scaledNumForeperiod = scale(numForeperiod)[,1],
         squaredScaledNumForeperiod = scaledNumForeperiod^2,
         scaledNumOneBackFP = scale(numOneBackFP)[,1])

#write_csv(data, "./Analysis/data.csv")
#write_csv(data2, "./Analysis/data2.csv")
#write_csv(summaryData, "./Analysis/summaryData.csv")
#write_csv(summaryData2, "./Analysis/summaryData2.csv")
