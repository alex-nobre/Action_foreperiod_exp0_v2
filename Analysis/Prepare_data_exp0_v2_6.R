
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

source("./Analysis/helper_functions.R")

# Read data
data <- read_csv("./Analysis/dataActionFPAll.csv")

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
  ungroup()

# Create column for difference between current and previous FP as a numeric variable:
data <- data %>%
  group_by(ID, block) %>%
  mutate(numOneBackFPDiff = c(NA, diff(numForeperiod))) %>%
  mutate(oneBackFPDiff = as.factor(numOneBackFPDiff), # to summarise RTs according to value of FP n-1
         squaredNumOneBackFPDiff = numOneBackFPDiff^2) %>% # quadratic term for difference between FP and FP n-1
  ungroup()

# Create column for previous orientation and for comparison of current and previous orientations
data <- data %>%
  mutate(seqOri = ifelse(lag(orientation)==orientation, 'same', 'different'),
         prevOri = lag(orientation)) %>%
  mutate(seqOri = as.factor(seqOri),
         prevOri = as.factor(prevOri))
  

# Remove trials without n-1 FP values (i.e., first of each block)
data <- data %>%
  filter(!is.na(oneBackFP), !is.na(twoBackFP))
#filter(!is.na(oneBackFP))

# Save data with error trials to assess accuracy
dataAll <- data

# Keep only trials with correct responses to analyze RT
data <- data %>%
  filter(!is.na(RT), Acc == 1)

# Create log10 of continuous indenpendent variables
data$numLogFP <- log10(data$numForeperiod)
data$logFP <- as.factor(data$numLogFP)
data$logOneBackFP <- log10(data$numOneBackFP)

dataAll$numLogFP <- log10(dataAll$numForeperiod)
dataAll$logFP <- as.factor(dataAll$numLogFP)
dataAll$logOneBackFP <- log10(dataAll$numOneBackFP)

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

# Create scaled predictors
data$scaledNumForeperiod <- scale(data$numForeperiod)[,1]
data$squaredScaledNumForeperiod <- data$scaledNumForeperiod^2
data$scaledNumOneBackFP <- scale(data$numOneBackFP)[,1]
data <- data %>%
  group_by(ID, block) %>%
  mutate(scaledNumOneBackFPDiff = scale(numOneBackFPDiff)[,1]) %>%
  ungroup()
data$squaredScaledNumOneBackFPDiff = data$scaledNumOneBackFPDiff^2

data2$scaledNumForeperiod <- scale(data2$numForeperiod)[,1]
data2$squaredScaledNumForeperiod <- data2$scaledNumForeperiod^2
data2$scaledNumOneBackFP <- scale(data2$numOneBackFP)[,1]
data2 <- data2 %>%
  group_by(ID, block) %>%
  mutate(scaledNumOneBackFPDiff = scale(numOneBackFPDiff)[,1]) %>%
  ungroup()
data2$squaredScaledNumOneBackFPDiff = data2$scaledNumOneBackFPDiff^2

dataAll$scaledNumForeperiod <- scale(dataAll$numForeperiod)[,1]
dataAll$squaredScaledNumForeperiod <- dataAll$scaledNumForeperiod^2
dataAll$scaledNumOneBackFP <- scale(dataAll$numOneBackFP)[,1]
dataAll <- dataAll %>%
  group_by(ID, block) %>%
  mutate(scaledNumOneBackFPDiff = scale(numOneBackFPDiff)[,1]) %>%
  ungroup()
dataAll$squaredScaledNumOneBackFPDiff = dataAll$scaledNumOneBackFPDiff^2


# Average data
summaryData <- data %>%
  group_by(ID,foreperiod,logFP,condition,
           orientation,prevOri,seqOri,
           oneBackFP,twoBackFP,oneBackFPDiff,
           block,counterbalance, testPos) %>%
  summarise(meanRT = mean(RT),
            meanLogRT = mean(logRT),
            meanRTzscore = mean(RTzscore),
            meanInvRT = mean(invRT),
            meanSeqEff = mean(oneBackEffect)) %>%
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
           oneBackFP,twoBackFP,oneBackFPDiff,
           block,counterbalance, testPos) %>%
  summarise(meanRT = mean(RT),
            meanLogRT = mean(logRT),
            meanRTzscore = mean(RTzscore),
            meanInvRT = mean(invRT),
            meanSeqEff = mean(oneBackEffect)) %>%
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


summaryDataAll <- dataAll %>%
  group_by(ID,foreperiod,condition,
           oneBackFP) %>%
  summarise(meanRT = mean(RT),
            meanAcc = mean(Acc),
            meanSeqEff = mean(oneBackEffect)) %>%
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
