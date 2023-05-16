# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 15:23:42 2022

@author: alpno
"""

import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt

subFileName = input('Enter File Name:\n')
subFile = './Data/' + subFileName

# Read 
subData = pd.read_csv(subFile)

# Remove unnecessary columns
subData = subData[['participant', 'date', 'Response.corr', 'blockCondition', 'block', 'orientation', 'foreperiod', 'corrAns', 'Response.rt', 'action_trigger.rt', 'Response.keys', 'counterbalance', 'extFixationDuration']]

# Rename condition columns for clarity
subData=subData.rename(columns={'blockCondition':'condition'})
subData=subData.rename(columns={'Response.rt':'RT'})

# Create columns for n-1 and n-2 foreperiods by block
subData['oneBackFP'] = subData.groupby(['block'])['foreperiod'].shift(1)
subData['twoBackFP'] = subData.groupby(['block'])['foreperiod'].shift(2)

# Create columns for n-1 and n-2 trial types
subData['prevOri'] = subData.groupby(['block'])['orientation'].shift(1)
subData['prevOri'] = subData.groupby(['block'])['orientation'].shift(2)

# Remove practice trials and lines with nan values
subData = subData[(subData['condition'] != 'practice') & (subData['condition'].notnull())]

# Check n and % of errors
print(len(subData[subData['Response.corr'] == 0]))
print((len(subData[subData['Response.corr'] == 0])/len(subData)) * 100)

# Errors by condition
print(len(subData[(subData['orientation']=='left') & (subData['Response.corr']==0)])/
      len(subData[subData['orientation']=='left'])*100)

print(len(subData[(subData['orientation']=='right') & (subData['Response.corr']==0)])/
      len(subData[subData['orientation']=='right'])*100)


# Keep only go trials with correct responses to analyze RT
subData = subData[(subData['RT'].notnull()) & (subData['Response.corr'] == 1)]

# Remove trials without n-1 FP values (i.e., first of each block)
subData = subData[subData['oneBackFP'].notnull()]

# Remove outliers
print('notclean: ' + str(len(subData)))
subData = subData[(subData['RT'] < 1) & (subData['RT'] > 0.1)]
print('clean: ' + str(len(subData)))


summaryData=subData.groupby(['foreperiod','condition','oneBackFP','prevOri','block','counterbalance'],
                           as_index=False)[['RT','Response.corr']].mean()

summaryPlot=sns.pointplot(x="foreperiod", y="RT",
                          data=summaryData)
                          
plt.show()
