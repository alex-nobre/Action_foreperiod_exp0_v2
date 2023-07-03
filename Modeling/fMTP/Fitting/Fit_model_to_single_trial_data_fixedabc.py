# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 11:21:35 2023

@author: alpno

Changes:
    - changed fmtp code to pass go/no-go trial sequence as argument to FPgonogo function
    - fixed issue of nan values in computation of z-scores

"""

import numpy as np
import pandas as pd
import os
from scipy import optimize as op
from scipy import stats
from sklearn.metrics import mean_squared_error

# Import for plotting
import matplotlib.pyplot as plt

# Import classes
import sys
sys.path.insert(0, './Modeling/fMTP/Fitting') # path to where classes are stored

from fmtp_single_trial import fMTP, FPexp
from hazard import fMTPhz
from fit import sort_fit, show_fit, get_fit 

#=================================== functions =============================

# def sse (parm, emp, sim):
#     SSE = sum(pow(emp.RT - simRTs, 2))
#     return SSE
# 
# def rmse (parm, emp, sim):
#     simRTs = temp2RT(sim.prep, parm[0], parm[1])
#     simzRTs = stats.zscore(simRTs)
#     RMSE = np.sqrt(mean_squared_error(emp.zRT, simzRTs))
#     #RMSE = np.sqrt(((emp.zRT - simzRTs)**2).mean())
#     #RMSE = np.sqrt(((emp.RT - simRTs)**2).mean())
#     return RMSE

def oneMinusCorr (parm, kValue, emp, conditionSimExp, conditionNaIdx):
    fmtp = fMTP(parm[0], parm[1], kValue)
    
    # Simulate experiment for condition
    simCondition, prep_fmtpCondition = conditionSimExp.run_exp(fmtp)
    
    # Average preparation at discrete FPs
    simCondition = simCondition.iloc[1:, :] # first trial has no preparation
    
    simCondition['distrib'] = distr
    
    # Remove nan's
    simCondition = simCondition.drop(index=conditionNaIdx)
    
    # Compute 1 - correlation between RT and prep
    negRT = (-1) * emp.RT
    prepCorr = np.corrcoef(negRT, simCondition.prep)[0][1]
    oneMinusCorr = 1 - prepCorr
    return oneMinusCorr, 
    

#============================================================================#
#========================= 1. All free parameters  ===========================    
#============================================================================#

empData = pd.read_csv('./Analysis/dataActionFPAll.csv')

IDList=np.unique(empData['ID']).tolist()

empData.RT=empData.RT*1000

# List of trial indices for single trial
rowsByCond = int((len(empData)/len(empData.ID.unique()))/2)

empDataAction = empData.loc[empData.condition=='action'].reset_index()
trialIndexAction = [[i for i in range(1,rowsByCond + 1)] for sub in range(int(len(empDataAction)/rowsByCond))]
trialIndexAction = [index for ID in trialIndexAction for index in ID]
empDataAction['trial']=trialIndexAction

empDataExternal = empData.loc[empData.condition=='external'].reset_index()
trialIndexExternal = [[i for i in range(1,rowsByCond + 1)] for sub in range(int(len(empDataExternal)/rowsByCond))]
trialIndexExternal = [index for ID in trialIndexExternal for index in ID]
empDataExternal['trial']=trialIndexExternal

# Generate experiment from model
FPs = np.arange(1.0, 3.0, 0.6)
distr = 'uni'

# List of k values (cannot be minimized in scipy.optimized since only integers are acceptable in model)
kList= [i for i in range(1,15)]

# Starting values for r and c
startVals=[-2.81, 0.0001]

# Lists to store results
hazardResults=[None]*(len(IDList) * len(kList))
seqEffResults=[None]*(len(IDList) * len(kList))


fitEl=0  
for i, iID in enumerate(IDList):
    print(i)
    
    # Action and external datasets for part
    IDDataAction=empDataAction.loc[empDataAction.ID==iID].reset_index()
    IDDataExternal=empDataExternal.loc[empDataExternal.ID==iID].reset_index()
    
    # FP lists for part
    IDFPAction=IDDataAction['foreperiod'].values.tolist()
    IDFPExternal=IDDataExternal['foreperiod'].values.tolist()
    
    # Remove first trial of each dataset for simulation
    IDDataAction = IDDataAction.iloc[1:,:]
    IDDataExternal = IDDataExternal.iloc[1:,:]
    
    # Remove nan's
    actionNaIdx = IDDataAction[IDDataAction['RT'].isna()].index.tolist()
    IDDataAction = IDDataAction.dropna(subset=['RT'])
    
    externalNaIdx = IDDataExternal[IDDataExternal['RT'].isna()].index.tolist()
    IDDataExternal = IDDataExternal.dropna(subset=['RT'])
    
    # Exps for action and external conditions for part
    actionExp = FPexp(FPs = FPs, FPList = IDFPAction, distribution = distr, tr_per_block = rowsByCond)
    externalExp = FPexp(FPs = FPs, FPList = IDFPExternal, distribution = distr, tr_per_block = rowsByCond)


    # Minimize error measure
    for ik in range(len(kList)):
        
        # Discrete k value
        k = kList[ik]
        
        #=========== Action ============#
        # Minimize r and c within this k value
        oneMinusCorrAction = op.minimize(oneMinusCorr, startVals, args = (k, IDDataAction, actionExp, actionNaIdx), bounds = ([-10, -0.1], [0, 0.0021]), method = 'L-BFGS-B')
        
        # Compute R2 from sim using minimized parameters
        fmtp = fMTP(oneMinusCorrAction.x[0], oneMinusCorrAction.x[1], k)
        simAction, prep_fmtpAction = actionExp.run_exp(fmtp)
        simAction = simAction.iloc[1:, :] # first trial has no preparation
        
        simAction['distrib'] = distr
        
        # Remove nan's
        simAction = simAction.drop(index=actionNaIdx)    
        RoneMinusCorrAction = np.corrcoef(simAction.prep, IDDataAction.RT)[0][1]**2
        
        #============= External ============#
        # Minimize r and c within this k value
        oneMinusCorrExternal = op.minimize(oneMinusCorr, startVals, args = (k, IDDataExternal, externalExp, externalNaIdx), bounds = ([-10, -0.1], [0, 0.0021]), method = 'L-BFGS-B')
        
        # Compute R2 from sim using minimized parameters
        fmtp = fMTP(oneMinusCorrExternal.x[0], oneMinusCorrExternal.x[1], k)
        simExternal, prep_fmtpExternal = externalExp.run_exp(fmtp)
        simExternal = simExternal.iloc[1:, :] # first trial has no preparation
        
        simExternal['distrib'] = distr
        
        # Remove nan's
        simExternal = simExternal.drop(index = externalNaIdx)
        RoneMinusCorrExternal = np.corrcoef(simExternal.prep, IDDataExternal.RT)[0][1]**2
        
        # Save to list
        hazardResults[fitEl]=(iID, k,
                           # sseAction.x[0], sseAction.x[1], sseAction.fun, RsseAction, 
                           # rmseAction.x[0], rmseAction.x[1], rmseAction.fun, RrmseAction,
                           oneMinusCorrAction.x[0], oneMinusCorrAction.x[1], oneMinusCorrAction.fun, RoneMinusCorrAction,
                           # sseExternal.x[0], sseExternal.x[1], sseExternal.fun, RsseExternal, 
                           # rmseExternal.x[0], rmseExternal.x[1], rmseExternal.fun, RrmseExternal,
                           oneMinusCorrExternal.x[0], oneMinusCorrExternal.x[1], oneMinusCorrExternal.fun, RoneMinusCorrExternal)
        
        fitEl+=1


# Join in a data frame
hazardResultsDF = pd.DataFrame(hazardResults,
                            columns = ['ID', 'k', 'rAction', 'cAction', 'One_minus_corr_Action', 'R2_1-corr_Action',
                                       'rExternal', 'cExternal', 'One_minus_corr_External', 'R2_1-corr_External'])

hazardResultsDF.to_csv('./Modeling/fMTP/Fitting/'+'fittingResults_krc.csv')                            
                            
