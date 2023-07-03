# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 15:46:34 2023

@author: alpno

Plots parameter values obtained from fitting process to check individual fitting.

Parameters for sequential effects are plotted in separate panel for each foreperiod for better visualization
"""

import pandas as pd
import numpy as np
import itertools
import os
import seaborn as sns
from scipy import stats

# Import for plotting
import matplotlib.pyplot as plt

# Import classes
import sys
sys.path.insert(0, './Modeling/fMTP/Fitting') # path to where classes are stored

from fmtp import fMTP, FPexp

#========================================= Functions =============================
def temp2RT (prepValue, a, b):
    return (a*prepValue + b)
  
#============================================================================#
#========================= 1. All free parameters  ===========================    
#============================================================================#

hazardResultsDF = pd.read_csv('./Modeling/fMTP/Fitting/fittingResults.csv')

# Find min for SSE in each condition
stCoefsAction = hazardResultsDF.loc[hazardResultsDF.groupby('ID').SSEAction.idxmin()].reset_index()
stCoefsAction = stCoefsAction[['ID', 'k', 'r', 'c', 'a SSE Action', 'b SSE Action', 'SSEAction']]
stCoefsExternal = hazardResultsDF.loc[hazardResultsDF.groupby('ID').SSEExternal.idxmin()].reset_index()
stCoefsExternal = stCoefsExternal[['k', 'r', 'c', 'a SSE External', 'b SSE External', 'SSEExternal']]

stCoefs = pd.DataFrame(pd.concat([stCoefsAction.rename(columns = {'ID':'ID', 'k':'kAction', 'r':'rAction', 'c':'cAction', 
                                                                  'a SSE Action':'a SSE Action', 'b SSE Action':'b SSE Action', 
                                                                  'SSEAction':'SSE_Action'}), 
                                  stCoefsExternal.rename(columns = {'k':'kExternal', 'r':'rExternal', 'c':'cExternal', 
                                                                    'a SSE External':'a SSE External', 
                                                                    'b SSE External':'b SSE External', 
                                                                    'SSEExternal': 'SSE_External'})], 
                                 axis = 1), columns = ['ID', 'kAction', 'rAction', 'cAction', 'a SSE Action', 'b SSE Action', 'SSE_Action',
                                                       'kExternal', 'rExternal', 'cExternal', 'a SSE External', 'b SSE External', 'SSE_External'])

# Plot SSE values
plt.scatter(stCoefs.index, stCoefs.SSE_Action)
plt.show()

# Check if they match (prob not)
fig = plt.figure()
ax = fig.subplots(3,2, sharey=True)
ax = ax.ravel()

titles = ['kAction', 'rAction', 'cAction', 'kExternal', 'rExternal', 'cExternal']

for i, panel in enumerate(ax):
    panel.hist(stCoefs.iloc[:,i+1])
    panel.set_title(titles[i])

plt.show()    

# Build plots with individual data and predictions
FPs = np.arange(1.0, 3.0, 0.6)
distr = 'uni'
exp = FPexp(FPs = FPs, distribution = distr, tr_per_block = 150)


# Plot by foreperiod only
empData = pd.read_csv('./Analysis/dataActionFPAll.csv')
empData.RT = empData.RT * 1000 # transform RT to ms
empData = empData.groupby(['ID', 'foreperiod', 'condition'], as_index=False)[['RT']].mean()

fig = plt.figure(constrained_layout = True, figsize = (15.96, 12))
ax = fig.subplots(5,6,sharey=False)
ax = ax.ravel()

titles = stCoefs.ID.tolist()

for iSub, sub in enumerate(titles):
    subData = stCoefs.iloc[iSub,:]#[stCoefs['ID']==sub]#
    
    # Subset participant's data by condition
    empSub = empData[empData['ID']==sub]
    empSubAction = empSub[empSub['condition']=='action']
    empSubExternal = empSub[empSub['condition']=='external']
    
    # Action condition
    kAction = subData['kAction']
    rAction = subData['rAction']
    cAction = subData['cAction']
    aAction = subData['a SSE Action']
    bAction = subData['b SSE Action']
    
    # Run model with current parameter values and simulate RTs for action condition
    fmtpAction=fMTP(rAction,cAction,kAction)
    state_discr_action, state_con_action = exp.run_exp(fmtpAction)
    state_discr_action=state_discr_action[1:]
    mean_state_action=state_discr_action.groupby(['FP']).mean().reset_index()
    mean_state_action['RT'] = temp2RT(mean_state_action.prep, aAction, bAction)
    
    # External condition
    kExternal = subData['kExternal']
    rExternal = subData['rExternal']
    cExternal = subData['cExternal']
    aExternal = subData['a SSE External']
    bExternal = subData['b SSE External']
    
    # Run model with current parameter values and simulate RTs for external condition
    fmtpExternal=fMTP(rExternal,cExternal,kExternal)
    state_discr_external, state_con_external = exp.run_exp(fmtpExternal)
    state_discr_external=state_discr_external[1:]
    mean_state_external=state_discr_external.groupby(['FP']).mean().reset_index()
    mean_state_external['RT'] = temp2RT(mean_state_external.prep, aExternal, bExternal)
    
    # Plot on single panel
    ax[iSub].plot(mean_state_action.FP, mean_state_action.RT,'.--', color = 'blue')
    ax[iSub].plot(empSubAction.foreperiod, empSubAction.RT, '.-', color = 'blue')
    ax[iSub].plot(mean_state_external.FP, mean_state_external.RT,'.--', color = 'orange')
    ax[iSub].plot(empSubExternal.foreperiod, empSubExternal.RT, '.-', color = 'orange')
    
fig.suptitle('ID fits')

figname='./Modeling/fMTP/Fitting/Plots/id_fit_SSE.png'
plt.savefig(figname,format='png')
plt.show()
    

#========================= Plot by sequential effects ===============================
empData = pd.read_csv('./Analysis/dataActionFPAll.csv')
empData.RT = empData.RT * 1000 # transform RT to ms
empData = empData.groupby(['ID', 'foreperiod', 'oneBackFP', 'condition'], as_index=False)[['RT']].mean()

fig = plt.figure(constrained_layout=True, figsize = (15.96,12))
subfigs=fig.subfigures(5,6)
subfigsCoords=list(itertools.product([0,1,2,3,4],[0,1,2,3,4,5]))

#subfigs = subfigs.ravel()
titles = stCoefs.ID.tolist()


for iSub in range(len(titles)):
    subData = stCoefs.iloc[iSub,:]
    empSub = empData[empData['ID']==titles[iSub]]
    empSubAction = empSub[empSub['condition']=='action']
    empSubExternal = empSub[empSub['condition']=='external']
    
    # Action condition
    kAction = subData['kAction']
    rAction = subData['rAction']
    cAction = subData['cAction']
    aAction = subData['a SSE Action']
    bAction = subData['b SSE Action']
    
    # Run model with current parameter values and simulate RTs for action condition
    fmtpAction=fMTP(rAction,cAction,kAction)
    state_discr_action, state_con_action = exp.run_exp(fmtpAction)
    state_discr_action=state_discr_action[1:]
    mean_state_action=state_discr_action.groupby(['FP', 'FPn_1']).mean().reset_index()
    mean_state_action['RT'] = temp2RT(mean_state_action.prep, aAction, bAction)
    
    #  External condition
    kExternal = subData['kExternal']
    rExternal = subData['rExternal']
    cExternal = subData['cExternal']
    aExternal = subData['a SSE External']
    bExternal = subData['b SSE External']
    
    # Run model with current parameter values and simulate RTs for external condition
    fmtpExternal=fMTP(rExternal,cExternal,kExternal)
    state_discr_external, state_con_external = exp.run_exp(fmtpExternal)
    state_discr_external=state_discr_external[1:]
    mean_state_external=state_discr_external.groupby(['FP', 'FPn_1']).mean().reset_index()
    mean_state_external['RT'] = temp2RT(mean_state_external.prep, aExternal, bExternal)
    
    # Plot on single panel
    ax=subfigs[subfigsCoords[iSub]].subplots(1,4,sharey=True)
    
    for idx, iFP in enumerate(np.unique(mean_state_action.FPn_1)):
        FP = round(iFP, 2)
        state_n1_action = mean_state_action[mean_state_action.FPn_1 == iFP]
        emp_n1_action = empSubAction[empSubAction.oneBackFP == FP]
        ax[idx].plot(state_n1_action.FP, state_n1_action.RT, '.--', color = 'blue')
        ax[idx].plot(emp_n1_action.foreperiod, emp_n1_action.RT, '.-', color = 'blue')
        
        state_n1_external = mean_state_external[mean_state_external.FPn_1 == iFP]
        emp_n1_external = empSubExternal[empSubExternal.oneBackFP == FP]
        ax[idx].plot(state_n1_external.FP, state_n1_external.RT, '.--', color = 'orange')
        ax[idx].plot(emp_n1_external.foreperiod, emp_n1_external.RT, '.-', color = 'orange')
        
        ax[idx].set_title(FP)
    
    subfigs[subfigsCoords[iSub]].suptitle(titles[iSub])

fig.suptitle('Sequential effects')

figname='./Modeling/fMTP/Fitting/Plots/seq_effects_id_fit_SSE.png'
plt.savefig(figname,format='png', bbox_inches = 'tight')

plt.show()


#=============== Compare parameters between conditions =========================

# Build dataset in long format for plotting
stCoefsAction = hazardResultsDF.loc[hazardResultsDF.groupby('ID').SSEAction.idxmin()].reset_index()
stCoefsAction = stCoefsAction[['ID', 'k', 'r', 'c', 'a SSE Action', 'b SSE Action', 'SSEAction']]
stCoefsAction['condition'] = pd.Series(['action'] * len(stCoefsAction))
stCoefsExternal = hazardResultsDF.loc[hazardResultsDF.groupby('ID').SSEExternal.idxmin()].reset_index()
stCoefsExternal = stCoefsExternal[['ID', 'k', 'r', 'c', 'a SSE External', 'b SSE External', 'SSEExternal']]
stCoefsExternal['condition'] = pd.Series(['external'] * len(stCoefsExternal))

stCoefsLong = pd.DataFrame(pd.concat([stCoefsAction.rename(columns = {'ID':'ID', 'condition':'condition', 'k':'k', 'r':'r', 'c':'c', 
                                                                   'a SSE Action':'aSSE', 'b SSE Action':'bSSE', 
                                                                   'SSEAction':'SSE'}), 
                                   stCoefsExternal.rename(columns = {'ID':'ID', 'condition':'condition', 'k':'k', 'r':'r', 'c':'c', 
                                                                     'a SSE External':'aSSE', 
                                                                     'b SSE External':'bSSE', 
                                                                     'SSEExternal': 'SSE'})], 
                                  axis = 0), columns = ['ID', 'condition', 'k', 'r', 'c', 'aSSE', 'bSSE', 'SSE'])
                                  

paramFig = plt.figure()
ax = paramFig.subplots(2, 2)

# Point plots
sns.stripplot(x = 'condition', y = 'k', data = stCoefsLong, jitter = 0.1, orient = 'v', ax = ax[0,0])
sns.pointplot(x = 'condition', y = 'k', data = stCoefsLong, ax = ax[0,0])

sns.stripplot(x = 'condition', y = 'r', data = stCoefsLong, jitter = 0.1, orient = 'v', ax = ax[0,1])
sns.pointplot(x = 'condition', y = 'r', data = stCoefsLong, ax = ax[0,1])

sns.stripplot(x = 'condition', y = 'c', data = stCoefsLong, jitter = 0.1, orient = 'v', ax = ax[1,0])
sns.pointplot(x = 'condition', y = 'c', data = stCoefsLong, ax = ax[1,0])


plt.show()

# Histograms
paramFig = plt.figure()
ax = paramFig.subplots(3, 2)

# Point plots
sns.histplot(x = "k", data = stCoefsLong[stCoefsLong.condition == "action"], ax = ax[0, 0])
sns.histplot(x = "k", data = stCoefsLong[stCoefsLong.condition == "external"], ax = ax[0, 1])

sns.histplot(x = "r", data = stCoefsLong[stCoefsLong.condition == "action"], ax = ax[1, 0])
sns.histplot(x = "r", data = stCoefsLong[stCoefsLong.condition == "external"], ax = ax[1, 1])

sns.histplot(x = "c", data = stCoefsLong[stCoefsLong.condition == "action"], ax = ax[2, 0])
sns.histplot(x = "c", data = stCoefsLong[stCoefsLong.condition == "external"], ax = ax[2, 1])


plt.show()

# Paired samples t-test
stats.wilcoxon(stCoefs.kAction, stCoefs.kExternal)
stats.wilcoxon(stCoefs.rAction, stCoefs.rExternal)
stats.wilcoxon(stCoefs.cAction, stCoefs.cExternal)

#===============================================================================#
#============================ 2. c and r fixed ==================================
#===============================================================================#
hazardResultsDF = pd.read_csv('./Modeling/fMTP/Fitting/fittingResults_k.csv')

# Find min for SSE in each condition
stCoefsAction = hazardResultsDF.loc[hazardResultsDF.groupby('ID').SSEAction.idxmin()].reset_index()
stCoefsAction = stCoefsAction[['ID', 'k', 'r', 'c', 'a SSE Action', 'b SSE Action', 'SSEAction']]
stCoefsExternal = hazardResultsDF.loc[hazardResultsDF.groupby('ID').SSEExternal.idxmin()].reset_index()
stCoefsExternal = stCoefsExternal[['k', 'r', 'c', 'a SSE External', 'b SSE External', 'SSEExternal']]

stCoefs = pd.DataFrame(pd.concat([stCoefsAction.rename(columns = {'ID':'ID', 'k':'kAction', 'r':'rAction', 'c':'cAction', 
                                                                  'a SSE Action':'a SSE Action', 'b SSE Action':'b SSE Action', 
                                                                  'SSEAction':'SSE_Action'}), 
                                  stCoefsExternal.rename(columns = {'k':'kExternal', 'r':'rExternal', 'c':'cExternal', 
                                                                    'a SSE External':'a SSE External', 
                                                                    'b SSE External':'b SSE External', 
                                                                    'SSEExternal': 'SSE_External'})], 
                                 axis = 1), columns = ['ID', 'kAction', 'rAction', 'cAction', 'a SSE Action', 'b SSE Action', 'SSE_Action',
                                                       'kExternal', 'rExternal', 'cExternal', 'a SSE External', 'b SSE External', 'SSE_External'])

  

# Build plots with individual data and predictions
FPs = np.arange(1.0, 3.0, 0.6)
distr = 'uni'
exp = FPexp(FPs = FPs, distribution = distr, tr_per_block = 150)


# Plot by foreperiod only
empData = pd.read_csv('./Analysis/dataActionFPAll.csv')
empData.RT = empData.RT * 1000 # transform RT to ms
empData = empData.groupby(['ID', 'foreperiod', 'condition'], as_index=False)[['RT']].mean()

fig = plt.figure(constrained_layout = True, figsize = (15.96, 12))
ax = fig.subplots(5,6,sharey=False)
ax = ax.ravel()

titles = stCoefs.ID.tolist()

for iSub, sub in enumerate(titles):
    subData = stCoefs.iloc[iSub,:]#[stCoefs['ID']==sub]#
    
    # Subset participant's data by condition
    empSub = empData[empData['ID']==sub]
    empSubAction = empSub[empSub['condition']=='action']
    empSubExternal = empSub[empSub['condition']=='external']
    
    # Action condition
    kAction = subData['kAction']
    rAction = subData['rAction']
    cAction = subData['cAction']
    aAction = subData['a SSE Action']
    bAction = subData['b SSE Action']
    
    # Run model with current parameter values and simulate RTs for action condition
    fmtpAction=fMTP(rAction,cAction,kAction)
    state_discr_action, state_con_action = exp.run_exp(fmtpAction)
    state_discr_action=state_discr_action[1:]
    mean_state_action=state_discr_action.groupby(['FP']).mean().reset_index()
    mean_state_action['RT'] = temp2RT(mean_state_action.prep, aAction, bAction)
    
    # External condition
    kExternal = subData['kExternal']
    rExternal = subData['rExternal']
    cExternal = subData['cExternal']
    aExternal = subData['a SSE External']
    bExternal = subData['b SSE External']
    
    # Run model with current parameter values and simulate RTs for external condition
    fmtpExternal=fMTP(rExternal,cExternal,kExternal)
    state_discr_external, state_con_external = exp.run_exp(fmtpExternal)
    state_discr_external=state_discr_external[1:]
    mean_state_external=state_discr_external.groupby(['FP']).mean().reset_index()
    mean_state_external['RT'] = temp2RT(mean_state_external.prep, aExternal, bExternal)
    
    # Plot on single panel
    ax[iSub].plot(mean_state_action.FP, mean_state_action.RT,'.--', color = 'blue')
    ax[iSub].plot(empSubAction.foreperiod, empSubAction.RT, '.-', color = 'blue')
    ax[iSub].plot(mean_state_external.FP, mean_state_external.RT,'.--', color = 'orange')
    ax[iSub].plot(empSubExternal.foreperiod, empSubExternal.RT, '.-', color = 'orange')
    
fig.suptitle('ID fits')

figname='./Modeling/fMTP/Fitting/Plots/id_fit_SSE_k.png'
plt.savefig(figname,format='png')
plt.show()
    

#========================= Plot by sequential effects ===============================
empData = pd.read_csv('./Analysis/dataActionFPAll.csv')
empData.RT = empData.RT * 1000 # transform RT to ms
empData = empData.groupby(['ID', 'foreperiod', 'oneBackFP', 'condition'], as_index=False)[['RT']].mean()

fig = plt.figure(constrained_layout=True, figsize = (15.96,12))
subfigs=fig.subfigures(5,6)
subfigsCoords=list(itertools.product([0,1,2,3,4],[0,1,2,3,4,5]))

#subfigs = subfigs.ravel()
titles = stCoefs.ID.tolist()


for iSub in range(len(titles)):
    subData = stCoefs.iloc[iSub,:]
    empSub = empData[empData['ID']==titles[iSub]]
    empSubAction = empSub[empSub['condition']=='action']
    empSubExternal = empSub[empSub['condition']=='external']
    
    # Action condition
    kAction = subData['kAction']
    rAction = subData['rAction']
    cAction = subData['cAction']
    aAction = subData['a SSE Action']
    bAction = subData['b SSE Action']
    
    # Run model with current parameter values and simulate RTs for action condition
    fmtpAction=fMTP(rAction,cAction,kAction)
    state_discr_action, state_con_action = exp.run_exp(fmtpAction)
    state_discr_action=state_discr_action[1:]
    mean_state_action=state_discr_action.groupby(['FP', 'FPn_1']).mean().reset_index()
    mean_state_action['RT'] = temp2RT(mean_state_action.prep, aAction, bAction)
    
    #  External condition
    kExternal = subData['kExternal']
    rExternal = subData['rExternal']
    cExternal = subData['cExternal']
    aExternal = subData['a SSE External']
    bExternal = subData['b SSE External']
    
    # Run model with current parameter values and simulate RTs for external condition
    fmtpExternal=fMTP(rExternal,cExternal,kExternal)
    state_discr_external, state_con_external = exp.run_exp(fmtpExternal)
    state_discr_external=state_discr_external[1:]
    mean_state_external=state_discr_external.groupby(['FP', 'FPn_1']).mean().reset_index()
    mean_state_external['RT'] = temp2RT(mean_state_external.prep, aExternal, bExternal)
    
    # Plot on single panel
    ax=subfigs[subfigsCoords[iSub]].subplots(1,4,sharey=True)
    
    for idx, iFP in enumerate(np.unique(mean_state_action.FPn_1)):
        FP = round(iFP, 2)
        state_n1_action = mean_state_action[mean_state_action.FPn_1 == iFP]
        emp_n1_action = empSubAction[empSubAction.oneBackFP == FP]
        ax[idx].plot(state_n1_action.FP, state_n1_action.RT, '.--', color = 'blue')
        ax[idx].plot(emp_n1_action.foreperiod, emp_n1_action.RT, '.-', color = 'blue')
        
        state_n1_external = mean_state_external[mean_state_external.FPn_1 == iFP]
        emp_n1_external = empSubExternal[empSubExternal.oneBackFP == FP]
        ax[idx].plot(state_n1_external.FP, state_n1_external.RT, '.--', color = 'orange')
        ax[idx].plot(emp_n1_external.foreperiod, emp_n1_external.RT, '.-', color = 'orange')
        
        ax[idx].set_title(FP)
    
    subfigs[subfigsCoords[iSub]].suptitle(titles[iSub])

fig.suptitle('Sequential effects')

figname='./Modeling/fMTP/Fitting/Plots/seq_effects_id_fit_SSE_k.png'
plt.savefig(figname,format='png', bbox_inches = 'tight')

plt.show()


#=============== Compare parameters between conditions =========================

# Build dataset in long format for plotting
stCoefsAction = hazardResultsDF.loc[hazardResultsDF.groupby('ID').SSEAction.idxmin()].reset_index()
stCoefsAction = stCoefsAction[['ID', 'k', 'r', 'c', 'a SSE Action', 'b SSE Action', 'SSEAction']]
stCoefsAction['condition'] = pd.Series(['action'] * len(stCoefsAction))
stCoefsExternal = hazardResultsDF.loc[hazardResultsDF.groupby('ID').SSEExternal.idxmin()].reset_index()
stCoefsExternal = stCoefsExternal[['ID', 'k', 'r', 'c', 'a SSE External', 'b SSE External', 'SSEExternal']]
stCoefsExternal['condition'] = pd.Series(['external'] * len(stCoefsExternal))

stCoefsLong = pd.DataFrame(pd.concat([stCoefsAction.rename(columns = {'ID':'ID', 'condition':'condition', 'k':'k', 'r':'r', 'c':'c', 
                                                                   'a SSE Action':'aSSE', 'b SSE Action':'bSSE', 
                                                                   'SSEAction':'SSE'}), 
                                   stCoefsExternal.rename(columns = {'ID':'ID', 'condition':'condition', 'k':'k', 'r':'r', 'c':'c', 
                                                                     'a SSE External':'aSSE', 
                                                                     'b SSE External':'bSSE', 
                                                                     'SSEExternal': 'SSE'})], 
                                  axis = 0), columns = ['ID', 'condition', 'k', 'r', 'c', 'aSSE', 'bSSE', 'SSE'])
                                  

paramFig = plt.figure()
ax = paramFig.subplots(1, 2)

# Point plots
sns.stripplot(x = 'condition', y = 'k', data = stCoefsLong, jitter = 0.1, orient = 'v', ax = ax[0])
sns.pointplot(x = 'condition', y = 'k', data = stCoefsLong, ax = ax[1])

plt.show()

# Histograms
paramFig = plt.figure()
ax = paramFig.subplots(1, 2)

# Point plots
sns.histplot(x = "k", data = stCoefsLong[stCoefsLong.condition == "action"], ax = ax[0])
sns.histplot(x = "k", data = stCoefsLong[stCoefsLong.condition == "external"], ax = ax[1])


plt.show()

# Paired samples t-test
stats.wilcoxon(stCoefs.kAction, stCoefs.kExternal)

#===============================================================================#
#============================ 3. k and c fixed ==================================
#===============================================================================#
hazardResultsDF = pd.read_csv('./Modeling/fMTP/Fitting/fittingResults_r.csv')

# Find min for SSE in each condition
stCoefsAction = hazardResultsDF.loc[hazardResultsDF.groupby('ID').SSEAction.idxmin()].reset_index()
stCoefsAction = stCoefsAction[['ID', 'k', 'r', 'c', 'a SSE Action', 'b SSE Action', 'SSEAction']]
stCoefsExternal = hazardResultsDF.loc[hazardResultsDF.groupby('ID').SSEExternal.idxmin()].reset_index()
stCoefsExternal = stCoefsExternal[['k', 'r', 'c', 'a SSE External', 'b SSE External', 'SSEExternal']]

stCoefs = pd.DataFrame(pd.concat([stCoefsAction.rename(columns = {'ID':'ID', 'k':'kAction', 'r':'rAction', 'c':'cAction', 
                                                                  'a SSE Action':'a SSE Action', 'b SSE Action':'b SSE Action', 
                                                                  'SSEAction':'SSE_Action'}), 
                                  stCoefsExternal.rename(columns = {'k':'kExternal', 'r':'rExternal', 'c':'cExternal', 
                                                                    'a SSE External':'a SSE External', 
                                                                    'b SSE External':'b SSE External', 
                                                                    'SSEExternal': 'SSE_External'})], 
                                 axis = 1), columns = ['ID', 'kAction', 'rAction', 'cAction', 'a SSE Action', 'b SSE Action', 'SSE_Action',
                                                       'kExternal', 'rExternal', 'cExternal', 'a SSE External', 'b SSE External', 'SSE_External'])

  

# Build plots with individual data and predictions
FPs = np.arange(1.0, 3.0, 0.6)
distr = 'uni'
exp = FPexp(FPs = FPs, distribution = distr, tr_per_block = 150)


# Plot by foreperiod only
empData = pd.read_csv('./Analysis/dataActionFPAll.csv')
empData.RT = empData.RT * 1000 # transform RT to ms
empData = empData.groupby(['ID', 'foreperiod', 'condition'], as_index=False)[['RT']].mean()

fig = plt.figure(constrained_layout = True, figsize = (15.96, 12))
ax = fig.subplots(5,6,sharey=False)
ax = ax.ravel()

titles = stCoefs.ID.tolist()

for iSub, sub in enumerate(titles):
    subData = stCoefs.iloc[iSub,:]#[stCoefs['ID']==sub]#
    
    # Subset participant's data by condition
    empSub = empData[empData['ID']==sub]
    empSubAction = empSub[empSub['condition']=='action']
    empSubExternal = empSub[empSub['condition']=='external']
    
    # Action condition
    kAction = subData['kAction']
    rAction = subData['rAction']
    cAction = subData['cAction']
    aAction = subData['a SSE Action']
    bAction = subData['b SSE Action']
    
    # Run model with current parameter values and simulate RTs for action condition
    fmtpAction=fMTP(rAction,cAction,kAction)
    state_discr_action, state_con_action = exp.run_exp(fmtpAction)
    state_discr_action=state_discr_action[1:]
    mean_state_action=state_discr_action.groupby(['FP']).mean().reset_index()
    mean_state_action['RT'] = temp2RT(mean_state_action.prep, aAction, bAction)
    
    # External condition
    kExternal = subData['kExternal']
    rExternal = subData['rExternal']
    cExternal = subData['cExternal']
    aExternal = subData['a SSE External']
    bExternal = subData['b SSE External']
    
    # Run model with current parameter values and simulate RTs for external condition
    fmtpExternal=fMTP(rExternal,cExternal,kExternal)
    state_discr_external, state_con_external = exp.run_exp(fmtpExternal)
    state_discr_external=state_discr_external[1:]
    mean_state_external=state_discr_external.groupby(['FP']).mean().reset_index()
    mean_state_external['RT'] = temp2RT(mean_state_external.prep, aExternal, bExternal)
    
    # Plot on single panel
    ax[iSub].plot(mean_state_action.FP, mean_state_action.RT,'.--', color = 'blue')
    ax[iSub].plot(empSubAction.foreperiod, empSubAction.RT, '.-', color = 'blue')
    ax[iSub].plot(mean_state_external.FP, mean_state_external.RT,'.--', color = 'orange')
    ax[iSub].plot(empSubExternal.foreperiod, empSubExternal.RT, '.-', color = 'orange')
    
fig.suptitle('ID fits')

figname='./Modeling/fMTP/Fitting/Plots/id_fit_SSE_r.png'
plt.savefig(figname,format='png')
plt.show()
    

#========================= Plot by sequential effects ===============================
empData = pd.read_csv('./Analysis/dataActionFPAll.csv')
empData.RT = empData.RT * 1000 # transform RT to ms
empData = empData.groupby(['ID', 'foreperiod', 'oneBackFP', 'condition'], as_index=False)[['RT']].mean()

fig = plt.figure(constrained_layout=True, figsize = (15.96,12))
subfigs=fig.subfigures(5,6)
subfigsCoords=list(itertools.product([0,1,2,3,4],[0,1,2,3,4,5]))

#subfigs = subfigs.ravel()
titles = stCoefs.ID.tolist()


for iSub in range(len(titles)):
    subData = stCoefs.iloc[iSub,:]
    empSub = empData[empData['ID']==titles[iSub]]
    empSubAction = empSub[empSub['condition']=='action']
    empSubExternal = empSub[empSub['condition']=='external']
    
    # Action condition
    kAction = subData['kAction']
    rAction = subData['rAction']
    cAction = subData['cAction']
    aAction = subData['a SSE Action']
    bAction = subData['b SSE Action']
    
    # Run model with current parameter values and simulate RTs for action condition
    fmtpAction=fMTP(rAction,cAction,kAction)
    state_discr_action, state_con_action = exp.run_exp(fmtpAction)
    state_discr_action=state_discr_action[1:]
    mean_state_action=state_discr_action.groupby(['FP', 'FPn_1']).mean().reset_index()
    mean_state_action['RT'] = temp2RT(mean_state_action.prep, aAction, bAction)
    
    #  External condition
    kExternal = subData['kExternal']
    rExternal = subData['rExternal']
    cExternal = subData['cExternal']
    aExternal = subData['a SSE External']
    bExternal = subData['b SSE External']
    
    # Run model with current parameter values and simulate RTs for external condition
    fmtpExternal=fMTP(rExternal,cExternal,kExternal)
    state_discr_external, state_con_external = exp.run_exp(fmtpExternal)
    state_discr_external=state_discr_external[1:]
    mean_state_external=state_discr_external.groupby(['FP', 'FPn_1']).mean().reset_index()
    mean_state_external['RT'] = temp2RT(mean_state_external.prep, aExternal, bExternal)
    
    # Plot on single panel
    ax=subfigs[subfigsCoords[iSub]].subplots(1,4,sharey=True)
    
    for idx, iFP in enumerate(np.unique(mean_state_action.FPn_1)):
        FP = round(iFP, 2)
        state_n1_action = mean_state_action[mean_state_action.FPn_1 == iFP]
        emp_n1_action = empSubAction[empSubAction.oneBackFP == FP]
        ax[idx].plot(state_n1_action.FP, state_n1_action.RT, '.--', color = 'blue')
        ax[idx].plot(emp_n1_action.foreperiod, emp_n1_action.RT, '.-', color = 'blue')
        
        state_n1_external = mean_state_external[mean_state_external.FPn_1 == iFP]
        emp_n1_external = empSubExternal[empSubExternal.oneBackFP == FP]
        ax[idx].plot(state_n1_external.FP, state_n1_external.RT, '.--', color = 'orange')
        ax[idx].plot(emp_n1_external.foreperiod, emp_n1_external.RT, '.-', color = 'orange')
        
        ax[idx].set_title(FP)
    
    subfigs[subfigsCoords[iSub]].suptitle(titles[iSub])

fig.suptitle('Sequential effects')

figname='./Modeling/fMTP/Fitting/Plots/seq_effects_id_fit_SSE_r.png'
plt.savefig(figname,format='png', bbox_inches = 'tight')

plt.show()

#=============== Compare parameters between conditions =========================

# Build dataset in long format for plotting
stCoefsAction = hazardResultsDF.loc[hazardResultsDF.groupby('ID').SSEAction.idxmin()].reset_index()
stCoefsAction = stCoefsAction[['ID', 'k', 'r', 'c', 'a SSE Action', 'b SSE Action', 'SSEAction']]
stCoefsAction['condition'] = pd.Series(['action'] * len(stCoefsAction))
stCoefsExternal = hazardResultsDF.loc[hazardResultsDF.groupby('ID').SSEExternal.idxmin()].reset_index()
stCoefsExternal = stCoefsExternal[['ID', 'k', 'r', 'c', 'a SSE External', 'b SSE External', 'SSEExternal']]
stCoefsExternal['condition'] = pd.Series(['external'] * len(stCoefsExternal))

stCoefsLong = pd.DataFrame(pd.concat([stCoefsAction.rename(columns = {'ID':'ID', 'condition':'condition', 'k':'k', 'r':'r', 'c':'c', 
                                                                   'a SSE Action':'aSSE', 'b SSE Action':'bSSE', 
                                                                   'SSEAction':'SSE'}), 
                                   stCoefsExternal.rename(columns = {'ID':'ID', 'condition':'condition', 'k':'k', 'r':'r', 'c':'c', 
                                                                     'a SSE External':'aSSE', 
                                                                     'b SSE External':'bSSE', 
                                                                     'SSEExternal': 'SSE'})], 
                                  axis = 0), columns = ['ID', 'condition', 'k', 'r', 'c', 'aSSE', 'bSSE', 'SSE'])
                                  

paramFig = plt.figure()
ax = paramFig.subplots(1, 2)

# Point plots
sns.stripplot(x = 'condition', y = 'r', data = stCoefsLong, jitter = 0.1, orient = 'v', ax = ax[0])
sns.pointplot(x = 'condition', y = 'r', data = stCoefsLong, ax = ax[1])

plt.show()

# Histograms
paramFig = plt.figure()
ax = paramFig.subplots(1, 2)

# Point plots
sns.histplot(x = "r", data = stCoefsLong[stCoefsLong.condition == "action"], ax = ax[0])
sns.histplot(x = "r", data = stCoefsLong[stCoefsLong.condition == "external"], ax = ax[1])


plt.show()

# Paired samples t-test
stats.wilcoxon(stCoefs.rAction, stCoefs.rExternal)

#===============================================================================#
#============================ 4. k and r fixed ==================================
#===============================================================================#

hazardResultsDF = pd.read_csv('./Modeling/fMTP/Fitting/fittingResults_c.csv')

# Find min for SSE in each condition
stCoefsAction = hazardResultsDF.loc[hazardResultsDF.groupby('ID').SSEAction.idxmin()].reset_index()
stCoefsAction = stCoefsAction[['ID', 'k', 'r', 'c', 'a SSE Action', 'b SSE Action', 'SSEAction']]
stCoefsExternal = hazardResultsDF.loc[hazardResultsDF.groupby('ID').SSEExternal.idxmin()].reset_index()
stCoefsExternal = stCoefsExternal[['k', 'r', 'c', 'a SSE External', 'b SSE External', 'SSEExternal']]

stCoefs = pd.DataFrame(pd.concat([stCoefsAction.rename(columns = {'ID':'ID', 'k':'kAction', 'r':'rAction', 'c':'cAction', 
                                                                  'a SSE Action':'a SSE Action', 'b SSE Action':'b SSE Action', 
                                                                  'SSEAction':'SSE_Action'}), 
                                  stCoefsExternal.rename(columns = {'k':'kExternal', 'r':'rExternal', 'c':'cExternal', 
                                                                    'a SSE External':'a SSE External', 
                                                                    'b SSE External':'b SSE External', 
                                                                    'SSEExternal': 'SSE_External'})], 
                                 axis = 1), columns = ['ID', 'kAction', 'rAction', 'cAction', 'a SSE Action', 'b SSE Action', 'SSE_Action',
                                                       'kExternal', 'rExternal', 'cExternal', 'a SSE External', 'b SSE External', 'SSE_External'])

  

# Build plots with individual data and predictions
FPs = np.arange(1.0, 3.0, 0.6)
distr = 'uni'
exp = FPexp(FPs = FPs, distribution = distr, tr_per_block = 150)


# Plot by foreperiod only
empData = pd.read_csv('./Analysis/dataActionFPAll.csv')
empData.RT = empData.RT * 1000 # transform RT to ms
empData = empData.groupby(['ID', 'foreperiod', 'condition'], as_index=False)[['RT']].mean()

fig = plt.figure(constrained_layout = True, figsize = (15.96, 12))
ax = fig.subplots(5,6,sharey=False)
ax = ax.ravel()

titles = stCoefs.ID.tolist()

for iSub, sub in enumerate(titles):
    subData = stCoefs.iloc[iSub,:]#[stCoefs['ID']==sub]#
    
    # Subset participant's data by condition
    empSub = empData[empData['ID']==sub]
    empSubAction = empSub[empSub['condition']=='action']
    empSubExternal = empSub[empSub['condition']=='external']
    
    # Action condition
    kAction = subData['kAction']
    rAction = subData['rAction']
    cAction = subData['cAction']
    aAction = subData['a SSE Action']
    bAction = subData['b SSE Action']
    
    # Run model with current parameter values and simulate RTs for action condition
    fmtpAction=fMTP(rAction,cAction,kAction)
    state_discr_action, state_con_action = exp.run_exp(fmtpAction)
    state_discr_action=state_discr_action[1:]
    mean_state_action=state_discr_action.groupby(['FP']).mean().reset_index()
    mean_state_action['RT'] = temp2RT(mean_state_action.prep, aAction, bAction)
    
    # External condition
    kExternal = subData['kExternal']
    rExternal = subData['rExternal']
    cExternal = subData['cExternal']
    aExternal = subData['a SSE External']
    bExternal = subData['b SSE External']
    
    # Run model with current parameter values and simulate RTs for external condition
    fmtpExternal=fMTP(rExternal,cExternal,kExternal)
    state_discr_external, state_con_external = exp.run_exp(fmtpExternal)
    state_discr_external=state_discr_external[1:]
    mean_state_external=state_discr_external.groupby(['FP']).mean().reset_index()
    mean_state_external['RT'] = temp2RT(mean_state_external.prep, aExternal, bExternal)
    
    # Plot on single panel
    ax[iSub].plot(mean_state_action.FP, mean_state_action.RT,'.--', color = 'blue')
    ax[iSub].plot(empSubAction.foreperiod, empSubAction.RT, '.-', color = 'blue')
    ax[iSub].plot(mean_state_external.FP, mean_state_external.RT,'.--', color = 'orange')
    ax[iSub].plot(empSubExternal.foreperiod, empSubExternal.RT, '.-', color = 'orange')
    
fig.suptitle('ID fits')

figname='./Modeling/fMTP/Fitting/Plots/id_fit_SSE_c.png'
plt.savefig(figname,format='png')
plt.show()
    

#========================= Plot by sequential effects ===============================
empData = pd.read_csv('./Analysis/dataActionFPAll.csv')
empData.RT = empData.RT * 1000 # transform RT to ms
empData = empData.groupby(['ID', 'foreperiod', 'oneBackFP', 'condition'], as_index=False)[['RT']].mean()

fig = plt.figure(constrained_layout=True, figsize = (15.96,12))
subfigs=fig.subfigures(5,6)
subfigsCoords=list(itertools.product([0,1,2,3,4],[0,1,2,3,4,5]))

#subfigs = subfigs.ravel()
titles = stCoefs.ID.tolist()


for iSub in range(len(titles)):
    subData = stCoefs.iloc[iSub,:]
    empSub = empData[empData['ID']==titles[iSub]]
    empSubAction = empSub[empSub['condition']=='action']
    empSubExternal = empSub[empSub['condition']=='external']
    
    # Action condition
    kAction = subData['kAction']
    rAction = subData['rAction']
    cAction = subData['cAction']
    aAction = subData['a SSE Action']
    bAction = subData['b SSE Action']
    
    # Run model with current parameter values and simulate RTs for action condition
    fmtpAction=fMTP(rAction,cAction,kAction)
    state_discr_action, state_con_action = exp.run_exp(fmtpAction)
    state_discr_action=state_discr_action[1:]
    mean_state_action=state_discr_action.groupby(['FP', 'FPn_1']).mean().reset_index()
    mean_state_action['RT'] = temp2RT(mean_state_action.prep, aAction, bAction)
    
    #  External condition
    kExternal = subData['kExternal']
    rExternal = subData['rExternal']
    cExternal = subData['cExternal']
    aExternal = subData['a SSE External']
    bExternal = subData['b SSE External']
    
    # Run model with current parameter values and simulate RTs for external condition
    fmtpExternal=fMTP(rExternal,cExternal,kExternal)
    state_discr_external, state_con_external = exp.run_exp(fmtpExternal)
    state_discr_external=state_discr_external[1:]
    mean_state_external=state_discr_external.groupby(['FP', 'FPn_1']).mean().reset_index()
    mean_state_external['RT'] = temp2RT(mean_state_external.prep, aExternal, bExternal)
    
    # Plot on single panel
    ax=subfigs[subfigsCoords[iSub]].subplots(1,4,sharey=True)
    
    for idx, iFP in enumerate(np.unique(mean_state_action.FPn_1)):
        FP = round(iFP, 2)
        state_n1_action = mean_state_action[mean_state_action.FPn_1 == iFP]
        emp_n1_action = empSubAction[empSubAction.oneBackFP == FP]
        ax[idx].plot(state_n1_action.FP, state_n1_action.RT, '.--', color = 'blue')
        ax[idx].plot(emp_n1_action.foreperiod, emp_n1_action.RT, '.-', color = 'blue')
        
        state_n1_external = mean_state_external[mean_state_external.FPn_1 == iFP]
        emp_n1_external = empSubExternal[empSubExternal.oneBackFP == FP]
        ax[idx].plot(state_n1_external.FP, state_n1_external.RT, '.--', color = 'orange')
        ax[idx].plot(emp_n1_external.foreperiod, emp_n1_external.RT, '.-', color = 'orange')
        
        ax[idx].set_title(FP)
    
    subfigs[subfigsCoords[iSub]].suptitle(titles[iSub])

fig.suptitle('Sequential effects')

figname='./Modeling/fMTP/Fitting/Plots/seq_effects_id_fit_SSE_c.png'
plt.savefig(figname,format='png', bbox_inches = 'tight')

plt.show()


#=============== Compare parameters between conditions =========================

# Build dataset in long format for plotting
stCoefsAction = hazardResultsDF.loc[hazardResultsDF.groupby('ID').SSEAction.idxmin()].reset_index()
stCoefsAction = stCoefsAction[['ID', 'k', 'r', 'c', 'a SSE Action', 'b SSE Action', 'SSEAction']]
stCoefsAction['condition'] = pd.Series(['action'] * len(stCoefsAction))
stCoefsExternal = hazardResultsDF.loc[hazardResultsDF.groupby('ID').SSEExternal.idxmin()].reset_index()
stCoefsExternal = stCoefsExternal[['ID', 'k', 'r', 'c', 'a SSE External', 'b SSE External', 'SSEExternal']]
stCoefsExternal['condition'] = pd.Series(['external'] * len(stCoefsExternal))

stCoefsLong = pd.DataFrame(pd.concat([stCoefsAction.rename(columns = {'ID':'ID', 'condition':'condition', 'k':'k', 'r':'r', 'c':'c', 
                                                                   'a SSE Action':'aSSE', 'b SSE Action':'bSSE', 
                                                                   'SSEAction':'SSE'}), 
                                   stCoefsExternal.rename(columns = {'ID':'ID', 'condition':'condition', 'k':'k', 'r':'r', 'c':'c', 
                                                                     'a SSE External':'aSSE', 
                                                                     'b SSE External':'bSSE', 
                                                                     'SSEExternal': 'SSE'})], 
                                  axis = 0), columns = ['ID', 'condition', 'k', 'r', 'c', 'aSSE', 'bSSE', 'SSE'])
                                  

paramFig = plt.figure()
ax = paramFig.subplots(1, 2)

# Point plots
sns.stripplot(x = 'condition', y = 'c', data = stCoefsLong, jitter = 0.1, orient = 'v', ax = ax[0])
sns.pointplot(x = 'condition', y = 'c', data = stCoefsLong, ax = ax[1])

plt.show()

# Histograms
paramFig = plt.figure()
ax = paramFig.subplots(1, 2)

# Point plots
sns.histplot(x = "c", data = stCoefsLong[stCoefsLong.condition == "action"], ax = ax[0])
sns.histplot(x = "c", data = stCoefsLong[stCoefsLong.condition == "external"], ax = ax[1])


plt.show()

# Paired samples t-test
stats.wilcoxon(stCoefs.cAction, stCoefs.cExternal)
