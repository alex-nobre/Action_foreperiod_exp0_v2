#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 21:03:56 2023

@author: cravo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from os import listdir
import glob
from scipy import stats
from statsmodels.formula.api import ols



import sys, os

folder_with_csv='G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_0_v2/Analysis'
#folder_with_csv='/media/cravo/5a7377a6-fa3b-4305-be1d-c3950af2cace/Dropbox/Dropbox/UFABC/Projetos_Alunos/AlexandreNobre/RStuff'
os.chdir(folder_with_csv)


data = pd.read_csv('data2.csv')
data['FPDiff']=data['logFP']-data['logOneBackFP']

unique_subs=data['ID'].unique()

params_action=np.zeros([len(unique_subs),4])
params_external=np.zeros([len(unique_subs),4])

for iSub,Sub in enumerate(unique_subs):
    
    data_sub=data.loc[data['ID']==Sub]
    
    model_action = ols("logRT ~ logFP*FPDiff", data_sub.loc[data_sub['condition']=='action']).fit()
    model_external = ols("logRT ~ logFP*FPDiff", data_sub.loc[data_sub['condition']=='external']).fit()
    
    params_action[iSub,:]=model_action.params.values
    params_external[iSub,:]=model_external.params.values



np.mean(params_action,axis=0)
np.mean(params_external,axis=0)
   
stats.ttest_1samp(params_action[:,0],0)
stats.ttest_1samp(params_action[:,1],0)
stats.ttest_1samp(params_action[:,2],0)
stats.ttest_1samp(params_action[:,3],0)

stats.ttest_1samp(params_external[:,0],0)
stats.ttest_1samp(params_external[:,1],0)
stats.ttest_1samp(params_external[:,2],0)
stats.ttest_1samp(params_external[:,3],0)
