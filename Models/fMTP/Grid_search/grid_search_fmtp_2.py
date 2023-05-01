# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 14:10:32 2022
Here, we perform a grid search for combinations of values for model parameters of fMTP that change the slope of the hazard function without changing sequential effects. We do this by quantifying the hazard and sequential effects and comparing them to see how they correlate.

For now, we vary only the parameters investigated in the paper (k, c and r).
10 values by parameters, linearly spaced between the smallest and largest value shown in the paper (except for k, which we vary from 1 to 10 to avoid floating values for k, which the model does not accept).
For each combination of values, simulate curves for hazard and sequential effects.
Quantify hazard effect by the slope of the FP X RT curve
Quantify the sequential effect by subtracting the slope for the smallest FP n-1 value from the largest FP n-1 value
Build scatterplots of hazard and sequential effect indices;
Check for:
a line (hazard and sequentia effect magnitudes correlate for the values of all parameters)
completely scattered points (each combination generates an effect that is unrelated to the other - this is unlikely)
clusters (correlation for some values, absence of correlation for others)
Hypothesis: as this model is in part an elaboration of trace conditions, it is unlikely that it may generate hazard effects dissociated from sequential effects: in trace conditioning, the hazard effect is a product of sequential effects; in fMTP, both are the products of a common memory mechanism, and thus should correlate.
If that is the case, a different model, such as dual-process models - e.g., as implemented by Grabenhorst et al., 2019/2021 - or a bayesian model might explain the results better.
For now, use our current data, which will (or will not) be confirmed by a second experiment.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
import seaborn as sns

import os
import sys

sys.path.insert(0, './Models/fMTP/Grid_search')

from fmtp import fMTP, FPexp

#================================== Constants and parameters ================================#
# Initialize values for simulation
FP = np.arange(1.0, 2.8, 0.6)
distr = 'uni'

# Set-up experiment using "FPexp"
exp = FPexp(FPs = FP, distribution = distr, tr_per_block = 150)
                                              
# create list of parameter values
kList=[i for i in range(1,11)]  #np.linspace(2,8,num=10)
#rList=np.linspace(1,4,num=10)
rList=np.linspace(-4,-1,num=10)
cList=np.linspace(0,0.0003,num=10)

# Simulate experiments and store in lists
simResults=[None]*(len(kList)*len(rList)*len(cList))   
    
simEl=0   
for ik in range(len(kList)):
    for ir in range(len(rList)):
        for ic in range(len(cList)):
            k=kList[ik]
            r=rList[ir]
            c=cList[ic]
            fmtp=fMTP(r,c,k)
            state_discr, state_con = exp.run_exp(fmtp)
            state_discr=state_discr[1:]
            
            # Hazard slope
            hazardLm=smf.ols(formula='prep~FP', data=state_discr).fit()
            hazardSlope=round(hazardLm.params[1],3)
            
            # Sequential effects difference of slopes
            lowFPn_1=min(np.unique(state_discr.FPn_1))
            highFPn_1=max(np.unique(state_discr.FPn_1))
            # low FP n-1
            state_lowFPn_1=state_discr[state_discr.FPn_1==lowFPn_1]
            lowSeqLm=smf.ols(formula='prep~FP',data=state_lowFPn_1).fit()
            lowSeqSlope=lowSeqLm.params[1]
            # High FP n-1
            state_highFPn_1=state_discr[state_discr.FPn_1==highFPn_1]
            highSeqLm=smf.ols(formula='prep~FP',data=state_highFPn_1).fit()
            highSeqSlope=highSeqLm.params[1]
            
            seqSlopeDiff=round(highSeqSlope-lowSeqSlope,3)
            simResults[simEl]=(hazardSlope,seqSlopeDiff)
            # simResults[simEl]=(hazardSlope,
            #                    round(highSeqSlope,3),
            #                    round(lowSeqSlope,3),
            #                    seqSlopeDiff)
            print(simEl)
            simEl+=1

simResults2=simResults
            

#=================== plot indices on a scatterplot ======================#
plt.scatter(*zip(*simResults2))
plt.show()

# Histograms of sequential effects
sepSimResults=list(set(zip(*simResults)))

plt.figure()
sns.histplot(sepSimResults[1],
             binwidth=0.1)
plt.show()

positiveSeqEff=[ef for ef in sepSimResults[0] if ef >0] # None are positive

# Plot on same panel, but with colors and annotations to indicate value of k        
chunkedList=list()
chunkSize=int(len(simResults)/len(kList))

for i in range(0,len(simResults),chunkSize):
    chunkedList.append(simResults[i:i+chunkSize])

fig=plt.figure()

kSublist=0
for ik in range(len(chunkedList)):
    k=kList[kSublist]
    plotSublist=chunkedList[kSublist]
    plt.scatter(*zip(*plotSublist))
    plt.annotate(('k='+str(k)),xy=(plotSublist[-1][0],plotSublist[-1][1]+0.002))
    kSublist+=1
    
plt.show()

figname='./Models/fMTP/Grid_search/Plots/hazardSeqSlopeScatterplot.png'
plt.savefig(figname,format='png')


# Plot separate facets for each k value (each facet including all combinations of r and c)
fig,ax=plt.subplots(nrows=2,ncols=5)

kSublist=0
for row in ax:
    for col in row:
        k=kList[kSublist]
        plotSublist=chunkedList[kSublist]
        col.scatter(*zip(*plotSublist))
        col.set_title('k='+str(k))
        kSublist+=1

plt.show()
plt.savefig('./Models/fMTP/Grid_search/Plots/hazardSeqSlopeScatterplot_facet.png',
            format='png')


# Separate plots for each r value (within a given value of k)
fig,ax=plt.subplots(nrows=2,ncols=5)
        
chunkedChunkedList=list()
chunkSize=int(len(cList))
kChunk=chunkedList[0]

for i in range(0,len(kChunk),chunkSize):
    chunkedChunkedList.append(kChunk[i:i+chunkSize])

rSublist=0
for row in ax:
    for col in row:
        plotSublist=chunkedChunkedList[rSublist]
        col.scatter(*zip(*plotSublist))
        rSublist+=1

plt.show()
plt.savefig('./Models/fMTP/Grid_search/Plots/rxc_k1_facet.png',
            format='png',bbox_inches='tight')
        
        
# Compute correlation
# hazardSlopeList=[value[0] for value in simResults]
# seqSlopeDiffList=[value[1] for value in simResults]

# np.corrcoef(hazardSlopeList, seqSlopeDiffList)
np.corrcoef(*zip(*simResults2))

corrList=list()
for ik in range(len(chunkedList)):
    corrSublist=chunkedList[ik]
    corrList.append(np.corrcoef(*zip(*corrSublist)))
    
# Without c

# List to store values
simResultsNoC=[None]*(len(kList)*len(rList))

simEl=0   
for ik in range(len(kList)):
    for ir in range(len(rList)):
        k=kList[ik]
        r=rList[ir]
        fmtp=fMTP(r,c,k)
        state_discr, state_con = exp.run_exp(fmtp)
        state_discr=state_discr[1:]
        
        # Hazard slope
        hazardLm=smf.ols(formula='prep~FP', data=state_discr).fit()
        hazardSlope=round(hazardLm.params[1],3)
        
        # Sequential effects difference of slopes
        lowFPn_1=min(np.unique(state_discr.FPn_1))
        highFPn_1=max(np.unique(state_discr.FPn_1))
        # low FP n-1
        state_lowFPn_1=state_discr[state_discr.FPn_1==lowFPn_1]
        lowSeqLm=smf.ols(formula='prep~FP',data=state_lowFPn_1).fit()
        lowSeqSlope=lowSeqLm.params[1]
        # High FP n-1
        state_highFPn_1=state_discr[state_discr.FPn_1==highFPn_1]
        highSeqLm=smf.ols(formula='prep~FP',data=state_highFPn_1).fit()
        highSeqSlope=highSeqLm.params[1]
        
        seqSlopeDiff=round(highSeqSlope-lowSeqSlope,3)
        simResultsNoC[simEl]=(hazardSlope,seqSlopeDiff)
        
        print(simEl)
        simEl+=1

sepSimResultsNoC=list(set(zip(*simResultsNoC)))
hazardList=list(sepSimResultsNoC[0])
reshapedHazardList=np.reshape(hazardList, (10,10))

hazardData=pd.DataFrame(reshapedHazardList.transpose(),
                        columns=kList,
                        index=[round(item,3) for item in rList])

hazardData.to_csv('./Models/fMTP/Grid_search/hazardData.csv')

seqEffList=list(sepSimResultsNoC[1])
reshapedSeqEffList=np.reshape(seqEffList, (10,10))
seqEffData=pd.DataFrame(reshapedSeqEffList.transpose(),
                        columns=kList,
                        index=rList)

seqEffData.to_csv('./Models/fMTP/Grid_search/SeqEffData.csv')
