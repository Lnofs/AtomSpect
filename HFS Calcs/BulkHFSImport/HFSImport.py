# -*- coding: utf-8 -*-
"""
Created on 06302025
This script takes an input file and does automated calculation and comparison for hyperfine constant values.

@author: %Leo Nofs
"""


from HFS_Const_Calcs import Ahfs, pererror, read_HFSInput, powerset

import csv, os, sys
import numpy as np
from matplotlib import pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir, os.path.pardir)))

pathfil = os.getcwd()
# filname = '\HFSInputs_Coupled.csv'
filname = '\HFSInputs.csv'

fil = pathfil + filname
#This function creates a dictionary of the input parameters. Which should be exactly what the calculation function needs!
HFSInputs = read_HFSInput(fil)

#%%CalculateInputs
#Make it so that ds is included for s states,
# Couple = ['ss']
Corrections = 'FHdelepsds'

#dict_keys(['Atom', 'Label', 'AMU', 'n1', 'E1', 'Eionize', 's1', 'l1', 'j1', 'Z', 'Za', 'I', 'mu_I', 'Q(b)', 'A_Exp', 'B_Exp', 'CorrFlags', 'Coupling', 'n2', 'E2', 's2', 'l2', 'j2', 'S', 'L', 'J', 'A_Calc', 'B_Calc', '', 'Source', 'Notes'])
for HFScalcs in HFSInputs:
        # HFScalcs['Coupling'] = HFSInputs['Coupling']
        if 'Corrflags' not in HFScalcs:
            HFScalcs['Corrflags'] = Corrections
        else:
            HFScalcs['Corrflags'] = HFScalcs['Corrflags'] + 'Corrections'

        # HFScalcs['Corrflags'] = HFSInputs['Corrflags'] #Set the flag for the correction factor for fast comparison of multiple settings.

        # HFScalcs['Corrflags'] = Corrections #Set the flag for the correction factor for fast comparison of multiple settings.
        # HFScalcs['Coupling'] = Couple
        HFScalcs['LongLable'] = f'{HFScalcs['Label']}, {HFScalcs['Corrflags']}'#, Coupling={cc}'
        HFScalcs['A_Calc'] = Ahfs(HFScalcs)['A']
        HFScalcs['pererror'] = pererror(HFScalcs['A_Exp'],HFScalcs['A_Calc'])



verboseresults = True
if verboseresults:
    for HFScalcs in HFSInputs:
        if HFScalcs['pererror']>20: #Print out all the cases of error being more than 30%.

            print(f'Label: {HFScalcs['Label']}') 
            print(f'A_Exp: {HFScalcs['A_Exp']}') 
            print(f'A_Calc: {HFScalcs['A_Calc']}') 
            print(f'%Error: {HFScalcs['pererror']}') 
            print(f'%Corrs: {HFScalcs['Corrflags']}') 
    
            print('')
              
              
              
#%%PlotInputDirectly
              
atomlist = set([dictvals['Atom'] for dictvals in HFSInputs])
plt.close('all')

# for atoms in atomlist:

plt.figure()
plt.title('Calculated Hyperfine A Constants')
newcolors = plt.cm.rainbow(np.linspace(0, 1, len(HFSInputs)))
plt.scatter(HFSInputs[0]['Label'], HFSInputs[0]['A_Exp'], color='black',alpha=0.75,marker='x', label='Exp')

for t, things in enumerate(HFSInputs):
    # if things['Atom']==atoms:
        plt.scatter(things['Label'],things['A_Calc'], s=100,color=newcolors[t],  label=things['LongLable'])
        plt.scatter(things['Label'], things['A_Exp'], color='black',alpha=0.75,marker='x')#, label='Exp')
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left",borderaxespad=0)
        plt.xticks(rotation=90)


plt.figure()
plt.title('Percent Error for Input File')
newcolors = plt.cm.rainbow(np.linspace(0, 1, len(HFSInputs)))

for t, things in enumerate(HFSInputs):
    # if things['Atom']==atoms:

        plt.scatter(things['Label'],things['pererror'], s=100,color=newcolors[t],  label=things['LongLable'])
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", ncols=2,borderaxespad=0)
        plt.xticks(rotation=90)



#%%InputSpan
#This section spans some input configurations to look for trends in calculation.
#The goal was to try to determine where calculations might be going wrong, as well as what plays the most impact.
#
Corrflaglist = ['F','H','R','del','eps','ds']

# Couple = ['ss','sl1','sl2','jj','None']
Couple = ['None']
spancorrflags = powerset(Corrflaglist)
resultslist = []

doinputspan=True
if doinputspan==True:
    
    minerror = []
    #dict_keys(['Atom', 'Label', 'AMU', 'n1', 'E1', 'Eionize', 's1', 'l1', 'j1', 'Z', 'Za', 'I', 'mu_I', 'Q(b)', 'A_Exp', 'B_Exp', 'CorrFlags', 'Coupling', 'n2', 'E2', 's2', 'l2', 'j2', 'S', 'L', 'J', 'A_Calc', 'B_Calc', '', 'Source', 'Notes'])
    for HFScalcs in HFSInputs:
        tempminerror = {'Label' : HFScalcs['Label']      ,
                        'A_Exp' : HFScalcs['A_Exp']
                        }
        tmin = 100000
    
        for cc in Couple:
            HFScalcs['Coupling'] = cc
            for corrs in spancorrflags:
                results = {}
                HFScalcs['Corrflags'] = corrs #Set the flag for the correction factor for fast comparison of multiple settings.
                results['A_Exp'] = HFScalcs['A_Exp']
                # results['B_Exp'] = HFScalcs['B_Exp']
                results['A_Calc'] = Ahfs(HFScalcs)['A']
                results['Label'] = HFScalcs['Label']
                results['Atom'] = HFScalcs['Atom']
                results['Corrflags'] = HFScalcs['Corrflags']
                results['Coupling'] = cc
    
                ttmin = pererror(HFScalcs['A_Exp'],results['A_Calc'])
                results['pererror'] = ttmin
                if ttmin<tmin:
                    tmin = np.min([tmin,ttmin])
                    tempminerror['A_Calc'] = Ahfs(HFScalcs)['A']
                    tempminerror['minerror'] = tmin
                    tempminerror['Corrflags'] = corrs
                    tempminerror['Coupling'] = cc
                    tempminerror['LongLable'] = f'{HFScalcs['Label']}, {corrs}'#, Coupling={cc}'
                    
                resultslist.append(results)
        minerror.append(tempminerror)
    
    
    
    
    
    # plt.close('all')
    
    plt.figure()
    plt.title('Corrections corresponding to minimum error')
    newcolors = plt.cm.rainbow(np.linspace(0, 1, len(minerror)))

    for tt, thingsa in enumerate(minerror):
        plt.scatter(thingsa['Label'],thingsa['A_Calc'], s=100,color=newcolors[tt], marker = 'o', label=thingsa['LongLable'])

        plt.scatter(thingsa['Label'],thingsa['A_Exp'], color='black',alpha=0.75,marker='x')#, label='Exp')
        
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left",borderaxespad=0)
        plt.xticks(rotation=90)
        
    plt.figure()
    plt.title('Min Percent Error')
    newcolors = plt.cm.rainbow(np.linspace(0, 1, len(minerror)))
    
    for t, things in enumerate(minerror):
        plt.scatter(things['Label'],things['minerror'], s=100,color=newcolors[t],  label=things['LongLable'])
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", ncols=2,borderaxespad=0)
        plt.xticks(rotation=90)
        
        


#%Print out all the
#%%
# verboseresults = True
if verboseresults:
    for HFScalcs in minerror:
        if HFScalcs['minerror']>30: #Print out all the cases of error being more than 30%.

            print(f'Label: {HFScalcs['Label']}') 
            print(f'A_Exp: {HFScalcs['A_Exp']}') 
            print(f'A_Calc: {HFScalcs['A_Calc']}') 
            print(f'%Error: {HFScalcs['minerror']}') 
            print(f'%Corrs: {HFScalcs['Corrflags']}') 
    
            print('')

# # plt.twinx()
# plt.figure()
# plt.plot(Transposemin[0][np.argsort(Transposemin[1])] ,np.argsort(Transposemin[1]))
# plt.tight_layout()

#%%PlotAs
PlotAllInOne = False
if PlotAllInOne:
    #Let the x axis be the different states/labels. Each label should be unique. Then will use different markers/colors for the different correction factors.
    plt.close('all')
    plt.figure()
    newcolors = plt.cm.rainbow(np.linspace(0, 1, len(spancorrflags)+1))
    
    for i,res in enumerate(resultslist):
        if i<len(spancorrflags):
            # print(i)
            plt.scatter(res['Label'],res['A_Calc'], s=3,color=newcolors[i%len(spancorrflags)],  label=res['Corrflags'])
        else:
            plt.scatter(res['Label'],res['A_Calc'], s=3,color=newcolors[i%len(spancorrflags)])
    
    ii=0
    for HFSins in HFSInputs:
        if ii==0:
            plt.scatter(HFSins['Label'], HFSins['A_Exp'], color='black',alpha=0.75,marker='x', label='Exp')
            ii+=1
        else:
            plt.scatter(HFSins['Label'], HFSins['A_Exp'], color='black',alpha=0.75,marker='x')
    
    plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", ncols=2,borderaxespad=0)
    plt.tight_layout()

#%%PerLevel (LSJ)
#Each column will now have a zoomed in with all of the individual correction factors 
PlotSingleLevel = False
if PlotSingleLevel:
        
    plt.close('all')
    
    Ncorr = len(spancorrflags)
    newcolors = plt.cm.rainbow(np.linspace(0, 1, len(spancorrflags)+1))
    
    # for i,HFScalcs in enumerate(HFSInputs[:2]):
    # plt.figure()
    # plt.title(HFScalcs['Label'])
    
    for jj,res in enumerate(resultslist):
        if jj%Ncorr==0:
            plt.figure()
            plt.title(f'{res['Label']}, Coupling={res['Coupling']} , A={res['A_Exp']} MHz')
            plt.hlines(res['A_Exp'], spancorrflags[0],Ncorr , color = 'black',label='Exp')
            plt.scatter(res['Corrflags'],res['A_Calc'], color = newcolors[jj%Ncorr])
    
        else:
            plt.scatter(res['Corrflags'],res['A_Calc'], color = newcolors[jj%Ncorr])
        plt.xticks(rotation=90)
        plt.grid(True,axis='x')
        # print(res['Corrflags'],res['A_Calc'])

    # for jj,res in enumerate(resultslist):
    #     print(res['Corrflags'],res['A_Calc'])

# for j,corrs in enumerate(spancorrflags):
#     if j==0:
#         plt.title(resultslist[j]['Label'])

#     plt.scatter(corrs,resultslist[j]['A_Calc'], color = newcolors[j])
#     # plt.plot([spancorrflags[0],spancorrflags[-1]], [HFScalcs['A_Exp'],HFScalcs['A_Exp']] , marker='x', color = 'black',label='Exp')



