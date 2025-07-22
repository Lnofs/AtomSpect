# -*- coding: utf-8 -*-
"""
Created on Sat Jul 27 10:06:38 2024

@author: Leo Nofs
This code does a comparison of what thee expected results will be given some example doppler broadening and instrument functions.

As of updating packages on Mar 31, 2025, this file no longer functions as expected. :(
"""


import  os, sys
import numpy as np
from matplotlib import pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir, os.path.pardir)))
from SpectSplit  import Zeeman_Main, Normalize,PlotFunction,read_Spectra, plotZfan, Convol_Spect,PlotFunction
pathfil = os.getcwd()
filname = '/Examples/SpectPredict/'


dpival = 800





spec405 = Normalize(np.load('spec_405.npy'))
wav405 = np.load('wave_405.npy')
spec430 = Normalize(np.load('spec_430.npy'))
wav430 = np.load('wave_430.npy')

specdat1 = [wav405,spec405]
specdat2 = [wav430,spec430]
#%%
        
# plt.close('all')        
    
Bangle = 90
Temps = [0.025,1,3]
specsteps = [0.05,0.01,0.005,0.002]
# functions = ['Lorrentzian','Gaussian','GaussianInst']
B_mag = [0.6,1.,3]
specres= 1

savefig = False

Temp_ev = 1
functions=['Gaussian']
specstep = 0.002222 
for function in functions:
    for i,Bvals in enumerate(B_mag):
        # for k,Temp_ev in enumerate(Temps):
            InputdeckW = {'s_ground':3 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                          's_excited':3 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                          'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
                          'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                          'E_ground': np.array([2951.29 ]) ,#Lowest J first, in cm^-1
                          'E_excited': np.array([26229.77,27488.11,27889.68]) , #Lowest J first, in cm^-1
                          'Bmag': Bvals , #Magnetic Field  [T]
                          'b_angle': Bangle , #Angle between LoS and Bmax
                          'plottitle':r'W I []5d^5(^6S)6s (^7S_3) -> []5d^5(^6S)6p (^7P)' f'B={Bvals}, LoS = {Bangle}, Resolution = {specstep}' , #Title for plotting (optional)     
                          'plot_window' : [395,435], #Min and max for plot window range (nm)
                          'amu' : 183 ,  #Weight in AMU
                          'specstep' : specstep , #Resolution of spectrometer in nm                      
                          'fxnwindow': 0.05, #How far the convolution evaluates from the last stick in nm
                          'Convfxn' : function , #Optional: 'Gaussian', 'Skewed'.
                          'Temp': Temp_ev*11602,  #Temperature in K. Used for Gaussian convolution
                          'specres' : specres ,                      
                           }             
            SpecPred = Zeeman_Main(InputdeckW)
            plotvars=['green', 2.5, '-.',InputdeckW['plottitle'] ,'PolLabel']
            for k in range(3):  #This is the number of bin windows            
                if k==0 and i==0:                
                    multax = PlotFunction(SpecPred,plotvars,plotwind=SpecPred['bin_windows'][k],
                                                   SpectrumPlot =specdat1 ,
                                                   Shape=(len(B_mag),len(SpecPred['bin_windows'])) ,                                               
                                                   position=[i,k],NormalizeSig=True,makefig=True)
                if k==2: #430 is a different spectrometer measurement, so this plots that in the third column
                    PlotFunction(SpecPred,plotvars,Shape=(len(B_mag),len(SpecPred['bin_windows'])) ,
                              axsin = multax,plotwind=SpecPred['bin_windows'][k],
                              SpectrumPlot =specdat2 ,
                              position=[i,k],NormalizeSig=True,makefig=False)            
                else:    
                    PlotFunction(SpecPred,plotvars,Shape=(len(B_mag),len(SpecPred['bin_windows'])) ,
                              axsin = multax,plotwind=SpecPred['bin_windows'][k],
                              SpectrumPlot =specdat1 ,
                              position=[i,k],NormalizeSig=True,makefig=False)
    plt.suptitle(f'Testing W lines for required resolution with {function} convolution ')



