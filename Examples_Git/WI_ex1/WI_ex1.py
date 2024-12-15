# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 17:08:11 2024

@author: Leo Nofs
"""

import csv, os, sys
import numpy as np
from matplotlib import pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir, os.path.pardir)))
from SpectZHFS  import Zeeman_Main, Normalize,PlotFunction,read_Spectra, plotZfan



doHFS=True

pathfil = os.getcwd()
filname = '\Wspec_ex1.csv'
fil = pathfil + filname

WSpect = read_Spectra(fil, headercount=1)





InputdeckW = {'s_ground':3 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
              's_excited':3 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
              'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
              'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
              'E_ground': np.array([2951.29 ]) ,#Lowest J first, in cm^-1
              'E_excited': np.array([26229.77,27488.11,27889.68]) , #Lowest J first, in cm^-1
              'Bmag': 0.6 , #Magnetic Field  [T]
              'b_angle': 65 , #Angle between LoS and Bmax
              'plottitle':'W I []6s (^7S_3) -> []6p (^7P_2) B=0.6T' , #Title for plotting (optional)     
              'plot_window' : [425,435], #Min and max for plot window range (nm)
              'amu' : 183 ,  #Weight in AMU
              'specstep' : 0.002222 , #Resolution of spectrometer in nm
               }




InputdeckWHFS = {'s_ground':3 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
              's_excited':3 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
              'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
              'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
              'E_ground': np.array([2951.29 ]) ,#Lowest J first, in cm^-1
              'E_excited': np.array([26229.77,27488.11,27889.68]) , #Lowest J first, in cm^-1
              'Bmag': 0.6 , #Magnetic Field  [T]
              'b_angle': 65 , #Angle between LoS and Bmax
              'plottitle':'W I []6s (^7S_3) -> []6p (^7P_2) B=0.6T' , #Title for plotting (optional)     
              'plot_window' : [425,435], #Min and max for plot window range (nm)
              'amu' : 183 ,  #Weight in AMU
              'specstep' : 0.002222 , #Resolution of spectrometer in nm
              'HFS_G' : [505.5e6] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
              'HFS_E': [0, 496.2e6, 440.5e6] , #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
              'I_spin' : -0.5, #Optional: nuclear spin, I, including parity. Tabulated by Stone et al.
              'mu_I' :  .11778476 , #Optional: Nuclear dipole moment as tabulated by Stone.             
               }




InputdeckW0 = {'s_ground':3 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
              's_excited':3 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
              'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
              'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
              'E_ground': np.array([2951.29 ]) ,#Lowest J first, in cm^-1
              'E_excited': np.array([26229.77,27488.11,27889.68]) , #Lowest J first, in cm^-1
              'Bmag': 0 , #Magnetic Field  [T]
              'b_angle': 65 , #Angle between LoS and Bmax
              'plottitle':'W I []6s (^7S_3) -> []6p (^7P_2) B=0.6T' , #Title for plotting (optional)     
              'plot_window' : [425,435], #Min and max for plot window range (nm)
              'amu' : 183 ,  #Weight in AMU
              'specstep' : 0.002222 , #Resolution of spectrometer in nm
               }

#Do the splitting calculation for the three cases
W0 = Zeeman_Main(InputdeckW0)
W1 = Zeeman_Main(InputdeckW)
W_HFS = Zeeman_Main(InputdeckWHFS)


#Make plotvars for plotfunction
WPlotvars = [['green', 2.5, '-.',InputdeckW['plottitle']] ,['orange', 2.5, '--',InputdeckW0['plottitle']] ,  ['brown', 2.5, '--',InputdeckWHFS['plottitle']]]
#Make plot object
WPlotObject = [W1,W0]

#Add HFS if desired
if doHFS==True:
    WPlotObject.append(W_HFS)
    
#Plot it
PlotFunction(WPlotObject,WPlotvars,plotwind = [429.4,429.55],plottitle='WI Spectra with LoS angle 65 degrees', SpectrumPlot =[WSpect['wavelength'],WSpect['signals'][0]] , NormalizeSig = True )
plt.rcParams['figure.figsize'] = [14,8]
plt.vlines(W1['wave_air'], np.zeros_like(W1['signal']),W1['signal'],color='green')
plt.savefig('W_ex1_430.jpg',dpi = 1200)




#%%
#This is a standard zeeman fan plot, the Bmag should be a range of values now
Brange = np.linspace(0,5,500)
InputdeckW2 = {'s_ground':3 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
              's_excited':3 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
              'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
              'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
              'E_ground': np.array([2951.29 ]) ,#Lowest J first, in cm^-1
              'E_excited': np.array([26229.77,27488.11,27889.68]) , #Lowest J first, in cm^-1
              'Bmag': Brange , #Magnetic Field  [T]
              'b_angle': 65 , #Angle between LoS and Bmax
              'plottitle':'W I []6s (^7S_3) -> []6p (^7P_2) B=0.6T' , #Title for plotting (optional)     
              'plot_window' : [425,435], #Min and max for plot window range (nm)
              'amu' : 183 ,  #Weight in AMU
              'specstep' : 0.002222 , #Resolution of spectrometer in nm
               }

# plotZfan(InputdeckW2)
  

