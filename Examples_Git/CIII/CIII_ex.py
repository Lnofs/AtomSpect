# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 17:08:11 2024

@author: Leo Nofs
"""

import csv, os, sys
import numpy as np
from matplotlib import pyplot as plt


#This is the base folder of the main function. Can be excluded if data files are in the same location as the 
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir, os.path.pardir)))
from SpectZHFS  import Zeeman_Main, Normalize,PlotFunction,Vac_to_air

pathfil = os.getcwd()

#Add the path for the npy spectra objects

sys.path.append(pathfil)


CIIISpec = np.load('CIII_2.24T.npy')

#Convert the spectrometer data from vacuum to air wavelengths
CIII_wav = [Vac_to_air(x) for x in CIIISpec[0]]


#%%
Los_Angle = 68.6 #Theta angle of LoS in degrees
Temperature = 40 #Temp in eV

InputdeckCIII_4649 = {'s_ground':1 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                      's_excited':1 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                      'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
                      'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                      'E_ground': np.array([238213])  ,#Lowest J first, in cm^-1
                      'E_excited':  np.array([259705.55,259711.22,259724.30]) , #Lowest J first, in cm^-1
                      'Bmag': 2.24 , #Magnetic Field  [T]
                      'b_angle': Los_Angle , #Angle between LoS and Bmax
                      'plot_window' : [462,470], #Min and max for plot window range (nm)
                      'amu' : 12.011 ,  #Weight in AMU
                      'specstep' : 0.005 , #Resolution of spectrometer in nm
                      'Convfxn': 'Gaussian', #Optional: 'Gaussian', 'Skewed'.
                      'Temp' : Temperature*11602, #Temperature in K. Used for Gaussian convolution
                      'plottitle': f'CIII [1s^2]2s3s (^3S_1) -> [1s^2]2s3p (^3P_J): B=2.24T , LoS = {Los_Angle}', #Title for plotting (optional)     

                      'fxnwindow': 1 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                      
                   }


InputdeckCIII_4663 = {'s_ground':1 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                      's_excited':1 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                      'l_ground': 1 , #Orbital Angular Momentum of ground state, int or half int
                      'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                      'E_ground': np.array([308216.58,308248.91,308317.29]) ,#Lowest J first, in cm^-1
                      'E_excited':  np.array([329685.38,329706.47,329743.57]) , #Lowest J first, in cm^-1
                      'Bmag': 2.24 , #Magnetic Field  [T]
                      'b_angle': Los_Angle , #Angle between LoS and Bmax
                      'plottitle': f'CIII [1s^2]2p3s (^3S_1) -> [1s^2]2p3p (^3P_J): B=2.24T , LoS = {Los_Angle}' , #Title for plotting (optional)     
                      'plot_window' : [462,470], #Min and max for plot window range (nm)
                      'amu' : 12.011 ,  #Weight in AMU
                      'specstep' : 0.025 , #Resolution of spectrometer in nm
                      'Convfxn': 'Skewed', #Optional: 'Gaussian', 'Skewed'.
                      'Temp' : 300, #Temperature in K. Used for Gaussian convolution
                      'fxnwindow': 0.5 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                      
                   }






CIII_4649 = Zeeman_Main(InputdeckCIII_4649)
CIII_4663 = Zeeman_Main(InputdeckCIII_4663)
    
    
CIII_Plotvars = [['green', 2.5, '-.',InputdeckCIII_4649['plottitle']]  ,  ['brown', 2.5, '--',InputdeckCIII_4663['plottitle']]]
CIII_scalevars = [1,1]
PlotFunction(CIII_4649,CIII_Plotvars[0],plotwind = [464.5,465.5], SpectrumPlot =[CIII_wav, Normalize(CIIISpec[1])],plottitle="CIII B = 2.24T , 68.6 degrees, Temp = 40eV" , NormalizeSig = True, NormalizeScale = CIII_scalevars  )

plt.vlines(CIII_4649['wave_air'], np.zeros_like(CIII_4649['signal']),CIII_4649['signal'],color='brown', label='Leos code')

# plt.savefig('CIII4649_Gaussian.png',dpi=1200)



   #%%

    

  
