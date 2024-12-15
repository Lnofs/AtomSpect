# -*- coding: utf-8 -*-
"""
Created on Sat Jul 27 10:06:38 2024

@author: Leo Nofs
"""

import csv, os, sys
import numpy as np
from matplotlib import pyplot as plt

#This is the base folder of the main function. Can be excluded if data files are in the same location as the 
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir, os.path.pardir)))
from SpectZHFS  import Zeeman_Main, Normalize,PlotFunction

#Add current working dirrectory to system path to be able to load the spectrometer files.
pathfil = os.getcwd()
sys.path.append(pathfil)

spec405 = np.load('spec_405.npy')
wav405 = np.load('wave_405.npy')
spec430 = np.load('spec_430.npy')
wav430 = np.load('wave_430.npy')

#Note that these W spectra do not fully match and is an active area of research.

Bmag = 2.3768942922

Bangles = [35]

for Bangle in Bangles:

    #These levels give lines at 400, 407, and 429nm.

    InputdeckW = {'s_ground':3 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                  's_excited':3 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                  'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
                  'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                  'E_ground': np.array([2951.29 ]) ,#Lowest J first, in cm^-1
                  'E_excited': np.array([26229.77,27488.11,27889.68]) , #Lowest J first, in cm^-1
                  'Bmag': Bmag , #Magnetic Field  [T]
                  'b_angle': Bangle , #Angle between LoS and Bmax
                  'plottitle':f'W I []6s (^7S_3) -> []6p (^7P_2) B=2.38T, LoS = {Bangle}' , #Title for plotting (optional)     
                  'plot_window' : [395,435], #Min and max for plot window range (nm)
                  'amu' : 183 ,  #Weight in AMU
                  'specstep' : 0.012222 , #Resolution of spectrometer in nm
                  'Convfxn' : 'Skewed' , #Optional: 'Gaussian', 'Skewed'.
                  'Temp': 300,  #Temperature in K. Used for Gaussian convolution
                   }



        
    W_ex2 = Zeeman_Main(InputdeckW)
    WPlotvars = ['green', 2.5, '-.',InputdeckW['plottitle']]
   
 

    PlotFunction(W_ex2,WPlotvars,plotwind = [400,401], SpectrumPlot =[wav405, Normalize(spec405)],plottitle="WI B = 2.3T , 35 degrees" , NormalizeSig = True )
    # plt.savefig('W_ex2_400.jpg',dpi = 1200)
 
    PlotFunction(W_ex2,WPlotvars,plotwind = [406.5,408], SpectrumPlot =[wav405, Normalize(spec405)],plottitle="WI B = 2.3T , 35 degrees" , NormalizeSig = True)
    # plt.savefig('W_ex2_407.jpg',dpi = 1200)

    
    PlotFunction(W_ex2, WPlotvars,plotwind=[429,430.5],SpectrumPlot=[wav430,Normalize(spec430)],plottitle="WI B = 2.3T"  , NormalizeSig = True)
    # plt.savefig('W_ex2_430.jpg',dpi = 1200)




    