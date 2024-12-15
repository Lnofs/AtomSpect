# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 13:26:58 2024

@author: Leo Nofs
"""

import csv, os, sys
import numpy as np
from matplotlib import pyplot as plt


sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir, os.path.pardir)))
from SpectZHFS  import Zeeman_Main, Normalize,PlotFunction,read_Spectra, plotZfan,Convol_Sticks, make_stickbins



#Flags for doing polarized or not
doHePol = True
doHeNonPol=True


if doHePol==True:
    pathfil = os.getcwd()
    filname = '\He_Spectra_1.5T_Polarized.csv'
    fil = pathfil + filname
    
    He_Pol_Spec = read_Spectra(fil, headercount=2)

   
    #This corrects for a misalignment in the calibration by 3 pixels. 
    He_PolShift = [0.1*x+ 0.006666 for x in He_Pol_Spec['wavelength']]    
    
     

    polangles = [0,30,70,90]


    fig, axs = plt.subplots(4,sharex=True,sharey=True, figsize = (6,8))
    colorlist = ['r','g','b','brown','m','c','y','black','seagreen']


    for i,angles in enumerate(polangles):

        InputdeckHePol = {'s_ground': 1 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                      's_excited': 1 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                      'l_ground': 1 , #Orbital Angular Momentum of ground state, int or half int
                      'l_excited':2 ,#Orbital Angular Momentum of excited state, int or half int
                      'E_ground': np.array([169087.8308131,  169086.8428979 ,  169086.7664725  ]) ,#Lowest J first, in cm^-1
                      'E_excited': np.array([191444.49804915,191444.47952984,191444.47832914]) , #Lowest J first, in cm^-1
                      'Bmag': 1.5 , #Magnetic Field  [T]
                      'b_angle': 90 , #Angle between LoS and Bmax
                      'Pol_angle': angles ,#Angle polarizing filter makes with max linear transmission, Optional
                      'plottitle': 'He I 1s2p -> 1s4d: B=1.5T, Polarized'  , #Title for plotting (optional)     
                      'plot_window' : [447,447.2], #Min and max for plot window range (nm), convolutions will only take place within this range
                      'amu' : 2  #Weight in AMU
                       }  
 
        
        HePol = Zeeman_Main(InputdeckHePol)

    
        axs[i].plot(HePol['SpecOut'][0],Normalize(HePol['SpecOut'][1], scaling=np.max(He_Pol_Spec['signals'][0])) ,color =colorlist[i],label = f'Pol_angle = {angles}' )
        axs[i].plot(He_PolShift,Normalize(He_Pol_Spec['signals'][i], scaling=np.max(He_Pol_Spec['signals'][0])), color= 'black', linestyle=':',label = 'Spectrometer')
        axs[i].legend(loc='center left', bbox_to_anchor=(1, 0.5))
   

    axs[0].set_xlim(447.1,447.2)
    axs[0].set_ylim(-.1,1.1)
    plt.ylabel("Normalized Intensity (Arb units)")
    plt.xlabel('Wavelength (nm)')
    plt.tight_layout()
    # plt.savefig('He_Pol.png',dpi=800)
    plt.show()


#%%
if doHeNonPol == True:
    pathfil = os.getcwd()

    fil = pathfil + '\He_Spectra_NonPolarized.csv'
    He_NonPol_Spec = read_Spectra(fil,scalewave = 0.1, headercount=1)  #Need to scale the wavelengths from angstrom to nm.
    
    
    Bmags = [3.2,2.5,1.5]

    fig, axs = plt.subplots(3,sharex=True,sharey=True, figsize = (10,12))
    colorlist = ['r','g','b','brown','m','c','y','black','seagreen']


    for i,Bvals in enumerate(Bmags):
        if i==2: #Need to shift the 1.5T spectra by 3 pixels.
            He_NonPol_Spec['wavelength'] = [x + 0.00666 for x in He_NonPol_Spec['wavelength']]

          
        InputdeckHeNonPol = {'s_ground': 1 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                      's_excited': 1 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                      'l_ground': 1 , #Orbital Angular Momentum of ground state, int or half int
                      'l_excited':2 ,#Orbital Angular Momentum of excited state, int or half int
                      'E_ground': np.array([169087.8308131,  169086.8428979 ,  169086.7664725  ]) ,#Lowest J first, in cm^-1
                      'E_excited': np.array([191444.49804915,191444.47952984,191444.47832914]) , #Lowest J first, in cm^-1
                      'Bmag': Bvals , #Magnetic Field  [T]
                      'b_angle': 90 , #Angle between LoS and Bmax
                      'plottitle': 'He I 1s2p -> 1s4d: B=1.5T, UnPolarized'  , #Title for plotting (optional)     
                      'plot_window' : [447,447.2], #Min and max for plot window range (nm), convolutions will only take place within this range
                      'amu' : 2  #Weight in AMU
                       }  
        HeNonpol = Zeeman_Main(InputdeckHeNonPol)

    
        axs[i].plot(HeNonpol['SpecOut'][0],Normalize(HeNonpol['SpecOut'][1]) ,color ='blue',alpha = 0.6,label ='Calculated')
        axs[i].plot(He_NonPol_Spec['wavelength'],Normalize(He_NonPol_Spec['signals'][i]),  color= 'black', linestyle=':', label=f'He Spectrometer,B={Bvals} T')

        axs[i].vlines(HeNonpol['reduced_sticks'][0], np.zeros_like(HeNonpol['reduced_sticks'][0]) , HeNonpol['reduced_sticks'][1])
        axs[i].legend(loc='best')
        axs[0].set_xlim(447.1,447.2)
        axs[0].set_ylim(-.1,1.1)
        axs[i].set(ylabel="Normalized Intensity (Arb units)")


    plt.xlabel('Wavelength (nm)')
    plt.tight_layout()
    fig.suptitle(InputdeckHeNonPol['plottitle'])

    # plt.savefig(f'He_NonPolarized.png',dpi=1200)





