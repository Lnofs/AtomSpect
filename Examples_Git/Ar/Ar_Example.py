# -*- coding: utf-8 -*-
"""

@author: Leo Nofs
"""

import csv, os, sys
import numpy as np
from matplotlib import pyplot as plt

#This adds the directy two folders up from the example file so taht the main module can be imported.
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir, os.path.pardir)))
from SpectZHFS import Zeeman_Main, Normalize,PlotFunction,read_Spectra

#Add path to the spectral data
pathfil = os.getcwd()
filname = '\Ar_specs.csv' #There are multiple magnetic field values in this file which will be split and parsed using read_Spectra utility.
fil = pathfil + filname

#Import and process the spectral data
Spec_Raw = read_Spectra(fil, headercount=2, keepheaders=True)
Ar_Spec = Spec_Raw[0]
Ar_B = [float(x) for x in  Spec_Raw[1][0][1:]]
   
#This corrects for a misalignment in the calibration by 3 pixels. 

shiftedAr = [x + (5-.3 -34*.002222) for x in Ar_Spec['wavelength'] ]  #Shifting because I think the raw data wasn't calibrated right


#Iterate over all of the Bfield values in the imported data.
for i,Bvals in enumerate(Ar_B):
    fig, axs = plt.subplots(1)
  
    InputdeckAr1 = {'s_ground': 0 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                  's_excited': 0 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                  'l_ground': 1 , #Orbital Angular Momentum of ground state, int or half int
                  'l_excited': 0 ,#Orbital Angular Momentum of excited state, int or half int
                  'E_ground': np.array([95399.8276]) ,#Lowest J first, in cm^-1
                  'E_excited': np.array([108722.6194]) , #Lowest J first, in cm^-1
                  'Bmag': Bvals , #Magnetic Field  [T]
                  'b_angle': 90 , #Angle between LoS and Bmax
                  'plottitle': f'Ar I 3p5 4s (1P) -> 3p5 4p (1S) B={Bvals}' , #Title for plotting (optional)     
                  'plot_window' : [745,755], #Min and max for plot window range (nm), convolutions will only take place within this range
                  'amu' : 40  #Weight in AMU
                   }
    
    #The energies here are slightly shifted to better match the observed spectra.
    InputdeckAr2= {'s_ground': 1 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                  's_excited': 1 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                  'l_ground' : 1 , #Orbital Angular Momentum of ground state, int or half int
                  'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                  'E_ground' : np.array([94541.62, 93750.5978,93131.89]) , #Lowest J first, in cm^-1
                  'E_excited': np.array([107054.6, 107482.7, 106224]) , #Lowest J first, in cm^-1 
                  'Bmag': Bvals , #Magnetic Field  [T]
                  'b_angle': 90 , #Angle between LoS and Bmax
                  'plottitle': f'Ar I 3p5 4s (3P) -> 3p5 4p (3P)  B={Bvals}' ,
                  'plot_window' : [745,755], #Min and max for plot window range (nm)
                  'amu' : 40  #Weight in AMU
                   }

           
    
    Ar1 = Zeeman_Main(InputdeckAr1)    
    axs.plot(shiftedAr,Ar_Spec['signals'][i] ,color ='black' , label='Spectrometer')
    axs.plot(Ar1['SpecOut'][0],Normalize(Ar1['SpecOut'][1]) , label = InputdeckAr1['plottitle'])# ,color =colorlist[i])
    
    #This spectra has two different Argon transitions in the range. So we fit the second peak to the local maximum in the second half by slicing.
    Ar2 = Zeeman_Main(InputdeckAr2)            
    axs.plot(Ar2['SpecOut'][0],Normalize(Ar2['SpecOut'][1] , scaling = np.max(Ar_Spec['signals'][i][690:])) ,  label=InputdeckAr2['plottitle'])# ,color =colorlist[i])
      
    axs.set_xlim(750,752)
    axs.set_ylim(-.1,1.1)
    plt.ylabel("Normalized Intensity (Arb units)")
    plt.xlabel('Wavelength (nm)')
    plt.tight_layout()
    plt.legend(loc='best')
    
    
    # plt.savefig(f'ArI B={Bvals}.png',dpi=800)

        
        

