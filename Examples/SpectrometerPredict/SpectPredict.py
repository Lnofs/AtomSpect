# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 10:06:38 2024

@author: Leo Nofs
This code demonstrates using sliders to predict what the spectrometer might see.

Also showcases using the slider and different data formats.
"""

import csv, os, sys
import numpy as np
from matplotlib import pyplot as plt


sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir, os.path.pardir)))
from AtomSpect  import AtomSpect_Main, Normalize,PlotFunction,read_Spectra, plotZfan,Convol_Sticks, make_stickbins,MakeSlider, ZeemanFan,Vac_to_air


pathfil = os.getcwd()





              
#%%


pathfil = os.getcwd()

fil = pathfil + '/He_Spectra_NonPolarized.csv'
He_NonPol_Spec = read_Spectra(fil,scalewave = 0.1, headercount=1)  #Need to scale the wavelengths from angstrom to nm.

#%%

InputHe = {'s_ground': 1 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
              's_excited': 1 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
              'l_ground': 1 , #Orbital Angular Momentum of ground state, int or half int
              'l_excited':2 ,#Orbital Angular Momentum of excited state, int or half int
              'E_ground': np.array([169087.8308131,  169086.8428979 ,  169086.7664725  ]) ,#Lowest J first, in cm^-1
              'E_excited': np.array([191444.49804915,191444.47952984,191444.47832914]) , #Lowest J first, in cm^-1
              'Bmag': 1.5 , #Magnetic Field  [T]
              'b_angle': 90 , #Angle between LoS and Bmax
              'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
              'specstep' : 0.002222 , #Resolution of spectrometer in nm
              # 'Skewness' : [0.62,0.99] , #How non-symmetric the Lorrentzian is
              # 'Pol_angle': 90 ,#Angle polarizing filter makes with max linear transmission, Optional
              'plottitle': 'He I 1s2p -> 1s4d: B=3.2T, UnPolarized'  , #Title for plotting (optional)     
              'plot_window' : [447.1,447.2], #Min and max for plot window range (nm), convolutions will only take place within this range
              'amu' : 2  ,#Weight in AMU
              'specres' : 10 , #How many steps per resolution are calculated. Higher makes a smoother curve.
              'DoLowSig' : 'Yes' , #If this exists in the input, lowfield signal strength will be calculated. 
              'DoHighSig': 'Yes' , #To use Strong field approximation, you MUST include EtermE and EtermG keys for term energy
 
              }  



InputdeckCIII_4649 = {'s_ground':1 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                      's_excited':1 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                      'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
                      'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                      'E_ground': np.array([238213])  ,#Lowest J first, in cm^-1
                      'E_excited':  np.array([259705.55,259711.22,259724.30]) , #Lowest J first, in cm^-1
                      'Bmag': 2.4, #Magnetic Field  [T]
                      'b_angle': 50 , #Angle between LoS and Bmax
                      'plot_window' : [464.5,465.5], #Min and max for plot window range (nm)
                      'amu' : 12.011 ,  #Weight in AMU
                      'specstep' : 0.0075 , #Resolution of spectrometer in nm
                      'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'.
                      'Temp' : 11602*7.9, #Temperature in K. Used for Gaussian convolution
                      'plottitle': r'CIII [1$s^2$]2s3s $(^3S_1)$ -> [1$s^2$]2s3p ($^3P_J$)', #Title for plotting (optional)     

                      # 'plottitle': f'CIII [1s^2]2s3s (^3S_1) -> [1s^2]2s3p (^3P_J): B=2.24T , LoS = {90}', #Title for plotting (optional)     
                      'ion_vel' : -1000, #Ion velocity in m/s
                      'fxnwindow': 1 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                      # 'Pol_angle': 50 ,#Angle polarizing filter makes with max linear transmission, Optional
                      # 'DoLowSig' : 'Yes' , #If this exists in the input, lowfield signal strength will be calculated. 
                      # 'DoHighSig' : 'Yes' , #If this exists in the input, lowfield signal strength will be calculated. 
                      
                   }
#
CIII_Z = AtomSpect_Main(InputdeckCIII_4649)


               

CIIISpecHighTemp = np.load('CIII_2.24T.npy') #This is the spectrometer data

CIIIwaveH = [Vac_to_air(x) for x in CIIISpecHighTemp[0]]
CIIISigH = Normalize(CIIISpecHighTemp[1])
CIII_PlotvarsH = [['green', 3, '-.','CIII']]                   


#%%
import xarray as xr
fil9 = f"{pathfil}/CIII_17eV.nc" #Low Temp from Supra dataset.

# sys.path.append(pathfil)


# fig1 = xr.open_dataset(fil, engine="h5netcdf")
fig9 = xr.open_dataset(fil9, engine="h5netcdf")

#Fiber44 of 20230215.015 available on Aurora dataset, publication states B = 2.51T, alpha = 71 deg, vi = 1km/s, ti = 17 eV
Fiber44sig = Normalize(fig9['spectrum_data'].values[1])
Fiber44wav = fig9['spectrum_data']['wavelength'].values[1]
AirFiber44 = [Vac_to_air(x) for x in Fiber44wav]

#%%
# plt.close('all')



# MakeSlider(InputHe,[He_NonPol_Spec['wavelength'],He_NonPol_Spec['signals'][0]],[447.1,447.2],do_sticks = 1,do_Polplot=0,polvary=0,banglevary=0)

# MakeSlider(InputdeckCIII_4649,[CIIIwaveH,CIIISigH],[464.5,465.5],do_sticks = 1,do_Polplot=1,polvary=0,banglevary=1)
# MakeSlider(InputdeckCIII_4649,[CIIIwaveL,CIIISigL],[464.5,465.5],do_sticks = 1,do_Polplot=1,polvary=0,banglevary=1)


MakeSlider(InputdeckCIII_4649,[AirFiber44,Fiber44sig],[464.5,465.5],do_sticks = 1,do_Polplot=1,polvary=0,banglevary=1)


