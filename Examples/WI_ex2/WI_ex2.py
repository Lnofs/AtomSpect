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
from SpectSplit  import Zeeman_Main, Normalize,PlotFunction

#Add current working dirrectory to system path to be able to load the spectrometer files.
pathfil = os.getcwd()
sys.path.append(pathfil)

spec405 = np.load('spec_405.npy')
wav405 = np.load('wave_405.npy')
spec430 = np.load('spec_430.npy')
wav430 = np.load('wave_430.npy')

doHFS=True
savefigs = False

#Note that these W spectra do not fully match and is an active area of research.
plt.close('all')
Bmag = 4.38

# Bangles = [35]
Bangles = [55]
# Polfilter = 57.5


Temp = 1
instrumres = 0.004

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
                  'specstep' : 0.018 , #Resolution of spectrometer in nm
                  'Convfxn' : 'Skewed' , #Optional: 'Gaussian', 'Skewed'.
                  'Skewness' : [0.62,0.99] , #How non-symmetric the Lorrentzian is

                  'Temp': 30*11602,  #Temperature in K. Used for Gaussian convolution
                  # 'Pol_angle' : Polfilter , 
                  'specres' : 2 ,
                   }

    # InputdeckW = {'s_ground':3 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
    #               's_excited':3 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
    #               'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
    #               'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
    #               'E_ground': np.array([2951.29 ]) ,#Lowest J first, in cm^-1
     
    #               'E_excited': np.array([26229.77,27488.11,27889.68]) , #Lowest J first, in cm^-1
    #               'Bmag': Bmag , #Magnetic Field  [T]
    #               'b_angle': Bangle , #Angle between LoS and Bmax
    #               # 'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
    #               # 'Temp' : Temp*11604,
    #               'Convfxn': 'Skewed', #Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
    #               'Skewness' : [0.62,0.99] , #How non-symmetric the Lorrentzian is
    #               'specstep' :instrumres , #Resolution of spectrometer in nm
    #               'plottitle':'W I []6s (^7S_3) -> []6p (^7P_2) B=0.6T' , #Title for plotting (optional)     
    #               'plot_window' : [425,435], #Min and max for plot window range (nm)
    #               'amu' : 183 ,  #Weight in AMU
    #                }
    
    
    
    
    InputdeckWHFS = {'s_ground':3 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                  's_excited':3 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                  'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
                  'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                  'E_ground': np.array([2951.29 ]) ,#Lowest J first, in cm^-1
                  'E_excited': np.array([26229.77,27488.11,27889.68]) , #Lowest J first, in cm^-1
                  'Bmag': 0.6 , #Magnetic Field  [T]
                  'b_angle': Bangle , #Angle between LoS and Bmax
                  'plottitle':'W I []6s (^7S_3) -> []6p (^7P_2) B=0.6T with HFS (PsiHFS=True)' , #Title for plotting (optional)     
                  'plot_window' : [425,435], #Min and max for plot window range (nm)
                  'amu' : 183 ,  #Weight in AMU
                  # 'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
                  # 'Temp' : Temp*11604,
                  'specstep' :instrumres , #Resolution of spectrometer in nm
                  # 'Skewness' : [0.62,0.99] , #How non-symmetric the Lorrentzian is
                  
                  'Convfxn': 'Skewed', #Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
                  'Skewness' : [0.62,0.99] , #How non-symmetric the Lorrentzian is
                  # 'Pol_angle': 65 ,#Angle polarizing filter makes with max linear transmission, Optional
    
                  'HFS_G' : [505e6] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                  'HFS_E': [0, 496.2e6, 440.5e6] , #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
    
                  # 'HFS_G' : [AHFS_G] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                  # 'HFS_E': AHFSE , #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                  'I_spin' : 0.5, #Optional: nuclear spin, I, including parity. Tabulated by Stone et al.
                  'mu_I' :  .11778476 , #Optional: Nuclear dipole moment as tabulated by Stone.    
                  'PsiHFS': True, 
                  # 'specres': 1 ,
                   }
    
    
    
    InputdeckW0 = {'s_ground':3 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                  's_excited':3 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                  'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
                  'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                  'E_ground': np.array([2951.29 ]) ,#Lowest J first, in cm^-1
                  'E_excited': np.array([26229.77,27488.11,27889.68]) , #Lowest J first, in cm^-1
                  'Bmag': 0 , #Magnetic Field  [T]
                  'b_angle': Bangle , #Angle between LoS and Bmax
                  'plottitle':'W I []6s (^7S_3) -> []6p (^7P_2) B=0T' , #Title for plotting (optional)     
                  'plot_window' : [425,435], #Min and max for plot window range (nm)
                  'amu' : 183 ,  #Weight in AMU
                  'Convfxn': 'Skewed', #Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
                  'Skewness' : [0.62,0.99] , #How non-symmetric the Lorrentzian is
                  # 'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
                  # 'Temp' : Temp*11604,
                  'specstep' :instrumres , #Resolution of spectrometer in nm
                   }
    
    InputdeckWHFS0 = {'s_ground':3 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                  's_excited':3 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                  'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
                  'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                  'E_ground': np.array([2951.29 ]) ,#Lowest J first, in cm^-1
                  'E_excited': np.array([26229.77,27488.11,27889.68]) , #Lowest J first, in cm^-1
                  'Bmag': 0 , #Magnetic Field  [T]
                  'b_angle': Bangle , #Angle between LoS and Bmax
                  'plottitle':'W I []6s (^7S_3) -> []6p (^7P_2) B=0T with HFS(PsiHFS=True)' , #Title for plotting (optional)     
                  'plot_window' : [425,435], #Min and max for plot window range (nm)
                  'amu' : 183 ,  #Weight in AMU
                  # 'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
                  # 'Temp' : Temp*11604,
                  'specstep' : instrumres , #Resolution of spectrometer in nm
                  # 'Skewness' : [0.62,0.99] , #How non-symmetric the Lorrentzian is
                  
                  'Convfxn': 'Skewed', #Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
                  'Skewness' : [0.62,0.99] , #How non-symmetric the Lorrentzian is
                  # 'Pol_angle': 65 ,#Angle polarizing filter makes with max linear transmission, Optional
    
                  'HFS_G' : [505e6] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                  'HFS_E': [0, 496.2e6, 440.5e6] , #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
    
                  # 'HFS_G' : [AHFS_G] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                  # 'HFS_E': AHFSE , #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                  'I_spin' : 0.5, #Optional: nuclear spin, I, including parity. Tabulated by Stone et al.
                  'mu_I' :  .11778476 , #Optional: Nuclear dipole moment as tabulated by Stone.    
                  'PsiHFS': True, 
    
                   }
        
specdat1 = [wav405,spec405]
specdat2 = [wav430,spec430]  
Bangle = 90
Temps = [0.025,1,3]
specsteps = [0.05,0.01,0.005,0.002]
functions = ['Lorrentzian','Gaussian','GaussianInst']
B_mag = [0.6,1.,3]
specres= 1

savefig = False

Temp_ev = 1
function='Gaussian'
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
            InputdeckWHFS = {'s_ground':3 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
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
                          'HFS_G' : [505e6] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                          'HFS_E': [0, 496.2e6, 440.5e6] , #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                          # 'HFS_G' : [AHFS_G] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                          # 'HFS_E': AHFSE , #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                          'I_spin' : 0.5, #Optional: nuclear spin, I, including parity. Tabulated by Stone et al.
                          'mu_I' :  .11778476 , #Optional: Nuclear dipole moment as tabulated by Stone.    
                          'PsiHFS': True, 
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

    