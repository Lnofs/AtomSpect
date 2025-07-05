# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 17:08:11 2024

@author: Leo Nofs
"""

import csv, os, sys
import numpy as np
from matplotlib import pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir, os.path.pardir)))
from SpectSplit  import Zeeman_Main, PlotFunction2,read_Spectra, plotZfan

from HFS_Const_Calcs import Ahfs, Ahfs_Calc




#%%Calculate HFSA for ground
W7Sin = {'Z' : 74   , #Atomic number
        'AMU': 183  , #Atomic mass in amu
        'E_ionize': 16.37 , #Ionization energy [eV]
        'E1': .40624084     , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
        'n1': 5    , #Princple quantum number.        
        'I' : 0.5    , #Nuclear spin angular momentum
        'muI': 0.11778476 ,    #Nuclear dipole moment (in units of mu_N)
        'l1'  :  0   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
        's1' : 2.5,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".
        'j1'  : 2.5   , #J = L +/- S coupled angular momentum 
        'Corrflags': 'deleps' , #Correction factors for the hyperfine calculation.
        #IF coupling is to be done, need J,L,S, j2,l2,s2 for the different possible schemes. It's possible that j1 may not be always known (especially outside of S levels).
        'Coupling': 'ss',
        'L' : 0 , #Coupled total L from the two electrons
        'S' : 3 , #Coupled total S of the two electrons
        'J' : 3 , #Coupled total J for the two electrons
        'A_Exp': 496, 

     }
W7SOutG = {'Z' : 74   , #Atomic number
        'AMU': 183  , #Atomic mass in amu
        'E_ionize': 7.86403 , #Ionization energy [eV]
        'E1': -.18667     , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
        'n1': 6     , #Princple quantum number.        
        'I' : 0.5    , #Nuclear spin angular momentum
        'muI': 0.11778476 ,    #Nuclear dipole moment (in units of mu_N)

        'l1'  :  0   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
        's1' :  .5 ,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".
        'j1'  : .5   , #J = L +/- S coupled angular momentum 
  
       'E2' :5.310495 ,
       'n2' : 7 ,
        # 'Corrflags': 'FRHdelepsds' , #Correction factors for the hyperfine calculation.


     }
try: 
    AHFS_G = 1e6*Ahfs(W7Sin, Input2=W7SOutG)['Acoupled']
except:
    AHFS_G = 1e6*Ahfs(W7Sin)['A']


# AHFS_G = Ahfs(W7Sin)['A']*1e6
print(AHFS_G/1e6)

#%%Calculate HFSA for Excited


W7SOut = {'Z' : 74   , #Atomic number
        'AMU': 183  , #Atomic mass in amu
        'E_ionize': 7.86403 , #Ionization energy [eV]
        'E1': -.18667     , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
        'n1': 6     , #Princple quantum number.        
        'I' : 0.5    , #Nuclear spin angular momentum
        'muI': 0.11778476 ,    #Nuclear dipole moment (in units of mu_N)

        'l1'  :  1   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
        's1' :  -.5 ,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".
        'j1'  : .5   , #J = L +/- S coupled angular momentum 
  
       'E2' :5.310495 ,
       'n2' : 7 ,
        'Corrflags': 'FRdelepsds' , #Correction factors for the hyperfine calculation.


     }

Jvals = [2,3,4]
AHFSE = []
for Js in Jvals:
    W7PJ = {'Z' : 74   , #Atomic number
            'AMU': 183  , #Atomic mass in amu
            'E_ionize': 16.37 , #Ionization energy [eV]
            'E1': .40624084     , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
            'n1': 5    , #Princple quantum number.        
            'I' : 0.5    , #Nuclear spin angular momentum
            'muI': 0.11778476 ,    #Nuclear dipole moment (in units of mu_N)
            'l1'  :  0   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
            's1' : 2.5,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".
            'j1'  : 2.5   , #J = L +/- S coupled angular momentum 
           'Corrflags': 'FRHdeleps' , #Correction factors for the hyperfine calculation.
            #IF coupling is to be done, need J,L,S, j2,l2,s2 for the different possible schemes. It's possible that j1 may not be always known (especially outside of S levels).
            'Coupling': 'sl2',
            'L' : 1 , #Coupled total L from the two electrons
            'S' : 3 , #Coupled total S of the two electrons
            'J' : Js , #Coupled total J for the two electrons
            'A_Exp': [0 ,496e6,440e6], #I get [604,415,39]
    
         }
    
    W7PinA = Ahfs(W7PJ)
    
    # print(Ahfs(W7PJ)['A']*1e6)
    # AHFSE.append(AHFS_E)

    AHFS_E = Ahfs(W7PJ, Input2=W7SOut)
    try:
        AHFSE.append(1e6*AHFS_E['Acoupled'])
    except:
        AHFSE.append(1e6*AHFS_E['A'])
 


    # print(AHFS_E['Acoupled'])
print(np.array(AHFSE)/1e6)

#%%
doHFS=True
savefigs = False
pathfil = os.getcwd()
filname = '\Wspec_ex1.csv'
fil = pathfil + filname

WSpect = read_Spectra(fil, headercount=1)


# plt.close('all')

Bangle = 65
Temp = 1
instrumres = 0.002
InputdeckW = {'s_ground':3 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
              's_excited':3 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
              'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
              'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
              'E_ground': np.array([2951.29 ]) ,#Lowest J first, in cm^-1
 
              'E_excited': np.array([26229.77,27488.11,27889.68]) , #Lowest J first, in cm^-1
              'Bmag': 0.6 , #Magnetic Field  [T]
              'b_angle': Bangle , #Angle between LoS and Bmax
              # 'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
              # 'Temp' : Temp*11604,
              'Convfxn': 'Skewed', #Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
              'Skewness' : [0.62,0.99] , #How non-symmetric the Lorrentzian is
              'specstep' :instrumres , #Resolution of spectrometer in nm
              'plottitle':'W I []6s (^7S_3) -> []6p (^7P_2) B=0.6T' , #Title for plotting (optional)     
              'plot_window' : [425,435], #Min and max for plot window range (nm)
              'amu' : 183 ,  #Weight in AMU
               }




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

              # 'HFS_G' : [505e6] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
              # 'HFS_E': [0, 496.2e6, 440.5e6] , #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                 'HFS_G' : [AHFS_G] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                 'HFS_E': AHFSE , #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
              'I_spin' :0.5, #Optional: nuclear spin, I, including parity. Tabulated by Stone et al.
              'mu_I' :  .11778476 , #Optional: Nuclear dipole moment as tabulated by Stone.    
              'PsiHFS': True, 
              # 'sortE' : True,
              # 'specres': 1 ,
               }


InputdeckWHFSB = {'s_ground':3 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
              's_excited':3 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
              'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
              'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
              'E_ground': np.array([2951.29 ]) ,#Lowest J first, in cm^-1
              'E_excited': np.array([26229.77,27488.11,27889.68]) , #Lowest J first, in cm^-1
              'Bmag': 0.6 , #Magnetic Field  [T]
              'b_angle': Bangle , #Angle between LoS and Bmax
              'plottitle':'W I []6s (^7S_3) -> []6p (^7P_2) B=0.6T with HFS' , #Title for plotting (optional)     
              'plot_window' : [425,435], #Min and max for plot window range (nm)
              'amu' : 183 ,  #Weight in AMU
              # 'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
              # 'Temp' : Temp*11604,
              'specstep' :instrumres , #Resolution of spectrometer in nm
              # 'Skewness' : [0.62,0.99] , #How non-symmetric the Lorrentzian is
              
              'Convfxn': 'Skewed', #Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
              'Skewness' : [0.62,0.99] , #How non-symmetric the Lorrentzian is

              # 'HFS_G' : [505e6] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
              # 'HFS_E': [1e6,1e9, 1e9] , #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]

              'HFS_G' : [AHFS_G] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
              'HFS_E': AHFSE , #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
              'I_spin' : 0.5, #Optional: nuclear spin, I, including parity. Tabulated by Stone et al.
              'mu_I' :  .11778476 , #Optional: Nuclear dipole moment as tabulated by Stone.    
              # 'PsiHFS': True, 
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

              # 'HFS_G' : [505e6] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
              # 'HFS_E': [0, 496.2e6, 440.5e6] , #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]

              'HFS_G' : [AHFS_G] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
              'HFS_E': AHFSE , #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
              'I_spin' : 0.5, #Optional: nuclear spin, I, including parity. Tabulated by Stone et al.
              'mu_I' :  .11778476 , #Optional: Nuclear dipole moment as tabulated by Stone.    
              'PsiHFS': True, 

               }

InputdeckWHFS0B = {'s_ground':3 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
              's_excited':3 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
              'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
              'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
              'E_ground': np.array([2951.29 ]) ,#Lowest J first, in cm^-1
              'E_excited': np.array([26229.77,27488.11,27889.68]) , #Lowest J first, in cm^-1
              'Bmag': 0 , #Magnetic Field  [T]
              'b_angle': Bangle , #Angle between LoS and Bmax
              'plottitle':'W I []6s (^7S_3) -> []6p (^7P_2) B=0T with HFS' , #Title for plotting (optional)     
              'plot_window' : [425,435], #Min and max for plot window range (nm)
              'amu' : 183 ,  #Weight in AMU
              # 'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
              # 'Temp' : Temp*11604,
              'specstep' : instrumres , #Resolution of spectrometer in nm

              'Convfxn': 'Skewed', #Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
              'Skewness' : [0.62,0.99] , #How non-symmetric the Lorrentzian is

              'HFS_G' : [AHFS_G] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
              'HFS_E': AHFSE , #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
              'I_spin' : 0.5, #Optional: nuclear spin, I, including parity. Tabulated by Stone et al.
              'mu_I' :  .11778476 , #Optional: Nuclear dipole moment as tabulated by Stone.    
              # 'PsiHFS': True, 

               }



#Do the splitting calculation for the three cases
W0 = Zeeman_Main(InputdeckW0)
W1 = Zeeman_Main(InputdeckW)


#Make plotvars for plotfunction
WPlotvars = [['green', 2.5, '-.',InputdeckW['plottitle']] ,['orange', 2.5, '--',InputdeckW0['plottitle']]         ]
#Make plot object
WPlotObject = [W1,W0]

# WPlotvars=[]
# WPlotObject=[]
if doHFS==True:
    W_HFS = Zeeman_Main(InputdeckWHFS)
    WPlotvars.append(['brown', 2.5, '--',InputdeckWHFS['plottitle']])    
    WPlotObject.append(W_HFS)

    W_HFSB = Zeeman_Main(InputdeckWHFSB)
    WPlotvars.append(['blue', 2.5, '--',InputdeckWHFSB['plottitle']])    
    WPlotObject.append(W_HFSB)
    
    # W_HFS0 = Zeeman_Main(InputdeckWHFS0)
    # WPlotvars.append(['magenta', 2.5, ':',InputdeckWHFS0['plottitle']])
    # WPlotObject.append(W_HFS0)
    
    # W_HFS0B = Zeeman_Main(InputdeckWHFS0B)
    # WPlotvars.append(['cyan', 2.5, ':',InputdeckWHFS0B['plottitle']])
    # WPlotObject.append(W_HFS0B)
    
title = f'Skewed Lorrentzian, Pixels = {instrumres} nm'

#Plot it
PlotFunction2(WPlotObject,WPlotvars,plotwind = [429.4,429.55],plottitle=title, SpectrumPlot =[WSpect['wavelength'],WSpect['signals']] , NormalizeSig = True )
plt.rcParams['figure.figsize'] = [14,8]
plt.vlines(W1['wave_air'], np.zeros_like(W1['signal']),W1['signal'],color='green')
plt.tight_layout()

if doHFS==True:
    plt.vlines(W_HFS['wave_air'], np.zeros_like(W_HFS['signal']),W_HFS['signal'],color='maroon')


# if savefigs:
#     plt.savefig('W_ex1_430.jpg',dpi = 1200)



#%%CompareHFSResults

# #Make plotvars for plotfunction
# WPlotvars = [['green', 2.5, '-.',InputdeckW['plottitle']] ,['orange', 2.5, '--',InputdeckW0['plottitle']]         ]
# #Make plot object
# WPlotObject = [W1,W0]

# WPlotvars=[]
# WPlotObject=[]
# if doHFS==True:
#     # W_HFS = Zeeman_Main(InputdeckWHFS)
#     # WPlotvars.append(['brown', 2.5, '--',InputdeckWHFS['plottitle']])    
#     # WPlotObject.append(W_HFS)

#     # W_HFSB = Zeeman_Main(InputdeckWHFSB)
#     # WPlotvars.append(['blue', 2.5, '--',InputdeckWHFSB['plottitle']])    
#     # WPlotObject.append(W_HFSB)
    
#     W_HFS0 = Zeeman_Main(InputdeckWHFS0)
#     WPlotvars.append(['purple', 2.5, ':',InputdeckWHFS0['plottitle']])
#     WPlotObject.append(W_HFS0)
    
#     W_HFS0B = Zeeman_Main(InputdeckWHFS0B)
#     WPlotvars.append(['blue', 2.5, ':',InputdeckWHFS0B['plottitle']])
#     WPlotObject.append(W_HFS0B)
    
# title = 'Comparison of H_HFS with PsiHFS flag'

# #Plot it
# PlotFunction(WPlotObject,WPlotvars,plotwind = [429.4,429.55],plottitle=title, SpectrumPlot =[WSpect['wavelength'],WSpect['signals'][0]] , NormalizeSig = True )
# plt.rcParams['figure.figsize'] = [14,8]
# plt.vlines(W1['wave_air'], np.zeros_like(W1['signal']),W1['signal'],color='green')
# plt.tight_layout()


# if doHFS==True:
#     plt.vlines(W_HFS['wave_air'], np.zeros_like(W_HFS['signal']),W_HFS['signal'],color='maroon')


# if savefigs:
#     plt.savefig('W_ex1_430_PsiHFSB=0.jpg',dpi = 1200)





