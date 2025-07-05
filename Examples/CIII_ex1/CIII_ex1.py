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
from SpectSplit  import Zeeman_Main, Normalize,PlotFunction,Vac_to_air, Convol_Spect, MultiSpec

pathfil = os.getcwd()

#Add the path for the npy spectra objects

sys.path.append(pathfil)


CIIISpec = np.load('CIII_2.24T.npy') #This is the spectrometer data

savefigs = False
doCombSpec = False
dopolfilterplot=True

#Convert the spectrometer data from vacuum to air wavelengths
CIII_wav = [Vac_to_air(x) for x in CIIISpec[0]]
Los_Angle = 68.6#Theta angle of LoS in degrees

Temperature = 28.1#Temp in eV
Polangle = 55.75


# plt.close('all')

InputdeckCIII_4649 = {'s_ground':1 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                      's_excited':1 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                      'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
                      'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                      'E_ground': np.array([238213])  ,#Lowest J first, in cm^-1
                      'E_excited':  np.array([259705.55,259711.22,259724.30]) , #Lowest J first, in cm^-1
                      'Bmag': 2.24, #Magnetic Field  [T]
                      'b_angle': Los_Angle , #Angle between LoS and Bmax
                      'plot_window' : [462,470], #Min and max for plot window range (nm)
                      'amu' : 12.011 ,  #Weight in AMU
                      'specstep' : 0.035 , #Resolution of spectrometer in nm
                      'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'.
                      'Temp' : Temperature*11604, #Temperature in K. Used for Gaussian convolution
                      'plottitle': f'CIII [1s^2]2s3s (^3S_1) -> [1s^2]2s3p (^3P_J): B=2.24T , LoS = {Los_Angle}', #Title for plotting (optional)     
                      'ion_vel' : -2500, #Ion velocity in m/s
                      'fxnwindow': 1 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                      # 'Pol_angle': Polangle ,#Angle polarizing filter makes with max linear transmission, Optional

                      
                   }


InputdeckCIII_4663 = {'s_ground':1 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                      's_excited':1 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                      'l_ground': 1 , #Orbital Angular Momentum of ground state, int or half int
                      'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                      'E_ground': np.array([308216.58,308248.91,308317.29]) ,#Lowest J first, in cm^-1
                      'E_excited':  np.array([329685.38,329706.47,329743.57]) , #Lowest J first, in cm^-1
                      'Bmag': 2.24 , #Magnetic Field  [T]
                      'b_angle': Los_Angle , #Angle between LoS and Bmax
                      'plottitle': f'CIII : B=2.24T , LoS = {Los_Angle}' , #Title for plotting (optional)     
                      'plot_window' : [462,470], #Min and max for plot window range (nm)
                      'amu' : 12.011 ,  #Weight in AMU
                      'specstep' : 0.035 , #Resolution of spectrometer in nm
                      'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'.
                      'Temp' : Temperature*11604, #Temperature in K. Used for Gaussian convolution
                      'fxnwindow': 1 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                      'ion_vel' : -2500, #Ion velocity in m/s
                      # 'Pol_angle': Polangle ,#Angle polarizing filter makes with max linear transmission, Optional

                   }





CIII_4649 = Zeeman_Main(InputdeckCIII_4649)

CIII_4663 = Zeeman_Main(InputdeckCIII_4663)
    
CIII_Plotvars = [['green', 3, '-.',r'Voigt + Instrum $\Delta\lambda = 0.035 nm$']  ,                   
                 ['brown', 2.5, '--',InputdeckCIII_4663['plottitle']]]

CIII_scalevars = [1,1]



CIII_Plotvars.append(['mediumblue', 1.5, ':','GuassianInstrum Broad'])
CIII_Plotvars.append(['crimson', 1.5, '--','Gaussian Broad (no shift)'])



#This plots my results and sticks, gaussian convolution, including the ion velocity.

PF = PlotFunction(CIII_4649,CIII_Plotvars[0],plotwind = [464.5,465.5],
             SpectrumPlot =[CIII_wav, Normalize(CIIISpec[1])],
             plottitle=r"CIII $[1s^2]2s3s (^3S_1) -> [1s^2]2s3p (^3P_J)$" f": B=2.24T , LoS = {Los_Angle}, T = {Temperature}, " ,
             NormalizeSig = True, NormalizeScale = CIII_scalevars  )

# plt.vlines(CIII_4649['wave_air'], np.zeros_like(CIII_4649['signal']),CIII_4649['signal'],linewidth = 1, alpha = 0.65,color='green')
  
dopolfilterplot=True
if dopolfilterplot:
    pilist = np.array([x[1] for x in CIII_4649['rad']])==0
    sigmaplus = np.array([x[1] for x in CIII_4649['rad']])==1
    sigmaminus = np.array([x[1] for x in CIII_4649['rad']])==-1

    pispec = [CIII_4649['wave_air'][pilist],CIII_4649['signal'][pilist]]
    sigmaplusspec = [CIII_4649['wave_air'][sigmaplus],CIII_4649['signal'][sigmaplus]]
    sigmaminusspec = [CIII_4649['wave_air'][sigmaminus],CIII_4649['signal'][sigmaminus]]
    
    plt.vlines(pispec[0], np.zeros_like(pispec[1]),pispec[1],linewidth = 2, alpha = 0.95,color='orange' , label='Pi component')
    plt.vlines(sigmaplusspec[0], np.zeros_like(sigmaplusspec[1]),sigmaplusspec[1],linewidth = 2, alpha = 0.95,color='purple' , label=r'$\sigma^+$ component')
    plt.vlines(sigmaminusspec[0], np.zeros_like(sigmaminusspec[1]),sigmaminusspec[1],linewidth = 2, alpha = 0.95,color='blue' , label=r'$\sigma^-$ component')

plt.legend(loc='best')
# handles, labels = taxs[0][0].get_legend_handles_labels()
# unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
# # taxs[0][0].legend(*zip(*unique),bbox_to_anchor=(3.5, 0.5), loc="outside center left")
# # plt.tight_layout()
# fig.legend(*zip(*unique),bbox_to_anchor=(.95,.5), loc="outside center right")


if savefigs:
    plt.savefig('CIII_4649.png',dpi=1200)

#%%

# plt.vlines(sigmaplusspec[0], np.zeros_like(sigmaplusspec[1]),sigmaplusspec[1],linewidth = 1, alpha = 0.65,color='green' , label='Pi component')
# plt.vlines(sigmaminusspec[0], np.zeros_like(sigmaminusspec[1]),sigmaminusspec[1],linewidth = 1, alpha = 0.65,color='green' , label='Pi component')

# plt.vlines(CIII_4649['wave_air'], np.zeros_like(CIII_4649['signal']),CIII_4649['signal'],linewidth = 1, alpha = 0.65,color='green')
# plt.vlines(CIII_4649['wave_air'], np.zeros_like(CIII_4649['signal']),CIII_4649['signal'],linewidth = 1, alpha = 0.65,color='green')


    

    # CIII_4649[CIII_4649==0] 
#%%ComboSpec
#This adds multiple different, overlapping spectra from different atoms, which will have different scalings needed.

if doCombSpec:
    thisTemp =28.1
    
   
    InputdeckOII =        {'s_ground':1.5 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                          's_excited':1.5 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                          'l_ground': 1 , #Orbital Angular Momentum of ground state, int or half int
                          'l_excited': 2 ,#Orbital Angular Momentum of excited state, int or half int
                          'E_ground': np.array([185235.281,185340.577,185499.124]) ,#Lowest J first, in cm^-1
                          'E_excited':  np.array([206730.762,206786.286,206877.865,207002.482]) , #Lowest J first, in cm^-1
                          'Bmag': 2.24 , #Magnetic Field  [T]
                          'b_angle': Los_Angle , #Angle between LoS and Bmax
                          'plottitle': f'OII [](^3P)3s (^4P_J) -> [](^3P)3p (^4D_J) : B=2.24T , LoS = {Los_Angle}' , #Title for plotting (optional)     
                          'plot_window' : [462,470], #Min and max for plot window range (nm)
                          'amu' : 16 ,  #Weight in AMU
                          'specstep' : 0.035 , #Resolution of spectrometer in nm
                          'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'.
                          'Temp' : thisTemp*11604, #Temperature in K. Used for Gaussian convolution
                          'fxnwindow': 1 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                          'ion_vel' : -2500, #Ion velocity in m/s
                          
                      }
    InputdeckOIII =        {'s_ground':1 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                          's_excited':1 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                          'l_ground': 1 , #Orbital Angular Momentum of ground state, int or half int
                          'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                          'E_ground': np.array([370329.18,370418.32,370526.49]) ,#Lowest J first, in cm^-1
                          'E_excited':  np.array([391830.76,391917.8,392209.53]) , #Lowest J first, in cm^-1
                          'Bmag': 2.24 , #Magnetic Field  [T]
                          'b_angle': Los_Angle , #Angle between LoS and Bmax
                          'plottitle': f'OIII 2s^2 2p (3P)4p -> 2s^22p(^2P)5s (^3P): B=2.24T , LoS = {Los_Angle}' , #Title for plotting (optional)     
                          'plot_window' : [462,470], #Min and max for plot window range (nm)
                          'amu' : 16 ,  #Weight in AMU
                          'specstep' : 0.035 , #Resolution of spectrometer in nm
                          'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'.
                          'Temp' : thisTemp*11604, #Temperature in K. Used for Gaussian convolution
                          'fxnwindow': 1 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                          'ion_vel' : -2500, #Ion velocity in m/s
                          
                      }
    

    OII_465 = Zeeman_Main(InputdeckOII)
    OIII_465 = Zeeman_Main(InputdeckOIII)

    
    
    
    #Trying a rough attempt at adding two different lines together before convolution (no padspec this time)
    
    # plt.close('all')

    
    multispecs = [CIII_4649,CIII_4663,OII_465,OIII_465]
    #This allows for scaling of the differeint possible contributions.
    scalers = [1,0.1,.02,0.02, 0]
    # scalers = [0,0,1,1]
    combspect = MultiSpec(multispecs,scalers)
    #Unfortunately, all included atoms will have the same temp, velocity, and weight currently.
    Combspec = Convol_Spect(combspect[0],combspect[1], [462,470], 12.011, 0.035,Temperature_in=11602*thisTemp,
                              functiontype = "GaussianInstrum", wind_size = 1, ionvel=-2500)
    


    CIII_Plotvars = [['green', 2.5, '-.','CIII 464.9']  ,  ['brown', 2.5, '--',InputdeckCIII_4663['plottitle']]]

    Normcombsig = Normalize(Combspec[0][1])

    plt.figure()
    plt.title( f'CIII  B=2.24T , LoS = {Los_Angle} 464nm line')
    plt.suptitle( r'$CIII: [1s^2]2p3s (^3P_J) -> [1s^2]2p3p (^3P_J)$ and $[1s^2]2s3s (^3S_1) -> [1s^2]2s3p (^3P_J)$')
    PlotFunction(CIII_4649,CIII_Plotvars[0],plotwind = [464.5,465.5], SpectrumPlot =[CIII_wav, Normalize(CIIISpec[1])],plottitle=f"B = 2.24T, Temp = {Temperature}eV, Los = {Los_Angle} deg," r"$v_i$ = -2500 m/s, $\Delta \lambda_{instrum}$ = 0.035 nm" , NormalizeSig = True, NormalizeScale = CIII_scalevars  ,makefig=False)
    
    plt.plot(CIII_wav, Normalize(CIIISpec[1]), color='black', label='Spectrometer Data')
    plt.plot(Combspec[0][0], Normcombsig, label='Combined CIII + OII lines')
    plt.twinx()
    plt.vlines(CIII_4663['wave_air'], np.zeros_like(CIII_4663['signal']),scalers[1]*CIII_4663['signal'],color='orange', label='CIII 4663, f=0.1')

    plt.vlines(OII_465['wave_air'], np.zeros_like(OII_465['signal']),scalers[2]*OII_465['signal'],color='maroon', label='OII, f=0.025 ')
    plt.vlines(OIII_465['wave_air'], np.zeros_like(OIII_465['signal']),scalers[3]*OIII_465['signal'],color='blue', label='OIII, f= 0.02')
    plt.vlines(CIII_4649['wave_air'], np.zeros_like(CIII_4649['signal']),CIII_4649['signal'],color='green', label='CIII 4649, f=1 ')

    # plt.plot(CIII_Dfitwav, Dfitnorm, color='magenta', label='DoroFit')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.xlim(464.5,465.5)
    plt.show()
    if savefigs:
        plt.savefig('CIII_Multispec.png',dpi=1200)

# if dopolfilterplot:
    

