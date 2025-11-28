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
from AtomSpect  import Zeeman_Main, Normalize,PlotFunction,Vac_to_air, Convol_Spect, MultiSpec

pathfil = os.getcwd()

#Add the path for the npy spectra objects

sys.path.append(pathfil)


# CIIISpec = np.load('CIII_2.21T_7.9ev.npy') #This is the spectrometer data

savedpi = 144
savefigs =1
doHighTemp = 1
doCombSpec = 0
doConvComp = 10
doLowTemp =0
do_title = 0
do_sticks=1
#Convert the spectrometer data from vacuum to air wavelengths

#7.9eV can match with LoS - 50, T=7.2eV, res = 0.0075 nm, 2.21T 

#7.9eV can match with LoS - 50.5, T=7.91eV, res = 0.0075 nm, 2.26T 
# Temperature = 7.9#Temp in eV


#%%28.1eV
if doHighTemp:
    plt.close('all')
    CIIISpec = np.load('CIII_2.24T.npy') #This is the spectrometer data
    CIII_wav = [Vac_to_air(x) for x in CIIISpec[0]]
    Los_Angle = 68#Theta angle of LoS in degrees
    Temperature = 28.1#Temp in eV
    # Polangle = 55.75
    Bmag = 2.24
    #Spetrometer resolution is defined in Gradic2022
    InputdeckCIII_4649 = {'s_ground':1 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                          's_excited':1 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                          'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
                          'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                          'E_ground': np.array([238213])  ,#Lowest J first, in cm^-1
                          'E_excited':  np.array([259705.55,259711.22,259724.30]) , #Lowest J first, in cm^-1
                          'Bmag': Bmag, #Magnetic Field  [T]
                          'b_angle': Los_Angle , #Angle between LoS and Bmax
                          'plot_window' : [464.5,465.3], #Min and max for plot window range (nm)
                          'spec_window' : [462,470], #Min and max convolution (nm)
                          'amu' : 12.011 ,  #Weight in AMU
                          'specstep' : 0.04 , #Resolution of spectrometer in nm
                          'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'.
                          'Temp' : Temperature*11604, #Temperature in K. Used for Gaussian convolution
                          'plottitle': r"CIII $[1s^2]2s3s (^3S_1) -> [1s^2]2s3p (^3P_J)$" f": B={Bmag}T , LoS = {Los_Angle}, T = {Temperature}, ", #Title for plotting (optional)     
                          'ion_vel' : -2500, #Ion velocity in m/s
                          'fxnwindow': 2 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                          # 'Pol_angle': Polangle ,#Angle polarizing filter makes with max linear transmission, Optional
    
                          
                       }
    InputdeckCIII_4649B = {'s_ground':1 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                          's_excited':1 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                          'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
                          'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                          'E_ground': np.array([238213])  ,#Lowest J first, in cm^-1
                          'E_excited':  np.array([259705.55,259711.22,259724.30]) , #Lowest J first, in cm^-1
                          'Bmag': Bmag, #Magnetic Field  [T]
                          'b_angle': Los_Angle , #Angle between LoS and Bmax
                          'plot_window' : [464.5,465.3], #Min and max for plot window range (nm)
                          'spec_window' : [462,470], #Min and max convolution (nm)
                          'amu' : 12.011 ,  #Weight in AMU
                          'specstep' : 0.04 , #Resolution of spectrometer in nm
                          'Convfxn': 'Gaussian', #Optional: 'Gaussian', 'Skewed'.
                          'Temp' : Temperature*11604, #Temperature in K. Used for Gaussian convolution
                          'plottitle': f'CIII [1s^2]2s3s (^3S_1) -> [1s^2]2s3p (^3P_J): B={Bmag}T LoS = {Los_Angle}', #Title for plotting (optional)     
                          # 'ion_vel' : -1700, #Ion velocity in m/s
                          'fxnwindow': 2 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                          # 'Pol_angle': Polangle ,#Angle polarizing filter makes with max linear transmission, Optional
    
                          
                       }
    
    
    
    
    CIII_4649 = Zeeman_Main(InputdeckCIII_4649) #Gaussian Instrum
    CIII_4649B = Zeeman_Main(InputdeckCIII_4649B) #Gaussian Only
    
        

    
    CIII_Plotvars = [['royalblue', 2, '--','Calculated' ,'D' ]    ]
    
    
    CIII_Plotvars.append(['crimson', 2, '--','Gaussian','o'])
    
    
    
    #This plots my results and sticks, gaussian convolution, including the ion velocity.

    thixas = PlotFunction(CIII_4649,CIII_Plotvars[0],plotwind = [464.5,465.3],
                 SpectrumPlot =[CIII_wav, Normalize(CIIISpec[1])],
              #   plottitle=r"CIII $[1s^2]2s3s (^3S_1) -> [1s^2]2s3p (^3P_J)$" f": B=2.24T , LoS = {Los_Angle}, T = {Temperature}, " ,
                 NormalizeSig = True,makefig=True ,fig_size=(10,8))

    if savefigs:
        plt.savefig('CIII_4649.pdf',dpi=savedpi, format='pdf')
        if do_title:      

            plt.suptitle(f'{InputdeckCIII_4649['plottitle']}')
            plt.savefig('CIII_4649B.pdf', dpi=savedpi, format='pdf')

    if do_sticks:
        # plt.close('all')
        thixas = PlotFunction([[0],[0]],CIII_Plotvars[0],plotwind = [464.5,465.3],
                     SpectrumPlot =[CIII_wav, Normalize(CIIISpec[1])],
                     makefig=True ,plotpol = 0)
        thixas[0][0].vlines(CIII_4649['wave_air'],np.zeros_like(CIII_4649['wave_air']), CIII_4649['signal'])
        tfig = thixas[0][0].get_figure()
        tfig.legends = []
        if savefigs:
            plt.savefig('CIII_4649Sticks.pdf',dpi=savedpi, format='pdf')



    if doConvComp:
        
        thixas = PlotFunction(CIII_4649,CIII_Plotvars[0],plotwind = [464.5,465.3],
                     SpectrumPlot =[CIII_wav, Normalize(CIIISpec[1])],
                  #   plottitle=r"CIII $[1s^2]2s3s (^3S_1) -> [1s^2]2s3p (^3P_J)$" f": B=2.24T , LoS = {Los_Angle}, T = {Temperature}, " ,
                     NormalizeSig = True,makefig=True ,plotpol = 0)
        PlotFunction(CIII_4649B,CIII_Plotvars[1],plotwind = [464.5,465.5],  axsin=thixas , 
                     NormalizeSig = True,makefig=0 ,plotpol = 0)

        if savefigs:
            plt.savefig('CIII_4649ConvComp.pdf',dpi=savedpi)
            if do_title:      
                # plt.tight_layout(rect=[0, 0, 0.75, .95])
                plt.suptitle(f'{InputdeckCIII_4649['plottitle']}')
                plt.savefig('CIII_4649ConvCompB.pdf', dpi=savedpi)
    
  

    

    #This adds multiple different, overlapping spectra from different atoms, which will have different scalings needed.
    
    if doCombSpec:
        plt.close('all')
        thisTemp =28.1
    
       
        InputdeckOII =        {'s_ground':1.5 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                              's_excited':1.5 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                              'l_ground': 1 , #Orbital Angular Momentum of ground state, int or half int
                              'l_excited': 2 ,#Orbital Angular Momentum of excited state, int or half int
                              'E_ground': np.array([185235.281,185340.577,185499.124]) ,#Lowest J first, in cm^-1
                              'E_excited':  np.array([206730.762,206786.286,206877.865,207002.482]) , #Lowest J first, in cm^-1
                              'Bmag': 2.22, #Magnetic Field  [T]
                              'b_angle': Los_Angle , #Angle between LoS and Bmax
                              'plottitle': f'OII [](^3P)3s (^4P_J) -> [](^3P)3p (^4D_J) : B=2.24T , LoS = {Los_Angle}' , #Title for plotting (optional)     
                              'plot_window' : [464.5,465.5], #Min and max for plot window range (nm)
                              'spec_window' : [462,470], #Min and max convolution (nm)
                              'amu' : 16 ,  #Weight in AMU
                              'specstep' : 0.04 , #Resolution of spectrometer in nm
                              'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'.
                              'Temp' : thisTemp*11604, #Temperature in K. Used for Gaussian convolution
                              'fxnwindow': 2 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                              # 'ion_vel' : -1700, #Ion velocity in m/s
                              
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
                              'plot_window' : [464.5,465.5], #Min and max for plot window range (nm)
                              'spec_window' : [462,470], #Min and max convolution (nm)
                              'amu' : 16 ,  #Weight in AMU
                              'specstep' : 0.02 , #Resolution of spectrometer in nm
                              'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'.
                              'Temp' : thisTemp*11604, #Temperature in K. Used for Gaussian convolution
                              'fxnwindow': 2 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                              # 'ion_vel' : -1700, #Ion velocity in m/s
                              
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
                          'plot_window' : [464.5,465.5], #Min and max for plot window range (nm)
                          'spec_window' : [462,470], #Min and max convolution (nm)
    
                          'amu' : 12.011 ,  #Weight in AMU
                          'specstep' : 0.02 , #Resolution of spectrometer in nm
                          'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'.
                          'Temp' : Temperature*11604, #Temperature in K. Used for Gaussian convolution
                          'fxnwindow': 1 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                          'ion_vel' : -2500, #Ion velocity in m/s
                          # 'Pol_angle': 90 ,#Angle polarizing filter makes with max linear transmission, Optional
    
                       }
    
    
        OII_465 = Zeeman_Main(InputdeckOII)
        OIII_465 = Zeeman_Main(InputdeckOIII)
    
        CIII_4663 = Zeeman_Main(InputdeckCIII_4663)
    
        
        
        #Trying a rough attempt at adding two different lines together before convolution (no padspec this time)
        
        # plt.close('all')
    
        
        multispecs = [CIII_4649,CIII_4663,OII_465,OIII_465]
        #This allows for scaling of the differeint possible contributions.
        scalers = [1,.1,.02,0.02, 0]
        # scalers = [0,0,1,1]
        combspect = MultiSpec(multispecs,scalers)
        #Unfortunately, all included atoms will have the same temp, velocity, and mass currently.
        Combspec = Convol_Spect(combspect[0],combspect[1], [462,470], 12.011, 0.035,Temperature_in=11602*thisTemp,
                                  functiontype = "GaussianInstrum", wind_size = 1, ionvel=-2500)
        
    
    
        CIII_Plotvars = [['green', 2.5, '-.','CIII 464.9']  ,  ['brown', 2.5, '--',InputdeckCIII_4663['plottitle']]]
    
        Normcombsig = Normalize(Combspec[0][1])
    
        fig, axs = plt.subplots(1 , figsize=(14,14)) 
        plt.xlabel('Wavelength (nm)', weight='semibold')
        plt.ylabel('Normalized Intensity (arb)', weight='semibold')

        axslines = axs.twinx()
        # axslines.vlines(CIII_4649['wave_air'], np.zeros_like(CIII_4649['signal']),CIII_4649['signal'],color='green', label='CIII 4649')
        # axslines.vlines(CIII_4663['wave_air'], np.zeros_like(CIII_4663['signal']),scalers[1]*CIII_4663['signal'],color='orange', label='CIII 4663, f=0.1')
    
        # axslines.vlines(OII_465['wave_air'], np.zeros_like(OII_465['signal']),scalers[2]*OII_465['signal'],color='maroon', label='OII, f=0.025 ')
        # axslines.vlines(OIII_465['wave_air'], np.zeros_like(OIII_465['signal']),scalers[3]*OIII_465['signal'],color='blue', label='OIII, f= 0.02')
        axslines.set_ylim(-0.1*np.max(CIII_4649['reduced_sticks'][1]),np.max(CIII_4649['reduced_sticks'][1]))
        axslines.stem(CIII_4649['wave_air'], CIII_4649['signal'], linefmt='C4', markerfmt='+', basefmt=" ",label='CIII 4649')
        axslines.stem(CIII_4663['wave_air'],scalers[1]*CIII_4663['signal'],linefmt='C1', markerfmt='D', basefmt=" ", label='CIII 4663, f=0.1')
    
        axslines.stem(OII_465['wave_air'], scalers[2]*OII_465['signal'],linefmt='C2', markerfmt='x', basefmt=" ", label='OII, f=0.025 ')
        axslines.stem(OIII_465['wave_air'], scalers[3]*OIII_465['signal'], linefmt='C3', markerfmt='o', basefmt=" ",label='OIII, f= 0.02')
        axslines.set_ylim(-0.1*np.max(CIII_4649['reduced_sticks'][1]),np.max(CIII_4649['reduced_sticks'][1]))
    
        axs.plot(Combspec[0][0], Normcombsig, label='Combined CIII + OII lines', marker='D',linewidth = 2, markersize = 5, markevery=2)
        combax = PlotFunction(CIII_4649,CIII_Plotvars[0],plotwind = [464.5,465.5], SpectrumPlot =[CIII_wav, Normalize(CIIISpec[1])],plottitle=f"B = 2.24T, Temp = {Temperature}eV, Los = {Los_Angle} deg," r"$v_i$ = -2500 m/s, $\Delta \lambda_{instrum}$ = 0.035 nm" ,
                     NormalizeSig=True, makefig=0,plotpol=0,axsin = axs)
        tfig=combax.get_figure()
        # tfig=combax[0][1].get_figure()
        tfig.legends = []

        #This block is to resulve some issues with legends. Normally, it is included in the PlotFunction, but as that was not used here,
        #It must be done like this.
        #Future work aims to find a better way to combine multiple spectra
        try:
            handles, labels = [x[0].get_legend_handles_labels() for x in combax]
            # pass
        except:
            try:
                
                handles, labels = combax[0][0].get_legend_handles_labels()
            except:
                handles, labels = combax.get_legend_handles_labels()
    
        handles2,labels2 = axslines.get_legend_handles_labels()
      
        handles.extend(handles2)
        labels.extend(labels2)
        # plt.rcParams['figure.figsize'] = [10,10]
 
        unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
        tfig.legend(*zip(*unique), bbox_to_anchor=(.5, .95), loc='outside center', ncol=4, fontsize=18)
        plt.tight_layout(rect=[0, 0, 1, 0.92])
        plt.subplots_adjust(hspace=0, wspace=0)
        if savefigs:
            plt.savefig('CIII_Multispec.pdf',format='pdf')
            if do_title:      
                # plt.tight_layout(rect=[0, 0, 0.75, .95])
                plt.suptitle(f'{InputdeckCIII_4649['plottitle']}')
                plt.savefig('CIII_MultispecB.pdf', dpi=savedpi)
    


#%%7.9eV
if doLowTemp:
    CIIISpec = np.load('CIII_2.21T_7.9ev.npy') #This is the spectrometer data
    CIII_wav = [Vac_to_air(x) for x in CIIISpec[0]]
    Los_Angle = 50.5#Theta angle of LoS in degrees
    Temperature = 7.01#Temp in eV
    # Polangle = 55.75
    Bmag = 2.26
    InputdeckCIII_4649L = {'s_ground':1 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                       's_excited':1 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                       'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
                       'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                       'E_ground': np.array([238213])  ,#Lowest J first, in cm^-1
                       'E_excited':  np.array([259705.55,259711.22,259724.30]) , #Lowest J first, in cm^-1
                       'Bmag': Bmag, #Magnetic Field  [T]
                       'b_angle': Los_Angle , #Angle between LoS and Bmax
                       'plot_window' : [464.5,465.5], #Min and max for plot window range (nm)
                       'spec_window' : [462,470], #Min and max convolution (nm)
                       'amu' : 12.011 ,  #Weight in AMU
                       'specstep' : 0.01, #Resolution of spectrometer in nm
                       'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'.
                       'Temp' : Temperature*11604, #Temperature in K. Used for Gaussian convolution
                       'plottitle': r"CIII $[1s^2]2s3s (^3S_1) -> [1s^2]2s3p (^3P_J)$" f": B={Bmag}T , LoS = {Los_Angle}, T = {Temperature}, ", #Title for plotting (optional)     
                       # 'specres' : 20,
                       'ion_vel' : -1700, #Ion velocity in m/s
                       'fxnwindow': 2 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                       # 'Pol_angle': Polangle ,#Angle polarizing filter makes with max linear transmission, Optional
 
                       
                    }
    CIII_4649L = Zeeman_Main(InputdeckCIII_4649L)

    CIII_Plotvars = ['royalblue', 2, '-.','Calculated' , 'D' ]
    thixas = PlotFunction(CIII_4649L,CIII_Plotvars,plotwind = [464.5,465.5],
                 SpectrumPlot =[CIII_wav, Normalize(CIIISpec[1])],
               # plottitle=r"CIII $[1s^2]2s3s (^3S_1) -> [1s^2]2s3p (^3P_J)$" f": B=2.24T , LoS = {Los_Angle}, T = {Temperature}, " ,
                 NormalizeSig = True,makefig=True )
    if savefigs:
        # plt.tight_layout()
        plt.savefig('CIII_4649LowTemp.pdf',dpi=savedpi)
        if do_title:      
            plt.tight_layout(rect=[0, 0, 1, .95])
            plt.suptitle(f'{InputdeckCIII_4649L['plottitle']}')
            plt.savefig('CIII_4649LowTempB.pdf', dpi=savedpi)
