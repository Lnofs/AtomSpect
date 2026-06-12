# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 13:26:58 2024
@author: Leo Nofs


This spectra was collected by Tommy Gonda, Matt Kriete, and David Ennis using the magnetic fields
of the Magnetized Dusty Plasma Experiment (MDPX) at Auburn University, in Auburn, AL

The Argon spectra was collected after reflection off a beam-splitting plate. This plate has a 
wavelength AND polarization mode dependent reflectivity coefficient. This modifies the relative
intensity ratio between linear (pi) and circular (sigma) components of the spectra.

For these Argon lines, the plate data sheet sees approximate values of (for 750 nm light)
TAlpha = np.sqrt(0.65) - For the pi component
TBeta = np.sqrt(0.35) - For the sigma component

This example showcases the need of the intermediate regime solution by looking at the 1s2p to 1s4d (447 nm) line for He I
with three magnetic field strengths (1.5T, 2.5T, 3.2T). This example also includes the inclusion of a linearly polarizing
filter in the spectrometer line of sight and the impact of the filter oreintation makes on the observed spectra.
    
The standard zeeman fan plots are additionally included here.
    


"""

import os
import sys
import numpy as np
from matplotlib import pyplot as plt


sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir, os.path.pardir)))

from AtomSpect import AtomSpect_Main, Normalize, PlotFunction, read_Spectra,  MakeSlider,  plotZfan


# Flags for doing polarized or not
doHePol = 1  # Shows all pols - inability to remove all of one pol
doHePol2 = 1 # Shows non-dipole component and inability to fully filter polarizations (Paschen Back)
doHePol3 = 1 # Same as hePol, but show the difference the film makes

doHeNonPol = 1# Shows multiple B fields with all three regimes
doHeNonPol2 = 1  # All 3 B fields, show sticks (dipole only)
doHeNonPol3 =  1# Single B - show non-dipole required to match

# savepdf = 0

savepdf = 1
savepng = 0
savedpi = 144
fanplot = 0 #Ground State Fan
fanplot2 = 0 #Excited State Fan
do_HeSlider = 0
do_Hepolslider = 0

closeplots = 0
# do_title=0 #Adds titles to the plots

if doHePol2 == True:

    plt.close('all')
    pathfil = os.getcwd()
    filname = '/He_Spectra_1.5T_Polarized.csv'
    fil = pathfil + filname

    He_Pol_Spec = read_Spectra(fil, headercount=2)

    # This corrects for a misalignment in the calibration by 3 pixels.
    He_PolShift = [0.1*x + 0.006666 for x in He_Pol_Spec['wavelength']]

    He_Pol_Spec = [He_Pol_Spec['signals'][0], He_Pol_Spec['signals'][-1]]

    polangles = [0, 90]

    colorlist = ['royalblue', 'g', 'r', 'brown', 'm', 'c', 'y', 'black', 'seagreen']

    HePol_PlotObj = []
    HePolCol = []
    for i, angles in enumerate(polangles):

        InputdeckHePol = {'s_ground': 1,  # Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                          's_excited': 1,  # Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                          'l_ground': 1,  # Orbital Angular Momentum of ground state, int or half int
                          'l_excited': 2,  # Orbital Angular Momentum of excited state, int or half int
                          'E_ground': np.array([169087.8308131,  169086.8428979,  169086.7664725]),  # Lowest J first, in cm^-1
                          'E_excited': np.array([191444.49804915, 191444.47952984, 191444.47832914]),  # Lowest J first, in cm^-1
                          'Bmag': 1.5,  # Magnetic Field  [T]
                          'b_angle': 90,  # Angle between LoS and Bmax
                          'Convfxn': 'Skewed',  # Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
                          'specstep': 0.002222,  # Resolution of spectrometer in nm
                          'Skewness': [0.62, 0.99],  # How non-symmetric the Lorrentzian is
                          'Pol_angle': angles,  # Angle polarizing filter makes with max linear transmission, Optional
                          'plottitle': 'He I 1s2p -> 1s4d: B=1.5T, Polarizing Filter',  # Title for plotting (optional)
                          'spec_window': [447.1, 447.2],  # Min and max for convolutions
                          'SpectrumData': [He_PolShift, He_Pol_Spec[i]],  # Spectrometer Data.
                          'plot_window': [447.1, 447.2],  # Min and max for plot window range (nm), convolutions will only take place within this range
                          'amu': 2,  # Weight in AMU
                          'specres': 10,  # How many steps per resolution are calculated. Higher makes a smoother curve.
                            'TAlpha' : np.sqrt(1-0.35) ,
                            'TBeta' : np.sqrt(1-0.65) ,
                          }

        HePol = AtomSpect_Main(InputdeckHePol)
        HePol_PlotObj = ['royalblue', 2.0, 'solid', 'Intermediate']

        if i == 0:
            taxs = PlotFunction(HePol, HePol_PlotObj, NormalizeSig=True,
                                Shape=[1,len(polangles)], makefig=True,
                                plotlabel=r'$\gamma$' + f'={angles}' + r'$^o$', plotnondip=True , fig_size = (19.5,9))
        else:
            PlotFunction(HePol, HePol_PlotObj,  NormalizeSig=True, axsin=taxs, position=[0, i],
                         makefig=False,  plotlabel=r'$\gamma$' + f'={angles}' + r'$^o$', plotnondip=True, legcols=6)

    # if do_title:      
    plt.tight_layout(rect=[0, 0, 1, .9])
    taxs[0][1].annotate('', xy=(447.174, .33), xytext=(447.174, 0.65), 
            arrowprops=dict(facecolor='black', shrink=0.05),
            )
    taxs[0][0].annotate('', xy=(447.153, .475), xytext=(447.1605, 0.65), 
            arrowprops=dict(facecolor='black', shrink=0.05),
            )
    

    if savepdf:

        plt.savefig('He_Pol2.pdf', dpi=savedpi)
    if savepng:    
            plt.tight_layout(rect=[0, 0, 1, .9])
            plt.rcParams['axes.titlepad'] = 315
            plt.title(f'{InputdeckHePol['plottitle']}')
            

            plt.savefig('He_Pol2.png', dpi=savedpi)

        
    plt.show()
# %%

if doHePol == True:
    plt.close('all')
    pathfil = os.getcwd()
    filname = '/He_Spectra_1.5T_Polarized.csv'
    fil = pathfil + filname

    He_Pol_Spec = read_Spectra(fil, headercount=2)

    # This corrects for a misalignment in the calibration by 3 pixels.
    He_PolShift = [0.1*x + 0.006666 for x in He_Pol_Spec['wavelength']]

    polangles = [0, 30, 70, 90]

    # fig, axs = plt.subplots(4,sharex=True,sharey=True, figsize = (16,10))
    colorlist = ['royalblue', 'g', 'r', 'brown', 'm', 'c', 'y', 'black', 'seagreen']

    HePol_PlotObj = []
    HePolCol = []
    for i, angles in enumerate(polangles):

        InputdeckHePol = {'s_ground': 1,  # Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                          's_excited': 1,  # Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                          'l_ground': 1,  # Orbital Angular Momentum of ground state, int or half int
                          'l_excited': 2,  # Orbital Angular Momentum of excited state, int or half int
                          'E_ground': np.array([169087.8308131,  169086.8428979,  169086.7664725]),  # Lowest J first, in cm^-1
                          'E_excited': np.array([191444.49804915, 191444.47952984, 191444.47832914]),  # Lowest J first, in cm^-1
                          'Bmag': 1.5,  # Magnetic Field  [T]
                          'b_angle': 90,  # Angle between LoS and Bmax
                          'Convfxn': 'Skewed',  # Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
                          'specstep': 0.002222,  # Resolution of spectrometer in nm
                          'Skewness': [0.62, 0.99],  # How non-symmetric the Lorrentzian is
                          'Pol_angle': angles,  # Angle polarizing filter makes with max linear transmission, Optional
                          'plottitle': 'He I 1s2p -> 1s4d: B=1.5T, Polarizing Filter',  # Title for plotting (optional)
                          'spec_window': [447, 447.2],  # Min and max for convolutions
                        
                          'SpectrumData': [He_PolShift, He_Pol_Spec['signals'][i]],  # Spectrometer Data.
                          'plot_window': [447.1, 447.2],  # Min and max for plot window range (nm), convolutions will only take place within this range
                          'amu': 2,  # Weight in AMU
                          'specres':3,  # How many steps per resolution are calculated. Higher makes a smoother curve.
                          # 'sortE' : True
                             'TAlpha' : np.sqrt(1-0.35) ,
                             'TBeta' : np.sqrt(1-0.65) ,

                          }

        HePol = AtomSpect_Main(InputdeckHePol)
        HePol_PlotObj = ['royalblue', 2.5, 'solid', 'Intermediate']

        if i == 0:  # In the future, automate this process of the plotting.
            taxs = PlotFunction(HePol, HePol_PlotObj, NormalizeSig=True,
                                Shape=[len(polangles), 1], makefig=True,
                                plotlabel=r'$\gamma$' + f'={angles}' + r'$^o$',marker_count=5)
        else:
            PlotFunction(HePol, HePol_PlotObj,  NormalizeSig=True, axsin=taxs[i][0], position=[0, i],
                         makefig=False,  plotlabel=r'$\gamma$' + f'={angles}' + r'$^o$',marker_count=5)
            
    # tfig = taxs[0][0].get_figure()   
    # plt.subplots_adjust(hspace=0, wspace=0)

    # savepdf = 1
    if savepdf:
        plt.savefig('He_Pol.pdf', dpi=savedpi)
    if savepng:      
        plt.tight_layout(rect=[0, 0, 1, .9])
        plt.rcParams['axes.titlepad'] = 450
        plt.title(f'{InputdeckHePol['plottitle']}')
        plt.savefig('He_PolB.png', dpi=savedpi)
        
        #%%

if doHePol3 == True:
    # plt.close('all')
    pathfil = os.getcwd()
    filname = '/He_Spectra_1.5T_Polarized.csv'
    fil = pathfil + filname

    He_Pol_Spec = read_Spectra(fil, headercount=2)

    # This corrects for a misalignment in the calibration by 3 pixels.
    He_PolShift = [0.1*x + 0.006666 for x in He_Pol_Spec['wavelength']]

    polangles = [0, 30, 70, 90]

    # fig, axs = plt.subplots(4,sharex=True,sharey=True, figsize = (16,10))
    colorlist = ['royalblue', 'g', 'r', 'brown', 'm', 'c', 'y', 'black', 'seagreen']

    HePol_PlotObj = []
    HePolCol = []
    for i, angles in enumerate(polangles):
        InputdeckHePol3 = {'s_ground': 1,  # Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                          's_excited': 1,  # Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                          'l_ground': 1,  # Orbital Angular Momentum of ground state, int or half int
                          'l_excited': 2,  # Orbital Angular Momentum of excited state, int or half int
                          'E_ground': np.array([169087.8308131,  169086.8428979,  169086.7664725]),  # Lowest J first, in cm^-1
                          'E_excited': np.array([191444.49804915, 191444.47952984, 191444.47832914]),  # Lowest J first, in cm^-1
                          'Bmag': 1.5,  # Magnetic Field  [T]
                          'b_angle': 90,  # Angle between LoS and Bmax
                          'Convfxn': 'Skewed',  # Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
                          'specstep': 0.002222,  # Resolution of spectrometer in nm
                          'Skewness': [0.62, 0.99],  # How non-symmetric the Lorrentzian is
                          'Pol_angle': angles,  # Angle polarizing filter makes with max linear transmission, Optional
                          'plottitle': 'He I 1s2p -> 1s4d: B=1.5T, Polarizing Filter',  # Title for plotting (optional)
                          'spec_window': [447, 447.2],  # Min and max for convolutions
                        
                          'plot_window': [447.1, 447.2],  # Min and max for plot window range (nm), convolutions will only take place within this range
                          'amu': 2,  # Weight in AMU
                          'specres':3,  # How many steps per resolution are calculated. Higher makes a smoother curve.
                          # 'sortE' : True
                          'SpectrumData': [He_PolShift, He_Pol_Spec['signals'][i]],  # Spectrometer Data.


                          }

        HePol3 = AtomSpect_Main(InputdeckHePol3)
        HePol_PlotVar3 = ['royalblue', 2.5, 'solid', 'Intermediate']

    
        InputdeckHeNonPol3b = {'s_ground': 1,  # Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                             's_excited': 1,  # Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                             'l_ground': 1,  # Orbital Angular Momentum of ground state, int or half int
                             'l_excited': 2,  # Orbital Angular Momentum of excited state, int or half int
                             'E_ground': np.array([169087.8308131,  169086.8428979,  169086.7664725]),  # Lowest J first, in cm^-1
                             'E_excited': np.array([191444.49804915, 191444.47952984, 191444.47832914]),  # Lowest J first, in cm^-1
                             'Bmag': 1.5,  # Magnetic Field  [T]
                             'b_angle': 90,  # Angle between LoS and Bmax
                             'Convfxn': 'Skewed',  # Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
                             'specstep': 0.002222,  # Resolution of spectrometer in nm
                             'Skewness': [0.62, 0.99],  # How non-symmetric the Lorrentzian is
                             'plottitle': 'He I 1s2p -> 1s4d: B=1.5T, Polarizing Filter ' ,#    Title for plotting (optional)
                             'plot_window': [447.1, 447.2],  # Min and max for plot window range (nm), convolutions will only take place within this range
                             'spec_window': [447, 447.2],  # Min and max for convolutions
                             'Pol_angle': angles,  # Angle polarizing filter makes with max linear transmission, Optional
                          'specres':3,  # How many steps per resolution are calculated. Higher makes a smoother curve.

                             'amu': 2,  # Weight in AMU
                          'SpectrumData': [He_PolShift, He_Pol_Spec['signals'][i]],  # Spectrometer Data.

                             # 'DoLowSig' : 'Y' , #If this exists in the input, lowfield signal strength will be calculated.
                             # 'DoHighSig' : 'Y' , #If this exists in the input, highfield signal strength will be calculated.
                             'sortE': True ,
                             # 'Pol_angle': 0,  # Angle polarizing filter makes with max linear transmission, Optional
            
                            'TAlpha' : np.sqrt(0.65) , #Beam Splitter Film Correction
                            'TBeta' : np.sqrt(0.35), #Beam Splitter Film Coating Correction
                             }
        
        HePol3b = AtomSpect_Main(InputdeckHeNonPol3b)

        
        colorlist = ['royalblue', 'g', 'r', 'brown', 'm', 'c', 'y', 'black', 'seagreen']

        HePol_PlotVar3b= ['red', 2.5, 'solid', 'Int - Corrected Pol','s']






        if i == 0:  # In the future, automate this process of the plotting.
            taxs = PlotFunction(HePol3b, HePol_PlotVar3b, NormalizeSig=True,
                                Shape=[len(polangles), 1], makefig=True,
                                plotlabel=r'$\gamma$' + f'={angles}' + r'$^o$',marker_count=15,plotpol=0)
            
            PlotFunction(HePol3, HePol_PlotVar3,  NormalizeSig=True, axsin=taxs[i][0], position=[0, i],
                         makefig=False,  plotlabel=r'$\gamma$' + f'={angles}' + r'$^o$',marker_count=11,plotpol=0)
        else:
            PlotFunction(HePol3, HePol_PlotVar3,  NormalizeSig=True, axsin=taxs[i][0], position=[0, i],
                         makefig=False,  plotlabel=r'$\gamma$' + f'={angles}' + r'$^o$',marker_count=11,plotpol=0)
            PlotFunction(HePol3b, HePol_PlotVar3b,  NormalizeSig=True, axsin=taxs[i][0], position=[0, i],
                         makefig=False,  plotlabel=r'$\gamma$' + f'={angles}' + r'$^o$',marker_count=15,plotpol=0, legcols=3)
            
    # tfig = taxs[0][0].get_figure()   
    # plt.subplots_adjust(hspace=0, wspace=0)

    # savepdf = 1
    if savepdf:
        plt.savefig('He_Pol3.pdf', dpi=savedpi)
    if savepng:      
        plt.tight_layout(rect=[0, 0, 1, .9])
        plt.rcParams['axes.titlepad'] = 450
        plt.title(f'{InputdeckHePol['plottitle']}')
        plt.savefig('He_Pol3B.png', dpi=savedpi)
# %%For the paper
if doHeNonPol == True: 
    pathfil = os.getcwd()
    plt.close('all')
    fil = pathfil + '/He_Spectra_NonPolarized.csv'
    He_NonPol_Spec = read_Spectra(fil, scalewave=0.1, headercount=1)  # Need to scale the wavelengths from angstrom to nm.

    Bmags = [3.2, 2.5, 1.5]

    for i, Bvals in enumerate(Bmags):
        if i == 2:  # Need to shift the 1.5T spectra by 3 pixels.
            He_NonPol_Spec['wavelength'] = [x + 0.00666 for x in He_NonPol_Spec['wavelength']]

        InputdeckHeNonPol = {'s_ground': 1,  # Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                             's_excited': 1,  # Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                             'l_ground': 1,  # Orbital Angular Momentum of ground state, int or half int
                             'l_excited': 2,  # Orbital Angular Momentum of excited state, int or half int
                             'E_ground': np.array([169087.8308131,  169086.8428979,  169086.7664725]),  # Lowest J first, in cm^-1
                             'E_excited': np.array([191444.49804915, 191444.47952984, 191444.47832914]),  # Lowest J first, in cm^-1
                             'Bmag': Bvals,  # Magnetic Field  [T]
                             'b_angle': 90,  # Angle between LoS and Bmax
                             'Convfxn': 'Skewed',  # Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
                             'specstep': 0.002222,  # Resolution of spectrometer in nm
                             'Skewness': [0.62, 0.99],  # How non-symmetric the Lorrentzian is
                             'plottitle': 'He I 1s2p -> 1s4d: B=1.5T, Low-Int-High',  # Title for plotting (optional)
                             'plot_window': [447.1, 447.2],  # Min and max for plot window range (nm), convolutions will only take place within this range
                             'spec_window': [447, 447.2],  # Min and max for convolutions

                             'amu': 2,  # Weight in AMU
                             'SpectrumData': [He_NonPol_Spec['wavelength'], He_NonPol_Spec['signals'][i]],  # Spectrometer Data.

                             'DoLowSig' : 'Y' , #If this exists in the input, lowfield signal strength will be calculated.
                             'DoHighSig' : 'Y' , #If this exists in the input, highfield signal strength will be calculated.
                             'EtermG': 169086.7664725, #Term energy of lower state in cm^-1

                             # 'EtermG': 169086.9076, #Term energy of lower state in cm^-1
                             'EtermE' : 191444.47832914, #Term energy of upper state in cm^-1
                             # 'EtermE' : 191444.4827, #Term energy of upper state in cm^-1
                             'TAlpha' : np.sqrt(1-0.35) ,
                             'TBeta' : np.sqrt(1-0.65) ,

                             }
        HeNonPol = AtomSpect_Main(InputdeckHeNonPol)

        colorlist = ['royalblue', 'g', 'r', 'brown', 'm', 'c', 'y', 'black', 'seagreen']
        HeNonPol_PlotVar = [] 
        HeNonPol_PlotVar.append(['royalblue', 2.5, 'solid', 'Intermediate','D'])
        HeNonPol_PlotVar.append( [colorlist[1], 1.5, '-.', 'Weak-field','X'])
        HeNonPol_PlotVar.append( [colorlist[2], 1.5, '--', 'Strong-field','o'])

        if i == 0:  # In the future, automate this process of the plotting.
            taxs = PlotFunction(HeNonPol, HeNonPol_PlotVar[0], NormalizeSig=True,
                                Shape=[1,len(Bmags)], makefig=True,
                                plotlabel=f'B={Bvals}T', plotnondip=0, plotpol=0 , fig_size = (19.5,9))
            PlotFunction(HeNonPol['SpecOutLow'], HeNonPol_PlotVar[1],  NormalizeSig=True, axsin=taxs, position=[i,0],
                         makefig=False,plotwind = InputdeckHeNonPol['plot_window'],  plotnondip=0,plotpol=0)
            PlotFunction(HeNonPol['SpecOutHigh'], HeNonPol_PlotVar[2],  NormalizeSig=True, axsin=taxs, position=[i,0],
                         makefig=False,plotwind = InputdeckHeNonPol['plot_window'],  plotnondip=0,plotpol=0)
     
        else:
            PlotFunction(HeNonPol, HeNonPol_PlotVar[0],  NormalizeSig=True, axsin=taxs, position=[0, i],
                            makefig=False,  plotlabel=f'B={Bvals}T', plotnondip=0, plotpol=0)         
            PlotFunction(HeNonPol['SpecOutHigh'], HeNonPol_PlotVar[2],  NormalizeSig=True, axsin=taxs, position=[0,i],
                         makefig=False,plotwind = InputdeckHeNonPol['plot_window'],  plotnondip=0,plotpol=0)
            
            PlotFunction(HeNonPol['SpecOutLow'], HeNonPol_PlotVar[1],  NormalizeSig=True, axsin=taxs, position=[0,i],
                         makefig=False,plotwind = InputdeckHeNonPol['plot_window'],  plotnondip=0,plotpol=0 , legcols=4)
    plt.tight_layout(rect=[0, 0, 1, .9])

    # if do_title:      
    #     plt.tight_layout(rect=[0, 0, 0.75, .95])
    #     plt.suptitle(f'{InputdeckHeNonPol['plottitle']}')
    if savepdf:
        plt.savefig('He_NonPolRegime.pdf', format='pdf')
    if savepng:      
        plt.tight_layout(rect=[0, 0, 1, .9])
        plt.rcParams['axes.titlepad'] = 405
        plt.title(f'{InputdeckHeNonPol['plottitle']}')
        plt.savefig('He_NonPolRegimeB.png', dpi=savedpi)
#%%
if doHeNonPol2 == True: 
    pathfil = os.getcwd()
    plt.close('all')
    fil = pathfil + '/He_Spectra_NonPolarized.csv'
    He_NonPol_Spec = read_Spectra(fil, scalewave=0.1, headercount=1)  # Need to scale the wavelengths from angstrom to nm.

    Bmags = [3.2, 2.5, 1.5]

    for i, Bvals in enumerate(Bmags):
        if i == 2:  # Need to shift the 1.5T spectra by 3 pixels.
            He_NonPol_Spec['wavelength'] = [x + 0.00666 for x in He_NonPol_Spec['wavelength']]

        InputdeckHeNonPol2 = {'s_ground': 1,  # Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                             's_excited': 1,  # Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                             'l_ground': 1,  # Orbital Angular Momentum of ground state, int or half int
                             'l_excited': 2,  # Orbital Angular Momentum of excited state, int or half int
                             'E_ground': np.array([169087.8308131,  169086.8428979,  169086.7664725]),  # Lowest J first, in cm^-1
                             'E_excited': np.array([191444.49804915, 191444.47952984, 191444.47832914]),  # Lowest J first, in cm^-1
                             'Bmag': Bvals,  # Magnetic Field  [T]
                             'b_angle': 90,  # Angle between LoS and Bmax
                             'Convfxn': 'Skewed',  # Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
                             'specstep': 0.002222,  # Resolution of spectrometer in nm
                             'Skewness': [0.62, 0.99],  # How non-symmetric the Lorrentzian is
                             'plottitle': 'He I 1s2p -> 1s4d: B=1.5T, Polarization Tracking',  # Title for plotting (optional)
                             'plot_window': [447.1, 447.2],  # Min and max for plot window range (nm), convolutions will only take place within this range
                             'spec_window': [447, 447.2],  # Min and max for convolutions

                             'amu': 2,  # Weight in AMU
                             'SpectrumData': [He_NonPol_Spec['wavelength'], He_NonPol_Spec['signals'][i]],  # Spectrometer Data.

                             # 'DoLowSig' : 'Y' , #If this exists in the input, lowfield signal strength will be calculated.
                             # 'DoHighSig' : 'Y' , #If this exists in the input, highfield signal strength will be calculated.
                             'sortE': True ,
                             # 'TAlpha' : np.sqrt(1-0.35) ,
                             # 'TBeta' : np.sqrt(1-0.65) ,
                             }
        HeNonPol2 = AtomSpect_Main(InputdeckHeNonPol2)
        
        InputdeckHeNonPol2b = {'s_ground': 1,  # Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                             's_excited': 1,  # Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                             'l_ground': 1,  # Orbital Angular Momentum of ground state, int or half int
                             'l_excited': 2,  # Orbital Angular Momentum of excited state, int or half int
                             'E_ground': np.array([169087.8308131,  169086.8428979,  169086.7664725]),  # Lowest J first, in cm^-1
                             'E_excited': np.array([191444.49804915, 191444.47952984, 191444.47832914]),  # Lowest J first, in cm^-1
                             'Bmag': Bvals,  # Magnetic Field  [T]
                             'b_angle': 90,  # Angle between LoS and Bmax
                             'Convfxn': 'Skewed',  # Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
                             'specstep': 0.002222,  # Resolution of spectrometer in nm
                             'Skewness': [0.62, 0.99],  # How non-symmetric the Lorrentzian is
                             'plottitle': 'He I 1s2p -> 1s4d: B=1.5T, Polarization Tracking',  # Title for plotting (optional)
                             'plot_window': [447.1, 447.2],  # Min and max for plot window range (nm), convolutions will only take place within this range
                             'spec_window': [447, 447.2],  # Min and max for convolutions

                             'amu': 2,  # Weight in AMU
                             'SpectrumData': [He_NonPol_Spec['wavelength'], He_NonPol_Spec['signals'][i]],  # Spectrometer Data.

                             # 'DoLowSig' : 'Y' , #If this exists in the input, lowfield signal strength will be calculated.
                             # 'DoHighSig' : 'Y' , #If this exists in the input, highfield signal strength will be calculated.
                             'sortE': True ,
                             # 'Pol_angle': 0,  # Angle polarizing filter makes with max linear transmission, Optional
            
                            'TAlpha' : np.sqrt(0.65) , #Beam Splitter Film Correction
                            'TBeta' : np.sqrt(0.35), #Beam Splitter Film Coating Correction
                             }
        
        HeNonPol2b = AtomSpect_Main(InputdeckHeNonPol2b)

        
        colorlist = ['royalblue', 'g', 'r', 'brown', 'm', 'c', 'y', 'black', 'seagreen']
        HeNonPol_PlotVar= ['royalblue', 2.5, 'solid', 'Int - Expected','s']
        HeNonPol_PlotVar2= ['red', 2.5, 'solid', 'Int -Corrected Pol']


        if i == 0:  # In the future, automate this process of the plotting.
            taxs = PlotFunction(HeNonPol2, HeNonPol_PlotVar, NormalizeSig=True,
                                Shape=[len(Bmags), 1], makefig=True,
                                plotlabel=f'B={Bvals}T', plotnondip=0, plotpol=0)
            PlotFunction(HeNonPol2b, HeNonPol_PlotVar2,  NormalizeSig=True, axsin=taxs[i][0], position=[0, i],
                            makefig=False,  plotlabel=f'B={Bvals}T', plotnondip=0, plotpol=0)   
 
        else:
            PlotFunction(HeNonPol2, HeNonPol_PlotVar,  NormalizeSig=True, axsin=taxs[i][0], position=[0, i],
                            makefig=False,  plotlabel=f'B={Bvals}T', plotnondip=0, plotpol=0)   
            PlotFunction(HeNonPol2b, HeNonPol_PlotVar2,  NormalizeSig=True, axsin=taxs[i][0], position=[0, i],
                            makefig=False,  plotlabel=f'B={Bvals}T', plotnondip=0, plotpol=0, legcols=3)  
    # if do_title:      
    #     plt.tight_layout(rect=[0, 0, 0.75, .95])
    #     plt.suptitle(f'{InputdeckHeNonPol2['plottitle']}')
    if savepdf:
        plt.savefig('He_NonPol2.pdf', format='pdf')
    if savepng:      
        plt.tight_layout(rect=[0, 0, 1, .9])
        plt.rcParams['axes.titlepad'] = 405
        plt.title(f'{InputdeckHeNonPol2['plottitle']}')
        plt.savefig('He_NonPol2B.png', dpi=savedpi)   

#%%
if doHeNonPol3 == True:
    pathfil = os.getcwd()
    plt.close('all')
    fil = pathfil + '/He_Spectra_NonPolarized.csv'
    He_NonPol_Spec = read_Spectra(fil, scalewave=0.1, headercount=1)  # Need to scale the wavelengths from angstrom to nm.

    Bmags = 3.2

    # He_NonPol_Spec['wavelength'] = [x + 0.00666 for x in He_NonPol_Spec['wavelength']]


    InputdeckHeNonPol3 = {'s_ground': 1,  # Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                      's_excited': 1,  # Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                      'l_ground': 1,  # Orbital Angular Momentum of ground state, int or half int
                      'l_excited': 2,  # Orbital Angular Momentum of excited state, int or half int
                      'E_ground': np.array([169087.8308131,  169086.8428979,  169086.7664725]),  # Lowest J first, in cm^-1
                      'E_excited': np.array([191444.49804915, 191444.47952984, 191444.47832914]),  # Lowest J first, in cm^-1
                      'Bmag': 3.2,  # Magnetic Field  [T]
                      'b_angle': 90,  # Angle between LoS and Bmax
                      'Convfxn': 'Skewed',  # Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
                      'specstep': 0.002222,  # Resolution of spectrometer in nm
                      'Skewness': [0.62, 0.99],  # How non-symmetric the Lorrentzian is
                      'plottitle': 'He I 1s2p -> 1s4d: B=3.2T, Non-Dipole',  # Title for plotting (optional)
                      'plot_window': [447.1, 447.2],  # Min and max for plot window range (nm), convolutions will only take place within this range
                      'amu': 2,  # Weight in AMU
                      'SpectrumData': [He_NonPol_Spec['wavelength'], He_NonPol_Spec['signals'][0]],  # Spectrometer Data.
                      'spec_window': [447, 447.2],  # Min and max for convolutions
                        'TAlpha' : np.sqrt(1-0.35) ,
                        'TBeta' : np.sqrt(1-0.65) ,
                      # 'DoLowSig': 'Y',  # If this exists in the input, lowfield signal strength will be calculated.
                      # 'DoHighSig': 'Y',  # If this exists in the input, highfield signal strength will be calculated.

                      }
    HeNonPol3 = AtomSpect_Main(InputdeckHeNonPol3)
    
    colorlist = ['royalblue', 'g', 'r', 'brown', 'm', 'c', 'y', 'black', 'seagreen']
    
    HeNonPol_PlotObj = [colorlist[0], 2.5, 'solid', 'Intermediate']
    
    taxs = PlotFunction(HeNonPol3, HeNonPol_PlotObj, NormalizeSig=True,
                        makefig=True,  plotlabel=f'B={Bmags}T', plotpol=1, plotnondip=1)


    # if do_title:      
    #     plt.tight_layout(rect=[0, 0, 0.75, .95])
    #     plt.suptitle(f'{InputdeckHeNonPol3['plottitle']}')
    if savepdf:
        plt.savefig('He_NonPol3.pdf', format='pdf')
    if savepng:      
        plt.tight_layout(rect=[0, 0, 1, 0.9])
        plt.rcParams['axes.titlepad'] = 315
        plt.title(f'{InputdeckHeNonPol3['plottitle']}')
        plt.savefig('He_NonPol3B.png', dpi=savedpi)

    # %%
if fanplot == True:
    InputdeckHeNonPolFan = {'s_ground': 1,  # Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                         's_excited': 1,  # Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                         'l_ground': 1,  # Orbital Angular Momentum of ground state, int or half int
                         'l_excited': 2,  # Orbital Angular Momentum of excited state, int or half int
                         'E_ground': np.array([169087.8308131,  169086.8428979,  169086.7664725]),  # Lowest J first, in cm^-1
                         'E_excited': np.array([191444.49804915, 191444.47952984, 191444.47832914]),  # Lowest J first, in cm^-1
                         # 'Bmag': 0.1,  # Magnetic Field  [T]
                         'b_angle': 90,  # Angle between LoS and Bmax
                         'Convfxn': 'Skewed',  # Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
                         'specstep': 0.002222,  # Resolution of spectrometer in nm
                         'Skewness': [0.62, 0.99],  # How non-symmetric the Lorrentzian is
                         'plottitle': 'He I 1s2p -> 1s4d',  # Title for plotting (optional)
                         'plot_window': [447, 447.2],  # Min and max for plot window range (nm), convolutions will only take place within this range
                         'amu': 2,  # Weight in AMU
                         'spec_window': [447, 447.2],  # Min and max for convolutions

                         'DoLowSig': 'Y',  # If this exists in the input, lowfield signal strength will be calculated.
                         'DoHighSig': 'Y',  # If this exists in the input, lowfield signal strength will be calculated.
                         'EtermG': 169086.7664725, #Term energy of lower state in cm^-1

                         # 'EtermG': 169086.9076, #Term energy of lower state in cm^-1
                         'EtermE': 191444.4827, #Term energy of upper state in cm^-1

                         # 'Pol_angle': 90 ,#Angle polarizing filter makes with max linear transmission, Optional
                         # 'sortE' : True

                         }
    
    # HeB0 = AtomSpect_Main(InputdeckHeNonPolFan)
    InputdeckHeNonPolFan['Bmag'] = np.linspace(0, 5, 1000)
    # Brange = np.linspace(0, 5, 100)
    flags = 'HGL'
    Zfantest = plotZfan(InputdeckHeNonPolFan, flags=flags,markercount = 35)
    # plt.ylim(169083.4, 169090)
    plt.show()
    # plt.tight_layout()
    if savepdf:
        plt.savefig('He_Fan.pdf', format='pdf')
    if savepng:
        plt.savefig('He_Fan.png', dpi=savedpi)
    flags_ = 'G'
    Zfantest = plotZfan(InputdeckHeNonPolFan, flags=flags_,markercount = 0)
    # plt.xlim(0,1)
    plt.show()
    # plt.tight_layout()
    if savepdf:
        plt.savefig('He_Fan_Int.pdf', format='pdf')
    if savepng:
        plt.savefig('He_Fan_Int.png', dpi=savedpi)
        
        #%%
# fanplot2 = 1
if fanplot2 == True:
    # plt.close('all')
    InputdeckHeNonPolFan = {'s_ground': 1,  # Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                         's_excited': 1,  # Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                         'l_ground': 1,  # Orbital Angular Momentum of ground state, int or half int 
                         'l_excited': 2,  # Orbital Angular Momentum of excited state, int or half int
                         'E_ground': np.array([169087.8308131,  169086.8428979,  169086.7664725]),  # Lowest J first, in cm^-1
                         'E_excited': np.array([191444.49804915, 191444.47952984, 191444.47832914]),  # Lowest J first, in cm^-1
                         'Bmag': 2.5,  # Magnetic Field  [T]
                         'b_angle': 90,  # Angle between LoS and Bmax
                         'Convfxn': 'Skewed',  # Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
                         'specstep': 0.002222,  # Resolution of spectrometer in nm
                         'Skewness': [0.62, 0.99],  # How non-symmetric the Lorrentzian is
                         'plottitle': 'He I 1s2p -> 1s4d',  # Title for plotting (optional)
                         'plot_window': [447, 447.2],  # Min and max for plot window range (nm), convolutions will only take place within this range
                         'amu': 2,  # Weight in AMU
                         'spec_window': [447, 447.2],  # Min and max for convolutions

                         'DoLowSig': 'Y',  # If this exists in the input, lowfield signal strength will be calculated.
                         'DoHighSig': 'Y',  # If this exists in the input, lowfield signal strength will be calculated.
                         'EtermG': 169086.7664725, #Term energy of lower state in cm^-1
                         'EtermE': 191444.47832914, #Term energy of upper state in cm^-1

                         # 'Pol_angle': 90 ,#Angle polarizing filter makes with max linear transmission, Optional
                         # 'sortE' : True

                         }
    InputdeckHeNonPolFan['Bmag'] = np.linspace(0, .075, 1000)
    # Brange = np.linspace(0, .25, 100)
    flags = 'ELH'
    # flags = 'E'

    Zfantest = plotZfan(InputdeckHeNonPolFan, flags=flags, markercount = 50)
    # plt.ylim(169083.4, 169090)
    savepdf=1
    if savepdf:
        plt.savefig('He_FanE.pdf',format='pdf')
    if savepng:
        plt.savefig('He_FanE.png', dpi=savedpi)

#%%

if do_HeSlider:
    pathfil = os.getcwd()
    # plt.close('all')
    fil = pathfil + '/He_Spectra_NonPolarized.csv'
    He_NonPol_Spec = read_Spectra(fil, scalewave=0.1, headercount=1)  # Need to scale the wavelengths from angstrom to nm.

    SliderHe = {'s_ground': 1,  # Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                's_excited': 1,  # Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                'l_ground': 1,  # Orbital Angular Momentum of ground state, int or half int
                'l_excited': 2,  # Orbital Angular Momentum of excited state, int or half int
                'E_ground': np.array([169087.8308131,  169086.8428979,  169086.7664725]),  # Lowest J first, in cm^-1
                'E_excited': np.array([191444.49804915, 191444.47952984, 191444.47832914]),  # Lowest J first, in cm^-1
                'Bmag': 1.5,  # Magnetic Field
                'b_angle': 90,  # Angle between LoS and Bmax
                'Convfxn': 'GaussianInstrum',  # Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
                'specstep': 0.002222,  # Resolution of spectrometer in nm
                # 'Skewness' : [0.62,0.99] , #How non-symmetric the Lorrentzian is
                # 'Pol_angle': 90 ,#Angle polarizing filter makes with max linear transmission, Optional
                'plottitle': 'He I 1s2p -> 1s4d: B=1.5T, Polarized',  # Title for plotting (optional)
                'spec_window': [447, 447.2],  # Min and max for convolutions

                'plot_window': [447.1, 447.2],  # Min and max for plot window range (nm), convolutions will only take place within this range
                'amu': 2,  # Weight in AMU
                'specres': 10,  # How many steps per resolution are calculated. Higher makes a smoother curve.
                'DoLowSig': 'Yes',  # If this exists in the input, lowfield signal strength will be calculated.
                'DoHighSig': 'Yes',  # If this exists in the input, lowfield signal strength will be calculated.
                'EtermG': 169086.7664725, #Term energy of lower state in cm^-1
                'EtermE': 191444.4827, #Term energy of upper state in cm^-1

                'fxnwindow': 0.25, #How far the convolution evaluates from the last stick in nm

                }
    MakeSlider(SliderHe, [He_NonPol_Spec['wavelength'], Normalize(He_NonPol_Spec['signals'][1])], [447.1, 447.2], banglevary=1, polvary=0)


if do_Hepolslider:
    pathfil = os.getcwd()
    filname = '/He_Spectra_1.5T_Polarized.csv'
    fil = pathfil + filname

    He_Pol_Spec = read_Spectra(fil, headercount=2)

    # This corrects for a misalignment in the calibration by 3 pixels.
    He_PolShift = [0.1*x + 0.006666 for x in He_Pol_Spec['wavelength']]

    polangles = [0, 30, 70, 90]

    # fig, axs = plt.subplots(4,sharex=True,sharey=True, figsize = (12,8))
    colorlist = ['royalblue', 'g', 'r', 'brown', 'm', 'c', 'y', 'black', 'seagreen']

    SliderHePol = {'s_ground': 1,  # Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                   's_excited': 1,  # Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                   'l_ground': 1,  # Orbital Angular Momentum of ground state, int or half int
                   'l_excited': 2,  # Orbital Angular Momentum of excited state, int or half int
                   'E_ground': np.array([169087.8308131,  169086.8428979,  169086.7664725]),  # Lowest J first, in cm^-1
                   'E_excited': np.array([191444.49804915, 191444.47952984, 191444.47832914]),  # Lowest J first, in cm^-1
                   'Bmag': 1.5,  # Magnetic Field
                   'b_angle': 90,  # Angle between LoS and Bmax
                   'Convfxn': 'GaussianInstrum',  # Optional: 'Gaussian', 'Skewed'. 'Lorrentzian' , 'GaussinInstrum'
                   'specstep': 0.002222,  # Resolution of spectrometer in nm
                   # 'Skewness' : [0.62,0.99] , #How non-symmetric the Lorrentzian is
                   'spec_window': [447, 447.2],  # Min and max for convolutions
                   'plottitle': 'He I 1s2p -> 1s4d: B=1.5T, Polarized',  # Title for plotting (optional)
                   'plot_window': [447.1, 447.2],  # Min and max for plot window range (nm), convolutions will only take place within this range
                   'amu': 2,  # Weight in AMU
                   'specres': 10,  # How many steps per resolution are calculated. Higher makes a smoother curve.
                   'fxnwindow': 0.25, #How far the convolution evaluates from the last stick in nm

                   # 'DoLowSig' : 'Yes' , #If this exists in the input, lowfield signal strength will be calculated.
                   }
    MakeSlider(SliderHePol, [He_PolShift, Normalize(He_Pol_Spec['signals'][-1])], [447.1, 447.2],
               banglevary=1, polvary=1, do_other=1)

if closeplots:
    plt.close('all')

#%% Save/load testing
#Testing the save and load functions and ensuring the loaded values work as intended.
# He_Pol_Specwav = [0.1*x +0.006666 for x in He_Pol_Spec['wavelength']]
# He_Pol_Spec2 = [Normalize(x) for x in He_Pol_Spec['signals']]
# Savenpy_Dict('He_asdf2',SliderHePol,spectra=[He_Pol_Specwav, He_Pol_Spec2],Normalize_sig = 1)
# #%%

# testHeDict,testHeLspec = Loadnpy_Dict('He_asdf2')
# #%%
# MakeSlider(testHeDict, [testHeLspec['wavelength'] ,Normalize(testHeLspec['signals'][-1]) ], [447.1, 447.2],
#                banglevary=1, polvary=1, do_other=1)


