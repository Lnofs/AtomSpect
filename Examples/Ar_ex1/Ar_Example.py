# -*- coding: utf-8 -*-
"""

@author: Leo Nofs
"""

import  os, sys
import numpy as np
from matplotlib import pyplot as plt
from brokenaxes import brokenaxes #This is only for adding a break in the axes to remove "dead space"
#This adds the directy two folders up from the example file so taht the main module can be imported.
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir, os.path.pardir)))
from AtomSpect import AtomSpect_Main,PlotFunction,read_Spectra , plotZfan ,plotarray



#Add path to the spectral data
pathfil = os.getcwd()
filname = '/Ar_specs.csv' #There are multiple magnetic field values in this file which will be split and parsed using read_Spectra utility.
fil = pathfil + filname

#Import and process the spectral data
Spec_Raw = read_Spectra(fil, headercount=2, keepheaders=True)
Ar_Spec = Spec_Raw[0]
Ar_B = [float(x) for x in  Spec_Raw[1][0][1:]]


#This corrects for a misalignment in the calibration by 3 pixels. 
shiftedAr = [x + (5-.3 -34*.002222) for x in Ar_Spec['wavelength'] ]  #Shifting because I think the raw data wasn't calibrated right

savefigs = 1
doAllB = 0 #Just the 750.4 nm line

doAllB2 = 0 #Both the 750.4 and 751.4 nm lines






if doAllB:
    plt.close('all')
    #Iterate over all of the Bfield values in the imported data.
    for i,Bvals in enumerate(Ar_B):
        # fig, axs = plt.subplots(1, figsize = (10,12))
        # plt.close('all')
    
        InputdeckAr1 = {'s_ground': 0 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                      's_excited': 0 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                      'l_ground': 1 , #Orbital Angular Momentum of ground state, int or half int
                      'l_excited': 0 ,#Orbital Angular Momentum of excited state, int or half int
                      'E_ground': np.array([95399.8276]) ,#Lowest J first, in cm^-1
                      'E_excited': np.array([108722.6194]) , #Lowest J first, in cm^-1
                      'Bmag': Bvals , #Magnetic Field  [T]
                      'b_angle': 90 , #Angle between LoS and Bmax
                      'plottitle': f'Ar I 3p5 4s (1P) -> 3p5 4p (1S)' , #Title for plotting (optional)     
                      'spec_window' : [745,755], #Min and max convolutions 
                      'plot_window' : [750.25,750.5], #Min and max for plot window range (nm)
                      'amu' : 40  , #Weight in AMU
                      'specstep' : 0.004 , #Pixel resoltion of spectrometer in nm
                      'SpectrumData': [shiftedAr, Ar_Spec['signals'][i]] , 
                      # 'DoHighSig' : 'Y' , #If this exists in the input, lowfield signal strength will be calculated. 
                      # 'DoLowSig' : 'Y' , #If this exists in the input, highfield signal strength will be calculated. 
                       }
        
    
        
        Ar1 = AtomSpect_Main(InputdeckAr1)    
        
        plotvars = ['royalblue', 1.5, 'solid', 'Exact']
        taxs = PlotFunction(Ar1, plotvars, NormalizeSig=True, makefig=True, plotlabel=f'B={Bvals}T')
                    
        
        if savefigs:
            # plt.savefig(f'ArI_B={Bvals}.pdf',dpi=1200,format='pdf')
            # plt.savefig(f'ArI_B={Bvals}.png',dpi=800)
    
            do_title=1
            if do_title:      
                plt.tight_layout(rect=[0, 0, 1, .90])
                plt.title(f'{InputdeckAr1['plottitle']}')
                # plt.savefig(f'ArI_B={Bvals}B.pdf',dpi=800)
                plt.savefig(f'ArI_B={Bvals}.png',dpi=144)



#AllB2 does both the 750.4 and 751.4 nm lines in each plot. These examples also include weak and strong field for each plot.
if doAllB2:
    plt.close('all')
    #Iterate over all of the Bfield values in the imported data.
    for i,Bvals in enumerate(Ar_B):
        # fig, axs = plt.subplots(1, figsize = (10,12))
        plt.close('all')
    
        InputdeckAr1 = {'s_ground': 0 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                      's_excited': 0 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                      'l_ground': 1 , #Orbital Angular Momentum of ground state, int or half int
                      'l_excited': 0 ,#Orbital Angular Momentum of excited state, int or half int
                      'E_ground': np.array([95399.8276]) ,#Lowest J first, in cm^-1
                      'E_excited': np.array([108722.6194]) , #Lowest J first, in cm^-1
                      'EtermG' : 95399.8276 , # Term energy for the ground level
                      'EtermE' : 108722.6194  ,# Term Energy for the excited level
                      'Bmag': Bvals , #Magnetic Field  [T]
                      'b_angle': 90 , #Angle between LoS and Bmax
                      'plottitle': f'Ar I 3p5 4s (1P) -> 3p5 4p (1S)' , #Title for plotting (optional)     
                      'spec_window' : [745,755], #Min and max convolutions 
                      'plot_window' : [750.25,751.6], #Min and max for plot window range (nm)
                      'amu' : 40  , #Weight in AMU
                      'specstep' : 0.004 , #Pixel resoltion of spectrometer in nm
                      'SpectrumData': [shiftedAr, Ar_Spec['signals'][i]] , 
                      'DoHighSig' : 'Y' , #If this exists in the input, lowfield signal strength will be calculated. 
                      'DoLowSig' : 'Y' , #If this exists in the input, highfield signal strength will be calculated. 
                       }
        
    
        

        InputdeckAr2 = {'s_ground': 1 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                      's_excited': 1 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                      'l_ground': 1 , #Orbital Angular Momentum of ground state, int or half int
                      'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                      # 'E_ground': np.array([94541.62, 93750.5978,93131.89]) , #Lowest J first, in cm^-1
                      # 'E_excited': np.array([107054.6, 107482.7, 106224]) , #Lowest J first, in cm^-1 
                      
                      'E_ground': np.array([94541.62, 93750.5978,93131.89]) , #Lowest J first, in cm^-1
                      'E_excited': np.array([107054.6, 107482.7, 106224]) , #Lowest J first, in cm^-1 
                      
                      #These term levels don't match exactly, otherwise it would be roughly 692 nm line.
                      'EtermG' : 93750.5978, # Term energy for the ground level
                      'EtermE' : 107054.6 , # Term Energy for the excited level
                      
                      'Bmag': Bvals ,  #Magnetic Field  [T]
                      'b_angle': 90 , #Angle between LoS and Bmax
                      'plottitle': f'Ar I 3p5 4s (3P) -> 3p5 4p (3P)' , #Title for plotting (optional)     
                      'spec_window' : [745,755], #Min and max convolutions 
                      'plot_window' : [750.25,751.6], #Min and max for plot window range (nm)
                      'amu' : 40  , #Weight in AMU
                      'specstep' : 0.004 , #Pixel resoltion of spectrometer in nm
                      'SpectrumData': [shiftedAr, Ar_Spec['signals'][i]] , 
                      'specres': 10,  # How many steps per resolution are calculated. Higher makes a smoother curve.
                      'DoHighSig' : 'Y' , #If this exists in the input, lowfield signal strength will be calculated. 
                      'DoLowSig' : 'Y' , #If this exists in the input, highfield signal strength will be calculated. 
                       }
       
    
        Ar1 = AtomSpect_Main(InputdeckAr1)    
        
        plotvars = ['royalblue', 1.5, 'solid', '750.4 nm - Int']
        plotvars1L = ['purple', 1, '-.', '750.4 nm - Weak', 'X']
        plotvars1H = ['forestgreen', 1, '-.', '750.4 nm - Strong', '>']

        taxs = PlotFunction(Ar1, plotvars, NormalizeSig=True, makefig=1,plotpol=0, plotlabel=f'B={Bvals}T', fig_size = (14,14),marker_count=20)
        
        PlotFunction(Ar1['SpecOutHigh'], plotvars1H, NormalizeSig=True, makefig=0,plotwind=[750.25,751.6],
                     axsin = taxs, plotpol=0)
        
        PlotFunction(Ar1['SpecOutLow'], plotvars1L, NormalizeSig=True, makefig=0, plotwind=[750.25,751.6],
                      axsin = taxs, plotpol=0)           
        Ar2 = AtomSpect_Main(InputdeckAr2)    
        
        plotvars2 = ['orange', 1, '-.', '751.4 nm - Int', 's']
        
        plotvars2L = ['magenta', 1, '-.', '751.4 nm - Weak', 'p']
        plotvars2H = ['darkturquoise', 1, '-.', '751.4 nm - Strong', 'd']


        PlotFunction(Ar2, plotvars2, NormalizeSig=True, makefig=0, NormalizeScale = np.max(Ar_Spec['signals'][i][690:]),
                     axsin = taxs, plotpol=0,legcols=3,marker_count=20)
        
        PlotFunction(Ar2['SpecOutHigh'], plotvars2H, NormalizeSig=True, makefig=0,plotwind=[750.25,751.6], NormalizeScale = np.max(Ar_Spec['signals'][i][690:]),
                     axsin = taxs, plotpol=0)
        
        PlotFunction(Ar2['SpecOutLow'], plotvars2L, NormalizeSig=True, makefig=0, plotwind=[750.25,751.6],NormalizeScale = np.max(Ar_Spec['signals'][i][690:]),
                     axsin = taxs, plotpol=0,legcols=3)
        
        # savefigs = 1
 
        if savefigs:
            plt.tight_layout(rect=[0, 0, 1, .90])
            plt.savefig(f'Ar2_B={Bvals}T.pdf',format='pdf')
            # plt.savefig(f'ArI_B={Bvals}.png',dpi=800)
    
            do_title=1
            if do_title:      
                plt.tight_layout(rect=[0, 0, 1, .90])
                plt.title(f'{InputdeckAr2['plottitle']}')
                # plt.savefig(f'ArI_B={Bvals}B.pdf',dpi=800)
                # plt.savefig(f'Ar2_B={Bvals}.png',format='pdf')


#%%TestingCombinedPlot for paper. The idea being to put both magnetic fields on the same plot, and have both transitions then.
doComboPlot =1
# savefigs=1
if doComboPlot:
    plt.close('all')
    B1 = 2.21
    B2 = 3.203
    Bcomb = [2.21,3.203]
    SpecLoc = [2,12]
    taxs = None
    for i in range(2):
        InputdeckAr1 = {'s_ground': 0 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                      's_excited': 0 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                      'l_ground': 1 , #Orbital Angular Momentum of ground state, int or half int
                      'l_excited': 0 ,#Orbital Angular Momentum of excited state, int or half int
                      'E_ground': np.array([95399.8276]) ,#Lowest J first, in cm^-1
                      'E_excited': np.array([108722.6194]) , #Lowest J first, in cm^-1
                      'EtermG' : 95399.8276 , # Term energy for the ground level
                      'EtermE' : 108722.6194  ,# Term Energy for the excited level
                      'Bmag': Bcomb[i] ,  #Magnetic Field  [T]
                      'b_angle': 90 , #Angle between LoS and Bmax
                      'plottitle': f'Ar I 3p5 4s (1P) -> 3p5 4p (1S)' , #Title for plotting (optional)     
                      # 'spec_window' : [745,755], #Min and max convolutions 
                      # 'plot_window' : [750.25,751.6], #Min and max for plot window range (nm)
                      'spec_window' : [745,755], #Min and max convolutions 
                      'plot_window' : [750.0,751.6], #Min and max for plot window range (nm)
                      'amu' : 40  , #Weight in AMU
                      'specstep' : 0.004 , #Pixel resoltion of spectrometer in nm
                      'SpectrumData': [shiftedAr, Ar_Spec['signals'][SpecLoc[i]] ] , 
                      'specres': 10,  # How many steps per resolution are calculated. Higher makes a smoother curve.

                      # 'DoHighSig' : 'Y' , #If this exists in the input, lowfield signal strength will be calculated. 
                      # 'DoLowSig' : 'Y' , #If this exists in the input, highfield signal strength will be calculated. 
                       }
    
        #This is delta of 13303.6742 cm-1
        InputdeckAr2 = {'s_ground': 1 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                      's_excited': 1 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                      'l_ground': 1 , #Orbital Angular Momentum of ground state, int or half int
                      'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                      # 'E_ground': np.array([94541.62, 93750.5978,93131.89]) , #Lowest J first, in cm^-1
                      # 'E_excited': np.array([107054.6, 107482.7, 106224]) , #Lowest J first, in cm^-1 
                      'E_ground': np.array([94541.62, 93750.5978,93131.89]) , #Lowest J first, in cm^-1
                      'E_excited': np.array([107054.6, 107482.7, 106224]) , #Lowest J first, in cm^-1 
                      'EtermG' : 93750.5978, # Term energy for the ground level
                      'EtermE' : 107054.6 , # Term Energy for the excited level
                      
                      'Bmag': Bcomb[i] ,  #Magnetic Field  [T]
                      'b_angle': 90 , #Angle between LoS and Bmax
                      'plottitle': f'Ar I 3p5 4s (3P) -> 3p5 4p (3P)' , #Title for plotting (optional)     
                      'spec_window' : [750,755] , #Min and max convolutions 
                      'plot_window' : [750.25,751.6], #Min and max for plot window range (nm)
                      'amu' : 40  , #Weight in AMU
                      'specstep' : 0.004 , #Pixel resoltion of spectrometer in nm
                      'SpectrumData': [shiftedAr, Ar_Spec['signals'][SpecLoc[i]] ] , 
                      'specres': 10,  # How many steps per resolution are calculated. Higher makes a smoother curve.
                      # 'DoHighSig' : 'Y' , #If this exists in the input, lowfield signal strength will be calculated. 
                      # 'DoLowSig' : 'Y' , #If this exists in the input, highfield signal strength will be calculated. 
                       }
       
        Ar1 = AtomSpect_Main(InputdeckAr1)    
        plotvars = ['royalblue', 1, 'solid', '750.4nm']
        
        if taxs is not None:
            PlotFunction(Ar1, plotvars, NormalizeSig=True, makefig=0,SpectrometerLS='dashed', axsin=taxs,position=[i,0], SpectrometerMarker = 'D', plotpol=0,marker_size=6)#, plotlabel=f'B={2.21}T')
            # PlotFunction(Ar1['SpecOutHigh'], ['red',1,'dotted', '750.4 nm - Strong','o'], NormalizeSig=0, makefig=0, axsin=taxs,position=[i,0], plotpol=0)#, plotlabel=f'B={2.21}T')

        else:
            taxs = PlotFunction(Ar1, plotvars, NormalizeSig=True, makefig=True,Shape = (2,1),SpectrometerLS='dashed', SpectrometerMarker = 'D', plotpol=0,fig_size=(14,12),marker_size=6)#, plotlabel=f'B={2.21}T')
            # PlotFunction(Ar1['SpecOutHigh'], ['red',1,'dotted', '750.4 nm - Strong','o'], NormalizeSig=0, makefig=0, axsin=taxs,position=[i,0], plotpol=0)#, plotlabel=f'B={2.21}T')

        
    
        Ar2 = AtomSpect_Main(InputdeckAr2)
        plotvars2 = ['red', 1, '-.', '751.4 nm', 's']
    
        PlotFunction(Ar2, plotvars2, NormalizeSig=True, makefig=0,SpectrometerLS='dotted', NormalizeScale = np.max(Ar_Spec['signals'][SpecLoc[i]][690:]),
                     axsin = taxs,position=[i,0], SpectrometerMarker = 'X', plotpol=0 ,  legcols=3,marker_size=6)
            
            
        # PlotFunction(Ar2['SpecOutLow'], ['magenta',1,'dotted', '751.4 nm - Low','o'], NormalizeSig=True, makefig=0,SpectrometerLS='dotted', NormalizeScale = np.max(Ar_Spec['signals'][SpecLoc[i]][690:]),
        #              axsin = taxs,position=[i,0], SpectrometerMarker = 'X', plotpol=0 , plotwind = InputdeckAr2['plot_window'], legcols=3)
            
        # PlotFunction(Ar2['SpecOutHigh'], ['green',1,'dotted', '751.4 nm - Strong','x'], NormalizeSig=True, makefig=0,SpectrometerLS='dotted', NormalizeScale = np.max(Ar_Spec['signals'][SpecLoc[i]][690:]),
        #              axsin = taxs,position=[i,0], SpectrometerMarker = 'X', plotpol=0 , plotwind = InputdeckAr2['plot_window'], legcols=3)
            
        taxs[i,0].text(.95, 0.95, f'B={Bcomb[i]}T', verticalalignment='top',
                      horizontalalignment='right', transform=taxs[i,0].transAxes)
        
        
        
        
    if savefigs:
        plt.savefig(f'ArI_Comb.pdf',format='pdf')
        # plt.savefig(f'ArI_B={Bvals}.png',dpi=800)

        do_title=0
        if do_title:      
            plt.tight_layout(rect=[0, 0, 1, .90])
            plt.title(f'{InputdeckAr1['plottitle']}')
            # plt.savefig(f'ArI_B={Bvals}B.pdf',dpi=800)
            plt.savefig(f'ArI_Comb.png',dpi=144)
#%%
ArFan = 0
if ArFan:
    plt.close('all')
    Brange = np.linspace(0, 5, 500)
    InputdeckAr2Fan = {'s_ground': 1 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
               's_excited': 1 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
               'l_ground': 1 , #Orbital Angular Momentum of ground state, int or half int
               'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
               # 'E_ground': np.array([94541.62, 93750.5978,93131.89]) , #Lowest J first, in cm^-1
               # 'E_excited': np.array([107054.6, 107482.7, 106224]) , #Lowest J first, in cm^-1 
               'E_ground': np.array([94541.62, 93750.5978,93131.89]) , #Lowest J first, in cm^-1
               'E_excited': np.array([107054.6, 107482.7, 106224]) , #Lowest J first, in cm^-1 
               'Bmag': Brange ,  #Magnetic Field  [T]
               'b_angle': 90 , #Angle between LoS and Bmax
               'plottitle': f'Ar I 3p5 4s (3P) -> 3p5 4p (3P)' , #Title for plotting (optional)     
               'spec_window' : [745,755], #Min and max convolutions 
               'plot_window' : [750.25,751.6], #Min and max for plot window range (nm)
               'amu' : 40  , #Weight in AMU
               'specstep' : 0.004 , #Pixel resoltion of spectrometer in nm
               # 'SpectrumData': [shiftedAr, Ar_Spec['signals'][SpecLoc[i]] ] , 
               'specres': 10,  # How many steps per resolution are calculated. Higher makes a smoother curve.
               # 'DoHighSig' : 'Y' , #If this exists in the input, lowfield signal strength will be calculated. 
               # 'DoLowSig' : 'Y' , #If this exists in the input, highfield signal strength will be calculated. 
                }
    flags_ = 'GHL'
    Zfantest = plotZfan(InputdeckAr2Fan, flags=flags_)

    flags_ = 'EHL'
    Zfantest = plotZfan(InputdeckAr2Fan, flags=flags_)


