# -*- coding: utf-8 -*-
"""

@author: Leo Nofs
"""

import  os, sys
import numpy as np
from matplotlib import pyplot as plt

#This adds the directy two folders up from the example file so taht the main module can be imported.
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir, os.path.pardir)))
from SpectSplit import Zeeman_Main, Normalize,Normalize,PlotFunction,read_Spectra

#Add path to the spectral data
pathfil = os.getcwd()
filname = '/Ar_specs.csv' #There are multiple magnetic field values in this file which will be split and parsed using read_Spectra utility.
fil = pathfil + filname

#Import and process the spectral data
Spec_Raw = read_Spectra(fil, headercount=2, keepheaders=True)
Ar_Spec = Spec_Raw[0]
Ar_B = [float(x) for x in  Spec_Raw[1][0][1:]]
savefigs = 1

#This corrects for a misalignment in the calibration by 3 pixels. 
shiftedAr = [x + (5-.3 -34*.002222) for x in Ar_Spec['wavelength'] ]  #Shifting because I think the raw data wasn't calibrated right

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
                  'Bmag': Bvals , #Magnetic Field  [T]
                  'b_angle': 90 , #Angle between LoS and Bmax
                  'plottitle': f'Ar I 3p5 4s (1P) -> 3p5 4p (1S)' , #Title for plotting (optional)     
                  'spec_window' : [745,755], #Min and max convolutions 
                  'plot_window' : [750.25,750.5], #Min and max for plot window range (nm)
                  'amu' : 40  , #Weight in AMU
                  'specstep' : 0.004 , #Pixel resoltion of spectrometer in nm
                  'SpectrumData': [shiftedAr, Ar_Spec['signals'][i]] , 
                  'DoHighSig' : 'Y' , #If this exists in the input, lowfield signal strength will be calculated. 
                  'DoLowSig' : 'Y' , #If this exists in the input, highfield signal strength will be calculated. 
                   }
    

    
    Ar1 = Zeeman_Main(InputdeckAr1)    
    
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

