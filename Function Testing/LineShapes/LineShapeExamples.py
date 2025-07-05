# -*- coding: utf-8 -*-
"""
Simple script to test the lineshape functions

@author: %Leo Nofs
"""


import os, sys

#This is the base folder of the main function. Required for importing the package. 
#This appends the path two folders down /../../ into the path to import the package.
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir, os.path.pardir)))
from SpectSplit import Gaussian, GaussianInstrum, Lorrentzian,Voigt,Normalize

import numpy as np
from matplotlib import pyplot as plt

savefigs =True #Flag for whether figures should be saved.

c_light = 299792458 #Speed of light in m/s
kboltz = 1.381e-23 #Boltzmann Constant
m_prot = 1.667e-27 #Mass of proton in kg



xRange = np.linspace(401.95e-9,402.05e-9, 1000)
xRange2 = [x*1e9 for x in xRange]


xpeak = 402e-9
peakplus = xpeak*(1 )
xpeak2 = xpeak*1e9

Temp_K = 3*11602
Atomic_Mass = 183 #Mass in amu
pixel_width = 0.002222 #Pixel width of spectrometer in nm
skewness = [1,.5]

ionv1 = 2500
ionv2 = 3500
ionv = [2500,3500]

peakplus = xpeak2*(1 + ionv1/3e8 )
peakminus = xpeak2*(1 - ionv2/3e8 )


Voigtv0 = Voigt(xRange,xpeak,Temp_K,Atomic_Mass,pixel_size=pixel_width*1e-9,v_ion = 0)
Voigtvplus = Voigt(xRange,xpeak,Temp_K,Atomic_Mass,pixel_size=pixel_width*1e-9,v_ion = ionv1)
Voigtvminus = Voigt(xRange,xpeak,Temp_K,Atomic_Mass,pixel_size=pixel_width*1e-9,v_ion = -ionv2 )
Voigtbidir = Voigt(xRange,xpeak,Temp_K,Atomic_Mass,pixel_size=pixel_width*1e-9,v_ion = ionv, bidir=True)

Gaussian0 = Gaussian(xRange,xpeak,Temp_K,Atomic_Mass,v_ion = 0)
Gaussianplus = Gaussian(xRange,xpeak,Temp_K,Atomic_Mass,v_ion = ionv1)
Gaussianminus = Gaussian(xRange,xpeak,Temp_K,Atomic_Mass,v_ion = -ionv2)
Gaussianbidir = Gaussian(xRange,xpeak,Temp_K,Atomic_Mass,v_ion = ionv, bidir=True)

#(x, x0, Temp_in, amu,pixel_size = 0,v_ion = 0):
GaussianInstrum0 = GaussianInstrum(xRange,xpeak,Temp_K,Atomic_Mass,pixel_size=pixel_width*1e-9,v_ion = 0)
GaussianInstrumplus = GaussianInstrum(xRange,xpeak,Temp_K,Atomic_Mass,pixel_size=pixel_width*1e-9,v_ion = ionv1)
GaussianInstrumminus = GaussianInstrum(xRange,xpeak,Temp_K,Atomic_Mass,pixel_size=pixel_width*1e-9,v_ion = -ionv2 )
GaussianInstrumbidir = GaussianInstrum(xRange,xpeak,Temp_K,Atomic_Mass,pixel_size=pixel_width*1e-9,v_ion =ionv, bidir=True)

#Lorrentzian(xrange, peak, pixel_size=0.002222 ,Skewness=[1,1], v_ion=0):
Lorrentzian0 = Lorrentzian(xRange,xpeak,pixel_size=pixel_width*1e-9,v_ion = 0)
Lorrentzianplus = Lorrentzian(xRange,xpeak,pixel_size=pixel_width*1e-9,v_ion = ionv1)
Lorrentzianminus = Lorrentzian(xRange,xpeak,pixel_size=pixel_width*1e-9,v_ion = -ionv2)
Lorrentzianbidir = Lorrentzian(xRange,xpeak,pixel_size=pixel_width*1e-9,v_ion =ionv, bidir=True)

#Skewed 1,.25
SkewLorrentzian0 = Lorrentzian(xRange,xpeak,pixel_size=pixel_width*1e-9,Skewness=skewness,v_ion = 0)
SkewLorrentzianplus = Lorrentzian(xRange,xpeak,pixel_size=pixel_width*1e-9,Skewness=skewness,v_ion = ionv1)
SkewLorrentzianminus = Lorrentzian(xRange,xpeak,pixel_size=pixel_width*1e-9,Skewness=skewness,v_ion = -ionv2)
SkewLorrentzianbidir = Lorrentzian(xRange,xpeak,pixel_size=pixel_width*1e-9,Skewness=skewness,v_ion =ionv, bidir=True)
#%%
plt.close('all')

filpath= 'Examples/LineShapes/'
#%%Voigt

plt.figure()
plt.title(f'Voigt: T={Temp_K/11602} eV, AMU = {Atomic_Mass}, Pixels = {pixel_width} nm')
plt.ylabel('Intensity (arb)')
plt.xlabel('Wavelength (nm)')
plt.plot(xRange2,Normalize(Voigtv0),label='Voigt, v=0')
plt.plot(xRange2,Normalize(Voigtvplus),label=f'Voigt, v={ionv1}')
plt.plot(xRange2,Normalize(Voigtvminus),label=f'Voigt , v=-{ionv2}')
plt.plot(xRange2,Normalize(Voigtbidir),linestyle = '--',label=f'Voigt, bidir {ionv}')
plt.vlines(xpeak2,0,1,linewidth = 2, color='black',label='Peak' )
plt.vlines(peakplus,0,1,linewidth = 2, color='black',linestyle = ':', label='Peak+' )
plt.vlines(peakminus,0,1,linewidth = 2, color='black',linestyle = '--', label='Peak-' )
plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left")
plt.rcParams['figure.figsize'] = [14,8]
plt.tight_layout()
if savefigs:
    plt.savefig(filpath+'VoigtExample.png',)

#%%Gaussian

plt.figure()
plt.title(f'Gaussian: T={Temp_K/11602} eV, AMU = {Atomic_Mass}')
plt.ylabel('Intensity (arb)')
plt.xlabel('Wavelength (nm)')
plt.plot(xRange2,Normalize(Gaussian0),label='Gaussian, v=0')
plt.plot(xRange2,Normalize(Gaussianplus),label=f'Gaussian, v={ionv1}')
plt.plot(xRange2,Normalize(Gaussianminus),label=f'Gaussian , v=-{ionv2}')
plt.plot(xRange2,Normalize(Gaussianbidir),linestyle = '--',label=f'Gaussian, bidir {ionv}')
plt.vlines(xpeak2,0,1,linewidth = 2, color='black',label='Peak' )
plt.vlines(peakplus,0,1,linewidth = 2, color='black',linestyle = ':', label='Peak+' )
plt.vlines(peakminus,0,1,linewidth = 2, color='black',linestyle = '--', label='Peak-' )
plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left")
plt.rcParams['figure.figsize'] = [14,8]
plt.tight_layout()
if savefigs:

    plt.savefig(filpath+'GaussianExample.png',)

#%%GaussianInstrum

plt.figure()
plt.title(f'GaussianInstrum: T={Temp_K/11602} eV, AMU = {Atomic_Mass}, Pixels = {pixel_width} nm')
plt.ylabel('Intensity (arb)')
plt.xlabel('Wavelength (nm)')
plt.plot(xRange2,Normalize(GaussianInstrum0),label='GuassianInstrum, v=0')
plt.plot(xRange2,Normalize(GaussianInstrumplus),label=f'GuassianInstrum, v={ionv1}')
plt.plot(xRange2,Normalize(GaussianInstrumminus),label=f'GuassianInstrum , v=-{ionv2}')
plt.plot(xRange2,Normalize(GaussianInstrumbidir),linestyle = '--',label=f'GuassianInstrum, bidir {ionv}')
plt.vlines(xpeak2,0,1,linewidth = 2, color='black',label='Peak' )
plt.vlines(peakplus,0,1,linewidth = 2, color='black',linestyle = ':', label='Peak+' )
plt.vlines(peakminus,0,1,linewidth = 2, color='black',linestyle = '--', label='Peak-' )
plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left")
plt.rcParams['figure.figsize'] = [14,8]
plt.tight_layout()
if savefigs:

    plt.savefig(filpath+'GaussianInsrumExample.png',)

#%%Lorrentzian

plt.figure()
plt.title(f'Lorrentzian: T={Temp_K/11602} eV, AMU = {Atomic_Mass}, Pixels = {pixel_width} nm')
plt.ylabel('Intensity (arb)')
plt.xlabel('Wavelength (nm)')
plt.plot(xRange2,Normalize(Lorrentzian0),label='Lorrentzian, v=0')
plt.plot(xRange2,Normalize(Lorrentzianplus),label=f'Lorrentzian, v={ionv1}')
plt.plot(xRange2,Normalize(Lorrentzianminus),label=f'Lorrentzian , v=-{ionv2}')
plt.plot(xRange2,Normalize(Lorrentzianbidir),linestyle = '--',label=f'Lorrentzian, bidir {ionv}')
plt.vlines(xpeak2,0,1,linewidth = 2, color='black',label='Peak' )
plt.vlines(peakplus,0,1,linewidth = 2, color='black',linestyle = ':', label='Peak+' )
plt.vlines(peakminus,0,1,linewidth = 2, color='black',linestyle = '--', label='Peak-' )
plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left")
plt.rcParams['figure.figsize'] = [14,8]
plt.tight_layout()
if savefigs:

    plt.savefig(filpath+'LorrentzianExample.png',)

#%%SkewedLorrentzian

plt.figure()
plt.title(f'Skewed Lorrentzian: T={Temp_K/11602} eV, AMU = {Atomic_Mass}, Pixels = {pixel_width} nm, Skewness = {skewness}')
plt.ylabel('Intensity (arb)')
plt.xlabel('Wavelength (nm)')
plt.plot(xRange2,Normalize(SkewLorrentzian0),label='Skewed Lorrentzian, v=0')
plt.plot(xRange2,Normalize(SkewLorrentzianplus),label=f'Skewed Lorrentzian, v={ionv1}')
plt.plot(xRange2,Normalize(SkewLorrentzianminus),label=f'Skewed Lorrentzian , v=-{ionv2}')
plt.plot(xRange2,Normalize(SkewLorrentzianbidir),linestyle = '--',label=f'Skewed Lorrentzian, bidir {ionv}')
plt.vlines(xpeak2,0,1,linewidth = 2, color='black',label='Peak' )
plt.vlines(peakplus,0,1,linewidth = 2, color='black',linestyle = ':', label='Peak+' )
plt.vlines(peakminus,0,1,linewidth = 2, color='black',linestyle = '--', label='Peak-' )
plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left")
plt.rcParams['figure.figsize'] = [14,8]
plt.tight_layout()
if savefigs:

    plt.savefig(filpath+'SkewedLorrentzianExample.png',)
#%%All at v=0



plt.figure()
plt.title('Comparing all Lineshapes without Doppler')
plt.ylabel('Intensity (arb)')
plt.xlabel('Wavelength (nm)')
plt.plot(xRange2,Normalize(Voigtv0),label='Voigt')
plt.plot(xRange2,Normalize(Gaussian0),label='Gaussian')
plt.plot(xRange2,Normalize(GaussianInstrum0),label='GuassianInstrum')
plt.plot(xRange2,Normalize(Lorrentzian0),label='Lorrentzian')
plt.plot(xRange2,Normalize(SkewLorrentzian0),label='Skewed Lorrentzian')
plt.vlines(xpeak2,0,1, color='black',label='Peak' )
plt.vlines(peakplus,0,1,linewidth = 2, color='black',linestyle = ':', label='Peak+' )
plt.vlines(peakminus,0,1,linewidth = 2, color='black',linestyle = '--', label='Peak-' )
plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left")
plt.rcParams['figure.figsize'] = [14,8]
plt.tight_layout()
if savefigs:

    plt.savefig(filpath+'All_NoDoppler.png',)
#%%All plus

plt.figure()
plt.title(f'Doppler+: T={Temp_K/11602} eV, AMU = {Atomic_Mass}, Pixels = {pixel_width} nm, Skewness = {skewness}, v={ionv1}')
plt.ylabel('Intensity (arb)')
plt.xlabel('Wavelength (nm)')
plt.plot(xRange2,Normalize(Voigtvplus),label='Voigt')
plt.plot(xRange2,Normalize(Gaussianplus),label='Gaussian')
plt.plot(xRange2,Normalize(GaussianInstrumplus),label='GuassianInstrum')
plt.plot(xRange2,Normalize(Lorrentzianplus),label='Lorrentzian')
plt.plot(xRange2,Normalize(SkewLorrentzianplus),label='Skewed Lorrentzian')
plt.vlines(xpeak2,0,1, color='black',label='Peak' )
plt.vlines(peakplus,0,1,linewidth = 2, color='black',linestyle = ':', label='Peak+' )
# plt.vlines(peakminus,0,1,linewidth = 2, color='black',linestyle = '--', label='Peak-' )
plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left")
plt.rcParams['figure.figsize'] = [14,8]
plt.tight_layout()
if savefigs:

    plt.savefig(filpath+'All_Vplus.png',)
#%%All minus

plt.figure()
plt.title(f'Doppler-: T={Temp_K/11602} eV, AMU = {Atomic_Mass}, Pixels = {pixel_width} nm, Skewness = {skewness}, v=-{ionv2}')
plt.ylabel('Intensity (arb)')
plt.xlabel('Wavelength (nm)')
plt.plot(xRange2,Normalize(Voigtvminus),label='Voigt')
plt.plot(xRange2,Normalize(Gaussianminus),label='Gaussian')
plt.plot(xRange2,Normalize(GaussianInstrumminus),label='GuassianInstrum')
plt.plot(xRange2,Normalize(Lorrentzianminus),label='Lorrentzian')
plt.plot(xRange2,Normalize(SkewLorrentzianminus),label='Skewed Lorrentzian')
plt.vlines(xpeak2,0,1, color='black',label='Peak' )
# plt.vlines(peakplus,0,1,linewidth = 2, color='black',linestyle = ':', label='Peak+' )
plt.vlines(peakminus,0,1,linewidth = 2, color='black',linestyle = '--', label='Peak-' )
plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left")
plt.rcParams['figure.figsize'] = [14,8]
plt.tight_layout()
if savefigs:

    plt.savefig(filpath+'All_Vminus.png',)
#%%All Bidir

plt.figure()
plt.title(f'Bidirectional: T={Temp_K/11602} eV, AMU = {Atomic_Mass}, Pixels = {pixel_width} nm, Skewness = {skewness}, v={ionv}')
plt.ylabel('Intensity (arb)')
plt.xlabel('Wavelength (nm)')
plt.plot(xRange2,Normalize(Voigtbidir),label='Voigt')
plt.plot(xRange2,Normalize(Gaussianbidir),label='Gaussian')
plt.plot(xRange2,Normalize(GaussianInstrumbidir),label='GuassianInstrum')
plt.plot(xRange2,Normalize(Lorrentzianbidir),label='Lorrentzian')
plt.plot(xRange2,Normalize(SkewLorrentzianbidir),label='Skewed Lorrentzian')
plt.vlines(xpeak2,0,1, color='black',label='Peak' )
plt.vlines(peakplus,0,1,linewidth = 2, color='black',linestyle = ':', label='Peak+' )
plt.vlines(peakminus,0,1,linewidth = 2, color='black',linestyle = '--', label='Peak-' )
plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left")
plt.rcParams['figure.figsize'] = [14,8]
plt.tight_layout()
if savefigs:

    plt.savefig(filpath+'All_Bidir.png',)



