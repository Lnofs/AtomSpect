#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 17:41:23 2025

@author: Leo Nofs

An example utilizing AtomSpect for Rb I, including HFS. 
"""
import  os, sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker
#This adds the directy two folders up from the example file so taht the main module can be imported.
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir, os.path.pardir)))
from AtomSpect  import AtomSpect_Main, Normalize,PlotFunction,read_Spectra, plotZfan,Convol_Sticks, make_stickbins,MakeSlider, ZeemanFan

doRBfans=0
savefigs=1
savdpi = 1200
InputRb = {'s_ground': .5 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                      's_excited': .5 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                      'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
                      'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                      'E_ground': np.array([0])  ,#Lowest J first, in cm^-1
                      'E_excited':  np.array([12578.950,12816.545]) , #Lowest J first, in cm^-1
                      'Bmag': 0, #Magnetic Field  [T]
                      'b_angle': 90 , #Angle between LoS and Bmax
                      'spec_window' : [779.75,780.25], #Min and max function convolution
                      'plot_window' :  [779.9,780.15], #Min and max for plot window range (nm)

                      'amu' : 85.011 ,  #Weight in AMU
                      'specstep' : 0.005, #Resolution of spectrometer in nm
                      'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'.
                      'Temp' : 300, #Temperature in K. Used for Gaussian convolution
                      'plottitle': f'Rb Testing', #Title for plotting (optional)     
                      'fxnwindow': .5 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                      # 'Pol_angle': Polangle ,#Angle polarizing filter makes with max linear transmission, Optional
                      'HFS_G' : [1011.910813e9] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                      'HFS_E': [120.527e6, 25e6], #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                      'HFSB_E': [0,25.88e6] ,
                      # 'PsiHFS': True, 

                      'I_spin' : 5/2, #Optional: nuclear spin, I, including parity. Tabulated by Stone et al.
                      'mu_I' :  1.35298, #Optional: Nuclear dipole moment as tabulated by Stone.   
                             }

InputRb_noHFS = {'s_ground': .5 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                      's_excited': .5 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                      'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
                      'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                      'E_ground': np.array([0])  ,#Lowest J first, in cm^-1
                      'E_excited':  np.array([12578.950,12816.545]) , #Lowest J first, in cm^-1
                      'Bmag': 0, #Magnetic Field  [T]
                      'b_angle': 90 , #Angle between LoS and Bmax
                      'spec_window' : [779.75,780.25], #Min and max function convolution
                      'plot_window' :  [779.9,780.15], #Min and max for plot window range (nm)

                      'amu' : 85.011 ,  #Weight in AMU
                      'specstep' : 0.005 , #Resolution of spectrometer in nm
                      
                      'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'.
                      'Temp' :300, #Temperature in K. Used for Gaussian convolution
                      'plottitle': f'Rb Testing', #Title for plotting (optional)     
                      'fxnwindow': .5 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                      # 'Pol_angle': Polangle ,#Angle polarizing filter makes with max linear transmission, Optional
                      # 'HFS_G' : [1011.910813e9] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                      # 'HFS_E': [120.527e6, 25e6], #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                      # 'HFSB_E': [0,25.88e6] ,

                      # 'I_spin' : 5/2, #Optional: nuclear spin, I, including parity. Tabulated by Stone et al.
                      # 'mu_I' :  1.35298, #Optional: Nuclear dipole moment as tabulated by Stone.   
                             }


InputRb2 = {'s_ground': .5 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                      's_excited': .5 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                      'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
                      'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                      'E_ground': np.array([0])  ,#Lowest J first, in cm^-1
                      'E_excited':  np.array([12578.950,12816.545]) , #Lowest J first, in cm^-1
                      'Bmag': 2, #Magnetic Field  [T]
                      'b_angle': 90 , #Angle between LoS and Bmax
                      'spec_window' : [777.75,782.25], #Min and max function convolution
                      'plot_window' :  [779.9,780.15], #Min and max for plot window range (nm)

                      'amu' : 85.011 ,  #Weight in AMU
                      'specstep' : 0.005 , #Resolution of spectrometer in nm
                      'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'.
                      'Temp' : 300, #Temperature in K. Used for Gaussian convolution
                      'plottitle': f'Rb Testing', #Title for plotting (optional)     
                      'fxnwindow': .05 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                      # 'Pol_angle': Polangle ,#Angle polarizing filter makes with max linear transmission, Optional
                      'HFS_G' : [1011.910813e9] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                      'HFS_E': [120.527e6, 25e6], #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                      'HFSB_E': [0,25.88e6] ,

                      'I_spin' : 5/2, #Optional: nuclear spin, I, including parity. Tabulated by Stone et al.
                      'mu_I' :  1.35298, #Optional: Nuclear dipole moment as tabulated by Stone.   
                             }

InputRb_noHFS2 = {'s_ground': .5 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                      's_excited': .5 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                      'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
                      'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                      'E_ground': np.array([0])  ,#Lowest J first, in cm^-1
                      'E_excited':  np.array([12578.950,12816.545]) , #Lowest J first, in cm^-1
                      'Bmag': 2, #Magnetic Field  [T]
                      'b_angle': 90 , #Angle between LoS and Bmax
                      'spec_window' : [777.75,782.25], #Min and max function convolution
                      'plot_window' :  [779.9,780.15], #Min and max for plot window range (nm)

                      'amu' : 85.011 ,  #Weight in AMU
                      'specstep' : 0.005 , #Resolution of spectrometer in nm
                      'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'.
                      'Temp' : 300, #Temperature in K. Used for Gaussian convolution
                      'plottitle': f'Rb Testing', #Title for plotting (optional)     
                      'fxnwindow': .05 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                      # 'Pol_angle': Polangle ,#Angle polarizing filter makes with max linear transmission, Optional
                      # 'HFS_G' : [1011.910813e9] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                      # 'HFS_E': [120.527e6, 25e6], #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                      # 'HFSB_E': [0,25.88e6] ,

                      # 'I_spin' : 5/2, #Optional: nuclear spin, I, including parity. Tabulated by Stone et al.
                      # 'mu_I' :  1.35298, #Optional: Nuclear dipole moment as tabulated by Stone.   
                             }
# MakeSlider(InputRb,[0,0],[779.75,780.5],do_sticks = 1,do_Polplot=1,polvary=0,banglevary=0)


plt.close('all')
Rb_PlotVars = [['green', 1.5, '-.','B = 0T - HFS' ,'x']]
Rb_PlotVars.append(['red', 1.5, ':',f'B = 0T - No HFS' ])


Rb_PlotVars2=[['blue', 1.5, '-.',f'B = 2T - HFS' ,'P']]
Rb_PlotVars2.append(['maroon', 1.5, ':',f'B = 2T - No HFS'  , 's'])


Rb_I = AtomSpect_Main(InputRb)
Rb_noHFS = AtomSpect_Main(InputRb_noHFS)
Rb_I2 = AtomSpect_Main(InputRb2)
Rb_noHFS2 = AtomSpect_Main(InputRb_noHFS2)


#This plots my results and sticks, gaussian convolution, including the ion velocity.

thixas = PlotFunction(Rb_I,Rb_PlotVars[0],plotwind = InputRb['plot_window'],
             # plottitle="Rb I",
             NormalizeSig = True ,makefig=True , plotpol = False,marker_count=10)
thixas[0][0].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.3f}"))
PlotFunction(Rb_noHFS,Rb_PlotVars[1],axsin = thixas,plotwind=InputRb['plot_window'] ,NormalizeSig=True,makefig=0,plotpol=False
             ,marker_count = 20)  
# thixas[0][0].text(0.05, 0.96, r'W$_i$ = 0.005 nm, B = 0T', verticalalignment='top',
#                       horizontalalignment='left', transform=thixas[0][0].transAxes)

# if savefigs:
#     plt.savefig('RbSpec.pdf', dpi = savdpi)

#%%
# thixas = PlotFunction(Rb_I2,Rb_PlotVars2[0],plotwind = InputRb2['plot_window'],
#              # plottitle="Rb I",
#              NormalizeSig = True ,makefig=True , plotpol = False, marker_count=5)
PlotFunction(Rb_I2,Rb_PlotVars2[0], axsin = thixas[0][0],plotwind =InputRb2['plot_window'] ,NormalizeSig=1,
             makefig=0,plotpol=False,marker_count=15)  
PlotFunction(Rb_noHFS2,Rb_PlotVars2[1], axsin = thixas[0][0],plotwind =InputRb2['plot_window'] ,NormalizeSig=1,
             makefig=0,plotpol=False,marker_count=15 ,legcols =2)  
plt.tight_layout(rect=[0, 0, 1, 0.9])

# thixas[0][0].text(0.05, 0.96, r'W$_i$ = 0.005 nm, B = 2T', verticalalignment='top',
#                       horizontalalignment='left', transform=thixas[0][0].transAxes)
if savefigs:
    plt.savefig('RbSpec0T2T.pdf', dpi = savdpi)
#%%   
if doRBfans:
    Bfan2 = np.linspace(0,.025,15000)

    
    Bfan = np.linspace(0,1,5000)
    InputdeckRbHFS_fan = {'s_ground':.5 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                      's_excited':.5 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                      'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
                      'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                      'E_ground': np.array([0]) ,#Lowest J first, in cm^-1
                      'E_excited': np.array([12578.95,12816.545])  , #Lowest J first, in cm^-1
                      'Bmag': Bfan, #Magnetic Field  [T]
                      'b_angle': 90 , #Angle between LoS and Bmax
                      'plottitle':'Rb85, I = 5/2, S = 1/2, L = 1' , #Title for plotting (optional)     
                      'plot_window' : [775,799] , #Min and max for plot window range (nm)
                      'amu' : 85 ,  #Weight in AMU
                      'Temp' : 300,
                      'HFS_G' : [1011e6]  , #Optional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                      'HFS_E': [120.527e6, 25.0354e6], #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                      'HFSB_E': [0,25.898e6] ,
                      'I_spin' : 5/2, #Optional: nuclear spin, I, including parity. Tabulated by Stone et al.
                      'mu_I' :  1.35298, #Optional: Nuclear dipole moment as tabulated by Stone.             
                       }

    # InputRb['Bmag'] = Bfan
    RBFanG = plotZfan(InputdeckRbHFS_fan,flags='G', markercount = 0)
    plt.xlim(0,.5)
    plt.ylim(-.4,.4)
    plt.tight_layout()
    plt.show()
    if savefigs:
        plt.savefig('RbFan_Ground.pdf', dpi=savdpi)
  


  
    RBFanE = plotZfan(InputdeckRbHFS_fan,flags='E',markercount = 0)
    plt.ylim(12816.52,12816.57)
    plt.xlim(0,0.015)
    plt.tight_layout()
    if savefigs:
        plt.savefig('RbFan_Excited.pdf', dpi=savdpi)
    plt.show()

#%%
doRbSlider = 0
if doRbSlider:
    
    InputRb2 = {'s_ground': .5 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                          's_excited': .5 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                          'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
                          'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                          'E_ground': np.array([0])  ,#Lowest J first, in cm^-1
                          'E_excited':  np.array([12578.950,12816.545]) , #Lowest J first, in cm^-1
                          'Bmag': 2, #Magnetic Field  [T]
                          'b_angle': 45 , #Angle between LoS and Bmax
                          'spec_window' : [779.82,780.535], #Min and max for plot window range (nm)
                          'plot_window' : [779.82,780.535], #Min and max for plot window range (nm)

                          'amu' : 85.011 ,  #Weight in AMU
                          'specstep' : 0.00001 , #Resolution of spectrometer in nm
                          'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'.
                          'Temp' : .1, #Temperature in K. Used for Gaussian convolution
                          'plottitle': f'Rb Testing With HFS', #Title for plotting (optional)     
                          'fxnwindow': .05 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                          # 'Pol_angle': 90 ,#Angle polarizing filter makes with max linear transmission, Optional
                          'HFS_G' : [1011.910813e9] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                          'HFS_E': [120.527e6, 25e6], #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                          'HFSB_E': [0,25.88e6] ,

                          'I_spin' : 5/2, #Optional: nuclear spin, I, including parity. Tabulated by Stone et al.
                          'mu_I' :  1.35298, #Optional: Nuclear dipole moment as tabulated by Stone.   
                                 }

    InputRb_noHFS2 = {'s_ground': .5 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                          's_excited': .5 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                          'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
                          'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                          'E_ground': np.array([0])  ,#Lowest J first, in cm^-1
                          'E_excited':  np.array([12578.950,12816.545]) , #Lowest J first, in cm^-1
                          'Bmag': 2, #Magnetic Field  [T]
                          'b_angle': 90 , #Angle between LoS and Bmax
                          'spec_window' : [779.82,780.535], #Min and max for plot window range (nm)
                          'plot_window' : [779.82,780.535], #Min and max for plot window range (nm)

                          'amu' : 85.011 ,  #Weight in AMU
                          'specstep' : 0.00001 , #Resolution of spectrometer in nm
                          'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'.
                          'Temp' : .1, #Temperature in K. Used for Gaussian convolution
                          'plottitle': f'Rb Testing - No HFS', #Title for plotting (optional)     
                          'fxnwindow': .05 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                          'Pol_angle': 90 ,#Angle polarizing filter makes with max linear transmission, Optional
                          # 'HFS_G' : [1011.910813e9] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                          # 'HFS_E': [120.527e6, 25e6], #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                          # 'HFSB_E': [0,25.88e6] ,

                          # 'I_spin' : 5/2, #Optional: nuclear spin, I, including parity. Tabulated by Stone et al.
                          # 'mu_I' :  1.35298, #Optional: Nuclear dipole moment as tabulated by Stone.   
                                 }
    
    # plt.close('all')
    MakeSlider(InputRb2,[0,0] , [779.02,781.135],do_sticks=0,do_Polplot=0,banglevary=0,polvary=0)
    MakeSlider(InputRb_noHFS2,[0,0] , [779.02,781.135],do_sticks=0,do_Polplot=0,banglevary=0,polvary=0)
