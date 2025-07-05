# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 18:11:16 2025

@author: Leo Nofs
"""

import csv, os, sys
import numpy as np
from matplotlib import pyplot as plt

#Goes up two folder levels
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir, os.path.pardir)))

#Goes up three folder levels and adds it to the path of the current file.
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir, os.path.pardir, os.path.pardir))) #You can keep adding os.path.pardir here to go up more folders.


from SpectSplit  import Zeeman_Main, Normalize,PlotFunction,read_Spectra, plotZfan,ZeemanFan
import SpectSplit as ss

def Raising(J,mj,states=[]):
    #I want to feed it a matrix to do this.
    ndim=3
    if mj<=J-1.:
        tjp = np.sqrt(J*(J+1) - mj*(mj+1))
        tmat = np.eye(ndim, k=1)
    else:
        tjp = 0
    # print(tjp)
    return tjp


def Lowering(J,mj,states=[]):
    #I want to feed it a matrix to do this.
    ndim=3
    if mj<=J+1.:
        tjp = np.sqrt(J*(J+1) - mj*(mj-1))
        tmat = np.eye(ndim, k=-1)
    else:
        tjp = 0
    # print(tjp)
    return tjp


def CouplingJJ(J1,mj1,J2,mj2):
    tjback = Raising(J1,mj1)*Lowering(J2,mj2) + Lowering(J1,mj1)*Raising(J2,mj2) + 2*mj1*mj2
    return tjback

def CouplingJJ2(J1,mj1,J2,mj2,J3,mj3):
    
    tjback = (Raising(J1,mj1) +Raising(J2,mj2) )*Lowering(J3,mj3) +( Lowering(J1,mj1) +Lowering(J2,mj2))*Raising(J3,mj3) + 2*(mj1+mj2)*mj3
    return tjback
#%%Quickly testing that the form of the HF interaction should match that of eq 5.14 in Townsend
HFS = {}
HFS['HFSA'] = [1000e6]
# HFS['HFSB'] = [1]
E_level = np.array([1])
I_state =.5
L_state = 0
S_state = .5
Jvals = np.arange(np.abs(L_state-S_state), (L_state+S_state+.1))



# HFSConst = [list(HFS['HFSA'])*int(2*Jvals[x] +1)*int(2*I_state + 1) for x in range(len(E_level))]
# HFS_Flat = [item for row in HFSConst for item in row]

# HFSConstB = [list(HFS['HFSB'])*int(2*Jvals[x] +1)*int(2*I_state + 1) for x in range(len(E_level))]
# HFS_FlatB = [item for row in HFSConstB for item in row]




SG = ss.stateprocess_HFS(L_state,S_state,I_state,HFS,1,E_level,1)

SZ = ss.HZeeman_HFS(SG,1)
#%%
Level_in=SG
Ldim=len(Level_in[-1]) #Dimension of the atomic level: Length of the # of possible states.
Psi=Level_in[0]           
PsiR=Level_in[0].T           
PsiL = np.zeros_like(Psi)
PsiS = np.zeros_like(Psi)
PsiS2 = np.zeros_like(Psi)

PsiI = np.zeros_like(Psi)
PsiH = np.zeros_like(Psi)

E0=np.eye(Ldim)  
#Create diagonal elements for the m values for each qm.     
tmLMat = np.diagflat([x[7] for x in Level_in[-1]])
tmSMat = np.diagflat([x[9] for x in Level_in[-1]])
tmIMat = np.diagflat([x[5] for x in Level_in[-1]])


for i in range(Ldim):
    tPsiL = Psi.T[i]            
    for j  in range(Ldim):
        # PsiL[j,i] = Psi[j](tmLMat@PsiR[j])     
        PsiS[i,j] = Psi[j]@tmSMat@PsiR[i] 
        PsiS2[i,j] =  np.linalg.multi_dot([tPsiL, tmSMat, Psi[j]]) 
        PsiI[i,j] = tPsiL@(tmIMat@PsiR[j])
#H_HFS = AIJcos(IJ) = AIJ*(F(F+1) - I(I+1) - J(J+1))/(IJ)

# tH_HFS = np.diagflat([x for x in Level_in[3]])

tH_HFS = np.diagflat([2*x for x in Level_in[2]])
# tH_HFS = np.diagflat([x for x in Level_in[3]])


muBEV = 5.7883818e-5 #bohr magneton in eV/Tesla
muBcm = muBEV/(1.23981e-4) #bohr magneton in cm^-1

#This is equivalent to the <n|Jz|m> calculation to make the pertubation matrix of zeeman energies    
#The += doesn't seem to change anything, but is kept as my original code used it
#HFS should be treated in the same way as Zeeman, with the matrix multiplication of the wave fxns.
for i in range(Ldim):

    for j  in range(Ldim):      
        PsiH[j,i] = PsiR[i]@(tH_HFS@PsiR[j]) 
        # PsiH[i,j] = np.linalg.multi_dot([PsiL[i], PsiHFS, Psi[j]])  

Htot = PsiL*1 + PsiS*2 - PsiI*1
Etot = np.asmatrix(E0 +PsiH + 1*Htot*(muBcm) ) #
# np.linalg.eigvalsh(Etot)
np.linalg.eigvalsh(PsiH)

# print(PsiL)
# # print(PsiS1)
# # print(PsiI1)
# Psi=Level_in[0]          
# PsiR=Level_in[0] 
# PsiL2 = np.zeros_like(Psi)
# PsiS = np.zeros_like(Psi)
# PsiI = np.zeros_like(Psi)
# PsiH = np.zeros_like(Psi)
# for i in range(Ldim):
                
#     for j  in range(Ldim):
#         PsiL2[i,j] = Psi[i]@(tmLMat@PsiR[j])     
#         PsiS[i,j] = Psi[i]@(tmSMat@PsiR[j])       
#         PsiI[i,j] = Psi[i]@(tmIMat@PsiR[j])       
#         PsiH[i,j] = Psi[i]@(tH_HFS@PsiR[j]) 
        
# print(PsiL2)

# print(PsiS)
# print(PsiI)

#%%
Ldim=len(SG[-1]) #Dimension of the atomic level: Length of the # of possible states.
SL = SG[-1]
Psi = SG[0]
PsiL = SG[0].T
# tHHFS = np.diagflat([x for x in SG[3] ])
Atest = [100e6]
tHHFS = np.diagflat([2*x for x in np.array(SG[2])])
PsiHFS = np.zeros_like(Psi)
PsiHFS2 = np.zeros_like(Psi)
PsiHFS3 = np.zeros_like(Psi)

hhev = 4.135668e-15 #Plank's constant in Ev/Hz

    #States are from bigsortedlist = [F,mF,J,mJ,I,mI,L,mL,S,mS]
                                    #[0,1 ,2,3 ,4,5 ,6,7 ,8,9 ]

for i in range(Ldim):
    tPsiL = Psi[i].T            
    for j in range(Ldim):
        # PsiHFS2[j,i] = i + 2*j
        
        # print(SL[j],i,j)
        # T1 =  CouplingJJ2(SL[j][6],SL[j][7],SL[j][8],SL[j][9],SL[i][4],SL[i][5])*Atest*(hhev/1.24e-4) 
        # print(SL[j][6],SL[j][7],SL[j][8],SL[j][9],SL[i][4],SL[i][5],i,j)
        
        # PsiHFS3[j,i] = CouplingJJ2(SL[j][6],SL[j][7],SL[j][8],SL[j][9],SL[i][4],SL[i][5])*Atest*(hhev/1.24e-4) 

        # PsiHFS[i,j] = CouplingJJ(SL[i][2],SL[i][3],SL[j][4],SL[j][5])*Atest*(hhev/1.24e-4) 
        # PsiHFS2[j,i] =PsiL[i]@(PsiHFS3@Psi[j]) 
        # PsiHFS2[i,j] =PsiL[i]@(tHHFS@Psi[j]) 

        PsiHFS2[j,i] =np.linalg.multi_dot([tPsiL, tHHFS, Psi[j]])  
        
ZPH2 = PsiHFS2
ZPH1 = SZ['PsiHFS']   
ZH0 = np.diagflat(SG[2]     )

# print(ZPH2)

import scipy as scp
# scp.linalg.eigvalsh(ZPH2,driver='evx',lower=True)
np.linalg.eigvalsh(ZPH2)

#%%
testcases= [True,True,False,True]
closeplots = False
closeallplots = True
if closeallplots:
    plt.close('all')

#%%#This is to match Bransden fig 5.12 (or attempt to)    Hydrogen 2p12, 2p32

if testcases[0]:
    print('1')
    Bfan = np.linspace(0,5,500)

    if closeplots:
        plt.close('all')
    Inputdeck3 =           {'s_ground':.5,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                              's_excited':.5 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                              'l_ground': 1 , #Orbital Angular Momentum of ground state, int or half int
                              'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                              'E_ground':  np.array([100,100.1]) ,#Lowest J first, in cm^-1
                              'E_excited':  np.array([82258.919113,82259.2850014]) , #Lowest J first, in cm^-1
                              'Bmag': Bfan , #Magnetic Field  [T]
                              'b_angle': 90 , #Angle between LoS and Bmax
                              'plottitle': f'Non-HFS, s = .5, L = 1 - Bransden fig 5.12' , #Title for plotting (optional)     
                              'plot_window' : [490,510], #Min and max for plot window range (nm)
                              # 'amu' : 16 ,  #Weight in AMU
                              # 'specstep' : 0.035 , #Resolution of spectrometer in nm
                              # 'Convfxn': 'GaussianInstrum', #Optional: 'Gaussian', 'Skewed'.
                              # 'Temp' :5*11604, #Temperature in K. Used for Gaussian convolution
                              # 'fxnwindow': 1 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                              # 'ion_vel' : -2500, #Ion velocity in m/s
}
    #Interesting crossover shifting at higher field between low and intermediate approximation.
    Inputdeck3Z = plotZfan(Inputdeck3,flags='EL')
    #%%
#This is to match Bransden fig 8.6    Arb atom


if testcases[1]:
    if closeplots:
        plt.close('all')    # Bfan = [0]*50
    Bfan = np.linspace(0,.1,100)
    # Bfan = [0,0,0]
    Inputdeck2 = {'s_ground':.5 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                  's_excited':0 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                  'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
                  'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                  'E_ground': np.array([1000]) ,#Lowest J first, in cm^-1
                  'E_excited': np.array([100])  , #Lowest J first, in cm^-1
                  'Bmag': Bfan, #Magnetic Field  [T]
                  'b_angle': 90 , #Angle between LoS and Bmax
                  'plottitle':'HFS , I = 1/2, J=1, L = 1, S=0, Bransden fig 8.6 (Arb atom)' , #Title for plotting (optional)     
                  'plot_window' : [775,799] , #Min and max for plot window range (nm)
                  'amu' : 85 ,  #Weight in AMU
                  'specstep' : 0.002222 , #Resolution of spectrometer in nm
                  'Convfxn' : 'Skewed' ,
                  'Temp' : 300,
                  # 'HFS_G' : [0]  , #Optional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                  # 'HFS_E': [0,0], #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                  'HFS_G' : [1e9]  , #Optional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                  'HFS_E': [120.527e6], #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                  # 'HFSB_E': [0,25.88e6] ,
                  'I_spin' : .5, #Optional: nuclear spin, I, including parity. Tabulated by Stone et al.
                  'mu_I' :  1.35298, #Optional: Nuclear dipole moment as tabulated by Stone.     
                  'PsiHFS': 'Yes', 

                   }
    
    Inputdeck2Z = plotZfan(Inputdeck2,flags='EL')
    # Inputdeck2Z = plotZfan(Inputdeck2,flags='GL')

    # plt.ylim(12816.,12816.57)
    # plt.xlim(0,1)
    # plt.figure()
    # for things in Inputdeck2Z['E']:
        
    #     plt.plot([0,.5,5],things)
    #%%Fuller Figure 3 (pg 904)
if testcases[2]:
    if closeplots:
        plt.close('all')    # Bfan = [0]*50
    Bfan = np.linspace(0,1,500)
    Fuller = {'s_ground':.5 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                  's_excited':.5 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                  'l_ground': 0 , #Orbital Angular Momentum of ground state, int or half int
                  'l_excited': 1 ,#Orbital Angular Momentum of excited state, int or half int
                  'E_ground': np.array([1000]) ,#Lowest J first, in cm^-1
                  'E_excited': np.array([12578.95,12816.545])  , #Lowest J first, in cm^-1
                  'Bmag': Bfan, #Magnetic Field  [T]
                  'b_angle': 90 , #Angle between LoS and Bmax
                  'plottitle':'Fuller Fig3, I = 3/2, S = 1/2, L = 1' , #Title for plotting (optional)     
                  'plot_window' : [775,799] , #Min and max for plot window range (nm)
                  'amu' : 85 ,  #Weight in AMU
                  # 'specstep' : 0.002222 , #Resolution of spectrometer in nm
                  # 'Convfxn' : 'Skewed' ,
                  'Temp' : 300,
                  # 'HFS_G' : [0]  , #Optional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                  # 'HFS_E': [0,0], #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                  'HFS_G' : [1e9]  , #Optional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                  'HFS_E': [120.527e6, 25e6], #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                  # 'HFSB_E': [0,25.88e6] ,
                  'I_spin' : 3/2, #Optional: nuclear spin, I, including parity. Tabulated by Stone et al.
                  'mu_I' :  1.35298, #Optional: Nuclear dipole moment as tabulated by Stone.             
                   }
    
    # RBFan = plotZfan(InputdeckRbHFS_fan,flags='GL')
    FullerFan = plotZfan(Fuller,flags='E')

    plt.ylim(12816.52,12816.57)
    plt.xlim(0,0.01)
  #%%  Actually Rubidium
if testcases[3]:
    if closeplots:
        plt.close('all')    # Bfan = [0]*50
    Bfan = np.linspace(0,1,500)
    # Bfan = [0,0,0]

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
                  # 'specstep' : 0.002222 , #Resolution of spectrometer in nm
                  # 'Convfxn' : 'Skewed' ,
                  'Temp' : 300,
                  # 'HFS_G' : [0]  , #Optional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                  # 'HFS_E': [0,0], #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                  'HFS_G' : [1011e6]  , #Optional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                  'HFS_E': [120.527e6, 25e6], #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                  # 'HFSB_E': [0,25.88e6] ,
                  'I_spin' : 5/2, #Optional: nuclear spin, I, including parity. Tabulated by Stone et al.
                  'mu_I' :  1.35298, #Optional: Nuclear dipole moment as tabulated by Stone.             
                  # 'PsiHFS': 'Yes', 

                   }
    
    # TesttZ = ZeemanFan(InputdeckRbHFS_fan,flags='GL')
    RBFan = plotZfan(InputdeckRbHFS_fan,flags='EL')

    plt.ylim(12816.1,12816.87)
    # plt.ylim(12816.520470813763,12816.562848707925)
    ZZR=RBFan['EHZeem']['AHFS0']
    
    plt.xlim(0,0.35)
    # plt.figure()
    # plt.ylim(12816.1,12816.87)

    # for things in RBFan['E']:
        
    #     plt.plot([0,.5,5],things)
#%%This is to compare to Grumer per Jonsson for He(2020)

Bfanvals = np.linspace(0,1,200)
dofan= False
if dofan:
    plt.close('all')
    #This fan plot matches Semenov 2000 for He 1s2s
    InputdeckHeFan = {'s_ground': 1 ,#Spin multiplicity,s, for ground state ^(2s +1)L_J , int or half int
                  's_excited': 1 ,#Spin multiplicity,s, for excited state ^(2s +1)L_J, int or half int
                  'l_ground': 1 , #Orbital Angular Momentum of ground state, int or half int
                  'l_excited':2 ,#Orbital Angular Momentum of excited state, int or half int
                  'E_ground': np.array([169087.8308131,  169086.8428979 ,  169086.7664725  ]) ,#Lowest J first, in cm^-1
                  'E_excited': np.array([191444.49804915,191444.47952984,191444.47832914]) , #Lowest J first, in cm^-1
                  'Bmag': Bfanvals , #Magnetic Field  [T]
                  'b_angle': 90 , #Angle between LoS and Bmax
                  'plottitle': 'He I 1s2p  , UnPolarized'  , #Title for plotting (optional)     
                  'plot_window' : [447,447.2], #Min and max for plot window range (nm), convolutions will only take place within this range
                  'amu' : 3 , #Weight in AMU
                  # 'Convfxn': 'Gaussian', #Optional: 'Gaussian', 'Skewed'.
                  'specstep': 0.002222, #nm per pixel of camera
                  # 'Temp' : 300, #Temperature in K. Used for Gaussian convolution
                  'HFS_G' : [0,0,0]  , #Optional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                  'HFS_E': [0,0,0], #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                  # 'HFSB_E': [0,25.88e6] ,
                  'I_spin' : 0.5, #Optional: nuclear spin, I, including parity. Tabulated by Stone et al.
                  'mu_I' :-2.12749772, #Optional: Nuclear dipole moment as tabulated by Stone.          

                   }  
    
    #Zfan plots as figures# 200 and 300 or reasons.
    Hfan = plotZfan(InputdeckHeFan,dolow=True)