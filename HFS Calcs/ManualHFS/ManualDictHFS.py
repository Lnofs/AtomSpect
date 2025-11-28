# -*- coding: utf-8 -*-
"""
Created on June 30 2025
This file takes the manual dictionary tests for the hyperfine constant values and puts them in a separate script
@author: %Leo Nofs
"""


# from AtomSpect import *
from HFS_Const_Calcs import Ahfs
#%%ManualDictTests
domanualDicts = True
if domanualDicts:
    # Couple = ['ss','sl1','sl2','jj','None']
    
    # for cc in Couple:
    #     Template = {'Z' : 3   , #Atomic number
    #             'AMU': 7  , #Atomic mass in amu
    #             'E_ionize': 7.239561 , #Ionization energy [eV]
    #             'E1': 1.847846    , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
    #             'n1': 2     , #Princple quantum number.        
    #             'I' : 1.5    , #Nuclear spin angular momentum
    #             'muI': 3.256427 ,    #Nuclear dipole moment (in units of mu_N)
    #             'l1'  :  1   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
    #             's1' : 0.5 ,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".
    #             'j1'  : 1.5   , #J = L +/- S coupled angular momentum 
    #             'Corrflags': 'NoCorr' , #Correction factors for the hyperfine calculation.
    #             # IF coupling is to be done, need J,L,S, j2,l2,s2 for the different possible schemes. It's possible that j1 may not be always known (especially outside of S levels).
    #             # 'Coupling' : 'sl1',
    #             # 'l2'  :  0   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc for the secondary electron
    #             # 's2' : 1 ,  #Spin of the secondary electron
    #             # 'j2'  : 1   , #J = L +/- S coupled angular momentum for the secondary electron
                
    #             # 'L' : 1 , #Coupled total L from the two electrons
    #             # 'S' : .5 , #Coupled total S of the two electrons
    #             # 'J' : .5 , #Coupled total J for the two electrons
    #             'A_Exp': -3 , 
    #          }
    #     ManualTest = Ahfs(Template)
    #     print(ManualTest)
    
    #%%Li72s12
    
    Li72s12 = {'Z' : 3  , #Atomic number
            'AMU': 7  , #Atomic mass in amu
            'E_ionize': 5.391715 , #Ionization energy [eV]
            'E1':0 , #Energy of electron [eV]            
            'n1': 2     , #Princple quantum number        
            'I' : 1.5    , #Nuclear spin angular momentum
            'muI': 3.256427 ,    #Nuclear dipole moment (in units of mu_N)
            'l1'  : 0   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
            's1' : 0.5 ,  #Spin of the electron in question. 
            'j1'  : .5   , #J = L +/- S coupled angular momentum 
            'Corrflags': 'ds' , #Correction factors for the hyperfine calculation. Possible flags: 'F H R del eps ds'
            'E2':3.37313 , #Energy of the second level for dsdn correction (Only used for s states typically)
            'n2':3 ,  #Principle quantum number of the second level for dsdn correction (Only used for s states typically)
            'A_Exp': 401 , #Experimental A in MHZ for easy comparison            
         }
    
    Li72s12A = Ahfs(Li72s12)
    print(f'Li72s12 : Calc = {Li72s12A}, Exp = {Li72s12['A_Exp']}')
    #%%Li2p12
    Li72p12 = {'Z' : 3  , #Atomic number
            'AMU': 7  , #Atomic mass in amu
            'E_ionize': 5.391715 , #Ionization energy [eV]
            'E1':1.847846  , #Energy of electron [eV]            
            'n1': 2     , #Princple quantum number        
            'I' : 1.5    , #Nuclear spin angular momentum
            'muI': 3.256427 ,    #Nuclear dipole moment (in units of mu_N)
            'l1'  : 1   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
            's1' : -0.5 ,  #Spin of the electron in question. 
            'j1'  : .5   , #J = L +/- S coupled angular momentum 
            'Corrflags': 'FRHdelepsds' , #Correction factors for the hyperfine calculation. Possible flags: 'F H R del eps ds'
            'E2': 3.834258 , #Energy of the second level for dsdn correction (Only used for s states typically)
            'n2':3 ,  #Principle quantum number of the second level for dsdn correction (Only used for s states typically)
            'A_Exp': 45 , #Experimental A in MHZ for easy comparison            
         }
    
    Li72p12A = Ahfs(Li72p12)
    print(f'Li72p12 : Calc = {Li72p12A}, Exp = {Li72p12['A_Exp']}')
 

    #%%Li2p32
    Li72p32 = {'Z' : 3  , #Atomic number
            'AMU': 7  , #Atomic mass in amu
            'E_ionize': 5.391715 , #Ionization energy [eV]
            'E1':1.847846  , #Energy of electron [eV]            
            'n1': 2     , #Princple quantum number        
            'I' : 1.5    , #Nuclear spin angular momentum
            'muI': 3.256427 ,    #Nuclear dipole moment (in units of mu_N)
            'l1'  : 1   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
            's1' : 0.5 ,  #Spin of the electron in question. 
            'j1'  : 1.5   , #J = L +/- S coupled angular momentum 
            'Corrflags': 'FRHdelepsds' , #Correction factors for the hyperfine calculation. Possible flags: 'F H R del eps ds'
            'E2': 3.834258 , #Energy of the second level for dsdn correction (Only used for s states typically)
            'n2':3 ,  #Principle quantum number of the second level for dsdn correction (Only used for s states typically)
            'A_Exp': -3 , #Experimental A in MHZ for easy comparison            
         }
    
    Li72p32A = Ahfs(Li72p32)
    print(f'Li72p32 : Calc = {Li72p32A}, Exp = {Li72p32['A_Exp']}')
       
  
    
  # %%Li73p12


    Li73p12 = {'Z' : 3 , #Atomic number
            'AMU': 7  , #Atomic mass in amu
            'E_ionize':5.391715 , #Ionization energy [eV]
            'E1': 3.834258  , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
            'n1': 3     , #Princple quantum number.        
            'I' : 1.5    , #Nuclear spin angular momentum
            'muI': 3.256427 ,    #Nuclear dipole moment (in units of mu_N)
            'l1'  : 1   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
            's1' : -0.5 ,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".
            'j1'  :.5   , #J = L +/- S coupled angular momentum 
            'Corrflags': 'FRHdelepsds' , #Correction factors for the hyperfine calculation.
            'E2': 4.521648  ,#Energy of the second level for dsdn correction (Only used for s states typically)
            'n2':4 ,  #Principle quantum number of the second level for dsdn correction (Only used for s states typically)
            'A_Exp': 13.5,  #Experimental A in MHZ for easy comparison      
            
         }
    Li73p12A = Ahfs(Li73p12)
    print(f'Li73p12 : Calc = {Li73p12A}, Exp = {Li73p12['A_Exp']}')
# %%Li73p32



    Li73p32 = {'Z' : 3 , #Atomic number
            'AMU': 7  , #Atomic mass in amu
            'E_ionize':5.391715 , #Ionization energy [eV]
            'E1': 3.834258  , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
            'n1': 3     , #Princple quantum number.        
            'I' : 1.5    , #Nuclear spin angular momentum
            'muI': 3.256427 ,    #Nuclear dipole moment (in units of mu_N)
            'l1'  : 1   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
            's1' : 0.5 ,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".
            'j1'  :1.5   , #J = L +/- S coupled angular momentum 
            'Corrflags': 'FHdelepsds' , #Correction factors for the hyperfine calculation.
            'E2': 4.521648  ,#Energy of the second level for dsdn correction (Only used for s states typically)
            'n2':4 ,  #Principle quantum number of the second level for dsdn correction (Only used for s states typically)
            'A_Exp': -0.96,  #Experimental A in MHZ for easy comparison      
            
         }
    Li73p32A = Ahfs(Li73p32)
    print(f'Li73p32 : Calc = {Li73p32A}, Exp = {Li73p32['A_Exp']}')
  
    #%%Na3s12
    Na3s12 = {'Z' :  11, #Atomic number
            'AMU': 23  , #Atomic mass in amu
            'E_ionize':5.13907696 , #Ionization energy [eV]
    
            'E1': 0    , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
            'n1': 3     , #Princple quantum number.        
            'I' : 1.5    , #Nuclear spin angular momentum
            'muI':2.217522 ,    #Nuclear dipole moment (in units of mu_N)
            'l1'  :  0   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
            's1' : 0.5 ,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".
            'j1'  : .5   , #J = L +/- S coupled angular momentum 
            # 'Corrflags': '' , #Correction factors for the hyperfine calculation.
            # IF coupling is to be done, need J,L,S, j2,l2,s2 for the different possible schemes. It's possible that j1 may not be always known (especially outside of S levels).
            # 'Coupling' : 'sl1' ,
            'E2': 3.191353,
            # 'E2': 2.9504483,
    
    
            'n2':4 , 
    
            'A_Exp': 885,  #Came from Allegrini
            # 'dnstar' : True,
            # 'ds2' :  True

         }
    
    
    Na3s12A = Ahfs(Na3s12)
    print(f'Na3s12 : Calc = {Na3s12A}, Exp = {Na3s12['A_Exp']}')
    #%%Na3d52
    
    Na3d52 = {'Z' :  11, #Atomic number
            'AMU': 23  , #Atomic mass in amu
            # 'E_ionize': 7.239561 , #Ionization energy [eV]
            'E_ionize':5.13907696 , #Ionization energy [eV]
    
            'E1': 3.6169733    , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
            'n1': 3     , #Princple quantum number.        
            'I' : 1.5    , #Nuclear spin angular momentum
            'muI': 2.217522 ,    #Nuclear dipole moment (in units of mu_N)
            'l1'  :  2   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
            's1' : 0.5 ,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".
            'j1'  : 2.5   , #J = L +/- S coupled angular momentum 
            'Corrflags': 'RFHdelepsds' , #Correction factors for the hyperfine calculation.
            # IF coupling is to be done, need J,L,S, j2,l2,s2 for the different possible schemes. It's possible that j1 may not be always known (especially outside of S levels).
            # 'Coupling' : 'sl1' ,
            'E2': 4.283498,
            # 'E2': 2.9504483,
    
    
            'n2':4 , 
    
            'A_Exp': 0.108 ,  #Came from Allegrini
            # 'dnstar' : True,
            # 'ds2' :  True

         }
    
    
    Na3d52A = Ahfs(Na3d52)
    print(Na3d52A)
    #%%Na3d32
    
    Na3d32 = {'Z' :  11, #Atomic number
            'AMU': 23  , #Atomic mass in amu
            # 'E_ionize': 7.239561 , #Ionization energy [eV]
            'E_ionize':5.13907696 , #Ionization energy [eV]
    
            'E1': 3.6169733    , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
            'n1': 3     , #Princple quantum number.        
            'I' : 1.5    , #Nuclear spin angular momentum
            'muI': 2.217522 ,    #Nuclear dipole moment (in units of mu_N)
            'l1'  :  2   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
            's1' : -0.5 ,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".
            'j1'  : 1.5   , #J = L +/- S coupled angular momentum 
            'Corrflags': 'RFHdelepsds' , #Correction factors for the hyperfine calculation.
            # IF coupling is to be done, need J,L,S, j2,l2,s2 for the different possible schemes. It's possible that j1 may not be always known (especially outside of S levels).
            # 'Coupling' : 'sl1' ,
            'E2': 4.283498,
            # 'E2': 2.9504483,
    
    
            'n2':4 , 
    
            'A_Exp': 0.527, 
            # 'dnstar' : True,
            # 'ds2' :  True

         }
    
    # d = [[item for item in v if item] for k, v in Template.items()]
    
    Na3d32A = Ahfs(Na3d32)
    print(Na3d32A)

    
      #%%Na3p32
    
    Na3p32 = {'Z' : 11   , #Atomic number
            'AMU': 23  , #Atomic mass in amu
            # 'E_ionize': 7.239561 , #Ionization energy [eV]
            'E_ionize':5.13907696 , #Ionization energy [eV]
    
            'E1': 2.103719    , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
            'n1': 3     , #Princple quantum number.        
            'I' : 1.5    , #Nuclear spin angular momentum
            'muI': 2.217522 ,    #Nuclear dipole moment (in units of mu_N)
            'l1'  :  1   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
            's1' : 0.5 ,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".
            'j1'  : 1.5   , #J = L +/- S coupled angular momentum 
            'Corrflags': 'FHdelepsds' , #Correction factors for the hyperfine calculation.
            # IF coupling is to be done, need J,L,S, j2,l2,s2 for the different possible schemes. It's possible that j1 may not be always known (especially outside of S levels).
            # 'Coupling' : 'sl1' ,
            'E2': 3.753091,
            # 'E2': 2.9504483,
    
    
            'n2':4 , 
    
            'A_Exp': 18.8 , 
            # 'dnstar' : True,
            # 'ds2' :  True
    
         }
    
    # d = [[item for item in v if item] for k, v in Template.items()]
    
    Na3p32A = Ahfs(Na3p32)
    print(Na3p32A)
    
    
#%%Na3p12
    Na3p12 = {'Z' : 11   , #Atomic number
            'AMU': 23  , #Atomic mass in amu
            # 'E_ionize': 7.239561 , #Ionization energy [eV]
            'E_ionize':5.13907696 , #Ionization energy [eV]
    
            'E1': 2.103719    , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
            'n1': 3     , #Princple quantum number.        
            'I' : 1.5    , #Nuclear spin angular momentum
            'muI': 2.217522 ,    #Nuclear dipole moment (in units of mu_N)
            'l1'  :  1   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
            's1' : 0.5 ,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".
            'j1'  : .5   , #J = L +/- S coupled angular momentum 
            'Corrflags': 'RFHdelepsds' , #Correction factors for the hyperfine calculation.
            # IF coupling is to be done, need J,L,S, j2,l2,s2 for the different possible schemes. It's possible that j1 may not be always known (especially outside of S levels).
            # 'Coupling' : 'sl1' ,
            'E2': 3.753091,
            # 'E2': 2.9504483,
    
    
            'n2':4 , 
    
            'A_Exp': 94.44,  #If I calculate this the same way I do for the other l-1/2, I get a much different answer...
            
            # 'dnstar' : True,
            # 'ds2' :  True
    
         }
    
    # d = [[item for item in v if item] for k, v in Template.items()]
    
    Na3p12A = Ahfs(Na3p12)
    print(Na3p12A)



      
    #%%Rb85Corr
    #Looks like there may be something to 
    Rb5s12 = {'Z' : 37   , #Atomic number
            'AMU': 85  , #Atomic mass in amu
            'E_ionize': 4.177128 , #Ionization energy [eV]
            'E1': 0   , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
            'n1': 5     , #Princple quantum number.        
            'I' : 2.5    , #Nuclear spin angular momentum
            'muI': 1.35298 ,    #Nuclear dipole moment (in units of mu_N)

            # 'muI': 1.339964 ,    #Nuclear dipole moment (in units of mu_N)
            'l1'  :  0   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
            's1' : 0.5 ,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".
            'j1'  : .5   , #J = L +/- S coupled angular momentum 
            # 'E2' :  1.68 , #Optional. For s electrons, needed for qdefect correction factor

            'E2' :  2.496113 , #Optional. For s electrons, needed for qdefect correction factor
            'n2' :  6   , #Optional, n2>n1. Needed for qdefect correction factor of s orbitals.
            'Corrflags' : 'FHepsdelds' ,
            'A_Exp': 1011, 
            # 'dnstar' : True,
            # 'ds2' :  True
            
                }
    
    Rb5s12A = Ahfs(Rb5s12)
    print(Rb5s12A)
# %%
    

    Rb5p12 = {'Z' : 37   , #Atomic number
            'AMU': 85  , #Atomic mass in amu
            'E_ionize': 4.177128 , #Ionization energy [eV]
            # 'E1': 2.5979    , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.

            'E1': 1.57923     , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
            'n1': 5     , #Princple quantum number.        
            'I' : 2.5    , #Nuclear spin angular momentum
            'muI': 1.35298 ,    #Nuclear dipole moment (in units of mu_N)
            'l1'  : 1   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
            's1' : 0.5 ,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".    
            'j1'  : .5   , #J = L +/- S coupled angular momentum 
            'E2' :  2.946702 , #Optional. For s electrons, needed for qdefect correction factor
            # 'E2' :  1.367
         #Optional. For s electrons, needed for qdefect correction factor
            'n2' :  6    ,#Optional, n2>n1. Needed for qdefect correction factor of s orbitals.
            # 'Corrflags' : 'FHdelepsds' ,
            'A_Exp': 120, 
            # 'dnstar' : True,
            # 'ds2' :  True
                }
    
    Rb5p12A = Ahfs(Rb5p12)
    print(Rb5p12A)
    
# %%
    

    # #%%TemplateDict
    # #I would likely suggest that the hyperfine constants are calculated separately from the Zeeman calculations.
    # #For ease of adding additional coupling, I'm here taking that s1 is the "outer" electron as then s2 can be added and coupled to an already existing calculation.
    # #This prevents the need to re-label the input parameters if testing different coupling schemes.
    
    
    Rb5p32 = {'Z' : 37   , #Atomic number
            'AMU': 85  , #Atomic mass in amu
            'E_ionize': 4.177128 , #Ionization energy [eV]
            # 'E1': 2.5979    , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.

            'E1': 1.57923     , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
            'n1': 5     , #Princple quantum number.        
            'I' : 2.5    , #Nuclear spin angular momentum
            'muI': 1.35298 ,    #Nuclear dipole moment (in units of mu_N)
            'l1'  : 1   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
            's1' : 0.5 ,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".    
            'j1'  : 1.5   , #J = L +/- S coupled angular momentum 
            'E2' :  2.946702 , #Optional. For s electrons, needed for qdefect correction factor
            # 'E2' :  1.367  ,
            'n2' :  6    ,#Optional, n2>n1. Needed for qdefect correction factor of s orbitals.
            'Corrflags' : 'Hdelepsds' ,
            'A_Exp': 25.009, 
            # 'dnstar' : True,
            # 'ds2' :  True
                }
    RBp32 = Ahfs(Rb5p32)
    print(RBp32)
    # #%%
    # # #%%TestSr
    # Sr5s6s3S1= {'Z' : 38   , #Atomic number
    #         'AMU': 87  , #Atomic mass in amu
    #         'E_ionize': (5.694867) , #Ionization energy [eV]
    #         'E1': 3.600349   , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
    #         # 'E1': 0  , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
    #         'n2': 6   , #Princple quantum number.        
    #         'I' :   4.5  , #Nuclear spin angular momentum
    #         # 'muI': 0.11778476 ,    #Nuclear dipole moment (in units of mu_N)
    #         'muI': -1.0928 ,    #Nuclear dipole moment (in units of mu_N)
    
    #         'l1'  :  0   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
    #         's1' : .5 ,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".
    #         'j1'  : .5   , #J = L +/- S coupled angular momentum 
    #         'E2' :   4.6400683 , #Optional. For s electrons, needed for qdefect correction factor

    #         # 'E2' :   4.6400683 , #Optional. For s electrons, needed for qdefect correction factor
    #         # 'E2' :  5.310495, #Optional. For s electrons, needed for qdefect correction factor
    #         'n1' :  5   ,#Optional, n2>n1. Needed for qdefect correction factor of s orbitals.
    #         # 'E2':0 ,
    #         'Corrflags': 'FHdeldeps' , #Correction factors for the hyperfine calculation.
    #         #IF coupling is to be done, need J,L,S, j2,l2,s2 for the different possible schemes. It's possible that j1 may not be always known (especially outside of S levels).
    #         'Coupling': 'ss',
    #         'l2'  :  0   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc for the secondary electron
    #         's2' : .5 ,  #Spin of the secondary electron
    #         'j2'  : .5   , #J = L +/- S coupled angular momentum for the secondary electron
            
    #         'L' : 0 , #Coupled total L from the two electrons
    #         'S' : 1 , #Coupled total S of the two electrons
    #         'J' : 1 , #Coupled total J for the two electrons
    #         'A_Exp': -542, 
    #         # 'ds2' :  True

    #      }
    
    # Sr5s6s3S1Copule = Ahfs(Sr5s6s3S1)
    # print(Sr5s6s3S1Copule)
    
    # #%%TestSr3P1
    # TestSr= {'Z' : 38   , #Atomic number
    #         'AMU': 87  , #Atomic mass in amu
    #         'E_ionize': 5.694867 , #Ionization energy [eV]
    #         'E1': 1.822888   , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
    #         'n1': 5     , #Princple quantum number.        
    #         'I' :   4.5  , #Nuclear spin angular momentum
    #         # 'muI': 0.11778476 ,    #Nuclear dipole moment (in units of mu_N)
    #         'muI': -1.0928 ,    #Nuclear dipole moment (in units of mu_N)
    
    #         'l1'  :  1   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
    #         's1' : 0.5 ,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".
    #         'j1'  : .5   , #J = L +/- S coupled angular momentum 
      
    #         'E2' :  4.206147, #Optional. For s electrons, needed for qdefect correction factor
    #         # 'E2' :  5.310495, #Optional. For s electrons, needed for qdefect correction factor
    #         'n2' :  6    ,#Optional, n2>n1. Needed for qdefect correction factor of s orbitals.
            
    #         'Corrflags': 'FHdelepsds' , #Correction factors for the hyperfine calculation.
    #         #IF coupling is to be done, need J,L,S, j2,l2,s2 for the different possible schemes. It's possible that j1 may not be always known (especially outside of S levels).
    #         'Coupling': 'jj',
    #         'l2'  :  0   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc for the secondary electron
    #         's2' : .5 ,  #Spin of the secondary electron
    #         'j2'  : .5   , #J = L +/- S coupled angular momentum for the secondary electron
            
    #         'L' : 0 , #Coupled total L from the two electrons
    #         'S' : 1 , #Coupled total S of the two electrons
    #         'J' : 1 , #Coupled total J for the two electrons
    #         'A_Exp': -260,  #Gets -278
    #      }
    
    # TestSr = Ahfs_Corrected(TestSr)
    # print(TestSr)
    # # %%TestSr3P2
    # TestSr= {'Z' : 38   , #Atomic number
    #         'AMU': 87  , #Atomic mass in amu
    #         'E_ionize': 5.694867 , #Ionization energy [eV]
    #         'E1': 1.822888   , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
    #         'n1': 5     , #Princple quantum number.        
    #         'I' :   4.5  , #Nuclear spin angular momentum
    #         # 'muI': 0.11778476 ,    #Nuclear dipole moment (in units of mu_N)
    #         'muI': -1.0928 ,    #Nuclear dipole moment (in units of mu_N)
    
    #         'l1'  :  0   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
    #         's1' : 0.5 ,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".
    #         'j1'  : .5   , #J = L +/- S coupled angular momentum 
    #         # 'E2' :  0, #Optional. For s electrons, needed for qdefect correction factor
    
    #         'E2' :  4.206147, #Optional. For s electrons, needed for qdefect correction factor
    #         # 'E2' :  5.310495, #Optional. For s electrons, needed for qdefect correction factor
    #         'n2' :  6    ,#Optional, n2>n1. Needed for qdefect correction factor of s orbitals.
            
    #         'Corrflags': 'delepsHds' , #Correction factors for the hyperfine calculation.
    #         #IF coupling is to be done, need J,L,S, j2,l2,s2 for the different possible schemes. It's possible that j1 may not be always known (especially outside of S levels).
    #         'Coupling': 'jj',
    #         'l2'  :  1   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc for the secondary electron
    #         's2' : .5 ,  #Spin of the secondary electron
    #         'j2'  :  1.5   , #J = L +/- S coupled angular momentum for the secondary electron
            
    #         'L' : 0 , #Coupled total L from the two electrons
    #         'S' : 1 , #Coupled total S of the two electrons
    #         'J' : 2 , #Coupled total J for the two electrons
    #         'A_Exp': -212, # Gets --214 if s electron is done as first. Definitely need to improve the a-calcs for hyperfine.
    #      }
    
    # TestSr = Ahfs_Corrected(TestSr)
    # print(TestSr)
    # #%%
    # # #%%MoreWTest
    # # #Testing [](^6S)6s (7S_3) state. s1 is the 6s electron, s2 is the (^6S) coupled electron
    # W7S3 = {'Z' : 74   , #Atomic number
    #         'AMU': 183  , #Atomic mass in amu
    #         'E_ionize': 7.86403 , #Ionization energy [eV]
    #         'E1': -.186672     , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
    #         'n1': 6     , #Princple quantum number.        
    #         'I' : 0.5    , #Nuclear spin angular momentum
    #         # 'muI': 0.11778476 ,    #Nuclear dipole moment (in units of mu_N)
    #         'muI': 0.1095569 ,    #Nuclear dipole moment (in units of mu_N)
    
    #         'l1'  :  0   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
    #         's1' : 2.5 ,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".
    #         'j1'  : 2.5   , #J = L +/- S coupled angular momentum 
      
    #         # 'E2' :  5.1238, #Optional. For s electrons, needed for qdefect correction factor
    #         'E2' :  5.310495, #Optional. For s electrons, needed for qdefect correction factor
    #         'n2' :  7    ,#Optional, n2>n1. Needed for qdefect correction factor of s orbitals.
            
    #         'Corrflags': 'Fdelepsds' , #Correction factors for the hyperfine calculation.
    #         #IF coupling is to be done, need J,L,S, j2,l2,s2 for the different possible schemes. It's possible that j1 may not be always known (especially outside of S levels).
    #         'Coupling': 'ss',
    #         'l2'  :  0   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc for the secondary electron
    #         's2' : .5 ,  #Spin of the secondary electron
    #         'j2'  : .5   , #J = L +/- S coupled angular momentum for the secondary electron
            
    #         'L' : 0 , #Coupled total L from the two electrons
    #         'S' : 3 , #Coupled total S of the two electrons
    #         'J' : 3 , #Coupled total J for the two electrons
    #         'A_Exp': 500, 
    #      }
    #%%
    # W7Scouple = Ahfs_Corrected(W7S3)
    # print(W7Scouple)
    #%%
    W7SOut = {'Z' : 74   , #Atomic number
            'AMU': 183  , #Atomic mass in amu
            'E_ionize': 7.86403 , #Ionization energy [eV]
            'E1': -.18667     , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
            'n1': 6     , #Princple quantum number.        
            'I' : 0.5    , #Nuclear spin angular momentum
            # 'muI': 0.11778476 ,    #Nuclear dipole moment (in units of mu_N)
            'muI': 0.1095569 ,    #Nuclear dipole moment (in units of mu_N)
    
            'l1'  :  0   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
            's1' :  .5 ,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".
            'j1'  : .5   , #J = L +/- S coupled angular momentum 
      
           'E2' :5.310495 ,
           'n2' : 7 ,
            'Corrflags': 'Fdelepsds' , #Correction factors for the hyperfine calculation.
            #IF coupling is to be done, need J,L,S, j2,l2,s2 for the different possible schemes. It's possible that j1 may not be always known (especially outside of S levels).
            # 'Coupling': 'jj',
            # 'l2'  :  0   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc for the secondary electron
            # 's2' : 2.5 ,  #Spin of the secondary electron
            # 'j2'  : 2.5   , #J = L +/- S coupled angular momentum for the secondary electron
            
            # 'L' : 1 , #Coupled total L from the two electrons
            # 'S' : 3 , #Coupled total S of the two electrons
            # 'J' : 3 , #Coupled total J for the two electrons
            'A_Exp': 496, 
    
         }
    
    W7SOutA = Ahfs(W7SOut)
    print(W7SOutA)
   #%%
    W7Sin = {'Z' : 74   , #Atomic number
            'AMU': 183  , #Atomic mass in amu
            # 'E1': -.18667     , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
            # 'E_ionize': 7.86403 , #Ionization energy [eV]

            'E_ionize': 16.37 , #Ionization energy [eV]
            'E1': .40624084     , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
            'n1': 5    , #Princple quantum number.        
            'I' : 0.5    , #Nuclear spin angular momentum
            # 'muI': 0.11778476 ,    #Nuclear dipole moment (in units of mu_N)
            'muI': 0.1095569 ,    #Nuclear dipole moment (in units of mu_N)
    
            'l1'  :  0   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
            's1' : 2.5,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".
            'j1'  : 2.5   , #J = L +/- S coupled angular momentum 
      
          
            'Corrflags': 'FHRdeleps' , #Correction factors for the hyperfine calculation.
            #IF coupling is to be done, need J,L,S, j2,l2,s2 for the different possible schemes. It's possible that j1 may not be always known (especially outside of S levels).
            # 'Coupling': 'jj',
            # 'E2': -.18667     , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.

            # 'E2' :5.310495 ,
            # 'n2' : 7,

            # 'l2'  :  0   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc for the secondary electron
            # 's2' : 2.5 ,  #Spin of the secondary electron
            # 'j2'  : 2.5   , #J = L +/- S coupled angular momentum for the secondary electron
            
            # 'L' : 1 , #Coupled total L from the two electrons
            # 'S' : 3 , #Coupled total S of the two electrons
            # 'J' : 3 , #Coupled total J for the two electrons
            'A_Exp': 496, 
    
         }
    
    W7SinA = Ahfs(W7Sin)
    print(W7SinA)
    #%%
def GenCouple(A,B,C,D):
    #Couples (A(A+1) + B(B+1) - C(C+1)) /(2D(D+1)) since it is such a common calculation, esp for testing.
    return (A*(A+1) + B*(B+1) - C*(C+1)) /(2*D*(D+1))
C1 = GenCouple(3,2.5,.5,3)
C2 = GenCouple(3,.5,2.5,3)

print(C1, C1*W7SinA['A'])
print(C2, C2*W7SOutA['A'])

A2 = (C1*W7SinA['A'] + C2*W7SOutA['A'])
print(A2)
# C2 = (1-C1)
#%%
    #%%
    # W7P4 = {'Z' : 74   , #Atomic number
    #         'AMU': 183  , #Atomic mass in amu
    #         'E_ionize': 7.86403 , #Ionization energy [eV]
    #         # 'E2': 0     , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
    #         # 'n2': 7     , #Princple quantum number.        
    #         'I' : 0.5    , #Nuclear spin angular momentum
    #         'muI': 0.11778476 ,    #Nuclear dipole moment (in units of mu_N)
    #         # 'muI': 0.1095569 ,    #Nuclear dipole moment (in units of mu_N)
    
    #         'l1'  :  1   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc
    #         's1' : 0.5 ,  #Spin of the electron in question. For coupling two electrons, s1 is the "inner" and s2 is "outer/excited".
    #         'j1'  : 1.5   , #J = L +/- S coupled angular momentum 
    #         'E1': -2.839698     , #Term energy for calculation [eV]. For l=0, this has the lower principle quantum number. Used for qdefect calc and correction.
    #         'n1': 6
    #         , #Princple quantum number.        
    
          
    #         'Corrflags': 'Hdeleps' , #Correction factors for the hyperfine calculation.
    #         #IF coupling is to be done, need J,L,S, j2,l2,s2 for the different possible schemes. It's possible that j1 may not be always known (especially outside of S levels).
    #         'Coupling': 'jj',
    #         'l2'  :  0   , #Orbital Angular Momentum S=0, L=1, P=2, D=3, F=4...etc for the secondary electron
    #         's2' : 2.5 ,  #Spin of the secondary electron
    #         'j2'  : 2.5   , #J = L +/- S coupled angular momentum for the secondary electron
            
    #         'L' : 1 , #Coupled total L from the two electrons
    #         'S' : 3 , #Coupled total S of the two electrons
    #         'J' : 4 , #Coupled total J for the two electrons
    #         'A_Exp': 440, ##I get 248, but lack of dsdn might be hurting me?
    
    #      }
    
    # W7P4couple = Ahfs_Corrected(W7P4)
    # print(W7P4couple)
    
    
    
    
    #%%OlderTests
    
    
    # #%%
    # # qdefG = 3.195211
    
    # corrfactors = Fcorrj(37,.5)*(deltacorr(37-4,85,.5)*epscorr(37-4,85,.5))/H2
    # print(A_Bransden(thisZ-4,thisn,qdefE,thismu,thisI,85,1,0.5))
    
    
    # #%% Trying the new Kopf2 which just uses the deltaFS, seems pretty close without corrections even, but is it limited to doublets?
    # #def Ahfs_Kopf2(Z,deltaFS,muI,I,l_in,j):
    
    # print('Rb85, Z=37, Kopf2')
    # thisZ = 37
    # thisn = 5
    # thisI = 2.5
    # thismu = 1.35298
    # thisQ = 0.277
    # qdefG = 3.195211 #Radzig has 3.131, Theo has 3.1181, I get 3.195211 (Gets closest honestly)
    # qdefE = 2.711483
    # # dsig = 1.04 #From 1-dsigma/dn, where dsigma is the defect value
    # dsig = dsigcorr(3.195211, 3.155019)
    
    # dFSe = 1.589049 - 1.559591 # deltaFineStruct_
    # F1 = Fcorrj(thisZ,.5)
    # print('A(Ground l=0, j=1/2)')
    # print(dsig*F1*Ahfs_Kopf2(thisZ-4,dFSg,qdefG,thismu,thisI,0,0.5)) 
    
    # H2 = Hcorr(thisZ-4,1)
    # F2 = Fcorrj(thisZ-4,.5) #K2 is for j = l-1/2
    
    # # F2 = Fcorr(thisZ,1,kval='k2') #K2 is for j = l-1/2
    # epscorrect = 1-epscorr(thisZ,85,.5)
    # print('A(Excited l=1, j=1/2)')
    # print(epscorrect*F2*Ahfs_Kopf(thisZ-4,thisn,qdefE,thismu,thisI,1,0.5)/H2) 
    
    
    # print('A(Excited l=1, j=3/2)')
    # F3 = Fcorr(thisZ,1,kval='k1') #k1 is for j = l+1/2
    # F3 = Fcorrj(thisZ,1.5) #k1 is for j = l+1/2
    
    # T322 = F3*Ahfs_Kopf(thisZ-4,thisn,qdefE,thismu,thisI,1,1.5)/H2
    # print(T322)
    
    # print('B(Excited l=1, j=3/2)')
    
    # B322 = Bhfs_kopf(T322,1,1.5,thisI,thismu,thisQ) #Calulated using 30.13 - using a and Q to get b. 
    
    # print(B322*hbar)
    
    
    # Rcorr1 = Rcorr(37-4,1,.5)
    # # #This is uing the direct calculation for B instead - Works for the HighZ, doesn't for the lowerZ...
    # # print('Direct Calc B with 30.11')
    # print(Bhfsconst(37-4,5,1,1.5,.277,Rcorr1))#Should be about 25.8, gets 26.858
    
    
    
    
    # print('A(Excited l=1, j=1/2)')
    # T322b = epscorrect*F2*Ahfs_Kopf2(37-4, 1.589049-1.559591,1.35298,2.5,1,0.5)/H2
    # print(T322b)
    
    # print('A(Excited l=1, j=3/2)')
    # T322b = F3*Ahfs_Kopf2(37-4, 1.589049-1.559591,1.35298,2.5,1,1.5)/H2
    # print(T322b)
    
    # print('B(Excited l=1, j=3/2)')
    
    # B322b = Bhfs_kopf(T322b,1,1.5,thisI,thismu,thisQ) 
    # print(B322b*hbar)
    
    
    # #%%
    
    # print('129I II, Z=53, ionized')
    # thisZ = 53
    # thisn = 5
    # thisI = 2.5
    # thismu = 2.81327
    # thisQ = .72
    # qdefG = 2.47
    
    # dsig =-3.77#From 1-dsigma/dn, where dsigma is the defect value, this one uses the dn*
    
    # H2 = Hcorr(thisZ,2)
    # F1 = Fcorrj(thisZ,4)
    # print('A(Lower l=2, j=4')
    # AI_G = dsig*F1*Ahfs_Kopf(thisZ,thisn,qdefG,thismu,thisI,2,4,Za=2)/H2
    # AI_G2 = dsig*F1*Ahfs_Kopf(thisZ-11,thisn,qdefG,thismu,thisI,2,4,Za=2)/H2
    # print(AI_G) 
    # print(AI_G2) 
    
    # print('B(Lower l=2, j=4')
    # BI_G = Bhfs_kopf(AI_G,2,4,thisI,thismu,thisQ)
    # BI_G2 = Bhfs_kopf(AI_G2,2,4,thisI,thismu,thisQ)
    # print(BI_G*hbar)
    # print(BI_G2*hbar)
    # #%%
    
    
    # # dsig = -3.77#From 1-dsigma/dn, where dsigma is the defect value, this one uses the dn*
    # qdefE = 3.47
    # thisn = 6
    
    # dsig = -3.52#From 1-dsigma/dn, where dsigma is the defect value, this one uses the dn*
    
    
    # H2 = Hcorr(thisZ,1)
    # F1 = Fcorrj(thisZ,3)
    # print('A(Upper l=1, j=3')
    # AI_E = dsig*F1*Ahfs_Kopf(thisZ,thisn,qdefE,thismu,thisI,1,3,Za=2)/H2
    # AI_E2 = dsig*F1*Ahfs_Kopf(thisZ-4,thisn,qdefE,thismu,thisI,1,3,Za=2)/H2
    # print(AI_E) 
    # print(AI_E2) 
    
    # print('B(Upper l=1, j=3')
    # BI_E = Bhfs_kopf(AI_E,1,3,thisI,thismu,thisQ) 
    # BI_E2 = Bhfs_kopf(AI_E2,1,3,thisI,thismu,thisQ) 
    # print(BI_E*hbar)
    # print(BI_E2*hbar)
    # #%%
    # H2 = Hcorr(thisZ,1)
    # F2 = Fcorr(thisZ,1,kval='k2') #K2 is for j = l-1/2
    # print('A(Excited l=1, j=1/2)')
    # print(F2*Ahfs_Kopf(thisZ,thisn,qdefE,thismu,thisI,1,0.5)/H2) 
    # print(F2*Ahfs_Kopf(thisZ-4,thisn,qdefE,thismu,thisI,1,0.5)/H2) 
    
    # F3 = Fcorr(thisZ,1,kval='k1') #k1 is for j = l+1/2
    # T32 = F3*Ahfs_Kopf(thisZ,thisn,qdefE,thismu,thisI,1,1.5)/H2
    # T322 = F3*Ahfs_Kopf(thisZ-4,thisn,qdefE,thismu,thisI,1,1.5)/H2
    
    # print('A(Excited l=1, j=3/2)')
    
    # print(T32)
    # print(T322)
    
    # print('B(Excited l=1, j=3/2)')
    
    # B32 = Bhfs_kopf(T32,1,1.5,thisI,thismu,thisQ) 
    # B322 = Bhfs_kopf(T322,1,1.5,thisI,thismu,thisQ) 
    
    # print(B32*hbar)
    # print(B322*hbar)
    
    
    
    # #%%
    # #Try some more hydrogenic now
    # #Starting with Li6
    
    # thisZ = 3
    # thisn = 2
    # thisI = 1
    # thismu = 0.820473
    # thisQ = -0.000822
    # qdefG = 0.418
    # qdefE = 0.0538
    
    # print('Li6')
    
    # F1 = Fcorr(thisZ,0)
    # print('A(Ground l=0, j=1/2)')
    # print(F1*Ahfs_Kopf(thisZ,thisn,qdefG,thismu,thisI,0,0.5)) 
    
    # H2 = Hcorr(thisZ,1)
    # F2 = Fcorr(thisZ,1,kval='k2') #K2 is for j = l-1/2
    # print('A(Excited l=1, j=1/2)')
    # print(F2*Ahfs_Kopf(thisZ,thisn,qdefE,thismu,thisI,1,0.5)/H2) 
    # print(F2*Ahfs_Kopf(thisZ-4,thisn,qdefE,thismu,thisI,1,0.5)/H2) 
    
    # F3 = Fcorr(thisZ,1,kval='k1') #k1 is for j = l+1/2
    # T32 = F3*Ahfs_Kopf(thisZ,thisn,qdefE,thismu,thisI,1,1.5)/H2
    # T322 = F3*Ahfs_Kopf(thisZ-4,thisn,qdefE,thismu,thisI,1,1.5)/H2
    
    # print('A(Excited l=1, j=3/2)')
    
    # print(T32)
    # print(T322)
    
    # print('B(Excited l=1, j=3/2)')
    
    # B32 = Bhfs_kopf(T32,1,1.5,thisI,thismu,thisQ) 
    # B322 = Bhfs_kopf(T322,1,1.5,thisI,thismu,thisQ) 
    
    # print(B32*hbar)
    # print(B322*hbar)
    
    
    # #%%
    # #Li7
    # print('Li7')
    
    # thisZ = 3
    # thisn = 2
    # thisI = -1.5
    # thismu = 3.2564625
    # thisQ =  -.0406
    # qdefG = .418
    # qdefE = 0.0538
    
    
    # F1 = Fcorr(thisZ,0)
    # print('A(Ground l=0, j=1/2)')
    # print(F1*Ahfs_Kopf(thisZ,thisn,qdefG,thismu,thisI,0,0.5)) 
    
    # H2 = Hcorr(thisZ,1)
    # F2 = Fcorr(thisZ,1,kval='k2') #K2 is for j = l-1/2
    # print('A(Excited l=1, j=1/2)')
    # print(F2*Ahfs_Kopf(thisZ,thisn,qdefE,thismu,thisI,1,0.5)/H2) 
    # print(F2*Ahfs_Kopf(thisZ-4,thisn,qdefE,thismu,thisI,1,0.5)/H2) 
    
    # F3 = Fcorr(thisZ,1,kval='k1') #k1 is for j = l+1/2
    # T32 = F3*Ahfs_Kopf(thisZ,thisn,qdefE,thismu,thisI,1,1.5)/H2
    # T322 = F3*Ahfs_Kopf(thisZ-4,thisn,qdefE,thismu,thisI,1,1.5)/H2
    
    # print('A(Excited l=1, j=3/2)')
    
    # print(T32)
    # print(T322)
    
    # print('B(Excited l=1, j=3/2)')
    
    # B32 = Bhfs_kopf(T32,1,1.5,thisI,thismu,thisQ) 
    # B322 = Bhfs_kopf(T322,1,1.5,thisI,thismu,thisQ) 
    
    # print(B32*hbar)
    # print(B322*hbar)
    
    # #%%
    # #Na23
    # print('Na23, Z=11')
    # thisZ = 11
    # thisn = 3
    # thisI = 1.5
    # thismu = 2.217522
    # thisQ =  .1045
    # qdefG = 1.3595
    # qdefE = 0.867
    # dsig = 1.0159
    
    # F1 = Fcorr(thisZ,0)
    # print('A(Ground l=0, j=1/2)')
    # print(dsig*F1*Ahfs_Kopf(thisZ,thisn,qdefG,thismu,thisI,0,0.5)) 
    
    # H2 = Hcorr(thisZ,1)
    # F2 = Fcorr(thisZ,1,kval='k2') #K2 is for j = l-1/2
    # print('A(Excited l=1, j=1/2)')
    # print(F2*Ahfs_Kopf(thisZ,thisn,qdefE,thismu,thisI,1,0.5)/H2) 
    # print(F2*Ahfs_Kopf(thisZ-4,thisn,qdefE,thismu,thisI,1,0.5)/H2) 
    
    # F3 = Fcorr(thisZ,1,kval='k1') #k1 is for j = l+1/2
    # T32 = F3*Ahfs_Kopf(thisZ,thisn,qdefE,thismu,thisI,1,1.5)/H2
    # T322 = F3*Ahfs_Kopf(thisZ-4,thisn,qdefE,thismu,thisI,1,1.5)/H2
    
    # print('A(Excited l=1, j=3/2)')
    
    # print(T32)
    # print(T322)
    
    # print('B(Excited l=1, j=3/2)')
    
    # B32 = Bhfs_kopf(T32,1,1.5,thisI,thismu,thisQ) 
    # B322 = Bhfs_kopf(T322,1,1.5,thisI,thismu,thisQ) 
    
    # print(B32*hbar)
    # print(B322*hbar)
    
    # #%%
    # #K39 n = 4
    # print('K39, Z=19')
    # thisZ = 19
    # thisn = 4
    # thisI = 1.5
    # thismu = 0.39147
    # thisQ =  0.049 #There are also galues of Q = 0.049 and 0.06 given,
    
    # qdefG = 2.1779
    # qdefE = 1.7123
    # dsig = 1
    
    # F1 = Fcorr(thisZ,0)
    # print('A(Ground l=0, j=1/2)')
    # print(dsig*F1*Ahfs_Kopf(thisZ,thisn,qdefG,thismu,thisI,0,0.5)) 
    
    # H2 = Hcorr(thisZ,1)
    # F2 = Fcorr(thisZ,1,kval='k2') #K2 is for j = l-1/2
    # print('A(Excited l=1, j=1/2)')
    # print(F2*Ahfs_Kopf(thisZ,thisn,qdefE,thismu,thisI,1,0.5)/H2) 
    # print(F2*Ahfs_Kopf(thisZ-4,thisn,qdefE,thismu,thisI,1,0.5)/H2) 
    
    # F3 = Fcorr(thisZ,1,kval='k1') #k1 is for j = l+1/2
    # T32 = F3*Ahfs_Kopf(thisZ,thisn,qdefE,thismu,thisI,1,1.5)/H2
    # T322 = F3*Ahfs_Kopf(thisZ-4,thisn,qdefE,thismu,thisI,1,1.5)/H2
    
    # print('A(Excited l=1, j=3/2)')
    
    # print(T32)
    # print(T322)
    
    # print('B(Excited l=1, j=3/2)')
    
    # B32 = Bhfs_kopf(T32,1,1.5,thisI,thismu,thisQ) 
    # B322 = Bhfs_kopf(T322,1,1.5,thisI,thismu,thisQ) 
    
    # print(B32*hbar)
    # print(B322*hbar)
    
    
    
    # #%%Try 183W for 6s ^7S_3
    
    # #Screening potential. Each of same group (same principle qn) contribute 0.35. n-1 contributes 0.85 and 1 for n-2 and below. (for s,p)
    # #For d,f. just 1 for all to the left?
    
    # #For example, Z* of tungsten is 74 - 70.55 - 3.45 for Za? Unclear.
    # #Using Z* from "screening potential" instead of quantum defect seems to work?
    # #I get the "correct" answer for the Ahfs here if a "nucleus deformation " is included per table 68 in Kopf...factor is 6.5
    # #Also gets 
    # #There are many ways to get near the expected HFS value, but it gets muddier...if I just use the values as calculated in my code
    # #I get the better agreement of the observed splitting in the devices...
    # FW = Fcorr(74,0) 
    # print(FW*Ahfs_Kopf(74,6,4.66,.11778476,-0.5,0,3 ,Za=1))  #505 MHz, get 3583 MHz, gets 520MHz if the screening potential is used instead of defect
    
    # FW2 = Fcorr(74,1) 
    # HcorrW = Hcorr(74,1)
    # print(FW2*Ahfs_Kopf(74,6,4.37,.11778476,-0.5,1,2 )/HcorrW)  #Unknown
    # print(FW2*Ahfs_Kopf(74,6,4.37,.11778476,-0.5,1,3 )/HcorrW)  #496 MHz
    
    # print(FW2*Ahfs_Kopf(74,6,4.37,.11778476,-0.5,1,4 )/HcorrW)  #440 MHz
    
    
    
    # #%%
    # #This is for the 7D1, sb about 1GHz
    # FW = Fcorr(74,2, kval='k2') 
    # HcorrW = Hcorr(74,2)
    # dsig = 1.0577
    # print(dsig*FW*Ahfs_Kopf(74-11,6,4.2,.11778476,-0.5,2,1 ,Za=4.5)/HcorrW) #Gets 1125 using Za = screened potential
    
    # #%%7F1 and 7F1 s.b about 732 MHz
    # FW = Fcorr(74,3, kval='k2') 
    # HcorrW = Hcorr(74,3)
    # dsig = 1.0577
    # print(dsig*FW*Ahfs_Kopf(74-11,6,4.2,.11778476,-0.5,3,1 ,Za=4.5)/HcorrW) #Gets 766 with screened potential, likely need to decrease defect more for F
    # #%%7F1 and 7F3 s.b about 437 MHz
    # #Looks like I have some extra difficulties for j=l...
    # FW = Fcorrj(74,3, kval='k2') 
    # HcorrW = Hcorr(74,3)
    # dsig = 1.0577
    # print(dsig*FW*Ahfs_Kopf(74-11,6,4.2,.11778476,-0.5,3,3 ,Za=4.5)/HcorrW) #Gets 766 with screened potential, likely need to decrease defect more for F
    
    
    # #%% Zinc67
    # FZn = Fcorr(30,1,kval='k1')
    # FZn2 = Fcorr(30,0,kval='k1')
    # HcorrZn = Hcorr(30,1)
    # Zn4p =FZn*Ahfs_Kopf(30,4,2.06,0.875479,2.5,1,2 )/HcorrZn
    # Zn4s = FZn2*Ahfs_Kopf(30,4,2.6408,0.875479,2.5,0,2 )
    
    # print(Zn4p)
    # print(Zn4s)
    
    # #%% #Ar 39
    # #In LS coupling, both of these should be P states..
    # FAr = Fcorr(18,0,kval='k1')
    # Ar4s322 = FAr*Ahfs_Kopf(18,4,2.0541,-1.588,-3.5,0,2 ) #Should be -285 (Williams) gets -291
    # print(Ar4s322)
    
    # FAr = Fcorr(18,1)
    # Ar4p323 = FAr*Ahfs_Kopf(18-4,4,1.6151,-1.588,-3.5,0,3 ) #Should be -135, gets -154
    # print(Ar4p323)
    
    # ArBP323 = Bhfs_kopf(Ar4p323,1,3,-3.5,-1.588,-0.12)/1.5 ##Should be approx 113 MHz, I get exactly 2x that.
    # print(ArBP323*hbar)
    
    
    
    # #%%
    
    
    # # print(Bhfsconst(53-11,(5-2.47),2,4,.7,.355781,Z_a=2)*-.81)#s.b about-233.8, gets -254, 
    # # print(Bhfsconst(53-11,(5-2.47),2,4,.7,Rcorr1,Z_a=2)*-3.77379)#s.b about-233.8, gets -247
    
    # # print(Bhfsconst(53-4,(6-3.47),1,3,.7,Rcorr(53,1,2),Z_a=2)*-1.97887)#s.b about--430, gets 
    # # print(Bhfsconst(53-4,(6-3.47),1,3,.7,.251328,Z_a=2)*-1.97887)#s.b about--430, gets 
    
