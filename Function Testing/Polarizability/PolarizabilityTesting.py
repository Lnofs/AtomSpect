# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %Leo Nofs
"""






def ScalarPol(E1,E2,j1,dipoleval):
    #This expresion comes from the Udel atom webpage. The results are in atomic units, which would be likely Hartree or Rydbergs (Facotr 2 diff)
    #In terms of eV, 1 Hartree = 2 Rydberg = 27.211386245 eV = 6579.683920THz
    #Seems this should be term energy as well, so E_ionize - E
    alphatemp = 0
    #This is generally assuming that there are going to be multiple energies and j values to go with them. Form comes from udel.edu.atom
    #E = hnu
    for Es,dpv in zip(E2,dipoleval):
        alphatemp = alphatemp + (dpv**2 )/((E1-Es)/27.211386245981) #Convert from ev to Hartree for the calc
    return [1000*2.48832E-8 * alphatemp*(2/(3*(2*j1 +1))) , alphatemp*(2/(3*(2*j1 +1)))] #Calculation gives result in 1000[a.u] it seems, so convert to Hz.
    # return alphatemp*(2/(3*(2*j1 +1)))

def TensorPol(E1,E2,j1,j2,dipoleval):
    alphatemp = 0
    #This is generally assuming that there are going to be multiple energies and j values to go with them. Form comes from udel.edu.atom
    
    for js,Es in zip(j2,E2):
        alphatemp = alphatemp + (dipoleval/(E1-E2))**2
    return alphatemp*(2/(3*(2*j1 +1)))

#%%
# E_io = 4.177128
# E_p12 =1.55959103
# E_p32 =1.58904906 #Excited energy Rb
# Epterm = 1.57923
# Sscalarpol=ScalarPol(4.177128,[0],.5,[1]) #Steck: 0.0794, Atom: .31828

# D1scalarpol=ScalarPol(E_io - E_p12,[0],.5,[1]) #Steck: 0.122306 , Atom: .805
# D2scalarpol=ScalarPol(E_io - E_p32-E_p12,[0],1.5,[1]) #Steck: 0.13408 , Atom: .868
# print(Sscalarpol*100*100)

# print(D1scalarpol*100*100)
# print(D2scalarpol*100*100)
# #%%
# D1scalarpol=ScalarPol(E_io - Epterm ,[0],.5,[1]) #Steck: 0.122306 , Atom: .805
# D2scalarpol=ScalarPol(E_io - Epterm-E_p12 ,[0],1.5,[1]) #Steck: 0.13408 , Atom: .868

# print(Sscalarpol)
# print(D1scalarpol)
# print(D2scalarpol)
#%%
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
        }

#%%Cstest
#Ev
# E_io = 3.893905
# E0 = 2.2981127
# E_p12 =1.385928617528
# E_p32 =1.454620692074
# E_p = 1.431723333892
# Sscalarpol=ScalarPol(E_io,[E0],.5,[1]) #Safro 6.237 [1000] au = .000155196 Hz/(V/m)^2 , Atom: 6.232 (1000 a.u), so 6237 in au

# #cm-1
# # E_io = 31406.4677325
# # E0 = 18535.5286
# # E_p12 =1.55959103
# # E_p32 =1.58904906 #Excited energy Rb
# # Sscalarpol=ScalarPol(E_io,[E0],.5,[1]) #Safro 6.237 [1000] au , Atom: 6.232 (1000 a.u), so 6237 in au
# #Hartree
# # E_io = .143098396099
# # E0 = .084454082
# # E_p12 =1.55959103
# # E_p32 =1.58904906 #Excited energy Rb
# # Sscalarpol=ScalarPol(E_io,[E0],.5,[1]) #Safro 6.237 [1000] au , Atom: 6.232 (1000 a.u), so 6237 in au


# #alpha(SI) = 2.48832E-8 * alpha[a.u]
# print(Sscalarpol)

# #For 6p douple
# D1scalarpol=ScalarPol(E_io,[0,E_p12],.5,[.5,.5])  #Safro 1.339 [1000] au = .000155196 Hz/(V/m)^2 
# D2scalarpol=ScalarPol(E_io-E_p12,[0,E_p32],1.5,[.5,.5])  #Safro 1.651 [1000] au = .000155196 Hz/(V/m)^2 


# print(D1scalarpol)
# print(D2scalarpol)
# #%%Cs7p and 8s

# E_io = 3.893905
# E0 = 2.2981127
# E_p12 =2.6985592
# E_p32 =2.721006
# E_p = 2.713524
# #alpha(SI) = 2.48832E-8 * alpha[a.u]

# Sscalarpol=ScalarPol(E_io,[E0],.5,[1,1]) #Safro 6.237 [1000] au = .000155196 Hz/(V/m)^2 , Atom: 6.232 (1000 a.u), so 6237 in au

# print(Sscalarpol)

# #For 7p douplet
# D1scalarpol=ScalarPol(E_io-E_p,[0],.5,[1,.5])  #Safro 29.88 [1000] au = .000155196 Hz/(V/m)^2 
# D2scalarpol=ScalarPol(E_io-E_p,[0,E_io-E_p32],1.5,[1])  #Safro 37.51 [1000] au = .000155196 Hz/(V/m)^2 


# print(D1scalarpol[1]*4)
# print(D2scalarpol[1])
