# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27,2024

@author: Leo Nofs
"""
import numpy as np
from scipy.special import gamma
from SpectSplit import read_HFSInput


from functools import reduce
from itertools import chain,combinations


import csv, os, sys
import numpy as np
from matplotlib import pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir, os.path.pardir)))



#%%Constant Declarations
#Just a bunch of atomic constants. In the future, should make a class that is called which contains these.
m_prot = 1.667e-27 #Mass of proton in kg
m_elec = 9.11e-31
Rfreq = 3.28984196025e15 #Rydberg Freq
c_light = 299792458 #Speed of light in m/s
alpha = 0.0072973525643
a_o = 5.29e-11 #Bohr radii in m
qe = 1.602e-16
eps0 = 8.85e-12 #freespace permittivity
hbar = 1.0546e-34 #Plank's constant hbar in Angular
hh = 2.0*np.pi*hbar #Plank's constant in frequency
muB = 9.27401e-24 #Bohr magneton J/T
muP = 2.79284734 #proton magneton (in terms of nuclear magnetons, nm)
muN = 1.5210322e-3 #Nucleon magneton (not nuclear dipole moment), in terms of muB
muN2 = muB/1836
mu_Bhh = muB/hh
mu_B = muB/hbar #Common value of Bohr_Mag/hbar which ensures the computation doesn't drop small terms. [Hz/T]
g_L = 1 #Lande Orbital g-factor, assuming infinite mass, otherwise is g_L = 1 - me/mnuc
g_S = 2.002319 #Lande electron spin g-factor
a_o = 5.29e-11 #Bohr radii in m
mu_0 = 4*np.pi*1e-7 #freespace mag permeability
eps0 = 8.85e-12 #freespace permittivity
kboltz = 1.381e-23 #Boltzmann Constant
m_prot = 1.667e-27 #Mass of proton in kg
m_elec = 9.11e-31
q_charge = 1.602e-19
hhev = 4.135668e-15 #Plank's constant in Ev/Hz
hhc = 1240e-13 #eV m for hc.
hcRyd=  13.605693122990 #hcRyd_inf - Infinite mass rydberg energy (CODATA)


#%%Functions

def powerset(iterable):
    #This is a modified version of the python website that outputs a list where the strings are all conbined together per entry.
    "Subsequences of the iterable from shortest to longest."
    # powerset([1,2,3]) → () (1,) (2,) (3,) (12) (13) (23) (123)
    s = list(iterable)
    sb = chain.from_iterable(combinations(s, r) for r in range(len(s)+1))
    ssb = [list(x) for x in sb]
    ssbb = [reduce(lambda x,y:x+y,things) for things in ssb[1:]]
    ssbb.append('NoCorr')
    return ssbb

def pererror(v1,v2):
    #Assuming that v1 is the exact value
    return 100*np.abs((v1 - v2)/v1)

def perdiff(v1,v2):
    #Assuming that v1 is the exact value
    return 200*(v1 - v2)/(v1+v2)


def Hz_to_nm(Ediff):
    #Converts from Hz to nm
    return 2.99792458E17 / Ediff



def ev_to_Hz(Ediff):
    #Converts from eV to MHz
    return (1e-6)*Ediff*2.41804e14
    
def cm_to_Hz(E_diff):
    #Converts from energy in cm^-1 to energy in Hz (E = hf = hc/lambda)
    return 100*c_light*E_diff *1E-6

def ev_to_cm(E_diff):
    #Convert from eV to cm^-1
    return E_diff/(1.2398e-4)
def cm_to_nm(E_diff):
    #Converts from energy in cm^-1 to wavelength in nm. 
    return 1e7/E_diff
def cm_to_ev(E_diff):
    #Convert from 1/cm to eV
    return E_diff*(1.2398e-4)
def Ahfs_dnu(deltanu,I):
    #From Blinder1993 regarding doublets. Works, but may not be useful
    return 2*deltanu/(2*I +1)



def Ahfs_Kopf(Z,n,ndefect,muI,I,l_in,j, Za=1):
    #Eqyatuib 26.8 for hydrogenic and small quantum numbers
    #Za is ionization, neutral = 1.
    thisgI = (muI/I)*(m_elec/m_prot) #this is giprime*
    if j==0:
        #There is no splitting for J=0
        return 0 
    else:
        
        tempval = hcRyd*alpha**2 * Z* Za**2 * thisgI / (j*(j+1)*(l_in + 0.5)*((n - ndefect)**3))
            # tempval = 13.605693122990 *alpha**2 * Z* Za**2 * thisgI / ((j*j)*(l_in)*((n - ndefect)**3))

    return ev_to_Hz(tempval)

def A_ratio(A2,mu1,I1,mu2,I2):
    #Using the accuracy of the isotope ratio for nuclear moments to calculate the HFS constant if one is known.
    return A2*mu1*I2/(mu2*I1)

#%%Corrections


def reducedmass(m1,m2):
    return 1/((1/m1) + (1/m2))

def correctedE(Mass):
    return reducedmass(Mass*m_prot,m_elec)*hcRyd/m_elec

def epscorr(Z,M,j,r_ro=1):
    #Based on Kopf 26.23, where <r^2/(r_0)^2> has limiting cases of: 3/5 for uniform, 1 for surface distribution of nuclear moment.
    # Uses same def of r_o as Kuhn (pg 337) r_o = 2*(A)**(1/3) *10^-15 m

    rtmp = 2e-15 *((M)**(1/3))
    thisrho = np.sqrt((j+.5)**2 - (alpha*Z)**2)
    return 1-((Z*rtmp/a_o)*(a_o/(2*Z*rtmp))**(2*(1-thisrho)))*r_ro

def deltacorr(Z,M,j):
    #Based on  Kopf 26.22 with r_o = 2*(A)**(1/3) *10^-15 m - Pg 337 Kuhn. From assuming nucleus is uniformly charged sphere.
    #Works based on table VI,2 on Kuhn pg 339!
    rtmp = 2e-15 *((M)**(1/3))
    ytemp = 2*Z*rtmp/a_o
    
    thisrho = np.sqrt((j+.5)**2 - (alpha*Z)**2)
    
    return 1-(((2*(1-thisrho)*thisrho*(2*thisrho+1))/((2*thisrho -1)*(gamma(2*thisrho +1))**2))*ytemp**(2*thisrho -1))


def Fcorrj(Z,j):
    #Eq 26.18 Kopf, adjusted to j values only. Matches all compared values in Kopf and Kuhn.
    #Assuming s2 = 1/2 always, which I guess is true yeah.
    #k = l+1 for j = l+.5, k1
    #k = -l for j = l-1/2 ,k2
    #Adding in from table 8, based on 26.25 for nonhydrogenic (corrections)

    rho = np.sqrt((np.abs(j) + 0.5)**2 - (alpha*Z)**2)    
    tF = (4*j*(j+1)*(j+.5))/(rho*(4* rho**2 -1))
    return tF

def Rcorr(Z,l,s,kval):
    #This is from Kopf 30.11, tested and matches Kuhn pg 337 # Only used in the calculation of the Bhfs constant?.
    
    #k = l+1 for j = l+.5, k1
    #k = l-1 for j = l-1/2 ,k2  
    #
    
    
    if kval == 'l+':
       j = l + np.abs(s)
       k = l+1
    else:
        j = l - np.abs(s)
        k=-l
    rho = np.sqrt((j + 0.5)**2 - (alpha*Z)**2)
    # rho = np.sqrt(k**2 - (alpha*Z)**2)
    if l==0:
        return 1
    else:
        tF = ((3*k*(k+1) - rho**2 + 1)*(l*(l+1)*(2*l+1)))/(rho*(rho**2 - 1)*(4* rho**2 -1))
        return tF



def dsigcorr(qdef1,qdef2):
    
    #This is 1 - dsigma/dn, wherein dn=1, so it's just the difference between two quantum defects of the same term.
    #I seem to get better results if I use dsigma/dn* for different n levels.
    #Thie supposedly is only valid for s states? Unclear why.
    return 1 - (qdef2 - qdef1)

def Hcorr(Z,l):
    
    
    #Kopf 26.1 - Matches compared values from Kuhn pg 362 and others. - Unsure if this is supposed to be used all the time or just for doublets.
    rhop = np.sqrt((np.abs(l)+1)**2 - (alpha*Z)**2) #For doubet l+1/2
    rhopp = np.sqrt(l**2 - (alpha*Z)**2) #For doublet l-1/2

    THcorr= 2*np.abs(l)*(np.abs(l)+1)*(rhop - rhopp - 1)/((alpha*Z)**2)
    return THcorr

def Isotope(Eterm,mlight,mheavy):
    #Heaviest isotope apparentlly always lies lower in energy per Kopfermann 33.1
    #This shift seems to work fine with cm-1 units as input, or at least pretty close.
    return Eterm + m_elec*Eterm*(1/mlight - 1/mheavy)/m_prot

def Isotopediff(Eterm,mlight,mheavy):
    return 1e7/Isotope(Eterm,1,1) - 1e7/Isotope(Eterm,mlight,mheavy)
#%%Coupling

def sscouple(S,s_in,s_out):
    #This calculates the coupling factors for the case of ss coupling per Kopfermann 27.4. Typically going to have s2=1/2
    #Returns [s1fac,s2fac]
    s_infac = (S*(S+1) + s_in*(s_in +1) - s_out*(s_out+1))/(2*S*(S +1)) #This is the correction for the "lower" electron
    s_outfac = (S*(S+1) + s_out*(s_out +1) - s_in*(s_in+1))/(2*S*(S +1)) #This is the correction for the "excited" electron
    return[s_infac,s_outfac]

def sl1couple(J,S,s_in,s_out,l_out):
    #This calculates the coupling factors for the case of sl1 for spin-spin is large compared to spin-orbital (of electron 2). 

    #Returns [sl1fac]
    sl1fac = ((J*(J+1) + S*(S +1) - l_out*(l_out+1))/(2*J*(J +1)))*((S*(S+1) + s_in*(s_in +1) - s_out*(s_out+1))/(2*S*(S +1)))
    return sl1fac

def sl2couple(J,s_in,j_out):
    return ((J*(J+1) + s_in*(s_in +1) - j_out*(j_out+1))/(2*J*(J +1)))

def jjcouple(J,j_out,s_in):
    #Kopfermann 27.19 for the case of jj coupling being fine and including the outermost electron.
    #Uses as1 and aj2
    jj_in = ((J*(J+1) + s_in*(s_in +1) - j_out*(j_out+1))/(2*J*(J +1))) #jj_in coupling
    jj_out = ((J*(J+1) + j_out*(j_out +1) - s_in*(s_in+1))/(2*J*(J +1))) #jj_out coupling
    return[jj_in,jj_out]

def jjcouple2(J,j_out,j_in):
    #Kopfermann 27.19 for the case of jj coupling being fine and including the outermost electron.
    #Uses as1 and aj2
    jj_in = ((J*(J+1) + j_in*(j_in +1) - j_out*(j_out+1))/(2*J*(J +1))) #jj_in coupling
    jj_out = ((J*(J+1) + j_out*(j_out +1) - j_in*(j_in+1))/(2*J*(J +1))) #jj_out coupling
    return[jj_in,jj_out]



def Ahfs_Kopf2(Z,deltaFS,muI,I,l_in,j):
    #Equation 26.11 for non-hydrogenic using deltaFS which is energy sep between
    #Za is ionization, neutral = 1.
    #This doesn't work for l=0 ground state (No term energy)
    thisgI = (muI/I)*(m_elec/m_prot) #this is giprime, sign is the same as I?
    if j==0:
        return 0
    else:
        tempval = deltaFS*l_in*(l_in+1)*thisgI/((l_in+0.5)*j*(j+1)*Z)
    return ev_to_Hz(tempval)

def gI_kopf(I,L,S,muI,Mass):
    #Not sure what this is, doesn't seem to be gI actually?
    tgi = g_S*(I*(I+1) + S*(S+1) - L*(L+1))/(2*I*(I+1)) +g_L*(I*(I+1) - S*(S+1) + L*(L+1))/(2*I*(I+1))
    # extrafactors = I*(I+1)*muI/(Mass*m_prot/m_elec)
    print('Do not use this function, it does not seem to correspond to anything')
    return tgi



def Bhfs_kopf(A_in,l,j,I,mu_I,Qin):
    #Trying 30.13 to get b from Q and a.... generally seems to work well enough.

    InvB = eps0*mu_0*muB*muB * 4 * l*(l+1)*mu_I * m_elec /(A_in*qe*qe*Qin*j * (2*j- 1)  * m_prot * I)
    tempB = 1/InvB
    return hbar*tempB

def Bhfsconst(Z,n,l,J,Q,Rcorr,Z_a=1):
    #Eq 30.11 (pg 154) in Kopf -- Need to do unit analysis to get this working
    #This also is in Kuhn, and works, but is off by a factor of 10 it seems for p states?
    
    if l==0 or J==0 or J==0.5:
        Btmp = 0
    else:
        rexpected = (Z* Z_a**2)/ (a_o**3 * n**3 *(l+.5)*l*(l+1))
        Btmp = 0.5*Rcorr*rexpected*Q*(qe**2)*((2*J -1)/(2*J +2))
    return cm_to_Hz(Btmp) #Converts to MHz

def A_Bransden(Z,n,ndef,muI,I,MnucAMU,l,j,Za=1):
    nucmass = reducedmass(m_elec,MnucAMU)
    a_mu = a_o*m_elec/nucmass
    thisgi = muI/I
    #This equation is modified to remove a division by a factor of Pi, and doing so works, but including correction factors would make it not.
    tempA = mu_0*thisgi*muB*muP*(1/(j*(j+1)*(2*l+1)))*Z*Za**2 / (a_mu**3 * (n-ndef)**3)
    return tempA

#%%AHFSCorrFunction


def Qdefect(E_term,n1,Za = 1,E_ryd = hcRyd):
    #Solves E_b = Eo Za**2 /(n- sigma)**2 for sigma and n*, where n* = n-sigma, E_b = E_ionize - E_Term
    #THe relation between Ionization energy and use in quantum defects needs more work.
    n_star = np.sqrt(E_ryd*Za*Za/(E_term))
    qdef = n1 - n_star
    return [n_star,qdef]

def QdefCorr(E_ionize,E1,n1,E2,n2, Za=1,E_ryd=hcRyd):
    #When doing this calculation, it's term energies that are used only. These are the same for all J (only LS resolved)

    ET1 = E_ionize - E1
    ET2 = E_ionize - E2
    
    QDef1 = Qdefect(ET1,n1,Za=Za,E_ryd = E_ryd)
    QDef2 = Qdefect(ET2,n2,Za=Za)
    return 1 - ((QDef2[1]-QDef1[1]))

def Ahfs_Calc(Input,verbose = False,detailed=False):
    '''
    

    Parameters
    ----------
    Input : Dict
        DESCRIPTION.
    verbose : TYPE, optional
        DESCRIPTION. The default is False.
    Input2 : TYPE, optional
        DESCRIPTION. The default is {}.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    #dict_keys(['Atom', 'Label', 'AMU', 'n1', 'E1', 'Eionize', 's1', 'l1', 'j1', 'Z', 'Za', 'I', 'mu_I', 'Q(b)', 'A_Exp', 'B_Exp', 'CorrFlags', 'Coupling', 'n2', 'E2', 's2', 'l2', 'j2', 'S', 'L', 'J', 'A_Calc', 'B_Calc', '', 'Source', 'Notes'])

    #using an input dictionary to do all of the calculations in a single swoop. Including the Z_effective corrections based on the l value.
   
    #Equation 26.8 for hydrogenic and small quantum numbers, but it seems to be the same for all systems, so just going to use with corrections.
    #Also is Eq 31 and 32 from Kuhn on pg 339
    #Za is ionization, neutral = 1. 
    
    #Verbose(bool) enables or disables error messages.
    muI = Input['muI']
    I = Input['I']
    l1 = Input['l1']
    j1=Input['j1']
    s1=Input['s1']

    
    thisZ = Input['Z'] #Number of charges, screened charge (Casamir1933) has Zi = Z for s, Z-4 for p, Z-11 for d, and Z-35 for 4f.
    n = Input['n1'] #Principle quantum number
    #For the Rubidium test, it seems that the correction factors were only needed for the ground s level.
    
    if 'Za' in Input:
        Za = Input['Za'] #Ionization/effective charge. Za=1 for neutral, 2 for singly ionized, etc...
    else:
        Za = 1

    if 'AMU' in Input:
        qdef = Qdefect(Input['E_ionize']-Input['E1'], Input['n1'], Za=Za,E_ryd=correctedE(Input['AMU'])) #This calculates [n*,sigma], always needed. Defined/discussed by Casamir (1936ish)
    else:
        qdef = Qdefect(Input['E_ionize']-Input['E1'], Input['n1'], Za=Za) #This calculates [n*,sigma], always needed. Defined/discussed by Casamir (1936ish)

    ndefect = qdef[1]
    # print(qdef)
    #It looks like the screening constant is possibly the same as delta?, so gonna try that out.

    # Zi = thisZ - ndefect       
                    
    # Testing match/case instead of a bunch of if statements for handling of Zi.
    match np.abs(l1):
        case 0: #s
            # Zi = thisZ - qdef[1]
            Zi = thisZ

        case 1: #p        
            Zi = thisZ -4
        case 2: #d
            if thisZ>12:
                Zi = thisZ - 11 
            else:
                Zi = qdef[0]-2
        case 3: #f
            Zi = thisZ-35 - qdef[0]
 
    thisgI = (muI/I)*(m_elec/m_prot) #this is giprime, unsure if a further mass correction is needed somewhere for nuclear weight

    #Since you can check for parts of a string, the corrections can be toggled as part of a single flag string.
    #Possible inclusions are: F, H, eps, del, ds for the corrections as expected. Still working on hydrogenic vs non and what they might entail other than couplings.
    
    #Will attempt to figure out what types of corrections are typically applicable and when.
    corrfactor = 1
    tFcorr = 1
    tHcorr = 1
    tepscorr = 1
    tdelcorr = 1
    tdscorr = 1
        
    if 'Corrflags' in Input:
        #Update: make some default correction flags based on observation of results.

            corrfactors = {}
            if 'R' in Input['Corrflags']:
                #This takes care of the two different components. This work to correct for Lithium, but none of the heavier elements.
                # jj = l1 + s1
            
                if j1 == (l1+.5):
                    tRcorr = Rcorr(thisZ,l1,s1,'l+')
                    corrfactor = corrfactor*tRcorr
                    corrfactors['Rcorr'] = tRcorr

                if j1 == (l1-.5):
                    tRcorr = Rcorr(thisZ,l1,s1,'l-')
                    corrfactor = corrfactor*tRcorr
                    corrfactors['Rcorr'] = tRcorr


            if 'F' in Input['Corrflags']:
                tFcorr = Fcorrj(np.abs(thisZ),j1)
                corrfactor = corrfactor*tFcorr
                corrfactors['Fcorr'] = tFcorr

            if 'H' in Input['Corrflags']:
                if np.abs(l1)>0:
                    tHcorr = Hcorr(np.abs(thisZ),l1)
                    corrfactor = corrfactor/tHcorr
                    corrfactors['Hcorr'] = tHcorr

            if 'eps' in Input['Corrflags']:
                if 'r_ro' in Input:                    
                    tepscorr = epscorr(np.abs(thisZ),Input['AMU'],j1,r_ro=Input['r_ro'])
                else:
                    tepscorr = epscorr(np.abs(thisZ),Input['AMU'],j1)
                corrfactor = corrfactor*tepscorr
                corrfactors['epscorr'] = tepscorr

            if 'del' in Input['Corrflags']:
                tdelcorr = deltacorr(np.abs(thisZ),Input['AMU'],j1)
                corrfactor = corrfactor*tdelcorr
                corrfactors['deltacorr'] = tdelcorr

            if 'ds' in Input['Corrflags']:
                #Trying the inclusion of ds for l>0 because I don't understand why it isn't included.  
                try:
                    if 'n2' in Input:
                        tdscorr = QdefCorr(Input['E_ionize'], Input['E1'], Input['n1'],Input['E2'], Input['n2'],Za=Za)
                        corrfactor = corrfactor*tdscorr
                        corrfactors['dscorr'] = tdscorr

                        # print(tdscorr)
                    else:
                        corrfactor=corrfactor
                        if verbose:
                            print('Need E2 and n2 to calculate dsdn correction factor - No dsdn factor included')
                except:             
                    corrfactor=corrfactor
                    if verbose:
                        print('Need E2 and n2 to calculate dsdn correction factor - No dsdn factor included')
            corrfactors['TotalCorr'] = corrfactor
                

                
    if j1==0:
        #There is no splitting for J=0
        return 0
    else:
      
        tempA = corrfactor*correctedE(Input['AMU']) *alpha**2 * Zi* Za**2 * thisgI / (j1*(j1+1)*(l1 + 0.5)*((n - ndefect)**3))
        if detailed:
            return [ev_to_Hz(tempA),corrfactor]
        else:
            return ev_to_Hz(tempA)

def Ahfs(Input,verbose = False, Input2 = {}):
    tempA = Ahfs_Calc(Input,verbose=verbose)
    data = {}

    if ('Coupling' in Input) and Input2:
            tempA2 = Ahfs_Calc(Input2,verbose=verbose)
            S = Input['S']
            s1 = Input['s1']
            s2 = Input2['s1']
            L = Input['L']
            l1 = Input['l1']
            l2 = Input2['l1']
            J = Input['J']
            j1 = Input['j1']
            j2 = Input2['j1']
            match Input['Coupling']:
              
                #Per Kuhn pg 339, the correction factors are for the small values, so per electron!
                #This assumes that Input1 is the inner electron, and input2 is the outer electron.
                    case 'ss':
                        sscoup = sscouple(S,s1,s2)
                        
                        tempout = tempA*sscoup[0] + tempA2*sscoup[1]
                        # print('do ss coupling') #Kopf 27.4
                    case 'sl1':
                        sl1coup = sl1couple(J,S,s1,s2,l2)
                        tempout = tempA*sl1coup
                        # print('do sl coupling for big spin-spin') #27.9
                    case 'sl2':
                        sl2coup = sl2couple(J,s1,j2)
                        tempout = tempA*sl2coup
                        # print('do sl coupling for big spin-orbital of electron 2 (jj coupling)') #Eq 27.15
                    case 'jj':
                        jjcoup = jjcouple(J,j2,s1)
                        tempout = jjcoup[0]*tempA +jjcoup[1]*tempA2
                        # print('do jj coupling inclusion of both electrons') #27.19  
                    case 'jj2':
                        jjcoup = jjcouple2(J,j2,j1)
                        tempout = jjcoup[0]*tempA +jjcoup[1]*tempA2
                        # print('do jj coupling inclusion of both electrons') #27.19  
                    case 'None':
                        
                        tempout = tempA
    
                    case _:
                        tempout = tempA
                        print('Something went wrong. Check Input2 parameters. Only A1 returned.')
            data['Ahfsout'] = tempA2
            data['Acoupled'] = tempout
            
    data['A'] = tempA
    return data
                    
    
def GenCouple(A,B,C,D):
    #Couples (A(A+1) + B(B+1) - C(C+1)) /(2D(D+1)) since it is such a common calculation, esp for testing.
    return (A*(A+1) + B*(B+1) - C*(C+1)) /(2*D*(D+1))

