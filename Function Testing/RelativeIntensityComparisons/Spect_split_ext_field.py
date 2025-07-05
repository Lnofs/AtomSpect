# -*- coding: utf-8 -*-
"""
This script includes an easy comparison function for C. Johnson's split_ext_field, which is based off Isler's fortran zeeman code.

@author: Leo Nofs
"""

import numpy as np
import scipy as scp
import csv
from matplotlib import pyplot as plt

from sympy.physics.wigner import wigner_3j as w3j, wigner_6j as w6j


#defining a bunch of physics constants needed for the calculations.

c_light = 299792458 #Speed of light in m/s
hbar = 1.0546e-34 #Plank's constant hbar in Angular
hh = 2.0*np.pi*hbar #Plank's constant in frequency
muB = 9.27401e-24 #Bohr magneton J/T
muP = 2.79284734 #proton magneton (in terms of nuclear magnetons, nm)
muN = 1.5210322e-3 #Nucleon magneton (not nuclear dipole moment), in terms of muB
muN2 = muB/1836
mu_B = muB/hbar #Common value of Bohr_Mag/hbar which ensures the computation doesn't drop small terms. [Hz/T]
g_L = 1 #Lande Orbital g-factor, assuming infinite mass, otherwise is g_L = 1 - me/mnuc
g_S = 2.002319 #Lande electron spin g-factor
a_o = 5.29e-11 #Bohr radii in m
mu_0 = 4*np.pi*1e-7 #freespace mag permeability
kboltz = 1.381e-23 #Boltzmann Constant
m_prot = 1.667e-27 #Mass of proton in kg
m_elec = 9.11e-31
hhev = 4.135668e-15 #Plank's constant in Ev/Hz

muBEV = 5.7883818e-5 #bohr magneton in eV/Tesla
muBcm = muBEV/(1.24e-4) #bohr magneton in cm^-1

from SpectSplit import *




def sym1(a,b,bb,c,cc,s):
    return (s*(s + 1)*(s - 2*a - 1)*(s - 2*a))/\
        ((2*b - 1)*(2*b)*(2*b+1)*(2*c-1)*(2*c)*(2*c+1))

def sym2(a,b,bb,c,cc,s):
    return (2*(s + 1)*(s - 2*a)*(s - 2*b)*(s - 2*c+1))/\
        ((2*b)*(2*b + 1)*(2*b + 2)*(2*c - 1)*(2*c)*(2*c + 1))

def sym3(a,b,bb,c,cc,s):
    return((s - 2*b - 1)*(s - 2*b)*(s - 2*c + 1)*(s - 2*c + 2))/\
        ((2*b + 1)*(2*b + 2)*(2*b + 3)*(2*c - 1)*(2*c)*(2*c + 1))

def sym4(a,b,bb,c,cc,s):
    return 1/((2*b)*(2*b + 1)*(2*b + 2)*(2*c)*(2*c + 1)*(2*c + 2))

def wigner_3j(ja,jb,ma,mb,dummus):
    """ This function (adapted from Isler zeeman6.f) will calculated the wigner 3j 
        symbol used in the energy splitting as well as the relative intensities of 
        the split lines. For these cases the top right element is always 1 which is
         why this function only has 5 inputs see Cowan TASS 17.6 (shift) and 17.34 (intensities).

       ( ja    jb     1    )
       (                   )
       ( ma    mb   dummus )

       :param ja: 
       :type  ja: int array

       :param jb: 
       :type  jb: in array

       :param ma: 
       :type  ma: int array

       :param mb: 
       :type  mb: int array

       :param dummus: 
       :type  dummus: int array


       "returns: array of 3j symbol values

    """
    nph = np.zeros_like(ja)
    anum1 = jb - ma
    anum2 = jb + ma
    wig = np.zeros_like(jb)

    den0=2*jb
    den1=2*jb + 1
    den2=2*jb + 2
    den3=2*jb + 3
    den4=2*jb + 4
    den5=2*jb + 5

    con1_inds = np.where( (dummus == 1) & (ja > jb))
    wig[con1_inds] = np.sqrt((anum1[con1_inds]*(anum1[con1_inds]+1.))/\
                             (den3[con1_inds]*den2[con1_inds]*den1[con1_inds]))
    nph[con1_inds] = np.abs(anum1[con1_inds]-1)

    con2_inds = np.where( (dummus == 0) & (ja > jb))
    wig[con2_inds] = np.sqrt( ((anum2[con2_inds]+1.)*(anum1[con2_inds]+1.)*2)/\
                              (den3[con2_inds]*den2[con2_inds]*den1[con2_inds]))
    nph[con2_inds] = np.abs(anum1[con2_inds] -1)

    con3_inds = np.where( (dummus == 1) & (ja == jb))
    wig[con3_inds] = np.sqrt( (anum1[con3_inds]*(anum2[con3_inds]+1.)*2.)/\
                              (den2[con3_inds]*den1[con3_inds]*den0[con3_inds]))
    nph[con3_inds] = np.abs(anum1[con3_inds] )

    con4_inds = np.where( (dummus == 0) & (ja == jb))
    wig[con4_inds] = np.sqrt( 1./(den1[con4_inds]*(jb[con4_inds]+1.)*jb[con4_inds]))*ma[con4_inds]
    nph[con4_inds] = np.abs(anum1[con4_inds] )


    wig=wig*(-1)**(nph+1)
    wig[np.where(wig==np.inf)] = 0
    wig[ np.isnan(wig)] = 0
    return wig
#%%

# j0 = 2
# jprime = 2
# mj= -1
# mjprime = 2
# qdum = (mjprime-mj)

# jat =jprime*np.ones([3,3])
# jbt = j0*np.ones([3,3])
# mat = mjprime*np.ones([3,3])
# mbt = mj*np.ones([3,3])

# dummust = qdum*np.ones([3,3])
# w3test = wigner_3j(jat,jbt,mat,mbt,dummust)

#%%
def wigner_6j(order, a, b, bb, c, cc):
    """ This function (adapted from Isler zeeman6.f) will calculated the wigner 6j 
        symbol used in the energy splitting as well as the relative intensities of 
        the split. See Cowan TASS 17.6 (shift) and 17.34 (intensities).

       (  order  b     c    )
       (                    )
       (   a     bb    cc   )

       :param order: 
       :type  order: int array

       :param a: 
       :type  a: in array

       :param b: 
       :type  b: int array

       :param bb: 
       :type  bb: int array

       :param c: 
       :type  cc: int array


       returns: array of 6j symbol values

    """
    sym6j = np.zeros_like(a)#setup 6j matrix

    con1 = np.where( (b == bb+1) & (c == cc+1))

    sym6j[con1] = np.sqrt(sym1(a[con1], b[con1], bb[con1],
                               c[con1], cc[con1], a[con1]+b[con1]+c[con1]))*\
                               (-1)**(np.abs( a[con1] + b[con1] +c[con1])+1)

    
    con2 = np.where( (b == bb-1 ) & (c == cc - 1))
    sym6j[con2] = np.sqrt( sym1(a[con2], cc[con2], bb[con2],
                                bb[con2], cc[con2], a[con2]+b[con2]+c[con2]))
    sym6j[con2] = sym6j[con2]*(-1)**(np.abs(a[con2] + b[con2] + c[con2]) +1)

    con3 = np.where( (b == bb) & (c == cc+1))
    sym6j[con3] = np.sqrt( sym2(a[con3], b[con3], bb[con3],
                                c[con3], cc[con3], a[con3] + b[con3] + c[con3]))*\
                                (-1)**(np.abs(a[con3] + b[con3] + c[con3]) +1)

    con4 = np.where( (b == bb+1) & ( c == cc))
    sym6j[con4] = np.sqrt( sym2(a[con4], cc[con4], bb[con4],
                                bb[con4] +1,cc[con4], a[con4] + cc[con4] + (bb[con4]+1)))*\
                                (-1)**( np.abs(a[con4] + cc[con4] + (bb[con4]+1))+1)

    con5 = np.where( (b == bb-1) & ( c == cc))

    sym6j[con5] = np.sqrt( sym2(a[con5], cc[con5], bb[con5], bb[con5], cc[con5],
                                a[con5] + bb[con5] + cc[con5]))*\
                                (-1)*(np.abs(a[con5] + bb[con5] + cc[con5])+1)
    

    con6 = np.where( (b == bb) & ( c == cc-1))
    sym6j[con6] = np.sqrt( sym2(a[con6], bb[con6], bb[con6], cc[con6], cc[con6],
                                a[con6] + bb[con6] + cc[con6]))*\
                                (-1)**(np.abs(a[con6] + bb[con6] + cc[con6])+1)

    con7 = np.where( (b == bb-1) & (c == c+1))
    sym6j[con7] = np.sqrt( sym3(a[con7],b[con7],bb[con7],c[con7],cc[con7],
                                a[con7] + b[con7] + c[con7]))*\
                                (-1)**(np.abs(a[con7] + b[con7] + c[con7])+1)

    con8 = np.where( (b == bb+1) & ( c == c-1))
    
    sym6j[con8] = np.sqrt( sym3(a[con8], cc[con8] -1, bb[con8], bb[con8] +1, cc[con8],
                                a[con8] + (cc[con8] -1) + (bb[con8] -1)))*\
                                (-1)**(np.abs(a[con8] + (cc[con8] -1) + (bb[con8] -1))+1)

    con9 = np.where( (b == bb-1) & ( c == cc+1))
    sym6j[con9] = np.sqrt( sym3(a[con9],cc[con9],bb[con9], bb[con9], cc[con9],
                                a[con9] + cc[con9] + bb[con9]))*\
                                (-1)**(np.abs(a[con9] + cc[con9] + bb[con9])+1)

    con10 = np.where( (b == bb+1) & ( c == cc-1))
    sym6j[con10] = np.sqrt( sym3(a[con10],bb[con10],bb[con10],cc[con10],cc[con10],
                                 a[con10] + bb[con10] + cc[con10]))*\
                                 (-1)**(np.abs(a[con10] + bb[con10] + cc[con10])+1)

    con11 = np.where( (b == bb) & ( c ==cc))
    sym6j[con11] = np.sqrt( sym4(a[con11], b[con11], bb[con11], c[con11], cc[con11],
                           a[con11] + b[con11] + c[con11]))*2*(b[con11]*(b[con11]+1)+\
                                         c[con11]*(c[con11]+1.)-a[con11]*(a[con11]+1.))*\
                                         (-1)**(np.abs(a[con11] + b[con11] +c[con11] +1) +1)
    

    sym6j[np.where(sym6j==np.inf)] = 0#get rid of intinfites and nans
    sym6j[ np.isnan(sym6j)] = 0

    return sym6j



def split_ext_field_clean(lvl_eng, bfield,orb,spin, b_angle = 90, gamma = 0,low_field_limit_remove_off_diag = False,rel_j_intensities=np.array([]), Filter=False ):
    """
    The equations governing the perturbation due to a magnetic field (Zeeman/Paschen-Back) effect
    are relatively simple. The splitting of the levels is calculated with Cowan 17.6 or Isler Phys. Plasma 1997 appendix
    <SLJ''M | H' | SLJM>.

    Then relative values of the intensities between the newly split levels can then be calculated with Cownan 17.34
    or with Isler equation. The Cowan discuss is very general referencing equation 14.78 which is way to general.
    The expression from Cowan Journal of the Optical Society of America Volume 58. Number 6 (1968) gives the simplification
    to transitions of interest for magnetic fusion equation 40. 


    These two equations are relatively simple. Most of the code is constructing the various bra and ket matrices
    as well as the Wigner 3 and 6 j symbols.


    This code follows the Isler Zeeman code written in Fortran. This implimentation is a pure python version which will 
    hopefully be easier to maintain and modify from the orginal Fortran. This implementation takes advantage of matrix operations
    in python making the code easier to map to equations from the papers. 


    Expansions will include non-dipole transitions (import for Tungsten)


    Low field limit to get back to the lande g-factor, take only diag of bmtx Cowan 17.8, 17.9
    




    :param lvl_eng: energy of upper and lower levels in the term of the transition of interest [cm-1]
    :param lvl_eng: list of float arrays

    :param bfield: magnetic field [gauss]
    :type bfield: scalar

    :param orb:Orbital angular momentum "L quantum number" upper then lower
    :type orb: in array

    :param spin: spin of the upper and lower (for dipole must be the same) code can only handle dipole currently
    :param spin: int array

    :param b_angle: angle the magnetic field makes with the line of sight [deg]
    :param b_angle: scalar

    :param gamma: angle in the direction having a transimission T1 makes with the magnetic field if theta = 90 [deg]
    :type  gamma: scalar

    :param low_field_limit_remove_off_diag: get rid of the off diagnonal terms if requested
    :type  low_field_limit_remove_off_diag: bool

    :param rel_j_intensities: relative intensities of the j levels in the LS term. All equal if not specified
    :type  rel_j_intensities: float array length of upper lvl_eng

    """

    bohr=4.66858e-5 #bohr magniton in pseudo energy units [cm-1 /gauss]

    #create list for the upper and lower energies
    #The energy of a level is referenced to the lowest energy of the term
    energy = []
    energy.append( np.zeros_like(lvl_eng[0]))
    energy.append( np.zeros_like(lvl_eng[1]))

    for j in range(0,len(energy[0])):
        energy[0][j]=lvl_eng[0][j]-lvl_eng[0][-1]

    for j in range(0,len(energy[1])):
        energy[1][j]=lvl_eng[1][j]-lvl_eng[1][-1]

    gamma_angle = np.deg2rad(gamma)
    theta = np.deg2rad(b_angle)#convert deg to radians

    numup = int((2*orb[0]+1)*(2*spin[0]+1))#number of uper levels
    numlow = int((2*orb[1]+1)*(2*spin[0]+1))#number of lower levels

    num_levels = np.array([int((2*orb[0]+1)*(2*spin[0]+1)),
                           int((2*orb[1]+1)*(2*spin[0]+1))])

    nnu=(orb[0]+spin[0])-np.abs(orb[0]-spin[0])+1
    nnl=(orb[1]+spin[0])-np.abs(orb[1]-spin[0])+1

    lcomp = np.zeros(2,dtype='int')
    lcomp[0] = nnu
    lcomp[1] = nnl
    if Filter==False:
        pol1 = 1
        pol2 = 1
    if Filter==True:
        pol1 = 1
        pol2 = 0
       
    spol1=np.sqrt(pol1)
    spol2=np.sqrt(pol2)

    epsz=np.sin(theta)*(spol1*np.cos(gamma)+spol2*np.sin(gamma))
    epsx=-1.*(spol1*np.sin(gamma)-spol2*np.cos(gamma))/np.sqrt(2)
    epsy=-np.cos(theta)*(spol1*np.cos(gamma)+spol2*np.sin(gamma))/np.sqrt(2)

    j_ul = [np.arange(np.abs(orb[0]-spin[0]),np.abs(orb[0]+spin[0])+1.,1.)[::-1],#upper j values
           np.arange(np.abs(orb[1]-spin[0]),np.abs(orb[1]+spin[0])+1.,1.)[::-1]]#lower j values

    wig3 =  np.zeros((numup,numup))

    mjs = [] #holds the possible mj values list upper then lower
    js =  [] #holds the possible j values in the same shape as the mjs
    mjsl = []
    jsl = []
    #find the possible mj's and organize them so matrix can be constructed
    js_inds = []

    for j in range(0,2):#loop over the upper then lower term
        tot = 0
        j_arr = np.zeros((num_levels[j]))

        j_arr_inds = np.zeros((num_levels[j]))
        jl_arr = np.zeros((num_levels[j]))
        mj_arr = np.zeros((num_levels[j]))
        mjl_arr = np.zeros((num_levels[j]))    
        for i in range(0,len(j_ul[j])):

            tmp = np.arange(-j_ul[j][i],j_ul[j][i]+1.,1.)#find the possible mj values
            mj_arr[tot:tot+len(tmp)] =  tmp[::-1]
            mjl_arr[tot:tot+len(tmp)] =tmp[::-1]
            j_arr[tot:tot+len(tmp)] = np.ones(len(tmp))*j_ul[j][i]# the mjs correspond to one j
            j_arr_inds[tot:tot+len(tmp)] = i
            tot = tot+len(tmp)
        mjs.append(mj_arr) #append the mj values for the upper or lower
        mjsl.append(mjl_arr) #append the mj values for the upper or lower    
        js.append(j_arr)#append the j values for the upper or lower



    mbra = []#mbra list upper shape (mj_up x mj_up) then lower (mj_low xmj_low)
    mket = []#same shape as mbra just the ket
    jbra = []#same shape as mbra just the jbra
    jket = []#same shape as mbra just the jket
    alpha = []
    phase1 = []
    phase2 = []
    for i in range(0,2):

        mbra_arr = np.zeros((num_levels[i],num_levels[i]))
        mket_arr = np.zeros((num_levels[i],num_levels[i]))
        jbra_arr = np.zeros((num_levels[i],num_levels[i]))
        jket_arr = np.zeros((num_levels[i],num_levels[i]))

        for j in range(0,num_levels[i]):
            mbra_arr[j,:] = mjs[i][j]
            mket_arr[:,j] = mjsl[i][j]

            jbra_arr[j,:] = js[i][j]#j_arr
            jket_arr[:,j] = js[i][j]#j_arr

        mbra.append(mbra_arr)# holds mj bra arrays for upper and lower
        mket.append(mket_arr)# holds mj k arrays for upper and lower

        jbra.append(jbra_arr)# holds j bra arrays for upper and lower
        jket.append(jket_arr)# holds j ket arrays for upper and lower

        alpha.append(jbra_arr + jket_arr +1)
        phase1.append(jbra_arr - mbra_arr + alpha[i])
        phase2.append(np.zeros_like(jket_arr))

    wig3 = []#holds the 3j values for the upper and lower
    wig6 = []#holds the 6j values for the upper and lower

    bmtx = []
    #now that the bra and kets have been constructed need to calculate the
    #3 and 6 j symbols so the perturbations of the hamiltonian can be calculated
    #for both the upper and lower levels (see Cowan 17.6 or Isler appendix Phys. Plasmas 5 (2) 1997
    for i in range(0,2):
        wig3_arr = np.zeros((num_levels[i],num_levels[i]))
        wig6_arr = np.zeros((num_levels[i],num_levels[i]))

        ja = np.copy(jbra[i])
        jb = np.copy(jket[i])

        ma = -1*np.copy(mbra[i])
        mb = np.copy(mket[i])

        dummus = np.zeros_like(ja)

        swap_inds = np.where( jket[i] > jbra[i])
        if(swap_inds):
            ja[swap_inds] = jket[i][swap_inds]
            jb[swap_inds] = jbra[i][swap_inds]

            ma[swap_inds] = mket[i][swap_inds]
            mb[swap_inds] = -1*mbra[i][swap_inds]
            phase2[i][swap_inds] = alpha[i][swap_inds]

        wig3_arr  = wigner_3j(ja,jb,ma,mb,dummus)*\
                            (-1)**(np.abs(phase1[i] + phase2[i]) + 1)

        if(np.where( ma != -1*mb)):
            wig3_arr[np.where( ma != -1*mb)] = 0.

        if( np.where(np.abs(ja-jb) >1)):
            wig3_arr[np.where(np.abs(ja-jb) >1)] = 0.

        if( np.where( (ja == 0) & (jb == 0))):
            wig3_arr[np.where( (ja == 0) & (jb == 0))] = 0.

        if( spin[0] == 0):
            wig3_arr[:] = 0

        if( np.where( (ja ==jb) & (ma == 0) & (mb == 0))):
            wig3_arr[np.where( (ja ==jb) & (ma == 0) & (mb == 0))] = 0.

        phase3 = np.ones_like(jket[i])*orb[i] + np.ones_like(jket[i])*spin[0] +\
                               jket[i] + 1

        order_arr = np.ones_like(jb)
        a = np.ones_like(jb)*orb[i]
        b = np.ones_like(jb)*spin[0]
        bb = np.ones_like(jb)*spin[0]
        c = ja
        cc = jb

        wig6_arr = wigner_6j(order_arr, a, b, bb, c, cc)*(-1)**(np.abs(phase3)+1)

        wig3.append(wig3_arr)
        wig6.append(wig6_arr)

        #apply Wigner-Eckart to the magnetic field perturbation hamiltonian
        #see Cowan 17.6 or Isler 1997 appendix
        # second part of the equation
        #(g_s -1 =1.002319) anomalous mag dipole moment        
        bmtx.append( wig3_arr*wig6_arr*np.sqrt((2*jbra[i] +1)*(2*jket[i] +1))*\
               np.sqrt((2*np.ones_like(jket[i])*spin[0]+1)*\
        (np.ones_like(jket[i])*spin[0] +1)*np.ones_like(jket[i])*spin[0])*1.002319)


        #see Cowan 17.6 first part of equation
        bmtx[i][np.diag_indices_from(bmtx[i])] = bmtx[i][np.diag_indices_from(bmtx[i])] +\
                                                mbra[i][np.diag_indices_from(bmtx[i])]#(mj)(delt_jj) 

        bmtx[i] = bmtx[i]*bohr*bfield  #Multiply the two parts of the 17.6 equation by (beta)(B_0)

        j_energy = np.zeros_like(np.diag(jket[i]))
        for k in range(0,len(j_ul[i])):
            j_energy[np.where(np.diag(jket[i]) == j_ul[i][k])[0]] = energy[i][k]

        #energy of the level without the magnetic field perturbation (E_j)(delta_jj)
        bmtx[i][np.diag_indices_from(bmtx[i])] = bmtx[i][np.diag_indices_from(bmtx[i])] + \
                                                     j_energy

    #In the low field limit the off-diagonal elements can be ignored
    if( low_field_limit_remove_off_diag):
        bmtx_tmp = []
        bmtx_tmp.append( np.zeros_like(bmtx[0]))#create new zeros matrix cause only want diag
        bmtx_tmp[0][np.diag_indices_from(bmtx[0])] = np.diag(bmtx[0])#keep the diagnonal terms calculated

        bmtx_tmp.append( np.zeros_like(bmtx[1]))
        bmtx_tmp[1][np.diag_indices_from(bmtx[1])] = np.diag(bmtx[1])
        bmtx = bmtx_tmp#update the bmtx matrix to

    eig1,vec1 = np.linalg.eig(bmtx[0])#eigen vals and vecs of the upper term
    eig2,vec2 = np.linalg.eig(bmtx[1])#eigen vals and vecs of the lower term

    split_energy_up = lvl_eng[0][-1] + eig1#modification of upper levels referenced to the lowest level of the term
    split_energy_low = lvl_eng[1][-1] + eig2#modification of lower levels referenced to the lowest level of the term

    wavevac = 1/(split_energy_up[:,None] - split_energy_low[None,:])#Standard def of wavelength (1/(up - low)) [cm]

    #create new bra and kets for the transition intensity calculation
    jbra_inten = np.zeros((len(jbra[0])*len(jbra[1]),len(jbra[0])*len(jbra[1])))
    mbra_inten = np.zeros((len(jbra[0])*len(jbra[1]),len(jbra[0])*len(jbra[1])))
    jket_inten = np.zeros((len(jbra[0])*len(jbra[1]),len(jbra[0])*len(jbra[1])))
    mket_inten = np.zeros((len(jbra[0])*len(jbra[1]),len(jbra[0])*len(jbra[1])))

    #populate the jbra and kets from the upper and lower levels
    #the bra is the upper and the ket is the lower
    for i in range(0,len(jbra[1])):
        jbra_inten[:,i::len(jbra[1])] = jbra[0][:,0]

    jbra_inten = jbra_inten.T

    for i in range(0,len(jbra[0])):
        jket_inten[:,len(jbra[1])*i:len(jbra[1])*(i+1)] = jket[1][0,:]

    jket_inten = jket_inten.T

    for i in range(0,len(mbra[1])):
        mbra_inten[:,i::len(mbra[1])] = mbra[0][:,0]

    mbra_inten = mbra_inten.T

    for i in range(0,len(mbra[0])):
        mket_inten[:,len(jbra[1])*i:len(jbra[1])*(i+1)] = mket[1][0,:]

    mket_inten = mket_inten.T


    q = mbra_inten - mket_inten #Cowan 14.29
    am = -1*mbra_inten
    bm = np.copy(mket_inten)
    dummus = np.copy(q)
    phase2 = np.zeros_like(mket_inten)
    phase1 = np.zeros_like(mket_inten)

    phase1 = q + jbra_inten - mbra_inten + jbra_inten + jket_inten +1

    neg_q_inds = np.where( q < 0)

    phase2[neg_q_inds] = jbra_inten[neg_q_inds] + jket_inten[neg_q_inds] + 1
    am[neg_q_inds] = mbra_inten[neg_q_inds]
    bm[neg_q_inds] = -1*mket_inten[neg_q_inds]
    dummus[neg_q_inds] = -1*q[neg_q_inds]

    ja = np.copy(jket_inten)
    ma = np.copy(bm)
    jb = np.copy(jbra_inten)
    mb = np.copy(am)
    phase3 = jbra_inten + jket_inten + 1

    jbra_greater = np.where(jbra_inten >= jket_inten)
    ja[jbra_greater] = jbra_inten[jbra_greater]
    ma[jbra_greater] = am[jbra_greater]
    jb[jbra_greater] = jket_inten[jbra_greater]
    mb[jbra_greater] = bm[jbra_greater]
    phase3[jbra_greater] = 0

    wig3_arr_inten = wigner_3j(ja,jb,ma,mb,dummus)*\
                     (-1)**(np.abs(phase1+phase2+phase3)+1)

    wig6_arr_inten = wigner_6j(np.ones_like(jbra_inten),
                               np.ones_like(jbra_inten)*spin[0],
                               np.ones_like(jbra_inten)*orb[0],
                               np.ones_like(jbra_inten)*orb[1],
                               jbra_inten,
                               jket_inten)*(-1)**(orb[1] + spin[0]+ jbra_inten + 2)


    afcf = np.kron(vec1[:,::-1],vec2[:,::-1])

    deltam = np.abs(mbra_inten - mket_inten)

    rad = np.ones_like(wig3_arr_inten)
    for i in range(0,len(rad)):
        for j in range(0,len(rad)):

            if(deltam[i,j] ==1):
                rad_tmp = np.sqrt(epsx**2 + epsy**2)
            if(deltam[i,j] == 0):
                rad_tmp = epsz
            rad[i,j] = rad_tmp


    #Cowan equation 17.34 Note that this is different from Isler 1997 and I am not sure
    #what the resolution is but it seems to be the same. Iser has (-1)^(3J' -m_j' +q)
    #corresponding to 3Jbra - mket. Note this is different from Cowan Jbra - mbra.
    #However, when this gets squared the negative sign does not matter so it is unclear to me at this
    #point if there is case where this difference actually matters.
    #Note that 3Jbra and Jbra are equilvalent with (-1)^(x).
    signal = np.sum(np.sqrt((2*jbra_inten+1)*(2*jket_inten+1))*afcf*wigner_3j(ja,jb,ma,mb,dummus)*\
                     wigner_6j(np.ones_like(jbra_inten),
                               np.ones_like(jbra_inten)*spin[0],
                               np.ones_like(jbra_inten)*orb[0],
                               np.ones_like(jbra_inten)*orb[1],
                               jbra_inten,
                               jket_inten)*\
                    (-1)**(phase1+phase2+phase3 + orb[1] +spin[1] + jbra_inten +2)*rad,axis=0)**2

    
    #If the relative intentieis of the j levels in the term were specified, weight the result by these
    if(rel_j_intensities.size != 0):

        qt = np.zeros_like(js[0])
        for i in range(0,len(rel_j_intensities)):
            qt[np.where(js[0] == j_ul[0][i])[0]] = rel_j_intensities[i]#create rel inten array from the j array

        qq = np.sum(np.einsum('i,ij->ij',qt,vec1[:,::-1])**2,axis=0)#find components from the eigen vectors

        st = []
        for i in range(0,len(qq)):
            st.append(np.ones(len(jbra[1]))*qq[i])#has to be the same pattern as the kron
                                                  #a1b1 a1b2 a1bn....anb1 anb2 anb3

        weights = np.asarray(st).flatten()#weights in
        signal = signal*weights



    #find mj sublevel purity and mjs labels
    #You cant map a single eigenvalue to a mj sublevel
    #this calculates the purity of the state
    mjs_pur = [] #hold the purity for both upper and lower
    mjs_label = []#hold the label or leading mjs

    #start with the upper
    mj_pur = np.zeros_like(vec1)
    mj_label = np.zeros_like(mjs[0],dtype='int')
    for i in range(0,len(mjs[0])):#loop over the mjs
        rr = np.zeros_like(mjs[0])
        rr[i] = 1 #only look at a particular mj value
        mj_pur[i,:] = np.sum(np.einsum('i,ij->ij',rr,vec1)**2,axis=0)
        #Shows fraction of eigenstate that is made up of mjs sublevel        
        mj_label[i] = np.argmax(mj_pur[i])#put a label on the largest value and call that associated mj
    
    mjs_pur.append(mj_pur)
    mjs_label.append(mj_label)

    #now find the purity of the lower
    mj_pur = np.zeros_like(vec2)
    mj_label = np.zeros_like(mjs[1],dtype='int')
    for i in range(0,len(mjs[1])):
        rr = np.zeros_like(mjs[1])
        rr[i] = 1
        mj_pur[i,:] = np.sum(np.einsum('i,ij->ij',rr,vec2)**2,axis=0)
        mj_label[i] = np.argmax(mj_pur[i])

    mjs_pur.append(mj_pur)
    mjs_label.append(mj_label)

    
    wavevac = wavevac.flatten()[::-1]#flatten and reverse the array so its in the same order as signal
    #convert the vacuum wavelength to air
    air_index = 1.+2.94981e6/(146.e8-(1./(wavevac**2)))+\
                  2.554e4/(41.e8-(1./(wavevac**2)))+6.4328e-5
    wavelengths= wavevac / air_index
    wavelengths = wavelengths * 1E7 #unit change to nm

    data = {}
    data['wave_air'] = wavelengths #wavelength in air (nm)
    data['wave_vac'] = wavevac*1E7 #wavelength in vacuum (nm)
    data['signal'] = signal
    data['split_energy_up'] = split_energy_up
    data['split_energy_low'] = split_energy_low
    data['eig_vec1'] = vec1 #eigen vectors of the upper 
    data['eig_vec2'] = vec2#eigen vectors of the lower
    data['mjs_purity'] = mjs_pur #purity fraction
    data['mjs_label'] = mjs_label#label in same order as mjs and js
    data['mjs'] = mjs#mj values both upper and lower
    data['js'] = js# j values both upper and lower
    return data


def split_ext_field_mod(lvl_eng, bfield,orb,spin, b_angle = 90, gamma = 0,low_field_limit_remove_off_diag = False,rel_j_intensities=np.array([]),Filter=False, T_Alpha = 1, T_Beta = 0):
    """
    The equations governing the perturbation due to a magnetic field (Zeeman/Paschen-Back) effect
    are relatively simple. The splitting of the levels is calculated with Cowan 17.6 or Isler Phys. Plasma 1997 appendix
    <SLJ''M | H' | SLJM>.

    Then relative values of the intensities between the newly split levels can then be calculated with Cownan 17.34
    or with Isler equation. The Cowan discuss is very general referencing equation 14.78 which is way to general.
    The expression from Cowan Journal of the Optical Society of America Volume 58. Number 6 (1968) gives the simplification
    to transitions of interest for magnetic fusion equation 40. 


    These two equations are relatively simple. Most of the code is constructing the various bra and ket matrices
    as well as the Wigner 3 and 6 j symbols.


    This code follows the Isler Zeeman code written in Fortran. This implimentation is a pure python version which will 
    hopefully be easier to maintain and modify from the orginal Fortran. This implementation takes advantage of matrix operations
    in python making the code easier to map to equations from the papers. 


    Expansions will include non-dipole transitions (import for Tungsten)


    Low field limit to get back to the lande g-factor, take only diag of bmtx Cowan 17.8, 17.9
    




    :param lvl_eng: energy of upper and lower levels in the term of the transition of interest [cm-1]
    :param lvl_eng: list of float arrays

    :param bfield: magnetic field [gauss]
    :type bfield: scalar

    :param orb:Orbital angular momentum "L quantum number" upper then lower
    :type orb: in array

    :param spin: spin of the upper and lower (for dipole must be the same) code can only handle dipole currently
    :param spin: int array

    :param b_angle: angle the magnetic field makes with the line of sight [deg]
    :param b_angle: scalar

    :param gamma: angle in the direction having a transimission T1 makes with the magnetic field if theta = 90 [deg]
    :type  gamma: scalar

    :param low_field_limit_remove_off_diag: get rid of the off diagnonal terms if requested
    :type  low_field_limit_remove_off_diag: bool

    :param rel_j_intensities: relative intensities of the j levels in the LS term. All equal if not specified
    :type  rel_j_intensities: float array length of upper lvl_eng

    """

    bohr=4.66858e-5 #bohr magniton in pseudo energy units [cm-1 /gauss]

    #create list for the upper and lower energies
    #The energy of a level is referenced to the lowest energy of the term
    energy = []
    energy.append( np.zeros_like(lvl_eng[0]))
    energy.append( np.zeros_like(lvl_eng[1]))

    for j in range(0,len(energy[0])):
        energy[0][j]=lvl_eng[0][j]-lvl_eng[0][-1]

    for j in range(0,len(energy[1])):
        energy[1][j]=lvl_eng[1][j]-lvl_eng[1][-1]



    numup = int((2*orb[0]+1)*(2*spin[0]+1))#number of uper levels
    numlow = int((2*orb[1]+1)*(2*spin[0]+1))#number of lower levels

    num_levels = np.array([int((2*orb[0]+1)*(2*spin[0]+1)),
                           int((2*orb[1]+1)*(2*spin[0]+1))])

    nnu=(orb[0]+spin[0])-np.abs(orb[0]-spin[0])+1
    nnl=(orb[1]+spin[0])-np.abs(orb[1]-spin[0])+1

    lcomp = np.zeros(2,dtype='int')
    lcomp[0] = nnu
    lcomp[1] = nnl
   
    if Filter==False:
        pol1 = 1
        pol2 = 1
        ealph = 1
        ebeta = 1
    else:
        pol1=T_Alpha
        pol2=T_Beta
        gammaangle = np.deg2rad(gamma)
        ealph = pol1*np.cos(gammaangle) + pol2*np.sin(gammaangle)
        ebeta = pol1*np.sin(gammaangle) + pol2*np.cos(gammaangle)
    theta_rad = np.deg2rad(b_angle)#convert deg to radians
   
    
    epsz = ealph*np.sin(theta_rad)
    
    e0 = epsz
    
    eps00 = np.sqrt(e0*e0)
    eps11 = 0.5*np.sqrt((np.cos(theta_rad)**2 * ealph**2 + ebeta**2 ))
    
    
   
    j_ul = [np.arange(np.abs(orb[0]-spin[0]),np.abs(orb[0]+spin[0])+1.,1.)[::-1],#upper j values
           np.arange(np.abs(orb[1]-spin[1]),np.abs(orb[1]+spin[1])+1.,1.)[::-1]]#lower j values

    wig3 =  np.zeros((numup,numup))

    mjs = [] #holds the possible mj values list upper then lower
    js =  [] #holds the possible j values in the same shape as the mjs
    mjsl = []
    jsl = []
    #find the possible mj's and organize them so matrix can be constructed
    js_inds = []

    for j in range(0,2):#loop over the upper then lower term
        tot = 0
        j_arr = np.zeros((num_levels[j]))

        j_arr_inds = np.zeros((num_levels[j]))
        jl_arr = np.zeros((num_levels[j]))
        mj_arr = np.zeros((num_levels[j]))
        mjl_arr = np.zeros((num_levels[j]))    
        for i in range(0,len(j_ul[j])):

            tmp = np.arange(-j_ul[j][i],j_ul[j][i]+1.,1.)#find the possible mj values
            mj_arr[tot:tot+len(tmp)] =  tmp[::-1]
            mjl_arr[tot:tot+len(tmp)] =tmp[::-1]
            j_arr[tot:tot+len(tmp)] = np.ones(len(tmp))*j_ul[j][i]# the mjs correspond to one j
            j_arr_inds[tot:tot+len(tmp)] = i
            tot = tot+len(tmp)
        mjs.append(mj_arr) #append the mj values for the upper or lower
        mjsl.append(mjl_arr) #append the mj values for the upper or lower    
        js.append(j_arr)#append the j values for the upper or lower



    mbra = []#mbra list upper shape (mj_up x mj_up) then lower (mj_low xmj_low)
    mket = []#same shape as mbra just the ket
    jbra = []#same shape as mbra just the jbra
    jket = []#same shape as mbra just the jket
    alpha = []
    phase1 = []
    phase2 = []
    for i in range(0,2):

        mbra_arr = np.zeros((num_levels[i],num_levels[i]))
        mket_arr = np.zeros((num_levels[i],num_levels[i]))
        jbra_arr = np.zeros((num_levels[i],num_levels[i]))
        jket_arr = np.zeros((num_levels[i],num_levels[i]))

        for j in range(0,num_levels[i]):
            mbra_arr[j,:] = mjs[i][j]
            mket_arr[:,j] = mjsl[i][j]

            jbra_arr[j,:] = js[i][j]#j_arr
            jket_arr[:,j] = js[i][j]#j_arr

        mbra.append(mbra_arr)# holds mj bra arrays for upper and lower
        mket.append(mket_arr)# holds mj k arrays for upper and lower

        jbra.append(jbra_arr)# holds j bra arrays for upper and lower
        jket.append(jket_arr)# holds j ket arrays for upper and lower

        alpha.append(jbra_arr + jket_arr +1)
        phase1.append(jbra_arr - mbra_arr + alpha[i])
        phase2.append(np.zeros_like(jket_arr))

    wig3 = []#holds the 3j values for the upper and lower
    wig6 = []#holds the 6j values for the upper and lower

    bmtx = []
    #now that the bra and kets have been constructed need to calculate the
    #3 and 6 j symbols so the perturbations of the hamiltonian can be calculated
    #for both the upper and lower levels (see Cowan 17.6 or Isler appendix Phys. Plasmas 5 (2) 1997
    for i in range(0,2):
        wig3_arr = np.zeros((num_levels[i],num_levels[i]))
        wig6_arr = np.zeros((num_levels[i],num_levels[i]))

        ja = np.copy(jbra[i])
        jb = np.copy(jket[i])

        ma = -1*np.copy(mbra[i])
        mb = np.copy(mket[i])

        dummus = np.zeros_like(ja)

        swap_inds = np.where( jket[i] > jbra[i])
        if(swap_inds):
            ja[swap_inds] = jket[i][swap_inds]
            jb[swap_inds] = jbra[i][swap_inds]

            ma[swap_inds] = mket[i][swap_inds]
            mb[swap_inds] = -1*mbra[i][swap_inds]
            phase2[i][swap_inds] = alpha[i][swap_inds]

        wig3_arr  = wigner_3j(ja,jb,ma,mb,dummus)*\
                            (-1)**(np.abs(phase1[i] + phase2[i]) + 1)

        if(np.where( ma != -1*mb)):
            wig3_arr[np.where( ma != -1*mb)] = 0.

        if( np.where(np.abs(ja-jb) >1)):
            wig3_arr[np.where(np.abs(ja-jb) >1)] = 0.

        if( np.where( (ja == 0) & (jb == 0))):
            wig3_arr[np.where( (ja == 0) & (jb == 0))] = 0.

        if( spin[0] == 0):
            wig3_arr[:] = 0

        if( np.where( (ja ==jb) & (ma == 0) & (mb == 0))):
            wig3_arr[np.where( (ja ==jb) & (ma == 0) & (mb == 0))] = 0.

        phase3 = np.ones_like(jket[i])*orb[i] + np.ones_like(jket[i])*spin[0] +\
                               jket[i] + 1

        order_arr = np.ones_like(jb)
        a = np.ones_like(jb)*orb[i]
        b = np.ones_like(jb)*spin[0]
        bb = np.ones_like(jb)*spin[0]
        c = ja
        cc = jb

        wig6_arr = wigner_6j(order_arr, a, b, bb, c, cc)*(-1)**(np.abs(phase3)+1)

        wig3.append(wig3_arr)
        wig6.append(wig6_arr)

        #apply Wigner-Eckart to the magnetic field perturbation hamiltonian
        #see Cowan 17.6 or Isler 1997 appendix
        # second part of the equation
        #(g_s -1 =1.002319) anomalous mag dipole moment        
        bmtx.append( wig3_arr*wig6_arr*np.sqrt((2*jbra[i] +1)*(2*jket[i] +1))*\
                np.sqrt((2*np.ones_like(jket[i])*spin[0]+1)*\
        (np.ones_like(jket[i])*spin[0] +1)*np.ones_like(jket[i])*spin[0])*1.002319)

        # bmtx.append( wig3_arr*wig6_arr*np.sqrt((2*jbra[i] +1)*(2*jket[i] +1))*\
        #        np.sqrt((2*np.ones_like(jket[i])*spin[0]+1)*\
        # (np.ones_like(jket[i])*spin[0] +1)*np.ones_like(jket[i])*spin[0])*2.002319)



        #see Cowan 17.6 first part of equation
        bmtx[i][np.diag_indices_from(bmtx[i])] = bmtx[i][np.diag_indices_from(bmtx[i])] +\
                                                mbra[i][np.diag_indices_from(bmtx[i])]#(mj)(delt_jj) 

        bmtx[i] = bmtx[i]*bohr*bfield  #Multiply the two parts of the 17.6 equation by (beta)(B_0)

        j_energy = np.zeros_like(np.diag(jket[i]))
        for k in range(0,len(j_ul[i])):
            j_energy[np.where(np.diag(jket[i]) == j_ul[i][k])[0]] = energy[i][k]

        #energy of the level without the magnetic field perturbation (E_j)(delta_jj)
        bmtx[i][np.diag_indices_from(bmtx[i])] = bmtx[i][np.diag_indices_from(bmtx[i])] + \
                                                     j_energy

    #In the low field limit the off-diagonal elements can be ignored
    if( low_field_limit_remove_off_diag):
        bmtx_tmp = []
        bmtx_tmp.append( np.zeros_like(bmtx[0]))#create new zeros matrix cause only want diag
        bmtx_tmp[0][np.diag_indices_from(bmtx[0])] = np.diag(bmtx[0])#keep the diagnonal terms calculated

        bmtx_tmp.append( np.zeros_like(bmtx[1]))
        bmtx_tmp[1][np.diag_indices_from(bmtx[1])] = np.diag(bmtx[1])
        bmtx = bmtx_tmp#update the bmtx matrix to

    eig1,vec1 = np.linalg.eig(bmtx[0])#eigen vals and vecs of the upper term
    eig2,vec2 = np.linalg.eig(bmtx[1])#eigen vals and vecs of the lower term

    split_energy_up = lvl_eng[0][-1] + eig1#modification of upper levels referenced to the lowest level of the term
    split_energy_low = lvl_eng[1][-1] + eig2#modification of lower levels referenced to the lowest level of the term

    wavevac = 1/(split_energy_up[:,None] - split_energy_low[None,:])#Standard def of wavelength (1/(up - low)) [cm]

    #create new bra and kets for the transition intensity calculation
    jbra_inten = np.zeros((len(jbra[0])*len(jbra[1]),len(jbra[0])*len(jbra[1])))
    mbra_inten = np.zeros((len(jbra[0])*len(jbra[1]),len(jbra[0])*len(jbra[1])))
    jket_inten = np.zeros((len(jbra[0])*len(jbra[1]),len(jbra[0])*len(jbra[1])))
    mket_inten = np.zeros((len(jbra[0])*len(jbra[1]),len(jbra[0])*len(jbra[1])))

    #populate the jbra and kets from the upper and lower levels
    #the bra is the upper and the ket is the lower
    for i in range(0,len(jbra[1])):
        jbra_inten[:,i::len(jbra[1])] = jbra[0][:,0]

    jbra_inten = jbra_inten.T

    for i in range(0,len(jbra[0])):
        jket_inten[:,len(jbra[1])*i:len(jbra[1])*(i+1)] = jket[1][0,:]

    jket_inten = jket_inten.T

    for i in range(0,len(mbra[1])):
        mbra_inten[:,i::len(mbra[1])] = mbra[0][:,0]

    mbra_inten = mbra_inten.T

    for i in range(0,len(mbra[0])):
        mket_inten[:,len(jbra[1])*i:len(jbra[1])*(i+1)] = mket[1][0,:]

    mket_inten = mket_inten.T

    #The dummus and relation to q has some errors
    q = mbra_inten - mket_inten #Cowan 14.29
    am = -1*mbra_inten
    bm = np.copy(mket_inten)
    dummus = np.copy(q)
    phase2 = np.zeros_like(mket_inten)
    phase1 = np.zeros_like(mket_inten)

    phase1 = q + jbra_inten - mbra_inten + jbra_inten + jket_inten +1

    neg_q_inds = np.where( q < 0) 

    phase2[neg_q_inds] = jbra_inten[neg_q_inds] + jket_inten[neg_q_inds] + 1
    am[neg_q_inds] = mbra_inten[neg_q_inds]
    bm[neg_q_inds] = -1*mket_inten[neg_q_inds]
    dummus[neg_q_inds] = -1*q[neg_q_inds]

    ja = np.copy(jket_inten)
    ma = np.copy(bm)
    jb = np.copy(jbra_inten)
    mb = np.copy(am)
    phase3 = jbra_inten + jket_inten + 1

    jbra_greater = np.where(jbra_inten >= jket_inten)
    ja[jbra_greater] = jbra_inten[jbra_greater]
    ma[jbra_greater] = am[jbra_greater]
    jb[jbra_greater] = jket_inten[jbra_greater]
    mb[jbra_greater] = bm[jbra_greater]
    phase3[jbra_greater] = 0

    wig3_arr_inten = wigner_3j(ja,jb,ma,mb,dummus)*\
                     (-1)**(np.abs(phase1+phase2+phase3)+1)

    wig6_arr_inten = wigner_6j(np.ones_like(jbra_inten),
                               np.ones_like(jbra_inten)*spin[0],
                               np.ones_like(jbra_inten)*orb[0],
                               np.ones_like(jbra_inten)*orb[1],
                               jbra_inten,
                               jket_inten)*(-1)**(orb[1] + spin[0]+ jbra_inten + 2)


    afcf = np.kron(vec1[:,::-1],vec2[:,::-1])

    deltam = np.abs(mbra_inten - mket_inten)

    rad = np.zeros_like(wig3_arr_inten)
    for i in range(0,len(rad)):
        rad_tmp = 0
        for j in range(0,len(rad)):
            rad_tmp=0
            if(deltam[i,j] ==1):
                rad_tmp = eps11
            if(deltam[i,j] == 0):
                rad_tmp = eps00
            rad[i,j] = rad_tmp


    #Cowan equation 17.34 Note that this is different from Isler 1997 and I am not sure
    #what the resolution is but it seems to be the same. Iser has (-1)^(3J' -m_j' +q)
    #corresponding to 3Jbra - mket. Note this is different from Cowan Jbra - mbra.
    #However, when this gets squared the negative sign does not matter so it is unclear to me at this
    #point if there is case where this difference actually matters.
    #Note that 3Jbra and Jbra are equilvalent with (-1)^(x).
    signal = np.sum(np.sqrt((2*jbra_inten+1)*(2*jket_inten+1))*afcf*wigner_3j(ja,jb,ma,mb,dummus)*\
                     wigner_6j(np.ones_like(jbra_inten),
                               np.ones_like(jbra_inten)*spin[0],
                               np.ones_like(jbra_inten)*orb[0],
                               np.ones_like(jbra_inten)*orb[1],
                               jbra_inten,
                               jket_inten)*\
                    (-1)**(phase1+phase2+phase3 + orb[1] +spin[1] + jbra_inten +2)*rad,axis=0)**2
    w3c0 = wigner_3j(ja,jb,ma,mb,dummus)
    w6c0 = wigner_6j(np.ones_like(jbra_inten),
              np.ones_like(jbra_inten)*spin[0],
              np.ones_like(jbra_inten)*orb[0],
              np.ones_like(jbra_inten)*orb[1],
              jbra_inten,
              jket_inten)
    sigmat = np.sum((np.sqrt((2*jbra_inten+1)*(2*jket_inten+1))*afcf*wigner_3j(ja,jb,ma,mb,dummus)*\
                     wigner_6j(np.ones_like(jbra_inten),
                               np.ones_like(jbra_inten)*spin[0],
                               np.ones_like(jbra_inten)*orb[0],
                               np.ones_like(jbra_inten)*orb[1],
                               jbra_inten,
                               jket_inten)*\
                    (-1)**(phase1+phase2+phase3 + orb[1] +spin[1] + jbra_inten +2)*rad)    )
    #If the relative intentieis of the j levels in the term were specified, weight the result by these
    if(rel_j_intensities.size != 0):

        qt = np.zeros_like(js[0])
        for i in range(0,len(rel_j_intensities)):
            qt[np.where(js[0] == j_ul[0][i])[0]] = rel_j_intensities[i]#create rel inten array from the j array

        qq = np.sum(np.einsum('i,ij->ij',qt,vec1[:,::-1])**2,axis=0)#find components from the eigen vectors

        st = []
        for i in range(0,len(qq)):
            st.append(np.ones(len(jbra[1]))*qq[i])#has to be the same pattern as the kron
                                                  #a1b1 a1b2 a1bn....anb1 anb2 anb3

        weights = np.asarray(st).flatten()#weights in
        signal = signal*weights



    #find mj sublevel purity and mjs labels
    #You cant map a single eigenvalue to a mj sublevel
    #this calculates the purity of the state
    mjs_pur = [] #hold the purity for both upper and lower
    mjs_label = []#hold the label or leading mjs

    #start with the upper
    mj_pur = np.zeros_like(vec1)
    mj_label = np.zeros_like(mjs[0],dtype='int')
    for i in range(0,len(mjs[0])):#loop over the mjs
        rr = np.zeros_like(mjs[0])
        rr[i] = 1 #only look at a particular mj value
        mj_pur[i,:] = np.sum(np.einsum('i,ij->ij',rr,vec1)**2,axis=0)
        #Shows fraction of eigenstate that is made up of mjs sublevel        
        mj_label[i] = np.argmax(mj_pur[i])#put a label on the largest value and call that associated mj
    
    mjs_pur.append(mj_pur)
    mjs_label.append(mj_label)

    #now find the purity of the lower
    mj_pur = np.zeros_like(vec2)
    mj_label = np.zeros_like(mjs[1],dtype='int')
    for i in range(0,len(mjs[1])):
        rr = np.zeros_like(mjs[1])
        rr[i] = 1
        mj_pur[i,:] = np.sum(np.einsum('i,ij->ij',rr,vec2)**2,axis=0)
        mj_label[i] = np.argmax(mj_pur[i])

    mjs_pur.append(mj_pur)
    mjs_label.append(mj_label)

    
    wavevac = wavevac.flatten()[::-1]#flatten and reverse the array so its in the same order as signal
    #convert the vacuum wavelength to air
    air_index = 1.+2.94981e6/(146.e8-(1./(wavevac**2)))+\
                  2.554e4/(41.e8-(1./(wavevac**2)))+6.4328e-5
    wavelengths= wavevac / air_index
    wavelengths = wavelengths * 1E7 #unit change to nm



    # for w3,w6 in zip(w3c0,w6c0):
        
    
    data = {}
    data['wave_air'] = wavelengths #wavelength in air (nm)
    data['wave_vac'] = wavevac*1E7 #wavelength in vacuum (nm)
    data['signal'] = signal
    data['split_energy_up'] = split_energy_up
    data['split_energy_low'] = split_energy_low
    data['eig_vec1'] = vec1 #eigen vectors of the upper 
    data['eig_vec2'] = vec2#eigen vectors of the lower
    data['mjs_purity'] = mjs_pur #purity fraction
    data['mjs_label'] = mjs_label#label in same order as mjs and js
    data['mjs'] = mjs#mj values both upper and lower
    data['js'] = js# j values both upper and lower
    data['afcf'] = afcf
    data['sigmat']=sigmat
    data['rad'] = rad
    data['deltaam'] = deltam
    data['mket_inten'] = mket_inten
    data['mbra_inten'] = mbra_inten
    data['jket_inten'] = jket_inten
    data['jbra_inten'] = jbra_inten
    
    data['ja'] = ja
    data['jb'] = jb
    data['ma'] = ma
    data['mb'] = mb
    data['w3mat'] = w3c0
    data['w6mat'] = w6c0
    # data2['mket_inten'] = mket_inten
    # data2['mbra_inten'] = mbra_inten
    # data2['jket_inten'] = jket_inten
    # data2['jbra_inten'] = jbra_inten

    data['dummus'] = dummus
    return data



def split_ext_fieldInputdeck(Inputdeck):
    s_ground = Inputdeck['s_ground']
    s_excited = Inputdeck['s_excited']
    l_ground = Inputdeck['l_ground']
    l_excited = Inputdeck['l_excited']
    E_ground = Inputdeck['E_ground']
    E_excited = Inputdeck['E_excited']
    Bmag = Inputdeck['Bmag']
    b_angle = Inputdeck['b_angle']
    
    Polangle = 0
    gamma=0
    haspolfilter = False
    if 'Pol_angle' in Inputdeck:
        gamma = Inputdeck['Pol_angle']
        haspolfilter = True
    
    
 
    cbfield = 10000*Bmag
    c_energy = []
    if len(E_excited)==1:
        c_energy.append(E_excited)
    else:
        c_energy.append(E_excited[::-1])    #Need to flip the energy ordering to match Curt's code.

    if len(E_ground)==1:
        c_energy.append(E_ground)
    else:
        c_energy.append(E_ground[::-1])
    # print(cbfield, c_energy)    
    if Bmag==0:
        limitoffdiag = True
    else:
        limitoffdiag = False
    curt_orig = split_ext_field_clean(c_energy,cbfield,np.array([l_excited,l_ground]),np.array([s_excited,s_ground]),b_angle=b_angle, gamma=np.pi*gamma/180, low_field_limit_remove_off_diag=limitoffdiag,Filter=haspolfilter)
    curt_mod = split_ext_field_mod(c_energy,cbfield,np.array([l_excited,l_ground]),np.array([s_excited,s_ground]),b_angle=b_angle, gamma=gamma, low_field_limit_remove_off_diag=limitoffdiag,Filter=haspolfilter)
    
   
    #Zeeman Signal Object is [[wave],[signal]]
    # c_mod = [curt_mod['wave_air'], curt_mod['signal']]
    # c_orig = [curt_orig['wave_air'], curt_orig['signal']]
    
    if 'Convfxn' in Inputdeck:
        Convtype = Inputdeck['Convfxn']
        stepsize = Inputdeck['specstep']
        convtemp = Inputdeck['Temp']
       
        if Inputdeck['Convfxn'] =='Gaussian':
            convtemp = Inputdeck['Temp']
    else:
        #This is the default configuration using the skewed lorrentzian that the Epscore Spectrometer uses.
        Convtype = 'Skewed'
        stepsize = 0.002222
        convtemp = 300
        
    # if 'SpectrumIn' in Inputdeck:
    #     #This top one doesn't use the spectrum input as there might be issues with that and the convolution functions interaction with the padping function
    #     ConvolvedSpect, redsticks = Convol_Spect(data['wave_air'], data['signal'], Inputdeck['plot_window'], Inputdeck['amu'], stepsize,Temperature_in = convtemp, functiontype = Convtype)    

    #     # ConvolvedSpect, redsticks = Convol_Spect(data['wave_air'], data['signal'], [np.min(Inputdeck['SpectrumIn'][0]) , np.max(Inputdeck['SpectrumIn'][0])], Inputdeck['amu'], stepsize,Temperature_in = convtemp, functiontype = Convtype)    
    # else:
    ConvolvedSpect_Orig, redsticks_Orig = Convol_Spect(curt_orig['wave_air'], curt_orig['signal'], Inputdeck['plot_window'], Inputdeck['amu'], stepsize,Temperature_in = convtemp, functiontype = Convtype)    
  
    ConvolvedSpect_Mod, redsticks_Mod = Convol_Spect(curt_mod['wave_air'], curt_mod['signal'], Inputdeck['plot_window'], Inputdeck['amu'], stepsize,Temperature_in = convtemp, functiontype = Convtype)    

    
    return [ConvolvedSpect_Orig,ConvolvedSpect_Mod, [curt_orig['wave_air'] , curt_orig['signal']]] 




#%%Testing
testing = False
if testing:
    
    bfield = 32000 #magnet field [gauss]
    b_angle = 90# angle the LOS makes with the magnetic field
    gamma = 0#allows for polarization of optical setup to be accounted for I haven't tested this at all
    low_field_limit_remove_off_diag = False # If true it will remove the off diagnonal term in the hamiltonian (high field limit)
    
    lvl_he_eng =[]
    #With this energy ordering, my data matches Curts within 1e-6 nm at worst, but it doesn't match the Isler data
    lvl_he_eng.append( np.array([191444.480929, 191444.4821307, 191444.5006512])) # energy of upper levels HIGHEST J FIRST
    lvl_he_eng.append(np.array([169086.7664725,169086.8428979, 169087.8308131])) # energy of lower levels
    
    
    
    #%%
    orb = np.array([2,1])
    spin = np.array([1,1])
    
    hezero = split_ext_field_mod(lvl_he_eng, 0, orb, spin,
                    b_angle = 90, gamma = gamma,
                    low_field_limit_remove_off_diag=True)
    
    
    
    ja_curt = [x[0] for x in  hezero['ja']] #This is the upper js
    jb_curt = [x[0] for x in  hezero['jb']] #This is the lower J's
    ma_curt = [x[0] for x in  hezero['ma']] #This is the upper mjs
    mb_curt = [x[0] for x in  hezero['mb']] #This is the lower mjs
    dumq_curt = [x[0] for x in  hezero['dummus']] #This is dummus
    
    jbra_curt = [x[0] for x in  hezero['jbra_inten']]#THis is upper J's again, They don't exactly match
    jket_curt = [x[0] for x in  hezero['jket_inten']]
    mbra_curt = [x[0] for x in  hezero['mbra_inten']]
    mket_curt = [x[0] for x in  hezero['mket_inten']]

    curt_qs = [ma-mb for ma,mb in zip(mbra_curt,mket_curt)]
    dummus_q = [[dumq_curt[i],ma-mb] for i,(mb,ma) in enumerate(zip(mbra_curt,mket_curt))] 
    
    
    jprimdiff = np.array(jbra_curt) - np.array(ja_curt)# When you add these two together, they match. For some reason adds to one and subtracts from the other
    jjdiff = np.array(jket_curt) - np.array(jb_curt)#
    
    
    mjprimediff =np.array(mbra_curt) - np.array(ma_curt) 
    mjjdiff =np.array(mket_curt) - np.array(mb_curt) 
    

    w3test = [x[0] for x in hezero['w3mat']]
    w6test = [x[0] for x in hezero['w6mat']]
    
    wig3test = []
    wig6test = []
   
    for ja,jb,ma,mb in zip(jbra_curt,jket_curt,mbra_curt,mket_curt):
        jja = ja*np.ones([3,3])
        jjb = jb*np.ones([3,3])
        mja = ma*np.ones([3,3])
        mjb = mb*np.ones([3,3])
        thisq = ma - mb
        dq = thisq*np.ones([3,3])
        tw3 = wigner_3j(jja,jjb,mja,mjb,dq)
        cw6 = wigner_6j(1*np.ones([3,3]),spin[0]*np.ones([3,3]),orb[0]*np.ones([3,3]),orb[1]*np.ones([3,3]),jja,jjb)
        lw6 = float(w6j(spin[0],orb[1],jb,1,ja,orb[0]))        
        wig6test.append([cw6[0][0] , lw6])
  
        Lw3 = 0
        # for i in range(-3,3):
        i=0
        Lw3 = Lw3 +  float(w3j(ja,jb,1,-ma,mb,thisq + i))
        wig3test.append([tw3[0][0] , Lw3])
        
    
    
    wig3test2 = []
    
    for ja,jb,ma,mb,dummusq in zip(ja_curt,jb_curt,ma_curt,mb_curt,dumq_curt):
        jja = ja*np.ones([135,135])
        jjb = jb*np.ones([135,135])
        mja = ma*np.ones([135,135])
        mjb = mb*np.ones([135,135])
        dq = dummusq*np.ones([135,135])
        
        tw3 = wigner_3j(jja,jjb,mja,mjb,dq)
    
        Lw3 = 0
        # for i in range(-1,1):
        i=0
        Lw3 = Lw3 +  float(w3j(ja,jb,1,ma,mb,dummusq + i))
        wig3test2.append([tw3[0][0] , Lw3])
    
    jac = ja_curt[14]
    jbc = jb_curt[14]
    mac = ma_curt[14]
    mbc = mb_curt[14]
    qdum = dumq_curt[14]
    # qdum2 = -(mjprime-mj)
    
    jat =jac*np.ones([15,15])
    jbt = jbc*np.ones([15,15])
    mat = mac*np.ones([15,15])
    mbt = mbc*np.ones([15,15])
    
    dummust = qdum*np.ones([15,15])
    w3test2 = wigner_3j(jat,jbt,mat,mbt,dummust)
    #I found where there is a discrepency in how dummus is calculated that makes the w3 not match, but are there other issues in raw calculation?
    #%% First do integer ones
    j1 = np.arange(0,5.1, 1)
    
    mj12 = np.arange(-5,5.1,1)
    
    # allres = []
    w3bigdiffresults = []
    for jj1 in j1:
        for jj2 in j1:
            for mj1 in mj12:
                for mj2 in mj12:
                    this_q = mj1 - mj2
                    
                    curtw3test = wigner_3j(jj1*np.ones([3,3]) , jj2*np.ones([3,3]), -mj1*np.ones([3,3]), mj2*np.ones([3,3]), this_q*np.ones([3,3])    )
                    leow3test = float(w3j(jj1,jj2,1,-mj1,mj2,this_q))
                    lw3_2 = 0
                    for i in range(-3,3):
                        lw3_2 = lw3_2 + float(w3j(jj1,jj2,1,-mj1,mj2,this_q + i))
                    # allres.append([curtw3test,leow3test])
                    check1 = (np.abs(mj1)<=jj1)
                    check2 = (np.abs(mj2)<=jj2)
                    check3 = np.abs(curtw3test[0][0]) - np.abs(leow3test) >1e-3 
                    allchecks = check1 & check2 & check3 
                    if allchecks:
                        w3bigdiffresults.append([jj1,jj2,mj1,mj2, curtw3test[0][0], leow3test, lw3_2])    
    
    #%% Now do half integer ones
    j1 = np.arange(0.5,5, 1)
    
    mj12 = np.arange(-4.5,5,1)
    
    allres = []
    w3bigdiffresults2 = []
    for jj1 in j1:
        for jj2 in j1:
            for mj1 in mj12:
                for mj2 in mj12:
                    this_q = mj1 - mj2
                    
                    curtw3test = wigner_3j(jj1*np.ones([3,3]) , jj2*np.ones([3,3]), -mj1*np.ones([3,3]), mj2*np.ones([3,3]), this_q*np.ones([3,3])    )
                    leow3test = float(w3j(jj1,jj2,1,-mj1,mj2,this_q))
                    lw3_2 = 0
                    for i in range(-3,3):
                        lw3_2 = lw3_2 + float(w3j(jj1,jj2,1,-mj1,mj2,this_q + i))
                    allres.append([curtw3test,leow3test])
                    check1 = (np.abs(mj1)<=jj1)
                    check2 = (np.abs(mj2)<=jj2)
                    check3 = np.abs(curtw3test[0][0]) - np.abs(leow3test) >1e-3 
                    allchecks = check1 & check2 & check3 
                    if allchecks:
                        w3bigdiffresults2.append([jj1,jj2,mj1,mj2, curtw3test[0][0], leow3test , lw3_2])    



