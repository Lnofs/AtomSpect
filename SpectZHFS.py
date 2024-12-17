# -*- coding: utf-8 -*-
"""
Zeeman_Main
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

def Aconst(Z,gi,n,l):
    if l==0:
        Atemp = (2/3)*mu_0*muN2*muB*g_S*gi*(1/np.pi)*((Z/(a_o*n))**3)/hh
    else:
        Atemp = (2/3)*mu_0*muN2*muB*g_S*gi*(1/np.pi) *(1/(l*(l+0.5)*(l+1)))*(Z/(a_o*n))**3 /hh
    return Atemp
        

def lande_gi(I,I_par,u_I):
    '''
    Generates the lande g_i factor given the spin (I), parity, and nuclear moment.
    All of these factors can be obtained from Stone (2005)
    https://doi.org/10.1016/j.adt.2005.04.001
    
    Parameters
    ----------
    I : Nuclear Spin (Integer or Half Integer)
        Nuclear Spin of the atom. 
    I_par : TYPE
        DESCRIPTION.
    u_I : TYPE
        DESCRIPTION.

    Returns
    -------
    thisgi : g_i
        The lande g_i factor for use in the zeeman shift.

    '''
    #This expression was reverse calculated using Stone and known values for Rubidium, which I no longer am confident in.
    thisgi = u_I*muN/(muP*I*I_par)
    return thisgi

def lande_gj(J,L,S):
    #Lande gJ as derived in Griffiths and others. Comes from projection of angular momentum and averaging of L and S
    if J==0:
        return 0
    else:
        return ( g_L*(J*(J+1) - S*(S+1) + L*(L+1))/(2*J*(J+1)) + g_S*((J*(J+1) + S*(S+1) - L*(L+1))/( 2*J*(J+1))) )

def lande_gf(F,I,J,L,S,g_I):
    #Lande gF, comes from averaging I, L, and S similar to gJ.
    thisgj = lande_gj(J,L,S)
    thisgf = thisgj*(F*(F+1) - I*(I+1) + J*(J+1)) / (2*F*(F+1)) + g_I*(F*(F+1) + I*(I+1) -J*(J+1) ) / (2*F*(F+1))
    return thisgf


def Delta_E_HFS(I,J,F, A):
    #This simple function can only do A constants for the time being. A has units of Hz. Typically A constants are given as h*value [Hz]
    #First need to define a constant K that will be used in the calculations
    #The conversions are to get the output in cm^-1
    K = F*(F+1) - I*(I+1) -J*(J+1)
    E_HFS_Split = 0.5*A*K *hhev/(1.24e-4)    
    return E_HFS_Split



#For conversion factors, 1 = (cm^-1)/(1.2398e-4 eV) = (cm^-1)/(100*c*Hz) = (1.2398e-4 eV)/(100*c*Hz)

def cm_to_nm(E_diff):
    #Converts from energy in cm^-1 to wavelength in nm. 
    return 1e7/E_diff

def cm_to_Hz(E_diff):
    #Converts from energy in cm^-1 to energy in Hz (E = hf = hc/lambda)
    return 100*c_light*E_diff

def cm_to_ev(E_diff):
    #Convert from cm^-1 to eV
    return 1.2398e-4*E_diff


def Vac_to_air(wavein):
    '''
    Converts wavelength from vacuum to air

    Parameters
    ----------
    wavein : TYPE
        Wavelength in (assumed in nm)

    Returns
    -------
    TYPE
        Wavelength_out in nm

    '''
    wavevac = wavein/1E7
    air_index = 1.+2.94981e6/(146e8-(1./(wavevac**2)))+2.554e4/(41e8-(1./(wavevac**2)))+6.4328e-5 
    # n_air = 1 + (0.05792105/(238.0185 - wavevac**(-2))) + (0.00167917/(57.362 - wavevac**(-2) ))
    wavelengths= wavevac / air_index
    return wavelengths*1e7    


def PlotFunction(Plotobject,plot_vars,plottitle = '', plotwind = [],SpectrumPlot = [], Nrows = 1, NormalizeSig = False, NormalizeScale = 1, makefig =True):
    #Plotobjects can [[[waves],[signals]] or a list of multiple plot objects. [[[wave1],[sig1]] ,[[wave2],[sig2]] ] etc
    #plot_vars are [color, linewidth, linestyle, label]
    #verymuch want to add the ability to make multiple plots in the same window using this function and indexing the shape of the multiplot.     
    tPlotObj =[]
    #First check if the plotobject is inhomogogenous, which is to be the case should the user mix LZeeman output objects and other spectra.
    try:
        np.shape(Plotobject) #Shape will fail for a mix of dictionary and arrays due to being inhomogenous. It shouldn't fail for normal [[waves],[signal]]
        if type(Plotobject) ==dict:
            tPlotObj.append(Plotobject['SpecOut'])    
        else:
            tPlotObj = Plotobject
    except:
        for i,x in enumerate(Plotobject): #If it fails, go through all elements as only multiple element inputs should be able to fail
            if type(x) ==dict:
                tPlotObj.append(x['SpecOut']) #If it's a dict, try to remove the SpecOut keyword from LZeeman output object.
            else:
                tPlotObj.append(x) #Otherwise, just append the array.
            #The result here is an array of pairwise datas, but they might be of different lengths.
    try:
        np.shape(tPlotObj)
    except:
        
        #If the np.shape fails a second time, then it means the spectra are just of different lengths. This is fine, so we just iterate through them and plot
        if makefig ==True:
            plt.figure()
            if SpectrumPlot!=[]:
                #Hardcode spectrum to be thickness 1 and black solid line. Only if a new plot is created.
                plt.plot(SpectrumPlot[0], SpectrumPlot[1],color='black',linewidth = 1.5, label='Spectrometer Data' )
        if plotwind ==[]:
            plottingwind = [np.min(tPlotObj[0][0]), np.max(tPlotObj[0][0]) ]
            plotwind = plottingwind
        else:
            plottingwind = plotwind  
            
        for i,(tplotobj,plotvars) in enumerate(zip(tPlotObj,plot_vars)):
            if NormalizeSig == True:
                if len(np.shape([NormalizeScale])) == 1: #Condition checks if the user input a single value for the NormalizeScale or a list. If it's a list, then iterate over the entries to scale as requested.
                    plt.plot(tplotobj[0],Normalize(tplotobj[1],scaling=NormalizeScale),color=plotvars[0],linestyle=plotvars[2],linewidth = plotvars[1], label=plotvars[3], alpha = 0.95)
                else:
                    plt.plot(tplotobj[0],Normalize(tplotobj[1],scaling=NormalizeScale[i]),color=plotvars[0],linestyle=plotvars[2],linewidth = plotvars[1], label=plotvars[3], alpha = 0.95)
                  
            else:
                plt.plot(tplotobj[0],tplotobj[1],color=plotvars[0],linestyle=plotvars[2],linewidth = float(plotvars[1]), label=plotvars[3], alpha = 0.95)
                
        plt.xlim(np.min(plottingwind),np.max(plottingwind))
        # plt.xlim(447.1, 448)
        plt.ylim(-.1,1.1)
        plt.xlabel('Wavelength (nm)',weight='semibold')
        plt.ylabel('Normalized Intensity (arb)',weight='semibold')
        plt.title(plottitle)
        plt.legend(loc='best')
        plt.tight_layout()
        plt.show()
        return 
    
    
    if makefig ==True:
        plt.figure()
        if SpectrumPlot!=[]:
            #Hardcode spectrum to be thickness 1 and black solid line. Only if a new plot is created.
            plt.plot(SpectrumPlot[0], SpectrumPlot[1],color='black',linewidth = 1.5, label='Spectrometer Data' )
    
    #This check allows the user to input multiple spectra objects or not.
    # return tPlotObj, NormalizeScale
    if len(np.shape(tPlotObj)) == 2: #The first case is for a single spectra object or dictionary. A dictionary returns np.shape(dict) = (), which has a length of 0.
        if NormalizeSig == True: #Check if Normalize Sig flag is True, if so then use the Normalize function to scale the results.
            plt.plot(tPlotObj[0],Normalize(tPlotObj[1],scaling=NormalizeScale),color=plot_vars[0],linestyle=plot_vars[2],linewidth = plot_vars[1], label=plot_vars[3], alpha = 0.95)
        else:
            plt.plot(tPlotObj[0],tPlotObj[1],color=plot_vars[0],linestyle=plot_vars[2],linewidth = plot_vars[1], label=plot_vars[3], alpha = 0.95)
        
        if plotwind ==[]:
            plottingwind = [np.min(tPlotObj[0]), np.max(tPlotObj[0]) ]
            plotwind = plottingwind
        else:
            plottingwind = plotwind        
         
    else:
        if plotwind ==[]:
            
            try:            
                
                plottingwind = [np.min(tPlotObj[0][0]), np.max(tPlotObj[0][0]) ]
            except:
                plottingwind = [np.min(tPlotObj[0]['SpecOut'][0]), np.max(tPlotObj[0]['SpecOut'][0]) ]
                
            plotwind = plottingwind
        else:
            plottingwind = plotwind
        
        # plotwind = [min(Plotobject[0][0]), max(Plotobject[0][0]) ]
        for i,(plotobj,plotvars) in enumerate(zip(tPlotObj,plot_vars)):
            # print(plotobj)
            try:
                t_plotobj = plotobj['SpecOut']
            except:
                t_plotobj = plotobj
            if len(np.shape(plot_vars)) ==1:
                plotvars = plot_vars
            if NormalizeSig == True:
                if len(np.shape([NormalizeScale])) == 1:
                    # print("A")
                    # plotvars = plot_vars[0]
                    plt.plot(t_plotobj[0],Normalize(t_plotobj[1],scaling=NormalizeScale),color=plotvars[0],linestyle=plotvars[2],linewidth = plotvars[1], label=plotvars[3], alpha = 0.95)
                else:
                    plt.plot(t_plotobj[0],Normalize(t_plotobj[1],scaling=NormalizeScale[i]),color=plotvars[0],linestyle=plotvars[2],linewidth = plotvars[1], label=plotvars[3], alpha = 0.95)
                    # print("B")
            else:
                if len(np.shape([NormalizeScale])) == 1:
                    # print("A")
                    # plotvars = plot_vars
                    plt.plot(t_plotobj[0],t_plotobj[1],color=plotvars[0],linestyle=plotvars[2],linewidth = plotvars[1], label=plotvars[3], alpha = 0.95)
                else:
                    plt.plot(t_plotobj[0],t_plotobj[1],color=plotvars[0],linestyle=plotvars[2],linewidth = plotvars[1], label=plotvars[3], alpha = 0.95)
                    # print("B")
                # plt.plot(t_plotobj[0],t_plotobj[1],color=plotvars[0],linestyle=plotvars[2],linewidth = float(plotvars[1]), label=plotvars[3], alpha = 0.95)
                # print("C")
    plt.xlim(np.min(plottingwind),np.max(plottingwind))
    # plt.xlim(447.1, 448)
    plt.ylim(-.1,1.1)
    plt.xlabel('Wavelength (nm)',weight='semibold')
    plt.ylabel('Normalized Intensity (arb)',weight='semibold')
    plt.title(plottitle)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.show()
    return


def Normalize(list_in ,scaling=1.0):
    '''
    

    Parameters
    ----------
    list_in : TYPE
        DESCRIPTION.
    scaling : TYPE, optional
        DESCRIPTION. The default is 1.0.

    Returns
    -------
    normlist : TYPE
        Noramlizes the an input function. Uses np.max as that is faster than max

    '''
    if scaling==1.0:
        normlist = [x/np.max(list_in) for x in list_in]
    else:   
        normlist = [scaling*x/np.max(list_in) for x in list_in]
    return normlist


def dipolestr(LG,LE):
    '''
    Wigner3j = (J'    J    1  )
               (-m_j' m_j  q  )   
    Wigner6j = {S    L'  J'  }
               {1    J    L  }                   
    Per Isler paper appendix                
    Parameters
    ----------
    LG, LE are arrays of quantum numbers. Formatted as LG=[J,mJ,L,mL,S,mS]

    Returns
    -------
    tdp3 : float
        Returns the relative dipole strength factor as a multiple of the reduced dipole matrix. Metcalf 4.32

    '''    
    thisq = LE[1] - LG[1]    
    tw3 = float(w3j(LG[0],1,LE[0],LG[1],thisq,-LE[1])) 
    tdp3 = ( (-1)**(LE[2] + LG[4] - LE[1]) ) * np.sqrt((2*LE[0]+1)*(2*LG[0]+1))*tw3*float(w6j(LE[2],LE[0],LG[4],LG[0],LG[2],1))
    return tdp3 



def make_Jmj(J_am):
    #This function takes an input general angular momentum J and creates a list of [J,mJ] for all corresponding mJ's.
    #J_am can be an array itself and the output will be a flat list of all possible J,mJ values, starting with the lowest [J,mj] and
    #cycling through the mJ's for each J.
# =============================================================================
# 
#     'Starting with the lowest J is likely to cause some issues and there should, instead, be a generalized way for the user to input the J-values ordering'
#     'The issue would be because the eigval function sorts in ascending order of energies, so I should do likewise, as that may correct for potential issues'
#     'with energy ordering and states as tracking past diagonalization is not possible.'
# =============================================================================
    ###
    if type(J_am)==int or type(J_am)==float:
        Jtmp = [[float(J_am),float(x)] for x in np.linspace(-J_am,J_am, int(2*J_am+1))]
    else:
        Jtmp = [[float(jvs),float(x)] for jvs in J_am for x in np.linspace(-jvs,jvs, int(2*jvs+1))]        
    return Jtmp

def lo_sk(xrange, peak, pixel_size=0.002222):
    # pixel_size=0.002222 #nm per pixel
    l_w = 0.62*pixel_size
    r_w = 0.99*pixel_size
    sl_func = (((l_w)**2 / (((xrange-peak)**2 + (l_w)**2))) *np.heaviside(peak-xrange,1) + ((r_w)**2 / (((xrange-peak)**2 + (r_w)**2))) *np.heaviside(xrange-peak,1))
    return sl_func



def polarization(mupper, mlower, Btheta,gamma=0,Pol_Filter=False,T_Alpha = 1, T_Beta=0):
    '''
    This polarization calculation is generally based on Isler1997 http://dx.doi.org/10.1063/1.872095
    A correction is made with respect to a negative sign in the definition of epsilon plus/minus to be more in line with typical notation.
    This correction removes the issue of an extra negative sign for epsilon_+1-1 when calculated out and the author having the magnitude only.

    Parameters
    ----------
    mupper : magnetic quantum number of upper level. mJ'
    mlower : magnetic quantum number of lower level, mJ
    Btheta : Angle between LoS and maximum magnetic field
    gamma : Angle polarization filter makes with respect to maximizing linear (pi) light transmission.
    Pol_Filter : Boolean, whether a polarization filter was in place for data collection
    T_Alpha : Transmission along one polarization axis. 
    T_Beta : Transmission coefficient along perpendicular polarization axis.

    Returns
    -------
    rad_tmp : Projection of polarization signal amplitude for emmitted light.

    '''

    
    if Pol_Filter==False:
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
    theta_rad = np.deg2rad(Btheta)#convert deg to radians
   

    
    epsz = ealph*np.sin(theta_rad)
    
    e0 = epsz
    
    eps00 = np.sqrt(e0*e0)
    # eps11 = np.sqrt(0.5*np.abs(np.real(e1p*e1m)))
    #THe 0.5 is pulled outside of the sqrt since our function only checks for delta_AM, and we would otherwise doublecount.
    
    eps11 = 0.5*np.sqrt((np.cos(theta_rad)**2) * (ealph**2) + (ebeta**2 ))
    # print(eps00, eps11)
    delta_am =  np.abs(np.array(mupper) - np.array(mlower))
    rad_tmp = 0
    if(delta_am ==1):
        rad_tmp = eps11
        
    elif(delta_am == 0):
        rad_tmp = eps00
    else:
        rad_tmp = 0
    return rad_tmp



from sympy.physics.quantum.cg import CG

def CGCoeff(j1,m1,j2,m2,j3,m3):
    #Calculates Clebsch-Gordan Coefficient and outputs a float instead of sympy result
    thisCG = float(CG(j1, m1, j2, m2, j3, m3).doit())
    return thisCG


def Zeemanlow(gf, mf, B):
    #Lowfield Zeeman approximation. Valid for J or F total AM.
    return (muBcm*gf*mf*B)

def stateprocess(L_state,S_state,E_level,thisB):
    '''
    Takes input state values and magnetic field to expand the state space. 
    L and S are used to calculate possible J values. all possible mJ values are also created. 
    E_level is an np.array and must match the length of possible J values. Since the values are pulled from NIST, this will likely happen.
    

    Parameters
    ----------
    L_state : Int or half int
        Orbital angular quantum number. S=0, P=1, D=2, F=3, etc
    S_state : Int or half int
        Spin quantum number, n = 2s+1 for multiplicity, so ^2S_J has (2*s +1) =2, so s=0.5, L=0
    E_level : np.array
        Energy in cm^-1, Lowest J value first.
    thisB : 
        Magnetic field in T

    Returns
    -------
    list
       [w3mat, Jmjs, LS,E_low, E_cm,E_SO_cm_Mat,sortedstates]
    w3mat is the state array of w3 coefficients corresponding to CG coefficients.
    Jmjs are the possible Jmjs #Deprecated
    LS is the L,mL,S,mS states. #Deprecated
    E_low is a list of energies corresponding to the lowfield zeeman approx
    E_cm is the sorted energies in cm^-1
    E_SO_cm_Mat is the sorted (diagonal) matrix of energies with repeated values for degeneracy of mJ
    sortedstates is [J,mJ,L,mL,S,mS] states corresponding to the sorted energies in the matrix. Off diagonals are excluded to this point.
    
    '''
    
    #Create an array that contains [L,ml,S,ms] for the given L,S values. Starts with lowest ml,mS values and itereates through ms first.
    LS0 = [[float(L_state), float(x), float(S_state), float(y)] for x in np.linspace(L_state, -L_state, int(2*L_state+1),dtype=float) for y in np.linspace(S_state, -S_state, int(2*S_state+1),dtype=float) ]
    # LS0 = [[L_state, x, S_state, y] for x in np.linspace(L_state, -L_state, int(2*L_state+1),dtype=float) for y in np.linspace(S_state, -S_state, int(2*S_state+1),dtype=float) ]
    
    #Now create the possible J values
    Jvals = np.arange(np.abs(L_state-S_state), (L_state+S_state+.1))#,dtype=float)

    #Use those J values to make all possible mJ values
    Jmjs0 = make_Jmj(Jvals)
    
    #Combine both lists together to be a single large [J,mJ,L,mL,S,mS] list. The dimensions should match for LS0 and Jmjs0.
    bigstatelist = [[x[0],x[1],y[0],y[1],y[2],y[3]] for x,y in zip(Jmjs0,LS0)] 
    
    #Expand the input energies which were given for each J to correspond to the number of mJs
    E_cm0 = [[E_level[x]]*int(2*Jvals[x] +1) for x in range(len(E_level))]
    #Flatten the results out
    E_cm_flat = [item for row in E_cm0 for item in row]
    
    #This step uses a low field approximation in order to sort the state-list based on the given magnetic field such that the relative state function order is known.
    #This is used for signal strength calculation as otherwise we may not know what order the functions are after they get diagonalized. This step might not be needed.
    
    biglist = []
    for i,x in enumerate(E_cm_flat):
        thislandegJ = lande_gj(Jmjs0[i][0],L_state,S_state) #caluclate lande gj factor for a given state
        if thisB == 0:
            thisZeemanEnergy =  x #This is to prevent issues with attempting to sort for B=0
        else:
            thisZeemanEnergy = Zeemanlow(thislandegJ,Jmjs0[i][1],thisB) + x # This is the low zeeman approximation total energy
        biglist.append([thisZeemanEnergy,x,Jmjs0[i],LS0[i],bigstatelist[i]])
    
    sortedlist = sorted(biglist)
    E_low = [x[0] for x in sortedlist]
    sortedstates = [x[-1] for x in sortedlist]
    E_cm = [x[1] for x in sortedlist] #This is the sorted energies, and is in cm^-1
    Jmjs = [x[2] for x in sortedlist]
    # print(Jmjs)
    LS = [x[3] for x in sortedlist]
    # print(LS)
    E_SO_cm_Mat = np.diagflat([item for item in E_cm]) #Despite the name, these energies are actually in cm-1
    
    
    #We now use wigner 3j to get an equivalent decoupling as per how CG coefficients work.
    #The result is a 2d where the columns are the |J,mJ> states and the rows determine what the corresponding |L,mL,S,mS> state was. 
    #Eventually would like a way to convert this data to a more human readable form 
    dim = len(Jmjs)
    w3mat = np.zeros((dim,dim))
    
    for i in range(dim):
        for j in range(dim):
            #This uses the Metcalf 4.30 for relation between CG and w3j
            w3mat[j,i] = ((-1)**(-LS[i][0] + LS[i][2] - Jmjs[j][1]))*np.sqrt(2*Jmjs[j][0] +1)*float(w3j(LS[i][0],LS[i][2], Jmjs[j][0], LS[i][1], LS[i][3], -Jmjs[j][1] ))
            # w3mat[j,i] = ((-1)**(-LS[i][0] + LS[i][2] - Jmjs[j][1]))*np.sqrt(2*Jmjs[j][0] +1)*w3j(float(LS[i][0]),float(LS[i][2]), float(Jmjs[j][0]), float(LS[i][1]), float(LS[i][3]), float(-Jmjs[j][1] ))

            # w3mat[j,i] = CGw3(LS[i][0],LS[i][2],Jmjs[j][0],LS[i][1],LS[i][3],-Jmjs[j][1])
 
    return [w3mat, Jmjs, LS,E_low, E_cm,E_SO_cm_Mat,sortedstates]


#%%


def ZeemanhighHFS(state_in,Ahfs,g_I, Bmag):
    #This is the high field approximation for zeeman splitting (with HFS). E_high = A*m_i*m_j + uB(gjmj + gimi)B 
    #FOr high field, can't use J, and must use s and L. Keep limited to diagonal though.
    thisgj = lande_gj(state_in[2],state_in[6],state_in[8])
    EPB_temp = Ahfs*state_in[3]*state_in[5] + mu_B*(thisgj*state_in[3] + g_I*state_in[5])*Bmag
    return EPB_temp

def CGw3(j1,m1,j2,m2,j3,m3):
    #Clebsh Gordon Coeff using wigner3j
        return ((-1)**(-j1 + j2 - m3))*np.sqrt(2*j3 +1)*float(w3j(float(j1),float(j2),float(j3),float(m1),float(m2),-float(m3)))


def stateprocess_HFS(L_state,S_state,I_state, HFS_Const,g_I,E_level,thisB):
        
    #Create an array that contains [L,ml,S,ms] for the given L,S values. Starts with lowest ml,mS values and itereates through ms first.
    # LS0 = [[L_state, x, S_state, y] for x in np.linspace(-L_state, L_state, int(2*L_state+1)) for y in np.linspace(-S_state, S_state, int(2*S_state+1)) ]
    #Create an array that contains [L,ml,S,ms,I,mI] for the given L,S,I values. Starts with Highest mL,mS values and itereates through mI first, then mS, then mL
    #Extra type-casting is required due to new numpy.float64 objects interaction with sympy wigner functions
    LSI0 = [[float(L_state), float(x), float(S_state), float(y), float(I_state), float(z)] for x in np.linspace(-L_state, L_state, int(2*L_state+1)) for y in np.linspace(-S_state, S_state, int(2*S_state+1))for z in np.linspace(I_state, -I_state, int(2*I_state+1)) ]
    
    #Now create the V values, starting with the lowest J state first.
    
    #Make all Js corresponding to L,S relation
    Jvals = np.arange(np.abs(L_state-S_state), (L_state+S_state+.1))
    
    #Make all mJs for the Js
    Jmjs0 = make_Jmj(Jvals)
    
    #Make a list of J,mJ,I,mI states
    JI0 = [[float(x[0]),float(x[1]),float(I_state), float(y)] for x in Jmjs0 for y in np.linspace(-I_state, I_state, int(2*I_state+1)) ]
    
    #Combine J and Is to make possible Fs
    Fvals = [np.arange(np.abs(-x + I_state), (x+ I_state+.1)) for x in Jvals]
    #Flatten the list to be one big list
    Fvals = [float(item) for row in Fvals for item in row]
    #Make the possible mFs for all Fs
    Fmfs0 = make_Jmj(Fvals)



    #Duplicate energies based on degeneracy and flatten the results
    E_cm0 = [[E_level[x]]*int(2*Jvals[x] +1)*int(2*I_state + 1) for x in range(len(E_level))]
    E0_cm_flat = [item for row in E_cm0 for item in row]
    #Duplicate HFS constants to match F levels based on degeneracy and flatten the results
    HFSConst = [[HFS_Const[x]]*int(2*Jvals[x] +1)*int(2*I_state + 1) for x in range(len(E_level))]
    HFS_Flat = [item for row in HFSConst for item in row]
    #Flatten the shape of the Energy array

    #Make a large list of [F,mF,J,mJ,I,mI,L,mL,S,mS]
    bigstatelist = [[x[0],x[1],y[0],y[1],y[2],y[3],z[0],z[1],z[2],z[3]] for x,y,z, in zip(Fmfs0,JI0,LSI0)]
    
    #This step uses a low field approximation in order to sort the state-list based on the given magnetic field such that the relative state function order is known.
    #This is used for signal strength calculation as otherwise we may not know what order the functions are after they get diagonalized. This step might not be needed.
        
    biglist = []
    for i,x in enumerate(E0_cm_flat):
        thislandegJ = lande_gj(bigstatelist[i][2],L_state,S_state) #Calculate lande gj factor for the given transition
        if thisB == 0:
            thisZeemanEnergy =  x + Delta_E_HFS(I_state,bigstatelist[i][2],bigstatelist[i][0],HFS_Flat[i]) #Only add hfs splitting if B=0
        else: #Lowfield zeeman approximation with hfs, which should be mf and gf.
            thisZeemanEnergy = Zeemanlow(thislandegJ,bigstatelist[i][3],thisB) + muBcm*g_I*thisB*bigstatelist[i][5] + x + Delta_E_HFS(I_state,bigstatelist[i][2],bigstatelist[i][0],HFS_Flat[i])
        biglist.append([thisZeemanEnergy,x,bigstatelist[i] , HFS_Flat[i]])
    
    #Sort the results and pull sorted values
    sortedlist = sorted(biglist)
    sortedstates = [x[2] for x in sortedlist]
    sortedHFC = [x[3] for x in sortedlist] #Sorted hyperfine constants
    E_cm = [x[1] for x in sortedlist]
    # E_approx = [x[0] for x in sortedlist]
    
    
    #Calculate Hyperfine energies with sorted states and values.
    E_HFSFlat = [E_cm[i] + Delta_E_HFS(states1[4],states1[2],states1[0], sortedHFC[i])    for i,states1 in enumerate(sortedstates)] #Hyperfine + E0 (spin orbital, Nist)
    E_HFS0 = [Delta_E_HFS(states1[4],states1[2],states1[0], sortedHFC[i])    for i,states1 in enumerate(sortedstates)] #Just Hyperfine energy

    E0_cm_Mat = np.diagflat([item for item in E_HFSFlat]) #Make the diagonal energy matrix
 
  
 
    #We now use wigner 3j to get an equivalent decoupling as per how CG coefficients work.
    #There are two coefficients as we go from |F,mF> -> |J,mJ>|I,mI> then another set for |J,mJ>|I,mI> -> |L,mL>|S,mS>|I,mI>. Keeping track of constants is slightly tricky.
    #There can be issues with overla
    #Eventually would like a way to convert this data to a more human readable form 
    dimhfs = len(Fmfs0)
    w3cg1 = np.zeros((dimhfs,dimhfs))
    w3cg2 = np.zeros((dimhfs,dimhfs))
    #Done in a single set of iterations. The check for states1[2] != states2[2] is saying that J != J' to calculate, otherwise it's zero.
    for i,states1 in enumerate(sortedstates):
        
        for j,states2 in enumerate(sortedstates):
            #This conditional is saying that J != J', as otherwise the different J levels will cause duplicate CG's.
            if states1[2]!= states2[2]:
                CG1 = 0
            else:
                #The CG coefficient couples F,mF to J,mJ;I,mI for CG1
                CG1 = CGw3(states2[2],states2[3],states2[4],states2[5], states1[0],states1[1])
            #This conditional ensures we don't allow nuclear spin to flip or change between states. states[5] is mI.   
            if states1[5]!=states2[5]:
                CG2 = 0
            else:
                CG2 = CGw3(states2[6],states2[7],states2[8],states2[9],states1[2],states1[3])

            
                w3cg2[i,j] = CG2
            w3cg1[i,j] = CG1
    
    w3cgprod = np.matmul(w3cg1,w3cg2)        

    return [w3cgprod,E0_cm_Mat,E_HFS0,sortedstates]






def Zeeman_signal(LevelG,LZG,LevelE,LZE,Bangle=90,gamma = 0,Filter=False):
    
    #Need to make a list of all possible combinations such that each upper state is repeated by the dim of the lower level, and a second list where the opposite is true.
    #ie, each lower level is repeated by the count of upper level states. This makes two lower x upper length arrays that will be useful for wigner calcs.
    bjlist = [[x,y] for y in LevelE[-1] for x in LevelG[-1]] #This list will be [[bigstatesG], [bigstatesE]]
    
    
    #need to define the polarization effect from the filter.
    rad = np.zeros((len(bjlist)))
    for i in range(len(bjlist)):
        rad[i] = polarization(bjlist[i][1][1],bjlist[i][0][1],Bangle,gamma,Pol_Filter=Filter)
   
    #Calculate the signal strength coefficients as the product of radiation polarization term and dipole strength.
    #This is essentially a row vector.
    Signalcoeff = []
    for i in range(len(bjlist)):
        Signalcoeff.append(rad[i]*dipolestr(bjlist[i][0],bjlist[i][1]))
    Signalcoeffa = np.array(Signalcoeff)
    
    afcf = np.kron(np.array(LZE[1]),np.array(LZG[1])) #This uses the tensor product after the magnetic field is taken into account.

    #Calculate the signal which is the sum of the square of the product of afcf, which si reduced dipole signal (including field effects) and the signal coeff. Isler1997
    signal_out = []
    for i in range(len(bjlist)):
        signal_out.append((np.sum(afcf[:,i]*Signalcoeffa)**2))
    sigout = np.array(signal_out)
        
    #Conver the output to wavelength in nm then account for shift from vacuum to air.
    wave_vac = np.array( [1e7/(y-x) for y in LZE[0] for x in LZG[0]])
    wave_air = Vac_to_air(wave_vac)

    return [sigout, wave_air,rad,bjlist,afcf,Signalcoeffa,wave_vac]


def StarkShift(LevelG,LevelE,Efield=1):
    #This is a preiminarly attempt at including stark shift. Currently issues exist involving both the summation component and the components of the dipole matrix.
    bjlist = [[x,y] for y in LevelE[-1] for x in LevelG[-1]]
    bigElist = np.array( [[x,y] for y in LevelE[-3] for x in LevelG[-3]])
    bigElist2 =np.append(LevelG[-3],LevelE[-3])
    # print(len(LevelG[-3]))
    # print(len(LevelE[-3]))
    dipole_str = []
    for i in range(len(bjlist)):
        # thisq = jlistbig[i][1][1] - jlistbig[i][0][1]
        dipole_str.append(dipolestr(bjlist[i][0],bjlist[i][1]))

    #I think we need to include the reduced dipole matrix still, which would be the kroneker product of the original wavefunctions, not the re-diagonalized.    
    #For stark shift, which is second order perturbation, we need to sum over the dipolestr calculation. See writeup for the relation details, but 
    #Delta_E = -E^2 sum(((dipolestr)**2)/(Ek - En)) for k!=n 
    #For that purpose, we shall construct two arrays where the elements are the results of these summations then combine them to obtain the quadratic stark shift.
    #The results will be a diagonal matrix and would be a similar procedure to how the quadratic zeeman shift would be calculated.    
    tstarkscaleG = []
    for i in range(len(LevelG[-3])):

        tstark = 0
        for j in range(len(bjlist)):
            # if i!=j:
                if bigElist2[i]==bigElist[j][1]:
                    tstark = tstark + 0
                else:
                    tstark = tstark + ((dipole_str[j])**2)*(1/(LevelG[-3][i] - bigElist[j][1]))
        tstarkscaleG.append(-Efield*Efield*tstark)           
            
    tstarkscaleE = []
    
    for i in range(len(LevelE[-3])):

        tstark = 0
        for j in range(len(bjlist)):
            # if i!=j:
                if bigElist2[i]==bigElist[j][1]:
                    tstark = tstark + 0
                else:
                    tstark = tstark + ((dipole_str[j])**2)*(1/(LevelE[-3][i]  - bigElist[j][0]))
        tstarkscaleE.append(-Efield*Efield*tstark)       
    StarkG = np.diag(tstarkscaleG)
    StarkE = np.diag(tstarkscaleE)
    # print(StarkE)
    # return bigElist2
    return [StarkG,StarkE]

   

def dipolestr_HFS(LG,LE):
    '''
    Wigner3j =  (F     1    F'  )
                (-m_f  q   -m_F )   
    Wigner6j1 = {L'    J'  S  }
                {J    L    1  } 
    Wigner6j2 = {J'   F'   I  }
                {F    J    1  }                   
    Per Metcalf 4.33 
    
    Input state objects are [F,mF,J,mJ,I,mI,L,mL,S,mS]               
    Parameters
    ----------

    Returns
    -------
    tdp3 : TYPE
        DESCRIPTION.

    '''
    #States are from bigsortedlist = [F,mF,J,mJ,I,mI,L,mL,S,mS]
                                    #[0,1 ,2,3 ,4,5 ,6,7 ,8,9 ]
                                    
    #This first block is from dan steck and is a slightly different formulation for the relative strength, but seems equivalent.
    thisq = LE[1] - LG[1]
   
    
    #This is Metcalf 4.33
    #Calculate the wigner 3j and two wigner 6j. 
    tw3 =float(w3j(LG[0],1,LE[0],LG[1],thisq,-LE[1]))
    tw61 = float(w6j(LE[6],LE[2],LG[8],LG[2],LG[6],1))
    tw62 = float(w6j(LE[2],LE[0],LG[4],LG[0],LG[2],1))
    
    tphase = (-1)**(1 + LE[6] + LG[8] + LG[2] + LE[2] + LG[4] - LE[1])
    troot = np.sqrt((2*LG[0] +1)*(2*LE[0]+1) * (2*LG[2] +1)*(2*LE[2] +1))

    tdp3 = tphase* troot*tw3*tw61*tw62

    return tdp3    
   
    

def Zeeman_signal_HFS(LevelG,LZG,LevelE,LZE,Bangle=90,gamma = 0,Filter=False):

    #Need to make a list of all possible combinations such that each upper state is repeated by the dim of the lower level, and a second list where the opposite is true.
    #ie, each lower level is repeated by the count of upper level states. This makes two lower x upper length arrays that will be useful for wigner calcs.
    bflist = [[x,y] for y in LevelE[-1] for x in LevelG[-1]] #This list will be [[bigstatesG], [bigstatesE]]
    
    
    #need to define the polarization effect from the filter.
    rad = np.zeros((len(bflist)))
    for i in range(len(bflist)):
        rad[i] = polarization(bflist[i][1][1],bflist[i][0][1],Bangle,gamma,Pol_Filter=Filter)

    #Calculate signal strength coefficient as product of dipole strength and polarization component
    Signalcoeff = []
    for i in range(len(bflist)):
        Signalcoeff.append(rad[i]*dipolestr_HFS(bflist[i][0],bflist[i][1]))
    Signalcoeffa = np.array(Signalcoeff)
    
    
    afcf = np.kron(np.array(LZE[1]),np.array(LZG[1]))#This uses the tensor product after the magnetic field is taken into account.
    
    #Calculate the signal. This is like the dot product of a column of the afcf reduced dipole state matrix with the row vector comprised of dipole strengths(with polarization)
    signal_out = []
    for i in range(len(bflist)):
        signal_out.append((np.sum(afcf[:,i]*Signalcoeffa)**2))
    sigout = np.array(signal_out)
    
    
    #Convert the result from cm^-1 to nm and include shift to air from vacuum
    wave_vac = np.array( [1e7/(y-x) for y in LZE[0] for x in LZG[0]])
    wave_air = Vac_to_air(wave_vac)

    return [sigout, wave_air,rad,bflist,afcf,Signalcoeffa,wave_vac]



def Zeeman_func_HFS(Level_in,Binput,gI):
    '''
    Calcualte the Zeeman shift including hyperfine effects.

    Parameters
    ----------
    Level_in : TYPE
        DESCRIPTION.
    Binput : TYPE
        DESCRIPTION.
    gI : TYPE
        DESCRIPTION.

    Returns
    -------
    list
        DESCRIPTION.

    '''
    Ldim=len(Level_in[-1]) #Dimension of the atomic level: Length of the # of possible states.
    Psi=Level_in[0]                
            
    PsiL = np.zeros_like(Psi)
    PsiS = np.zeros_like(Psi)
    PsiI = np.zeros_like(Psi)
    
   
    #Create diagonal elements for the m values for each qm.     
    tmLMat = np.diagflat([x[7] for x in Level_in[-1]])
    tmSMat = np.diagflat([x[9] for x in Level_in[-1]])
    tmIMat = np.diagflat([x[5] for x in Level_in[-1]])
    
    
    #This is equivalent to the <n|Jz|m> calculation to make the pertubation matrix of zeeman energies    
    #The += doesn't seem to change anything, but is kept as my original code used it

    for i in range(Ldim):
        tPsiL = Psi[i]            
        for j in range(Ldim):
            PsiL[i,j] = PsiL[i,j] + np.linalg.multi_dot([tPsiL, tmLMat, Psi[j]])
            PsiS[i,j] = PsiS[i,j] + np.linalg.multi_dot([tPsiL, tmSMat, Psi[j]])  
            PsiI[i,j] = PsiI[i,j] + np.linalg.multi_dot([tPsiL, tmIMat, Psi[j]])            
    H_Zeeman = g_S*PsiS + g_L*PsiL + gI*PsiI #Combine all components
    
    
    #Combine the energy matrices
    scaledMat = np.asmatrix(Level_in[1] + Binput*H_Zeeman*(muBcm) )
    #Diagonalize and calculate eigenvalues and eigenvectors.
    eigZ,eigvecZ = np.linalg.eigh(scaledMat)#, driver = "evx")
    eigZ2 = eigZ

    return [eigZ2,eigvecZ]

    


#%%
def Zeeman_func(Level_in,Binput):
    # Level=[CGmat, Jmjs, LS, E_SO_Hz_Mat]
    Ldim=len(Level_in[-1]) #Dimension of the atomic level: Length of the # of possible states.
    #Spi is wavefunction basis matrix. COmprised of CG coefficients.
    Psi=Level_in[0]                
            
    PsiL = np.zeros_like(Psi)
    PsiS = np.zeros_like(Psi)
    
    #Create diagonal elements for the m values for each qm.     
    tmLMat = np.diagflat([Level_in[2][x][1] for x in range(Ldim)])
    tmSMat = np.diagflat([Level_in[2][x][3] for x in range(Ldim)])
    
    #This is equivalent to the <n|Jz|m> calculation to make the pertubation matrix of zeeman energies    
    #The += doesn't seem to change anything, but is kept as my original code used it
    for i in range(Ldim):
        tPsiL = Psi[i]            
        for j in range(Ldim):
            PsiL[i,j] = PsiL[i,j] + np.linalg.multi_dot([tPsiL, tmLMat, Psi[j]])
            PsiS[i,j] = PsiS[i,j] + np.linalg.multi_dot([tPsiL, tmSMat, Psi[j]])  
    #Combine the energy matricies for mL and mS
    H_Zeeman = g_S*PsiS + g_L*PsiL
    

    #Diagonalize and calculate the eigenvalues and eigenvectors
    scaledMat = np.asmatrix(Level_in[5] + Binput*H_Zeeman*(muBcm ))
    eigZ,eigvecZ = scp.linalg.eigh(scaledMat)#, driver = "evx") #Using driver="evx" seems to change some of the signs of the resulting matricies. In the 0-field case, the eigenvectors should all be +1 I believe.


    return [eigZ,eigvecZ]

def ZeemanStark_func(Level_in,StarkIn,Binput):
    #This function adds the stark shift energies into the matrix for diagonalization. Currently the stark shift matrix calculation function is incorrect.
    # Level=[CGmat, Jmjs, LS, E_SO_Hz_Mat]
    Ldim=len(Level_in[-1]) #Dimension of the atomic level: Length of the # of possible states.
    #Psi is wavefunction basis matrix. COmprised of CG coefficients.
    Psi=Level_in[0]                
            
    PsiL = np.zeros_like(Psi)
    PsiS = np.zeros_like(Psi)
    #Create diagonal elements for the m values for each qm.      
    tmLMat = np.diagflat([Level_in[2][x][1] for x in range(Ldim)])
    tmSMat = np.diagflat([Level_in[2][x][3] for x in range(Ldim)])
    
    #This is equivalent to the <n|Jz|m> calculation to make the pertubation matrix of zeeman energies    
    #The += doesn't seem to change anything, but is kept as my original code used it
    for i in range(Ldim):
        tPsiL = Psi[i]            
        for j in range(Ldim):
            PsiL[i,j] = PsiL[i,j] + np.linalg.multi_dot([tPsiL, tmLMat, Psi[j]])
            PsiS[i,j] = PsiS[i,j] + np.linalg.multi_dot([tPsiL, tmSMat, Psi[j]])  
    H_Zeeman = g_S*PsiS + g_L*PsiL
    

        # EZtemp = Binput[0]*H_Zeeman*(muB/hh)
        # Etot_temp = np.linalg.eigvalsh(Level[3] + Binput[0]*H_Zeeman*(muB/hh)  )
        #Uncomment in case of errors.
    scaledMat = np.asmatrix(Level_in[5] + StarkIn + Binput*H_Zeeman*(muBcm ))
    eigZ,eigvecZ = scp.linalg.eigh(scaledMat)#, driver = "evx")


    return [eigZ,eigvecZ]



    #%%
def ZeemanFan(Inputdeck):
    #Take the required values from the input deck.
    S_ground = Inputdeck['s_ground']
    S_excited = Inputdeck['s_excited']
    L_ground = Inputdeck['l_ground']
    L_excited = Inputdeck['l_excited']
    E_ground = Inputdeck['E_ground']
    E_excited = Inputdeck['E_excited']
    Bfield_ext = Inputdeck['Bmag']
    Bext_angle = Inputdeck['b_angle']
    
    #Generate empty lists for storing values in later. Ground, excited, low and exact. Look to add high approximation too.
    GZ_fan = []
    EZ_fan = []
    GZ_fan_Low = []
    EZ_fan_Low = []
     
    #Check if HFS is being included.
    if 'I_spin' in Inputdeck:
        G_fanE_High = np.zeros_like(Bfield_ext)
        E_fanE_High = np.zeros_like(Bfield_ext)         
        for Bvals in Bfield_ext:
            
            HFSConstsG = Inputdeck['HFS_G']
            HFSConstsE = Inputdeck['HFS_E']
            #Calculate the lande g_i factor from the nuclear factor
            I_spin = np.abs(Inputdeck['I_spin'])
            I_parity = np.sign(Inputdeck['I_spin'])
            mu_I = Inputdeck['mu_I']
            thisg_I = lande_gi(I_spin, I_parity, mu_I)
            
            #State process to generate the complte state list for the ground and excited levels. This includes degeneracy.
            Ground_Level_HFS = stateprocess_HFS(L_ground, S_ground, I_spin, HFSConstsG, thisg_I, E_ground, Bvals)
            Excited_Level_HFS = stateprocess_HFS(L_excited, S_excited, I_spin, HFSConstsE, thisg_I, E_excited, Bvals)
            
            #Take processed state objects and calculate zeeman splitting (with HFS)
            Ground_Zeeman_HFS= Zeeman_func_HFS(Ground_Level_HFS, Bvals, thisg_I)
            Excited_Zeeman_HFS = Zeeman_func_HFS(Excited_Level_HFS, Bvals, thisg_I)    
            #For each magnetic field, append the eigenvalues for the new energies.
            GZ_fan.append( Ground_Zeeman_HFS[0] )                               
            EZ_fan.append (Excited_Zeeman_HFS[0])
    else:
        
        for Bvals in Bfield_ext:
            #State process to generate the complte state list for the ground and excited levels. This includes degeneracy.
            Ground_Level = stateprocess(L_ground,S_ground,E_ground,Bvals)
            Excited_Level = stateprocess(L_excited, S_excited, E_excited, Bvals)
            
            #Take processed state objects and calculate zeeman splitting (without HFS)
            Ground_Zeeman= Zeeman_func(Ground_Level, Bvals)
            Excited_Zeeman = Zeeman_func(Excited_Level, Bvals) 
            #Append the eigenvalues for energy for each magnetic field value
            GZ_fan.append( Ground_Zeeman[0] )                               
            EZ_fan.append (Excited_Zeeman[0])
            #Do the same but using lowfield approximation, which is calculated in the stateprocess function
            GZ_fan_Low.append( Ground_Level[3])
            EZ_fan_Low.append( Excited_Level[3])

    return [GZ_fan,EZ_fan,GZ_fan_Low,EZ_fan_Low]



def plotZfan(Inputdeck,dolow=True,savefig = False, figname = 'placeholder'):
    '''
    Takes the input deck (which should have Bmag as a range of values), and calculates the typical zeeman fan.
    Both exact solution and lowfield are plotted on the same figure if dolow=True.

    Parameters
    ----------
    Inputdeck : dict
        DESCRIPTION.
    dolow : Boolean, optional
        DESCRIPTION. The default is True.
    savefig : Boolean, optional
        DESCRIPTION. The default is False.
    figname : String, optional
        DESCRIPTION. The default is 'placeholder'.

    Returns
    -------
    Returns exact and lowfield zeeman fan incase different plotting is required.
    [Gfan,Efan,Gfan_low,Efan_low]

    '''
    #Takes the inputdeck and a range of magnetic fields to calculate a typical zeeman fan.
    #Exact and lowfield are plotted on the same figure if that option is requested.
    
    Zfantest = ZeemanFan(Inputdeck) #Calculate
    #Pull energy arrays from Zfancalculation
    Gfan = np.array(Zfantest[0]).T
    Efan =np.array(Zfantest[1]).T
    Gfan_low = np.array(Zfantest[2]).T
    Efan_low = np.array(Zfantest[3]).T
    Brange=Inputdeck['Bmag']

    #Make a new plot for the ground state
    plt.figure()
    for i,things in enumerate(Gfan):
        plt.plot(Brange,things, color ='blue')#, label = "Ground")
        
    plt.plot(Brange[0],Gfan[0][0], color ='blue', label = "Ground") #This makes a single point for the purpose of the legend.

    if dolow==True:  
        for i,things in enumerate(Gfan_low):    
            plt.plot(Brange,things, color ='maroon')#, label = "Ground_Low")
        plt.plot(Brange[0],Gfan_low[0][0], color ='maroon', label = "Ground_Low") #This makes a single point for the purpose of the legend.
   
    #Setting plot parameters
    plt.xlabel('Magnetic Field (T)',weight='semibold')
    plt.ylabel('Energy (cm^-1)',weight='semibold')
    plt.legend(loc='best')
    plt.title(Inputdeck['plottitle'])
    plt.tight_layout() 
    if savefig==True:
        plt.savefig(f'{figname}_Ground.png',dpi=800) #Save the figure if requested

     #Make a new plot for the excited state
    plt.figure()
    for i,things in enumerate(Efan):
        plt.plot(Brange,things, color ='red')#, label = "Excited")
    plt.plot(Brange[0],Efan[0][0], color ='red', label = "Excited")#This makes a single point for the purpose of the legend.
    

    if dolow==True:  
        for i,things in enumerate(Efan_low):    
            plt.plot(Brange,things, color ='orange')#, label = "Excited_Low")
        plt.plot(Brange[0],Efan_low[0][0], color ='orange', label = "Excited_Low")#This makes a single point for the purpose of the legend.

    #Setting plot parameters
    plt.xlabel('Magnetic Field (T)',weight='semibold')
    plt.ylabel('Energy (cm^-1)',weight='semibold')
    plt.legend(loc='best')
    plt.title(Inputdeck['plottitle'])
    plt.tight_layout()
    if savefig==True:
        plt.savefig(f'{figname}_Excited.png',dpi=800) #Save the figure if requested
    return  [Gfan,Efan,Gfan_low,Efan_low]
    
        
  #%%  
    
    
def Zeeman_Main(Inputdeck):
    '''
    This is the main zeeman function. THe input is a dictionary with many optional keywords for more modular functionaly. These optional keywords are checked for to branch the function.
    

    Parameters
    ----------
    Inputdeck :     InputdeckExample = {
                     's_ground': 3  , #Spin multiplicity for ground state ^(2s +1)L_J , int or half int
                     's_excited': 3 , #Spin multiplicity for excited state ^(2s +1)L_J, int or half int
                     'l_ground': 0  , #Orbital Angular Momentum of ground state, int or half int
                     'l_excited': 1 , #Orbital Angular Momentum of excited state, int or half int
                     'E_ground': np.array([2951.29 ]) , #np.array of lower state, Lowest J first, cm^-1
                     'E_excited': np.array([26229.77,27488.11,27889.68]) ,  #np.array of upper state, Lowest J first, cm^-1
                     'Bmag': .6 , #Magnetic Field  [T]
                     'b_angle': 60 , #Angle between LoS and Bmax
                     'Pol_angle': 0 , #Angle polarizing filter makes with max linear transmission
                     'amu' : 183 , #Weight in AMU
                     'Convfxn': 'Gaussian' , #Optional: 'Gaussian', 'Skewed'. Future work will allow custom functions or modifying the skewed lorrentzian
                     'Temp' : 300, #Temperature in K. Used for Gaussian convolution
                     'specstep' : 0.002222 , #stepsize of linefunction [nm]
                     'fxnwindow': 0.5 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                     'plottitle': 'Excample Title' , #Title for plotting (optional)                     
                     'plot_window' : [425,435], #Min and max for plot window range (nm), convolutions will only take place within this range
                     'HFS_G' : [505.5e6] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                     'HFS_E': [0, 496.2e6, 440.5e6] , #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                     'I_spin' : -0.5, #Optional: nuclear spin, I, including parity. Tabulated by Stone et al.
                     'mu_I' :  .11778476 , #Optional: Nuclear dipole moment as tabulated by Stone.              
                     }

    Returns
    -------
    data : TYPE
        DESCRIPTION.

    '''
    

    #Define local variables from the inputdeck
    S_ground = Inputdeck['s_ground']
    S_excited = Inputdeck['s_excited']
    L_ground = Inputdeck['l_ground']
    L_excited = Inputdeck['l_excited']
    E_ground = Inputdeck['E_ground']
    E_excited = Inputdeck['E_excited']
    Bfield_ext = Inputdeck['Bmag']
    Bext_angle = Inputdeck['b_angle']
    
    
    #Set polarization filter angle to 0 and define no polarization filter.
    Polangle = 0
    Polfilter = False
    if 'Pol_angle' in Inputdeck: #Check if the inputdeck has a pol_angle and adjust parameters accordingly.
        Polangle = Inputdeck['Pol_angle']
        Polfilter = True

    if 'I_spin' in Inputdeck: #Check if the input deck has nuclear spin as to calculate HFS splitting and zeeman effect or not.
        HFSConstsG = Inputdeck['HFS_G']
        HFSConstsE = Inputdeck['HFS_E']

        #Use spin, spin parity, and nuclear dipole moment to calculate lande g_i nuclear factor. Might be off by a factor of 1836?
        I_spin = np.abs(Inputdeck['I_spin'])
        I_parity = np.sign(Inputdeck['I_spin'])
        mu_I = Inputdeck['mu_I']
        thisg_I = lande_gi(I_spin, I_parity, mu_I)
        
        #Process the wavefunction matrix and energies. Preliminary sorting is done so that the final state order should be roughly as expected.
        Ground_Level_HFS = stateprocess_HFS(L_ground, S_ground, I_spin, HFSConstsG, thisg_I, E_ground, Bfield_ext)
        Excited_Level_HFS = stateprocess_HFS(L_excited, S_excited, I_spin, HFSConstsE, thisg_I, E_excited, Bfield_ext)
        
        #Take the sorted states, the state matrix, and energy matrix and calculate the zeeman splitting (with HFS)
        Ground_Zeeman_HFS= Zeeman_func_HFS(Ground_Level_HFS, Bfield_ext, thisg_I)
        Excited_Zeeman_HFS = Zeeman_func_HFS(Excited_Level_HFS, Bfield_ext, thisg_I)
        #Use the eigenvectors and values from the diagonalization of the hamiltonian to calculate the relative signal strength. This function uses a combination of Isler1997 and Metcalf 4.32/4.33
        Z_sig_HFS = Zeeman_signal_HFS(Ground_Level_HFS, Ground_Zeeman_HFS ,Excited_Level_HFS,Excited_Zeeman_HFS,Bangle=Bext_angle,gamma=Polangle, Filter=Polfilter )  

        #Output a bunch of data
        data = {}
        data['signal'] = Z_sig_HFS[0]
        data['wave_vac'] = Z_sig_HFS[6]
        data['wave_air'] = Z_sig_HFS[1]
        data['HFSG_Energy'] = Ground_Level_HFS[2]
        data['HFSE_Energy'] = Excited_Level_HFS[2]
        data['statelist_Ground'] = Ground_Level_HFS[-1]
        data['statelist_Excited'] = Excited_Level_HFS[-1]
        data['RedDipoleOp'] = Z_sig_HFS[-3] #This is the tensor/kronecker product of the upper and lower state matrices after the field has been included. afcf in calculations. Might not be what is expected.
        data['w3coeffs_ground'] = Ground_Level_HFS[0]
        data['w3coeffs_excited'] = Excited_Level_HFS[0]


    #If there is no nuclear spin in the input deck, do a similar calculation but no hpyerfine splitting ins included nor nuclear zeeman splitting.
    else:
        #Process the wavefunction matrix and energies. Preliminary sorting is done so that the final state order should be roughly as expected.
        Ground_Level = stateprocess(L_ground,S_ground,E_ground,Bfield_ext)
        Excited_Level = stateprocess(L_excited,S_excited,E_excited,Bfield_ext)
        data = {}
        #Preliminary workflow for attempting stark shift calculation. Presently doesn't provide reasonable results.
        if 'Efield' in Inputdeck:
            L_Stark = StarkShift(Ground_Level, Excited_Level, Efield=Inputdeck['Efield'])

            data['Groundlvl'] = Ground_Level
            L_ZSG = ZeemanStark_func(Ground_Level, L_Stark[0], Bfield_ext)
            L_ZSE = ZeemanStark_func(Excited_Level, L_Stark[1], Bfield_ext)
            Z_sig = Zeeman_signal(Ground_Level, L_ZSG, Excited_Level,L_ZSE,Bangle=Bext_angle,gamma=Polangle, Filter=Polfilter ) 
            data['starktest'] = L_Stark   
        #IF no electrc field is included in input deck, calculate normal zeeman shift.
        else:
            #Use ground and excited states to calculate zeeman splitting and diangonalize the energy hamiltonian to obtain eigenvalues and vectors.
            Ground_Zeeman = Zeeman_func(Ground_Level, Bfield_ext)
            Excited_Zeeman = Zeeman_func(Excited_Level, Bfield_ext)
            #Use the new eigenvectors as the basis to calculate the relative signal strengths for the transitions
            Z_sig = Zeeman_signal(Ground_Level, Ground_Zeeman, Excited_Level,Excited_Zeeman,Bangle=Bext_angle,gamma=Polangle, Filter=Polfilter )              
            
            data['EigGround'] = Ground_Zeeman
            data['EigExcited'] =Excited_Zeeman

        data['signal'] = Z_sig[0]
        data['wave_vac'] = Z_sig[6]
        data['wave_air'] = Z_sig[1]
        data['statelist_Ground'] = Ground_Level[-1]
        data['statelist_Excited'] = Excited_Level[-1]
        
        #Energy level output with lowfield approximation of zeeman splitting output.
        data['E0_Ground'] = Ground_Level[3]
        data['E0_Excited'] = Excited_Level[3]

        data['RedDipoleOp'] = Z_sig[-3]        

    #This check is for Convfxn. There is currently Gaussian and skewed (which is a skewed lorrentzian based on an instrument function).
    #Would like to add the ability for the user to meodify the skewed function, add their own instrument function lineshape, and regular lorrentzian profile
    #Might re-think some of these conditionals for the possible inputs and how they are checked for. 
    if 'Convfxn' in Inputdeck:
        Convtype = Inputdeck['Convfxn']
        stepsize = Inputdeck['specstep']
        convtemp = Inputdeck['Temp']
       
        if Inputdeck['Convfxn'] =='Gaussian':
            convtemp = Inputdeck['Temp']
        if 'fxnwindow' in Inputdeck:            
            fxnwindow = Inputdeck['fxnwindow']
        else:
            fxnwindow = 0.1 #This is how far the bin window will extend beyond the edge of the window. May need to be extended if large broadening.
    else:
        #Setting the convolution window for the function

        #This is the default configuration using the skewed lorrentzian that the Epscore Spectrometer uses.
        Convtype = 'Skewed'
        convtemp = 300
        
        if 'fxnwindow' in Inputdeck:            
            fxnwindow = Inputdeck['fxnwindow']
        else:
            fxnwindow = 0.1 #This is how far the bin window will extend beyond the edge of the window. May need to be extended if large broadening.
        if 'specstep' in Inputdeck:
            
            stepsize = Inputdeck['specstep']
        else:
            stepsize = 0.002222 #This is the resolution of the spectrometer pixels in nm. For smooth curves, I typically make this size 1/10th that of the spectrometer.
        

    ConvolvedSpect, redsticks = Convol_Spect(data['wave_air'], data['signal'], Inputdeck['plot_window'], Inputdeck['amu'], stepsize,Temperature_in = convtemp, functiontype = Convtype , wind_size=fxnwindow)    
       
    data['SpecOut'] = ConvolvedSpect
    data['reduced_sticks'] = redsticks
    return data


def make_stickbins(sticks_in,signals_in, windsize=.5):
    #This function collects the signals into wavelength bins where a bin ends when the next signal stick is more than 1nm from the top of the bin.
    #
    # reduceme = np.array(signals_in)>0.01*np.max(signals_in)/len(signals_in) #This was an attempt to exclude anything too small, but it doesn't seem to always work right.
    # sigred = np.array(signals_in)[reduceme]
    # wavered = np.array(sticks_in)[reduceme]
    
    sigred = signals_in
    wavered = sticks_in
    #Set initial min and max of the wavelength window
    min0 = np.min(wavered)
    max0 = np.max(wavered)
    
    
    tlen = len(wavered) #Length of wavered
    
    #Set both tmin and tmax to the bottom of the overall wavelength window.
    tmax = min0
    tmin = min0
    
    wavewind = []
    
    while tmax!=max0: #While the temp upper end of the window is not the overall max, keep doing this
        
        if tlen==1: #The idea is that the top will keep going up until the window only has one stick in it, then the bottom is brought up to be the first stick after this point.
        #This window then goes again for the next bin.
            wavewind.append([tmin-windsize,tmax+windsize]) #append the previous bin that was found via the the else conditionals below.
            
            tred1 = wavered>tmax #since tmax is the top of the previous bin, find all sticks above this
            tmin = np.min(wavered[tred1]) #New bin min is the first stick of the next bin.
            
            tred2 = wavered<tmin+windsize
            tred3 = tred1 & tred2
            ttop = wavered[tred3]  
            tmax = np.max(ttop)
            tlen = len(ttop)
    
        else:
            tred1 = wavered>=tmax       #Where reduced wave is greater than or equal to tmax
            tred2 = wavered<=(tmax+windsize) #Where reducedwave is less than or equal to tmax+windsize
            tred3 = tred1 & tred2 #Combine both conditionals
            ttop = wavered[tred3] #ttop is now wavered[tmax<x<tmax+windsize]
            tmax = np.max(ttop)
            tlen = len(ttop)
    #Add the last bin.
    wavewind.append([tmin-windsize,tmax+windsize])
    wavebins = []
    signalbins = []
    for windows in wavewind:
        l1 = wavered>windows[0]
        l2 = wavered<windows[1]
        l3 = l1 & l2 #This lets you combine two if conditionals and use them as index slices for lists.
        wavebins.append(wavered[l3])
        signalbins.append(sigred[l3])
    data= {}
    data['binned_waves'] = wavebins
    data['binned_signals'] = signalbins
    data['bin_windows'] = wavewind
    return data


def DopplerFWHM(x0,Temp_in,amu):
    #This is the one that is actually used, it's doppler in wavelength
    return (1/x0)*np.sqrt(8*np.log(2))*np.sqrt(kboltz*Temp_in/(amu*m_prot))


def Gaussian(x, x0, Temp_in, amu, Intensity=1):
    #Gaussian in wavelength. Temperature is in units of K, not eV.
    tdop = DopplerFWHM(x0,Temp_in,amu)
    xx = c_light*(1/x-1/x0)/tdop
    return np.exp(-np.log(2)*xx**2)


def Convol_Sticks(x0,xwind,Temp_in, amu_in, stepsize_nm, function = "Skewed"):
    '''
    Assume that the inputs are all in units of nm.
    The difference in this function is that I'm assuming that this subroutine will be called after the sticks have already been limited to a larger window.
    This means that the full lineshape for a given convolution may not always be plotted, but then ensures that any given larger convolution will always end up the same length of array.
    '''
    # hstep = stepsize_nm # step size in nm!! (Want 10x the steps per pixe for better shape?)
    hstep = 0.1*stepsize_nm # step size in nm!! (Want 10x the steps per pixe for better shape?)

    xmin = np.min(xwind)
    xmax = np.max(xwind)

    numwavesteps = int(round((xmax-xmin)/hstep))
    waveout = np.linspace(xmin, xmax, numwavesteps)

    if function=="Skewed":    
        #This one just does convolution over the entire range, should be fine because the range is fed by the binned stick windows.
        tlorry = lo_sk(waveout,x0,pixel_size =stepsize_nm)
        specout = tlorry
    elif function=="Gaussian":
        tlorry = Gaussian(1e-9*waveout,1e-9*x0,Temp_in,amu_in)
        specout = tlorry

    return [waveout,np.array(specout)]

def Convol_Spect(waves_in, signal_in, wave_wind, atomic_weight_amu, resolution_nm,Temperature_in=300, functiontype = "Skewed", wind_size =.05):
    
    #This block reduces the waves and signals to be only values that are greater than 0.01*max/numsticks

    reducer1 =  np.array(signal_in)>1e-9
    
    reducer2 = np.logical_and(waves_in>np.min(wave_wind),waves_in<np.max(wave_wind))
    # reducer3 = reducer2
    reducer3 = reducer1 & reducer2
    twavered = waves_in[reducer3]
    tsigred = signal_in[reducer3]
    
    #This function makes "bins" of anytime there are sticks that are within +/-1 nm of any other sticks (including itself)
    #This vastly speeds up the calculations as anytime this isn't true will just be padded with 0 signal
    #0.5nm is like 50 fwhm at 100000K or something.
    tbins = make_stickbins(twavered, tsigred, windsize = wind_size)
       
    tempwavs = []
    tempsigs = []
    for sticks, sigs, winds in zip(tbins['binned_waves'], tbins['binned_signals'], tbins['bin_windows']):
        twav = np.linspace(winds[0], winds[1], int(round(10*(winds[1] - winds[0])/resolution_nm))) #Factor of 10 for better resolution
        # twav = np.linspace(winds[0], winds[1], int(round((winds[1] - winds[0])/resolution_nm))) #Factor of 10 for better resolution

        tsig = np.zeros_like(twav)
        for wavs, sigvals in zip(sticks,sigs):
            tspec = Convol_Sticks(wavs, [winds[0],winds[1]], Temperature_in, atomic_weight_amu,resolution_nm,function=functiontype)
            tsig = tsig + sigvals*tspec[1]
        tempwavs.extend(twav)
        tempsigs.extend(tsig)
    spectholder = [tempwavs,tempsigs]
    
    return     spectholder, [twavered, tsigred]





def read_Spectra(filename, scalewave = 1,headercount = 1,keepheaders = False):
    '''
    read_Spectra is a general spectra import that is flexible enough for multiple columns of spectrometer data (Say for multiple field strengths)
    The default assumption is that the wavelengths are in nm, but you can set scalewave to adjust as you input if needed.
    If your data file has headers with details you wish to keep correlated, the keepheaders flag will keep all headers as an output.

    Parameters
    ----------
    filename : string/path
        Destination/name of the file to be imported
    scalewave : float, optional
        Scaling for wavelengths in should it be desired to quickly convert from A to nm or etc.The default is 1.
    headercount : float, optional
        How many header rows are in the data.The default is 1.
    keepheaders : BOOL, optional
        Flag for whether the data within the header of the file is kept and returned for the user.The default is False.

    Returns
    -------
    specdata_out : TYPE
        If the header data is not kept:
            Output is [wavelengths,signals]
        If headers are kept, output is [[wavelengths,signals],[headers]]

    '''
    with open(f'{filename}', mode='r', newline='') as readme:
        reader = csv.reader(readme, delimiter=',', quoting=csv.QUOTE_NONE)
        headers = []
        for i in range(headercount):
            headers.append(next(reader))
        specdat = []
        
        #Read through the rows and scale the wavelength by the scalefactor.
        for rows in reader:       
            specdat.append([float(x) for x in rows])       
        
        specdat = np.array(specdat).T #Transpose so that the ith entry is all wavelengths or intensities.
        specwav = specdat[0]
        
        Spectraout = []
        [Spectraout.append(x) for i,x in enumerate(specdat[1:])]    
    
        #This small bit gives magnetic field values based on the headers info
        # Bvals = [float(x) for x in headers[0][1:]]    
    
        specdata_out = {}
        specdata_out['wavelength'] = [scalewave*x for x in specwav] #output wavelength, assumed in nm.
        # specdata_out['Bvals'] = Bvals
        specdata_out['signals'] = Spectraout
        
    if keepheaders==True:
        return [specdata_out , headers]
    else:
        return specdata_out  

