# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 16:07:14 2023

@author: Diash
"""
import numpy as np
from matplotlib import pyplot as plt
import os


def polarization0(mupper, mlower, Btheta,gamma=0,phi=0,Pol_Filter=False):
    #GAMMA is taken to be the angle of the filter with respect to the maximum of the signal. IE: 0 corresponds to the maximum signal for a given filter.
    #This is the function as exactly copied from Curt's original code. Pol_Filter doesn't matter but kept as a flag for easy of copy/paste
    #The Filter flag is still here, as the thought is the user would go in and change those values if you had polarization or not.
    #This function says gamma is in degrees, but is originally lacking a conversion, and as such polarization1 is assuming the user knew about needing rads and fixing spol values.
    
    
    if Pol_Filter==False:
        spol1 = 1
        spol2 = 1
    else:
        spol1=1
        spol2=0
    # gammaangle = np.deg2rad(gamma)
    gammaangle = gamma
    theta_rad = np.deg2rad(Btheta)#convert deg to radians
    

    
    # epsz = ealph*np.sin(theta_rad)
    # epsx = -ebeta*np.sin(phi_rad) - ealph*np.cos(theta_rad)*np.cos(phi_rad)
    # epsy = ebeta*np.cos(phi_rad) - ealph*np.cos(theta_rad)*np.sin(phi_rad)
    
    # e0 = epsz
    # e1p = (-1/np.sqrt(2))*(epsx + 1j*epsy)
    # e1m = (1/np.sqrt(2))*(epsx - 1j*epsy)
    
    # eps00 = np.sqrt(e0*e0)
    # eps11 = np.sqrt(0.5*np.abs(np.real(e1p*e1m)))
    
    epsz=np.sin(theta_rad)*(spol1*np.cos(gammaangle)+spol2*np.sin(gammaangle))
    epsx=-1.*(spol1*np.sin(gammaangle)-spol2*np.cos(gammaangle))/np.sqrt(2)
    epsy=-np.cos(theta_rad)*(spol1*np.cos(gammaangle)+spol2*np.sin(gammaangle))/np.sqrt(2)
    
    delta_am =  np.abs(np.array(mupper) - np.array(mlower))
    rad_tmp = 0
    if(delta_am ==1):
        rad_tmp = np.sqrt(epsx**2 + epsy**2)
    if(delta_am == 0):
        rad_tmp = epsz
    return rad_tmp

def polarization1(mupper, mlower, Btheta,gamma=0,phi=0,Pol_Filter=False):
    #GAMMA is taken to be the angle of the filter with respect to the maximum of the signal. IE: 0 corresponds to the maximum signal for a given filter.
    #This is the function as exactly copied from Curt's original code. Pol_Filter doesn't matter but kept as a flag for easy of copy/paste
    #This is the same as the original function, but it converts gamma into a radian function to avoid some of the "silly" mistake, but there are still many issues.
    
    
    if Pol_Filter==False:
        spol1 = 1
        spol2 = 1
    else:
        spol1=1
        spol2=0
    gammaangle = np.deg2rad(gamma)
    # gammaangle = gamma
    theta_rad = np.deg2rad(Btheta)#convert deg to radians
    

    
    # epsz = ealph*np.sin(theta_rad)
    # epsx = -ebeta*np.sin(phi_rad) - ealph*np.cos(theta_rad)*np.cos(phi_rad)
    # epsy = ebeta*np.cos(phi_rad) - ealph*np.cos(theta_rad)*np.sin(phi_rad)
    
    # e0 = epsz
    # e1p = (-1/np.sqrt(2))*(epsx + 1j*epsy)
    # e1m = (1/np.sqrt(2))*(epsx - 1j*epsy)
    
    # eps00 = np.sqrt(e0*e0)
    # eps11 = np.sqrt(0.5*np.abs(np.real(e1p*e1m)))
    
    epsz=np.sin(theta_rad)*(spol1*np.cos(gammaangle)+spol2*np.sin(gammaangle))
    epsx=-1.*(spol1*np.sin(gammaangle)-spol2*np.cos(gammaangle))/np.sqrt(2)
    epsy=-np.cos(theta_rad)*(spol1*np.cos(gammaangle)+spol2*np.sin(gammaangle))/np.sqrt(2)
    
    delta_am =  np.abs(np.array(mupper) - np.array(mlower))
    rad_tmp = 0
    if(delta_am ==1):
        rad_tmp = np.sqrt(epsx**2 + epsy**2)
    if(delta_am == 0):
        rad_tmp = epsz
    return rad_tmp
def polarization2(mupper, mlower, Btheta,gamma=0,phi=0,Pol_Filter=False,T_Alpha = 1, T_Beta=0):
    #GAMMA is taken to be the angle of the filter with respect to the maximum of the signal. IE: 0 corresponds to the maximum signal for a given filter.
    #I assume that this function is such that the signal doesn't change at all without a filter and that the original Polarization DOES. CORRECT
    
    
    
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
   

    phi_rad = np.deg2rad(phi)
    
    epsz = ealph*np.sin(theta_rad)
    epsx = -ebeta*np.sin(phi_rad) - ealph*np.cos(theta_rad)*np.cos(phi_rad)
    epsy = ebeta*np.cos(phi_rad) - ealph*np.cos(theta_rad)*np.sin(phi_rad)
    
    e0 = epsz
    e1p = (-1/np.sqrt(2))*(epsx + 1j*epsy)
    e1m = (1/np.sqrt(2))*(epsx - 1j*epsy)
    
    eps00 = np.sqrt(e0*e0)
    eps11 = np.sqrt(0.5*np.abs(np.real(e1p*e1m)))
    
    delta_am =  np.abs(np.array(mupper) - np.array(mlower))
    rad_tmp = 0
    if(delta_am ==1):
        rad_tmp = eps11
    if(delta_am == 0):
        rad_tmp = eps00
    return rad_tmp


def polarization3(mupper, mlower, Btheta,gamma=0,Pol_Filter=False,T_Alpha = 1, T_Beta=0, phi=0):
    #GAMMA is taken to be the angle of the filter with respect to the maximum of the signal. IE: 0 corresponds to the maximum signal for a given filter.

    #This function removes phi analytically, but keeps phi as a kwarg here for ease of testing things.
    
    
    
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
    eps11 = .5*np.sqrt((np.cos(theta_rad)**2) * (ealph**2) + (ebeta**2 ))
    
    delta_am =  np.abs(np.array(mupper) - np.array(mlower))
    rad_tmp = 0
    if(delta_am ==1):
        rad_tmp = eps11
    if(delta_am == 0):
        rad_tmp = eps00
    return rad_tmp


#Testing will be focused on pi and sigma transitions. So mlower=mupper for pi and mlower-mupper = 1 and letting the angles vary. The goal being to see how the signal will change.
#Realizing that I'll likely need to include the actual spectra as an example as the relative components may change too.



vary180 = np.arange(0,181,30)
# vary360 = np.arange(0,360,10)

# vary180=[0,30]
vary180b=[0, 30, 60, 90]
#First: Vary Btheta
Thetavary = np.linspace(0,180,360)
Gammavary = np.linspace(0,360,360)


# plt.ioff()
os.chdir(r"C:/LN_Auburn/Coding/Splitting Code/LN_SplitField/PolarizationFigs/")
thetavar2 = np.arange(0,181,10)
#%%
'''
1 - Original Function, No Filter, span all, unfixed gamma angle
'''
plt.close('all')
#Theta from 0 to 180, phi and gamma both from 0 to 360, showing no change in the output except for the theta dependence.
# plt.figure()
# fig
fig, (ax1, ax2) = plt.subplots(2)
# plt.title("Varying theta (b_angle) without Filter, for gamma,phi in [0,30,60,90,120,150,180] (Orig Pol Function)")
fig.suptitle("Original, No Filter, Unfixed Gamma: Shows gamma dependence")
colorlist = ['r','g','b','brown','m','c','y','black','seagreen']
for a,colors in zip(vary180,colorlist):
    for b in vary180:
        orig_pi_NF = [polarization0(0,0,x,gamma=a,phi=b) for x in Thetavary]
        ax1.plot(Thetavary, orig_pi_NF, color=colors)# , label= f'Pi Transition, gamma={a},phi={b}')
        ax1.set_title("Pi Transition")
        # ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
        orig_sigma_NF = [polarization0(1,0,x,gamma=a,phi=b) for x in Thetavary]
        ax1.plot([0],[0], color = colors, label = f'$\gamma$ ={a}, $\phi$ = {b}')
        ax2.plot(Thetavary, orig_sigma_NF, color=colors)# , label= f'Sigma Transition, gamma={a},phi={b}')
        ax2.set_title("Sigma Transitions")
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax1.label_outer()
ax2.label_outer()
plt.ylim(-.5,1.5)
plt.xlabel('Theta (Bfield angle to LOS (deg)')
plt.show()

#%%
'''
2 - Original Function, No Filter, span all, unfixed gamma angle, but separated by gamma, showing all phi are stacked, so fxn is independent of phi
'''
plt.close('all')

fig, axs = plt.subplots(7,4,sharex=True,sharey=True)

fig.suptitle("Original, No Filter, Unfixed Gamma vs Fixed gamma: Fixed still has intensity>1" )

for lb,b in enumerate(vary180):
    for a in vary180:
        pi_nofilter =  [polarization0(0,0,x,gamma=b,phi=a) for x in Thetavary]
        axs[lb,0].plot(Thetavary,pi_nofilter, color=colorlist[lb])
        axs[lb,0].set(ylabel=f'$\gamma$ = {b}')
        sigma_nofilter =  [polarization0(1,0,x,gamma=b,phi=a) for x in Thetavary]
        axs[lb,1].plot(Thetavary,sigma_nofilter, color=colorlist[lb])
        axs[lb,0].set_ylim(-1.25,1.25)
        # axs[lb,1].plot([0],[0],color = colorlist[lb],label = f'$\gamma$ = {b}, $\phi$ = {a}')
        
        pi_nofilter1 =  [polarization1(0,0,x,gamma=b,phi=a) for x in Thetavary]
        axs[lb,2].plot(Thetavary,pi_nofilter1, color=colorlist[lb])
        axs[lb,2].set(ylabel=f'$\gamma$ = {b}')
        sigma_nofilter1 =  [polarization1(1,0,x,gamma=b,phi=a) for x in Thetavary]
        axs[lb,3].plot(Thetavary,sigma_nofilter1, color=colorlist[lb])
        axs[lb,3].plot([0],[0],color = colorlist[lb],label = f'$\gamma$ = {b}, $\phi$ = {a}')
        axs[lb,3].legend(loc='center left', bbox_to_anchor=(1, 0.5))
        
               
        
        axs[lb,0].set_ylim(-1.25,1.25)
        
axs[0,0].set_title("Pi - Orig")
axs[0,1].set_title("Sigma - Orig")
axs[0,2].set_title("Pi - Gamma Deg Fix")
axs[0,3].set_title("Sigma Gamma Deg Fix")
# plt.ylim(-.5,1.5)
plt.xlabel('Theta (Bfield angle to LOS (deg)')
plt.show()



#%%
'''
3 - Orignal GammaFix vs Full Isler Polarization Derived (With Phi)
'''
plt.close('all')

fig, axs = plt.subplots(7,4,sharex=True,sharey=True)

fig.suptitle("Original, No Filter, Unfixed Gamma" )

for lb,b in enumerate(vary180):
    for a in vary180:
        pi_nofilter1 =  [polarization1(0,0,x,gamma=b,phi=a) for x in Thetavary]
        axs[lb,0].plot(Thetavary,pi_nofilter1, color=colorlist[lb])
        axs[lb,0].set(ylabel=f'$\gamma$ = {b}')
        sigma_nofilter1 =  [polarization1(1,0,x,gamma=b,phi=a) for x in Thetavary]
        axs[lb,2].plot(Thetavary,sigma_nofilter1, color=colorlist[lb])
        axs[lb,0].set_ylim(-.25,1)
        # axs[lb,1].plot([0],[0],color = colorlist[lb],label = f'$\gamma$ = {b}, $\phi$ = {a}')
        
        pi_nofilter2 =  [polarization2(0,0,x,gamma=b,phi=a) for x in Thetavary]
        axs[lb,1].plot(Thetavary,pi_nofilter2, color=colorlist[lb])
        # axs[lb,2].set(ylabel=f'$\gamma$ = {b}')
        sigma_nofilter2 =  [polarization2(1,0,x,gamma=b,phi=a) for x in Thetavary]
        axs[lb,3].plot(Thetavary,sigma_nofilter2, color=colorlist[lb])
        axs[lb,3].plot([0],[0],color = colorlist[lb],label = f'$\gamma$ = {b}, $\phi$ = {a}')
        axs[lb,3].legend(loc='center left', bbox_to_anchor=(1, 0.5))
        
        pi_withfilter1 =  [polarization1(0,0,x,gamma=b,phi=a,Pol_Filter=True) for x in Thetavary]
        sigma_withfilter1 =  [polarization1(1,0,x,gamma=b,phi=a,Pol_Filter=True) for x in Thetavary]
        
        axs[lb,0].plot(Thetavary,pi_withfilter1, color=colorlist[lb], linestyle = ':', linewidth = 5)        
        axs[lb,2].plot(Thetavary,sigma_withfilter1, color=colorlist[lb], linestyle = ':', linewidth = 5)

        pi_withfilter2 =  [polarization2(0,0,x,gamma=b,phi=a,Pol_Filter=True) for x in Thetavary]
        sigma_withfilter2 =  [polarization2(1,0,x,gamma=b,phi=a,Pol_Filter=True) for x in Thetavary]
        
        axs[lb,1].plot(Thetavary,pi_withfilter2, color=colorlist[lb], linestyle = ':', linewidth = 5)        
        axs[lb,3].plot(Thetavary,sigma_withfilter2, color=colorlist[lb], linestyle = ':', linewidth = 5)
        
axs[0,0].set_title("Pi - Orig")
axs[0,2].set_title("Sigma - Orig")
axs[0,1].set_title("Pi - Full")
axs[0,3].set_title("Sigma - Full")
# plt.ylim(-.5,1.5)
plt.xlabel('Theta (Bfield angle to LOS (deg)')
plt.show()

#%%
'''
3 - Full Isler  vs No Phi
'''
plt.close('all')

fig, axs = plt.subplots(7,2,sharex=True,sharey=True)

fig.suptitle("Full Isler vs NoPhi, With and Without Filter" )

for lb,b in enumerate(vary180):
    for a in vary180:
        pi_nofilter2 =  [polarization2(0,0,x,gamma=b,phi=a) for x in Thetavary]
        axs[lb,0].plot(Thetavary,pi_nofilter2, color=colorlist[lb])
        axs[lb,0].set(ylabel=f'$\gamma$ = {b}')
        sigma_nofilter2 =  [polarization2(2,0,x,gamma=b,phi=a) for x in Thetavary]
        axs[lb,1].plot(Thetavary,sigma_nofilter2, color=colorlist[lb])
        axs[lb,0].set_ylim(-.25,1)
        # axs[lb,1].plot([0],[0],color = colorlist[lb],label = f'$\gamma$ = {b}, $\phi$ = {a}')
        
        pi_nofilter3 =  [polarization3(0,0,x,gamma=b,phi=a) for x in Thetavary]
        axs[lb,0].plot(Thetavary,pi_nofilter3, color=colorlist[lb])
        # axs[lb,2].set(ylabel=f'$\gamma$ = {b}')
        sigma_nofilter3 =  [polarization2(3,0,x,gamma=b,phi=a) for x in Thetavary]
        axs[lb,1].plot(Thetavary,sigma_nofilter3, color=colorlist[lb])
        # axs[lb,1].plot([0],[0],color = colorlist[lb],label = f'$\gamma$ = {b}, $\phi$ = {a}')
        
        
        
        pi_withfilter2 =  [polarization2(0,0,x,gamma=b,phi=a,Pol_Filter=True) for x in Thetavary]
        sigma_withfilter2 =  [polarization2(1,0,x,gamma=b,phi=a,Pol_Filter=True) for x in Thetavary]
        
        axs[lb,0].plot(Thetavary,pi_withfilter2, color=colorlist[lb], linestyle = ':', linewidth = 5)        
        axs[lb,1].plot(Thetavary,sigma_withfilter2, color=colorlist[lb], linestyle = ':', linewidth = 5)

        pi_withfilter3 =  [polarization3(0,0,x,gamma=b,phi=a,Pol_Filter=True) for x in Thetavary]
        sigma_withfilter3 =  [polarization3(1,0,x,gamma=b,phi=a,Pol_Filter=True) for x in Thetavary]
        
        axs[lb,0].plot(Thetavary,pi_withfilter3, color=colorlist[lb], linestyle = ':', linewidth = 5)        
        axs[lb,1].plot(Thetavary,sigma_withfilter3, color=colorlist[lb], linestyle = ':', linewidth = 5)
        
    axs[lb,1].plot([0],[0],color = colorlist[lb],label = f'$\gamma$ = {b} - NoFilter')
    axs[lb,1].plot([0],[0],color = colorlist[lb],linestyle=':',label = f'$\gamma$ = {b} Filter')
    axs[lb,1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
   

# axs[0,1].plot([0],[0],color = colorlist[lb],label = f'$\gamma$ = {b}, $\phi$ = {a}')
        
axs[0,0].set_title("Pi - Both")
axs[0,1].set_title("Sigma - Both")
# axs[0,1].set_title("Pi - No-Phi")
# axs[0,3].set_title("Sigma - No-Phi")
# plt.ylim(-.5,1.5)
plt.xlabel('Theta (Bfield angle to LOS (deg)')
plt.show()

#%%
'''
3 - Original vs fixed gamma, no filter, pi only
'''
# plt.close('all')

fig, axs = plt.subplots(7,2,sharex=True,sharey=True)

fig.suptitle("Original vs Fixed Gamma, Filter On/Off, Pi Only")

for lb,b in enumerate(vary180):
    for a in vary180:
       
        pi_nofilter =  [polarization0(0,0,x,gamma=b,phi=a) for x in Thetavary]
        pi_withfilter =  [polarization0(0,0,x,gamma=b,phi=a,Pol_Filter=True) for x in Thetavary]
        
        axs[lb,0].plot(Thetavary,pi_withfilter, color=colorlist[lb], linestyle = ':', linewidth = 5)
        axs[lb,0].plot(Thetavary,pi_nofilter, color=colorlist[lb])
               
        pi_nofilter1 =  [polarization1(0,0,x,gamma=b,phi=a) for x in Thetavary]
        pi_withfilter1 =  [polarization1(0,0,x,gamma=b,phi=a,Pol_Filter=True) for x in Thetavary]
        
        axs[lb,1].plot(Thetavary,pi_withfilter1, color=colorlist[lb], linestyle = ':', linewidth = 5)
        axs[lb,1].plot(Thetavary,pi_nofilter1, color=colorlist[lb])
        
        
        axs[lb,0].set_ylim(-1.25,1.25)
        axs[lb,1].plot([0],[0],color = colorlist[lb],label = f'$\gamma$ = {b}, $\phi$ = {a}')
        axs[lb,1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
        
    
axs[0,0].set_title("Original")
axs[0,1].set_title("FixedGamma")
# plt.ylim(-.5,1.5)
plt.xlabel('Theta (Bfield angle to LOS (deg)')
plt.show()

#%%
'''
4 - Fixed Gamma vs Function, no filter, pi only
'''
# plt.close('all')

fig, axs = plt.subplots(7,2,sharex=True,sharey=True)

fig.suptitle("Original vs Fixed Gamma, Filter On/Off, Pi Only")

for lb,b in enumerate(vary180):
    for a in vary180:
       
        pi_nofilter =  [polarization1(0,0,x,gamma=b,phi=a) for x in Thetavary]
        pi_withfilter =  [polarization1(0,0,x,gamma=b,phi=a,Pol_Filter=True) for x in Thetavary]
        
        axs[lb,0].plot(Thetavary,pi_withfilter, color=colorlist[lb], linestyle = ':', linewidth = 5)
        axs[lb,0].plot(Thetavary,pi_nofilter, color=colorlist[lb])
               
        pi_nofilter1 =  [polarization2(0,0,x,gamma=b,phi=a) for x in Thetavary]
        pi_withfilter1 =  [polarization2(0,0,x,gamma=b,phi=a,Pol_Filter=True) for x in Thetavary]
        
        axs[lb,1].plot(Thetavary,pi_withfilter1, color=colorlist[lb], linestyle = ':', linewidth = 5)
        axs[lb,1].plot(Thetavary,pi_nofilter1, color=colorlist[lb])
        
        
        axs[lb,0].set_ylim(-1.25,1.25)
        axs[lb,1].plot([0],[0],color = colorlist[lb],label = f'$\gamma$ = {b}, $\phi$ = {a}')
        axs[lb,1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
axs[0,0].set_title("Original")
axs[0,1].set_title("FixedGamma")
# plt.ylim(-.5,1.5)
plt.xlabel('Theta (Bfield angle to LOS (deg)')
plt.show()
#%%
'''
4 - Original vs fixed gamma Function, No Filter, sigma only
'''
# plt.close('all')

fig, axs = plt.subplots(7,2,sharex=True,sharey=True)

fig.suptitle("Original vs Fixed Gamma, No Filter, Sigma Only")

for lb,b in enumerate(vary180):
    for a in vary180:
        sigma_nofilter =  [polarization0(1,0,x,gamma=b,phi=a) for x in Thetavary]
        axs[lb,0].plot(Thetavary,sigma_nofilter, color=colorlist[lb])
        axs[lb,0].set(ylabel=f'$\gamma$ = {b}')
        sigma_nofilter_fixedgamma =  [polarization1(1,0,x,gamma=b,phi=a) for x in Thetavary]
        axs[lb,1].plot(Thetavary,sigma_nofilter_fixedgamma, color=colorlist[lb])
        axs[lb,0].set_ylim(-1.25,1.25)
        axs[lb,1].plot([0],[0],color = colorlist[lb],label = f'$\gamma$ = {b}, $\phi$ = {a}')
        axs[lb,1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
axs[0,0].set_title("Original")
axs[0,1].set_title("FixedGamma")
# plt.ylim(-.5,1.5)
plt.xlabel('Theta (Bfield angle to LOS (deg)')
plt.show()
#%%
'''
5 - Fixed Gamma vs Fixed-Phi, No Filter, pi only
'''
# plt.close('all')

fig, axs = plt.subplots(7,2,sharex=True,sharey=True)

fig.suptitle("Fixed Gamma vs Fixed-Phi, No Filter, Pi Only")

for lb,b in enumerate(vary180):
    for a in vary180:
        pi_nofilter_fixedgamma =  [polarization1(0,0,x,gamma=b,phi=a) for x in Thetavary]
        axs[lb,0].plot(Thetavary,pi_nofilter_fixedgamma, color=colorlist[lb])
        axs[lb,0].set(ylabel=f'$\gamma$ = {b}')
        pi_nofilter_fixedphi =  [polarization2(0,0,x,gamma=b,phi=a) for x in Thetavary]
        axs[lb,1].plot(Thetavary,pi_nofilter_fixedphi, color=colorlist[lb])
        axs[lb,0].set_ylim(-1.25,1.25)
        axs[lb,1].plot([0],[0],color = colorlist[lb],label = f'$\gamma$ = {b}, $\phi$ = {a}')
        axs[lb,1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
axs[0,0].set_title("Fixed Gamma")
axs[0,1].set_title("Fixed-Phi")
# plt.ylim(-.5,1.5)
plt.xlabel('Theta (Bfield angle to LOS (deg)')
plt.show()
#%%
'''
6 - Fixed Gamma vs Fixed-Phi, No Filter, pi only
'''
# plt.close('all')

fig, axs = plt.subplots(7,2,sharex=True,sharey=True)

fig.suptitle("Fixed Gamma vs Fixed-Phi, No Filter, Sigma Only")

for lb,b in enumerate(vary180):
    for a in vary180:
        sigma_nofilter_fixedgamma =  [polarization1(1,0,x,gamma=b,phi=a) for x in Thetavary]
        axs[lb,0].plot(Thetavary,sigma_nofilter_fixedgamma, color=colorlist[lb])
        axs[lb,0].set(ylabel=f'$\gamma$ = {b}')
        sigma_nofilter_fixedphi =  [polarization2(1,0,x,gamma=b,phi=a) for x in Thetavary]
        axs[lb,1].plot(Thetavary,sigma_nofilter_fixedphi, color=colorlist[lb])
        axs[lb,0].set_ylim(0,1.25)
        axs[lb,1].plot([0],[0],color = colorlist[lb],label = f'$\gamma$ = {b}, $\phi$ = {a}')
        axs[lb,1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
axs[0,0].set_title("Fixed Gamma")
axs[0,1].set_title("Fixed-Phi")
# plt.ylim(-.5,1.5)
plt.xlabel('Theta (Bfield angle to LOS (deg)')
plt.show()
#%%
'''
7 - Fixed-Phi vs Fixed-NoPhi, No Filter+Filter, Pi+Sigma
'''
# plt.close('all')

fig, axs = plt.subplots(7,4,sharex=True,sharey=True)

fig.suptitle("Fixed Gamma vs Fixed-Phi, No Filter, Sigma Only")

for lb,b in enumerate(vary180):
    for a in vary180:
        axs[lb,0].set(ylabel=f'$\gamma$ = {b}')
        axs[lb,0].set_ylim(0,1.25)

        pi_nf_phi =  [polarization2(0,0,x,gamma=b,phi=a, Pol_Filter=False) for x in Thetavary]
        axs[lb,0].plot(Thetavary,pi_nf_phi, color=colorlist[lb])

        pi_nf_nophi =  [polarization3(0,0,x,gamma=b,phi=a, Pol_Filter=False) for x in Thetavary]
        axs[lb,0].plot(Thetavary,pi_nf_nophi, color=colorlist[lb])
        
        sigma_nf_phi =  [polarization2(1,0,x,gamma=b,phi=a, Pol_Filter=False) for x in Thetavary]
        axs[lb,1].plot(Thetavary,sigma_nf_phi, color=colorlist[lb])

        sigma_nf_nophi =  [polarization3(1,0,x,gamma=b,phi=a, Pol_Filter=False) for x in Thetavary]
        axs[lb,1].plot(Thetavary,sigma_nf_nophi, color=colorlist[lb])

        pi_wf_phi =  [polarization2(0,0,x,gamma=b,phi=a, Pol_Filter=True) for x in Thetavary]
        axs[lb,0].plot(Thetavary,pi_nf_phi, color=colorlist[lb])

        pi_wf_nophi =  [polarization3(0,0,x,gamma=b,phi=a, Pol_Filter=True) for x in Thetavary]
        axs[lb,0].plot(Thetavary,pi_nf_nophi, color=colorlist[lb])
        
        sigma_wf_phi =  [polarization2(1,0,x,gamma=b,phi=a, Pol_Filter=True) for x in Thetavary]
        axs[lb,1].plot(Thetavary,sigma_nf_phi, color=colorlist[lb])

        sigma_wf_nophi =  [polarization3(1,0,x,gamma=b,phi=a, Pol_Filter=True) for x in Thetavary]
        axs[lb,1].plot(Thetavary,sigma_nf_nophi, color=colorlist[lb])
    

        
        
axs[0,0].set_title("Fixed Gamma")
axs[0,1].set_title("Fixed-Phi")
# plt.ylim(-.5,1.5)
plt.xlabel('Theta (Bfield angle to LOS (deg)')
plt.show()