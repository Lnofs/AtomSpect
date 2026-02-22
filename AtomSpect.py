# -*- coding: utf-8 -*-
"""
AtomSpect main function definitions.
@author: Leo Nofs
"""

import numpy as np
import scipy as scp
import csv
from matplotlib import pyplot as plt
from matplotlib import ticker
import tkinter as tk
from tkinter import simpledialog

from matplotlib.widgets import Button, Slider

from sympy.physics.quantum.cg import CG
from sympy.physics.wigner import wigner_3j as w3j, wigner_6j as w6j


# %%Constants

c_light = 299792458  # Speed of light in m/s
hbar = 1.0546e-34  # Plank's constant hbar in Angular
hh = 2.0*np.pi*hbar  # Plank's constant in frequency
muB = 9.27401e-24  # Bohr magneton J/T
muP = 2.79284734  # proton magneton (in terms of nuclear magnetons, nm)
muN = 1.5210322e-3  # Nucleon magneton (not nuclear dipole moment), in terms of muB
muN2 = muB/1836
mu_B = muB/hbar  # Common value of Bohr_Mag/hbar which ensures the computation doesn't drop small terms. [Hz/T]
g_L = 1  # Lande Orbital g-factor, assuming infinite mass, otherwise is g_L = 1 - me/mnuc
g_S = 2.002319  # Lande electron spin g-factor
a_o = 5.29e-11  # Bohr radii in m
mu_0 = 4*np.pi*1e-7  # freespace mag permeability
eps0 = 8.85e-12  # freespace permittivity
kboltz = 1.381e-23  # Boltzmann Constant
m_prot = 1.667e-27  # Mass of proton in kg
m_elec = 9.11e-31
q_charge = 1.602e-19
hhev = 4.135668e-15  # Plank's constant in Ev/Hz

muBEV = 5.7883818e-5  # bohr magneton in eV/Tesla
muBcm = muBEV/(1.23981e-4)  # bohr magneton in cm^-1

# specres = 10 #Number of points evalulated per pixel or function stepsize. 10 gives smoother curves for synthetic spectra vs spectrometer data


# %%Lande Factors
def lande_gi(I, u_I):
    '''
    Generates the lande g_i factor given the spin (I)and nuclear moment.
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
    # Technically, this should be the reduced mass here.
    thisgi = (u_I/I)*(m_elec/m_prot)  # This is gi_prime and is used in zeeman calculations.
    return thisgi


def lande_gj2(J, L, S):
    # Lande gJ as derived in Griffiths and others. Comes from projection of angular momentum and averaging of L and S
    #This has one less step of calculation than lande_gj function but should be equivalent.
    if J == 0:
        return 0
    else:
        return (g_L*(J*(J+1) - S*(S+1) + L*(L+1))/(2*J*(J+1)) + g_S*((J*(J+1) + S*(S+1) - L*(L+1))/(2*J*(J+1))))


def lande_gj(J, L, S):
    # Lande gJ as derived in Griffiths and others. Comes from projection of angular momentum and averaging of L and S
    if J == 0: #There is no lande effect if J=0 due to a singularity (1/J)
        return 0
    else:
        return 1 + (J*(J+1) + S*(S+1) - L*(L+1))/(2*J*(J+1))


def lande_gf(F, I, J, L, S, g_I):
    # Lande gF, comes from averaging I, L, and S similar to gJ.
    thisgj = lande_gj(J, L, S)
    #There is no lande effect if F=0, due to a singularity (1/F)
    if F == 0:
        return 0
    else:
        thisgf = thisgj*(F*(F+1) - I*(I+1) + J*(J+1)) / (2*F*(F+1)) + g_I*(F*(F+1) + I*(I+1) - J*(J+1)) / (2*F*(F+1))
        return thisgf


def Delta_E_HFS(I, J, F, A, B):
    #Calculate the shift in energy due to the hyperfine effect, where A and B are the hyperfine constants.
    #These are in units of Hz typically. Output is in units of cm^-1
    #hhev/(1.24e-4) is a conversion factor that takes the Hz and results in cm^-1
    #See SteckRb for a good overview of hyperfine. Kopfermann is a very detailed source as well.
    
    # First need to define a constant K that will be used in the calculations
    K = F*(F+1) - I*(I+1) - J*(J+1)
    #Calculate A_hfs
    E_Ahfs = 0.5*A*K * hhev/(1.24e-4)
    #Check for the cases when the B_hfs contribution must be zero (I=0,1/2, J=0,1/2, or if B=0 )
    if I == np.abs(0.5) or J == 0.5 or I == 0 or J == 0 or B == 0:
        E_Bhfs = 0
    else:
        #If none of those cases are true, calculate the shift from the quadropole term.
        E_Bhfs = B*((0.75*K*(K+1) - I*(I+1)*J*(J+1))/(2*I*(2*I - 1)*J*(2*J - 1)))*(hhev/1.24e-4)
    return E_Ahfs + E_Bhfs


# For conversion factors, 1 = (cm^-1)/(1.2398e-4 eV) = (cm^-1)/(100*c*Hz) = (1.2398e-4 eV)/(100*c*Hz)

# %%Energy Conversion
def Hz_to_nm(Ediff):
    # Converts from Hz to nm
    return 2.99792458E17 / Ediff


def nm_to_Hz(Ediff):
    # Converts from nm to Hz
    return 2.99792458E17 / Ediff


def ev_to_Hz(Ediff):
    # Converts from eV to Mz
    return (1e-6)*Ediff*2.41804e14


def cm_to_nm(E_diff):
    # Converts from energy in cm^-1 to wavelength in nm.
    return 1e7/E_diff


def cm_to_Hz(E_diff, MHz=False):
    # Converts from energy in cm^-1 to energy in Hz (E = hf = hc/lambda)
    if MHz:
        return 100*c_light*E_diff/1e6  # Return result in MHz.
    else:
        return 100*c_light*E_diff


def cm_to_ev(E_diff):
    # Convert from cm^-1 to eV
    return 1.2398e-4*E_diff


def ev_to_cm(E_diff):
    # Convert from eV to cm^-1
    return E_diff/(1.2398e-4)
# %%Wavelength


def Vac_to_air(wavein):
    '''
    Converts wavelength from vacuum to air by using the index of refaction.
    
    The calculation is based on Morton2000, so is done for "dry" air and
    1 atmosphere of presure, 14% C and 0.045%CO_2 by volume.
    
    Parameters
    ----------
    wavein : TYPE
        Wavelength in (assumed in nm)

    Returns
    -------
    TYPE
        Wavelength_out in nm

    '''

    # Using Morton2000
    wavevac_angs = wavein*10 #Convert to angstrom to match the paper
    #Define a common factor
    s = (1e4)/wavevac_angs
    #Calculate n_air
    n_air = 1 + 8.34254e-5 + (2.406147e-2)/(130 - s**2) + (1.5998e-4)/(38.9 - s**2)
    #Calculate the new wavelength based on the index of refraction
    wavelengths = wavein / n_air

    return wavelengths


def air_to_vac(wavein):
    '''
    Converts wavelength from vacuum to air, doing the opposite of Vac_to_air
    The calculation is based on Morton2000, so is done for "dry" air and
    1 atmosphere of presure, 14% C and 0.045%CO_2 by volume.
    
    Parameters
    ----------
    wavein : TYPE
        Wavelength in (assumed in nm)
    
    Returns
    -------
    TYPE
        Wavelength_out in nm

    '''

    # Using Morton2000
    wavevac_angs = wavein*10 #Convert to angstrom to match the paper
    #Define a common factor
    s = (1e4)/wavevac_angs
    #Calculate n_air
    n_air = 1 + 8.34254e-5 + (2.406147e-2)/(130 - s**2) + (1.5998e-4)/(38.9 - s**2)
    #Calculate the new wavelength based on the index of refraction
    wavelengths = wavein * n_air

    return wavelengths

# %%A_Ein


def A_Ein(E_in, dipolesig):
    '''
    Calculate the Einstein A coefficient using the dipole moment.
    Bransden QM 
    A_ab = (womega^3 * e*2)/(3\pi eps0hbar*c**3) * |<a|r|b>|^2

    There is also the intensity expression from Bethe and Salpeter - 59.10
    
    I = (4 q^2 omega^4)/(3 c^3) |r|^2 - again pol and spacially averaged.
    
    Compare to 59.11 for Einstein A as
    A = (4 q^2 omega^3)/(3 hbarc^3) |r|^2 - again pol and spacially averaged.
    
    These are in erg/sec units.
    There is definitely interesting physics and math to discover in these relations
    As well as why the expected expressions do not precisely yield the A coefficients.
    What is the missing component, and what is my calulation from the dipolestr actually giving?
    
    
    The e
    Parameters
    ----------
    E_in : float
        Energy in wavelength (nm) for the given line
    dipolesig : float
        Scalar of the dipole sig strength for the given transition.

    Returns
    -------
    tAein : TYPE
        DESCRIPTION.

    '''
 
    # omega = (Ea-Eb)/hbar
    tomega = nm_to_Hz(E_in)/(2*np.pi) #Note also the lack of hbar here, just h.
    #This expression looks close to eq 4.119 in Bransden Atoms (Pg190/182)
    tAein = ((tomega**2 * q_charge**2)/(np.pi*eps0*hh*c_light**3)) * dipolesig  # Not exactly as given by papers, but works?
   
    return tAein

# %%Plot


def plotarray(shape, figure_size=(10, 8),fontsize=18):
    '''
    This is a general function which makes a blank array of plots with shape (rows,col)
    and returns the fig object and array of plots. This is used with PlotFunction to feed plot data
    towards a streamlined plotting experience.
    A generalized version of these functions would be ideal, but the scope here is for use
    in this code.

    Parameters
    ----------
    shape : tuple
        Dimensions of the plot (Rows, Colums)
    figure_size : tuple (Width, Height)
        Plot size in inches for the resulting figure. The default is (10, 8).
    fontsize : Int
        Size for the font in the figure. The default is 18.

    Returns
    -------
    fig : natplotlib figure object
        callable and assignable figure object.
    plotmat : array
        Array of plots where individual axes can be indexed for plotting on.

    '''
    
    fig = plt.figure(figsize=figure_size)  #Make the main figure
    plotmat = np.zeros(shape, dtype='O') #Make an empty np array of datatype 'Object"
    
    #Set the font size
    font = {
            'weight' : 'normal',
            'size'   : fontsize}
    plt.rc('font', **font)


    #1. Go through the array and define the plots in the subfigure.
    #Iterate through the range and set that the columns and rows have shared axes.
    
    for rows in range(shape[0]):
        for cols in range(shape[1]):
            #First row and column have unique instructions since they define what is shared.
            if rows == 0:
                if cols == 0:
                    plotmat[rows, cols] = plt.subplot2grid(shape, (rows, cols))
                else:
                    plotmat[rows, cols] = plt.subplot2grid(shape, (rows, cols), sharey=plotmat[0, 0])

            else:
                plotmat[rows, cols] = plt.subplot2grid(shape, (rows, cols), sharey=plotmat[0, 0], sharex=plotmat[0, cols])
   
    # Get rid fo the labels and ticks that are covered up for x axis.
    for rows in range(shape[0]-1):
        for cols in range(shape[1]):
            plotmat[rows, cols].tick_params('x', labelbottom=False)
    # Get rid of y axis labels, but only for the later columns.
    for rows in range(shape[0]):
        for cols in range(1, shape[1]):
            plotmat[rows, cols].tick_params('y', labelleft=False)
            
    #Set the x and y labels for standard use case of AtomSpect.
    fig.supxlabel('Wavelength (nm)', weight='semibold')
    fig.supylabel('Intensity (arb)', weight='semibold')
    #Adjust the spacing to be closer together
    plt.subplots_adjust(hspace=0, wspace=0)

    #Return a figure and plotmatrix objects that can be called for other functions.
    return fig, plotmat



def parse_plot_data(plot_object):
    '''
    Standardizes input into a list of [wave, intensity] arrays.
    This handles the Zeeman output dictionary, a single [wave,intensity] array,
    or a list containing any combination of the two.
    '''
    # If it is a dictionary, assume it is the output of Zeeman_func, and pull out the Spectra [cite: 221]
    if isinstance(plot_object, dict) and 'SpecOut' in plot_object:
        return [plot_object['SpecOut']]
    
    # If it's a single [wave, intensity] pair [cite: 221]
    if isinstance(plot_object, (list, np.ndarray)) and len(np.shape(plot_object)) == 2:
        return [plot_object]
        
    # If it's a mixed list of dictionaries and arrays [cite: 222, 223]
    parsed_data = []
    if isinstance(plot_object, list):
        for item in plot_object:
            if isinstance(item, dict) and 'SpecOut' in item:
                parsed_data.append(item['SpecOut']) # Pull out SpecOut keyword from LZeeman object [cite: 222]
            else:
                parsed_data.append(item) # Otherwise, just append the array [cite: 223]
    return parsed_data

def get_target_axis(taxs, position):
    '''
    Safely extracts a single target axis regardless of taxs array dimensions.
    '''
    # If taxs is already a single axis object (e.g. from a manual slice), return it
    if not isinstance(taxs, np.ndarray):
        return taxs
    
    # Standardize position to [row, col]. 
    # If only one index is provided, treat it as the column index.
    row = position[0] if len(position) > 0 else 0
    col = position[1] if len(position) > 1 else 0
    
    # Handle indexing based on the array dimensions returned by plotarray
    # plotarray typically returns a 2D array even for a single row.
    if taxs.ndim == 1:
        return taxs[col]
    else:
        return taxs[row, col]

def consolidate_legend(tfig, plotpol=True, plotnondip=False, legcols=None, fontsize=18):
    '''
    Gathers all handles and labels from EVERY axis in the figure (including twins)
    and generates a single, unique legend at the top.
    '''
    tfig.legends = [] # Clear existing legends
    all_handles = []
    all_labels = []

    # Iterate through all axes in the figure (twins are stored in tfig.axes as well)
    for ax in tfig.axes:
        h, l = ax.get_legend_handles_labels()
        all_handles.extend(h)
        all_labels.extend(l)

    # Filter for unique labels to prevent 'Spectrometer' etc. from repeating
    unique_map = {}
    for handle, label in zip(all_handles, all_labels):
        if label not in unique_map:
            unique_map[label] = handle
            
    if not unique_map:
        return tfig

    labels = list(unique_map.keys())
    handles = list(unique_map.values())

    # Determine layout adjustment and column count
    if legcols:
        cols = legcols
        rect_adj = [0, 0, 1, 0.95]
    elif plotnondip:
        cols = 3
        rect_adj = [0, 0, 1, 0.9]
    elif plotpol:
        cols = 5
        rect_adj = [0, 0, 1, 0.95]
    else:
        cols = 2
        rect_adj = [0, 0, 1, 0.95]

    plt.tight_layout(rect=rect_adj)
    tfig.legend(handles, labels, bbox_to_anchor=(0.5, 0.95), 
                loc='outside center', ncol=cols, fontsize=fontsize)
    return tfig

def PlotFunction(Plotobject, plot_vars, plottitle='', plotwind=None, SpectrumPlot=[], Shape=[1, 1], 
                 position=[0, 0], axsin=[], NormalizeSig=False, NormalizeScale=1, makefig=True, 
                 fig_size=(10, 8), plotpol=True, plotlabel='', plotnondip=False, dotwin=False, 
                 marker=None, marker_count=None, fontsize=None, ScaleSpectrometer=None, 
                 SpectrometerLabel='Spectrometer', SpectrometerMarker='d', SpectrometerColor='black',
                 SpectrometerLS='solid', legcols=None):
    
    '''
    Parameters
    ----------
    Plotobject : What is being plotted.
        This can be [x,y] data or a Zeeman Object dictionary output. It used to handle lists of possible inputs
        But I would suggest using the ability to define an existing axes or plotarray instead.
    plot_vars : list [color,linethickness, linestyle, label (optional), markerstyle (optional)]
        List of definable parameters for a given plot such as color, linestyle, label, line thickness, and marker style
    plottitle : str, optional
        String for the plot title
    plotwind : list or array, optional
        Sets the limits for upper and lower bound of the plot window. If not defined, will try to see if there is an 
        PlotObject['input'] which is the input dictionary for the main code which is contained in the output of the main function.
        This input can include the plot_window keyword.
        Otherwise, defaults to either the limits of the spectrometer data or the full range of the spectra convolution
    SpectrumPlot : Array, optional
        This is an array of spectrometer data [wavelength(nm),intensity]. Does not need to be normalized, can be using ScaleSpectrometer flag.
        This gets plotted in black, but can be changed.
    Shape : list or tuple, optional
        Shape of the plot array that is to be created if makefig=True is set. Ignored otherwise.
    position : index (list or int), optional [row,column]
        Position in the plot array for the given object to be plotted. Useful for making one plot then adding many overplots
    axsin : array of plot objects, optional
        This is the plotmatrix as output from PlotArray or returned as an otput of PlotFunction
        If this is fed into PlotFunction here, this is where the object will be plotted. position information is thus required.
    NormalizeSig : BOOL, optional
        Whether a pass is done on the data to normalize it (making the maximum height as 1).
        Note that the sums of the signal strenghts are already normalized to equal 1. This just uniformly scales the plot outputs.
        The default is False.
    NormalizeScale : INT, optional
        A given curve can be scaled so that the maximum height is different than 1.
        Useful for relative intensities of multiple species to match the observed spectra. The default is 1.
    makefig : BOOL, optional
        Whether a figure object should be created. Need Shape unless 1x1 plot is desired. The default is True.
        Set to FALSE 
    fig_size : tuple (Width, Height)
        Plot size in inches for the resulting figure. The default is (10, 8).
    plotpol : BOOL, optional
        Whether the polarization components are plotted as stems in addition to the spectra. The default is True.
    plotlabel : STR, optional
        Label for the plot. The default is ''.
    plotnondip : BOOL, optional
        Whether non-dipole transitions should be plotted as well (purple star stem plots)
        Despite only being a dipole based code, non-dipole components appear in the He results. The default is False.
    dotwin : BOOL, optional
        Whether a twin axis is plotted with a separate scale for the sticks than the main plot. The default is False.
    marker : INT, optional
        Index of markerlist. These are ['o','8','D','s','P','^','v','X']. Can also be 
    marker_count : INT, optional
        How many points between each marker on the plot. Allows for easier reading of the plot. The default is None.
    fontsize : TYPE, optional
       Size of the font
    ScaleSpectrometer : float, optional
       If provided, scales the spectrometer height so that the maximum value is that as provided. Default would be 1.
    SpectrometerLabel : STR, optional
        Optional label for spectrometer data. The default is 'Spectrometer'.
    SpectrometerMarker : STR, optional
        Ability to change the spectrometer data marker. The default is 'd'.
    SpectrometerColor : STR, optional
        Color for spectrometer data. The default is 'black'.
    SpectrometerLS : STR, optional
        Spectrometer Linestyle. The default is 'solid'.
    legcols : INT, optional
        How many columns the legend will have.. The default is None.
    
    Returns
    -------
    The plot array of axes objects that can be called and indexed for adding plot objects.
    After creation, this array can be fed in as axsin multiple times to keep adding data.
    From an index axes, it is possible to get additional information with axes operators.


    '''

    
    # --- 1. Initialization ---
    font_size = fontsize if fontsize else 18
    markcount = marker_count if marker_count else 20
    spectrometer_scale = ScaleSpectrometer if ScaleSpectrometer else 1
    
    if len(np.shape(plot_vars)) == 1:
        working_vars = [plot_vars]
    else:
        working_vars = plot_vars

    resolved_markers = []
    for v in working_vars:
        if len(v) > 4: 
            resolved_markers.append(v[4])
        elif marker: 
            resolved_markers.append(marker)
        else: 
            resolved_markers.append('o')

    # --- 2. Figure and Axis Setup ---
    if makefig:
        tfig, taxs = plotarray(Shape, figure_size=fig_size, fontsize=font_size)
    else:
        taxs = axsin
    
    if not isinstance(taxs, np.ndarray):
        target_axis = taxs
    else:
        try:
            target_axis = taxs[position[0], position[1]]
        except (IndexError, TypeError):
            try:
                target_axis = taxs[position[0]]
            except (IndexError, TypeError):
                target_axis = taxs.flatten()[0]
    
    tfig = target_axis.get_figure()

    # --- 3. Data Extraction and Spectrometer Input Logic ---
    tPlotObj = []
    main_dict = None
    
    # Identify the primary dictionary for metadata and sticks
    if isinstance(Plotobject, dict):
        main_dict = Plotobject
        tPlotObj.append(Plotobject.get('SpecOut', []))
    else:
        for x in Plotobject:
            if isinstance(x, dict):
                tPlotObj.append(x.get('SpecOut', []))
                if main_dict is None: main_dict = x
            else:
                tPlotObj.append(x)

    # Resolve Spectrometer Data: Check direct input first, then the Plotobject dictionary
    final_spec_data = SpectrumPlot
    if len(final_spec_data) == 0 and main_dict is not None:
        # Check for SpectrumData in the 'Input' sub-dictionary
        if 'Input' in main_dict and 'SpectrumData' in main_dict['Input']:
            final_spec_data = main_dict['Input']['SpectrumData']

    # Plot Spectrometer data if found in either location
    if len(final_spec_data) > 0:
        target_axis.plot(final_spec_data[0], Normalize(final_spec_data[1], scaling=spectrometer_scale), 
                         color=SpectrometerColor, marker=SpectrometerMarker, markevery=10, 
                         linewidth=1.5, label=SpectrometerLabel, linestyle=SpectrometerLS)

    # --- 4. Plotting Calculated Spectra ---
    for i, (tplot_data, p_vars) in enumerate(zip(tPlotObj, working_vars)):
        y_data = tplot_data[1]
        if NormalizeSig:
            scale_val = NormalizeScale[i] if isinstance(NormalizeScale, list) else NormalizeScale
            y_data = Normalize(tplot_data[1], scaling=scale_val)
            
        target_axis.plot(tplot_data[0], y_data, color=p_vars[0], linestyle=p_vars[2], 
                         linewidth=p_vars[1], label=p_vars[3], alpha=0.8,
                         marker=resolved_markers[i], markevery=markcount, markersize=9)

    # --- 5. Styling and Limits ---
    plottingwind = plotwind
    if main_dict and 'Input' in main_dict and not plotwind:
        plottingwind = main_dict['Input'].get('plot_window', None)

    if plottingwind:
        target_axis.set_xlim(np.min(plottingwind), np.max(plottingwind))
    target_axis.set_ylim(-0.1, 1.1)
    
    if plottitle: 
        target_axis.set_title(plottitle)
    if plotlabel:
        target_axis.text(0.05, 0.95, f'{plotlabel}', verticalalignment='top', 
                         transform=target_axis.transAxes, fontsize=font_size-2)

    # --- 6. Polarization Stems ---
    if plotpol and main_dict and 'pi_sticks' in main_dict:
        stem_ax = target_axis.twinx() if dotwin else target_axis
        
        if dotwin and 'reduced_sticks' in main_dict:
            stem_ax.set_ylim(-0.1 * np.max(main_dict['reduced_sticks'][1]), 
                             np.max(main_dict['reduced_sticks'][1]))
        
        stem_ax.stem(main_dict['pi_sticks'][0], main_dict['pi_sticks'][1], 
                     linefmt='C0', basefmt=" ", label='$\pi$')
        stem_ax.stem(main_dict['sigmaminus_sticks'][0], main_dict['sigmaminus_sticks'][1], 
                     linefmt='C2', markerfmt='_', basefmt=" ", label=r'$\sigma^-$')
        stem_ax.stem(main_dict['sigmaplus_sticks'][0], main_dict['sigmaplus_sticks'][1], 
                     linefmt='C3', markerfmt='x', basefmt=" ", label=r'$\sigma^+$')
        
        if plotnondip and 'other_sticks' in main_dict:
            m, _, _ = stem_ax.stem(main_dict['other_sticks'][0], main_dict['other_sticks'][1], 
                                   linefmt='C4', markerfmt='*', basefmt=" ", label=r'$\Delta m >1$')
            plt.setp(m, markersize=10, markerfacecolor='none')

    # --- 7. Legend and Final Layout ---
    consolidate_legend(tfig, plotpol=plotpol, plotnondip=plotnondip, legcols=legcols, fontsize=font_size)
    
    if makefig:
        plt.subplots_adjust(hspace=0, wspace=0)
        plt.show()

    return taxs
# %%Plot_Slider
# This section defines the plot_slider which creates interactve plots where parameters can be changed. This is a proof of concept, but is already looking very promising. Further UI work is needed.


def Zslider(B, thisTemp, thisRes, Input):
    '''
    This function is used in the Slider iteration. It takes the input values and then feeds the value into Zeeman_Main
    This is for the calculations each time the slider is changed.
    Bmag, spectrometer resolution, and Temperature (eV) are all allowed to change.
    '''
    Input['Bmag'] = B
    Input['specstep'] = thisRes
    Input['Temp'] = thisTemp*11602

    return Zeeman_Main(Input)


def Zslider_Pol(B, thisTemp, thisRes, Bangle, PolAngle, Input):
    '''
    Due to the nature of the polarization complexity in the computation, this function is separate from the main slider.
    This is used when polarization is included.
    '''
    if not PolAngle:
        try:
            del Input['PolAngle']
        except:
            pass
    else:
        Input['Pol_angle'] = PolAngle

    Input['Bmag'] = B
    Input['specstep'] = thisRes
    Input['Temp'] = thisTemp*11602
    Input['b_angle'] = Bangle

    return Zeeman_Main(Input)


def enable_manual_input(fig, slider):
    """Allow manual entry of slider values by double-clicking on the value text."""
    valtext = slider.valtext  # this is a matplotlib.text.Text object

    def on_double_click(event):
        if event.dblclick and valtext.contains(event)[0]:
            # popup Tk dialog for numeric input
            root = tk.Tk()
            root.withdraw()  # hide main window
            try:
                new_val = simpledialog.askfloat("Input", f"Enter {slider.label.get_text()} value:",
                                                initialvalue=slider.val)
                if new_val is not None:
                    slider.set_val(new_val)  # updates slider + fires callbacks
            finally:
                root.destroy()

    fig.canvas.mpl_connect("button_press_event", on_double_click)



def MakeSlider(Inputdict, Spectra, plotwind=None, do_sticks=True, do_Polplot=True, banglevary=False, polvary=False, do_other = False):
    '''
    Generate the plot with changeable sliders that allow for real-time interaction with the spectral results.
    All flags should be self-evident, with do_other referring to non-dipole transitions being shown.
    
    
    
    '''
    
    # Initial params
    init_B = 2.5  # Initial B field (T)
    init_Temp = 0.05  # Initial temp (300K)
    init_Res = 0.002  # Initial pixel size (nm)
    init_bang = 90  # Initially LoS is perpendicular to B
    init_polang = 90  # Initially block all circular light

    #Make the blank figure and set parameters
    figs, axs = plt.subplots(figsize=(16, 10))
    plt.rcParams.update({'font.size': 14})
    figs.subplots_adjust(left=0.2, bottom=0.2, right=0.8)

    #Plot the spectrometer data
    axs.plot(Spectra[0], Spectra[1], color="black", label="Spectrometer")

    #Check for polvary
    if polvary and not banglevary:
        Mycalc0 = Zslider_Pol(init_B, init_Temp, init_Res, init_bang, init_polang, Inputdict)
    elif banglevary and not polvary:
        Mycalc0 = Zslider_Pol(init_B, init_Temp, init_Res, init_bang, None, Inputdict)
    elif polvary and banglevary:
        Mycalc0 = Zslider_Pol(init_B, init_Temp, init_Res, init_bang, init_polang, Inputdict)

    else:
        Mycalc0 = Zslider(init_B, init_Temp, init_Res, Inputdict)
        Bang_Slider = None
        Polang_Slider = None
    if banglevary:
        axampBang = figs.add_axes([0.85, 0.25, 0.0225, 0.63])
        Bang_Slider = Slider(axampBang, "B_Angle to LoS", 0, 90, valstep=0.5, valinit=init_bang,
                             orientation="vertical")
    if polvary:
        axamppolang = figs.add_axes([0.93, 0.25, 0.0225, 0.63])
        Polang_Slider = Slider(axamppolang, "Pol Filter",  0, 90, valstep=0.5, valinit=init_polang,
                               orientation="vertical")
    # artists
    line, = axs.plot(Mycalc0['SpecOut'][0],
                     Normalize(Mycalc0['SpecOut'][1]), color="tab:brown",
                     lw=3.5, label="Exact", linestyle=':',marker='o', markevery=10,markersize=7, alpha=0.9)

    # line = axs.scatter(Mycalc0['SpecOut'][0],
    #                  Normalize(Mycalc0['SpecOut'][1]),color="tab:brown" ,
    #                  lw=1, label="Exact",marker = '.' , alpha = 0.75)
    numcols = 5
    if 'DoLowSig' in Inputdict:
        possibleanswers = ['Y', 'True', True, 1, 'YES', 'y', 'yes', 'Yes']
        if Inputdict['DoLowSig'] in possibleanswers:

            line2, = axs.plot(Mycalc0['SpecOutLow'][0],
                              Normalize(Mycalc0['SpecOutLow'][1]),
                              lw=2, alpha=.9, label="Lowfield", linestyle=(5, (5, 2, 1, 2)), marker='D', markevery=25, markersize=5)
            numcols = 4
    else:
        line2 = None

    if 'DoHighSig' in Inputdict:
        possibleanswers = ['Y', 'True', True, 1, 'YES', 'y', 'yes', 'Yes']
        if Inputdict['DoHighSig'] in possibleanswers:

            line3, = axs.plot(Mycalc0['SpecOutHigh'][0],
                              Normalize(Mycalc0['SpecOutHigh'][1]),
                              lw=2, alpha=.9, label="Highfield", linestyle=(5, (5, 2, 1, 2)), marker='x', markevery=25, markersize=5)
            numcols = 4
            if 'DoLowSig' in Inputdict:  # This changes the number of colums if BOTH high and low are in the inputflags.
                if Inputdict['DoLowSig'] in possibleanswers:
                    numcols = 3

    else:
        line3 = None

    if plotwind == None:
        try:
            axs.set_xlim(Inputdict['plot_window'][0], Inputdict['plot_window'][1])
        except:
            pass
    else:
        try:
            axs.set_xlim(plotwind[0], plotwind[1])
        except:
            print('Something went wrong, setting plotwind range as min/max of spectra')
            axs.set_xlim(np.min(Mycalc0['SpecOut'][0]), np.max(Mycalc0['SpecOut'][0]))

    # axs.legend(loc="best")
    vstem = None
    vstem = None
    # vl_sp = vl_pi = vl_sm = vl_o = None
    if do_Polplot:
        do_sticks = True
    if do_sticks and not do_Polplot:
        # vl = axs.vlines(Mycalc0['reduced_sticks'][0], 0, Mycalc0['reduced_sticks'][1])
        vmarks, vstem, vbase = axs.stem(Mycalc0['reduced_sticks'][0], Mycalc0['reduced_sticks'][1], linefmt='C0', basefmt=" ")

    elif do_sticks and do_Polplot:
        # vl_sp = axs.vlines(Mycalc0['sigmaplus_sticks'][0], np.zeros_like(Mycalc0['sigmaplus_sticks'][1]), Mycalc0['sigmaplus_sticks'][1], lw=2,color="tab:purple", label="σ+")
        # vl_pi = axs.vlines(Mycalc0['pi_sticks'][0], np.zeros_like(Mycalc0['pi_sticks'][1]), Mycalc0['pi_sticks'][1], color="tab:green",lw=2, label="π")
        # vl_sm = axs.vlines(Mycalc0['sigmaminus_sticks'][0], np.zeros_like(Mycalc0['sigmaminus_sticks'][1]), Mycalc0['sigmaminus_sticks'][1],lw=2, color="tab:red", label="σ−")
        # vl_o = axs.vlines(Mycalc0['other_sticks'][0], np.zeros_like(Mycalc0['other_sticks'][1]), Mycalc0['other_sticks'][1],lw=2, color="tab:cyan", label="??")

        vmarksp, vstemsp, vbasesp = axs.stem(Mycalc0['pi_sticks'][0], Mycalc0['pi_sticks'][1], linefmt='C0', markerfmt='o', basefmt=" ", label='$\pi$ component')
        vmarkpi, vstempi, vbasepi = axs.stem(Mycalc0['sigmaminus_sticks'][0], Mycalc0['sigmaminus_sticks'][1], linefmt='C2', markerfmt='_', basefmt=" ", label=r'$\sigma^-$ component')
        vmarksm, vstemsm, vbasesm = axs.stem(Mycalc0['sigmaplus_sticks'][0], Mycalc0['sigmaplus_sticks'][1], linefmt='C3', markerfmt='+', basefmt=" ",  label=r'$\sigma^+$ component')

        propspi = vmarkpi.properties()
        propssp = vmarksp.properties()
        propssm = vmarksm.properties()
        
        if do_other:
            vmarko, vstemo, vbaseo = axs.stem(Mycalc0['other_sticks'][0], Mycalc0['other_sticks'][1], linefmt='C4', markerfmt='D', basefmt=" ",  label='other component')
            propso = vmarko.properties()

        # axs.legend(loc="best")
    # sliders

    # figs.legend(bbox_to_anchor=(.5, 0.15),loc='outside center')
    figs.legend(bbox_to_anchor=(.5, .075), loc='outside center', ncol=numcols)  # Excluding highfield and lowfield.

    axampB = figs.add_axes([0.025, 0.25, 0.0225, 0.63])
    B_Slider = Slider(axampB, "B (T)", 0, 5, valinit=init_B,
                      orientation="vertical", valstep=0.01)

    axampT = figs.add_axes([0.075, 0.25, 0.0225, 0.63])
    Te_Slider = Slider(axampT, "T (eV)", 0.01, 30, valinit=init_Temp,
                       orientation="vertical")

    axampR = figs.add_axes([0.125, 0.25, 0.0225, 0.63])
    Res_Slider = Slider(axampR, r"$W_i$ (nm)", 0.002, 0.1, valstep=0.001, valinit=init_Res,
                        orientation="vertical")

    # reset button
    # resetax = figs.add_axes([0.8, 0.025, 0.1, 0.04])
    # button = Button(resetax, "Reset", hovercolor="0.975")
    plt.ion()

    artists = [line]
    if line2:
        artists.append(line2)

    if vstem:
        artists.append(vstem)

    def update(val):
        if polvary and not banglevary:
            Mycalc = Zslider_Pol(B_Slider.val, Te_Slider.val, Res_Slider.val, Inputdict['b_angle'], Polang_Slider.val, Inputdict)
        elif banglevary and not polvary:
            Mycalc = Zslider_Pol(B_Slider.val, Te_Slider.val, Res_Slider.val, Bang_Slider.val, None, Inputdict)
        elif polvary and banglevary:
            Mycalc = Zslider_Pol(B_Slider.val, Te_Slider.val, Res_Slider.val, Bang_Slider.val, Polang_Slider.val, Inputdict)
        else:
            Mycalc = Zslider(B_Slider.val, Te_Slider.val, Res_Slider.val, Inputdict)
        # update lines
        line.set_xdata(Mycalc['SpecOut'][0])
        line.set_ydata(Normalize(Mycalc['SpecOut'][1]))
        if line2:
            line2.set_xdata(Mycalc['SpecOutLow'][0])
            line2.set_ydata(Normalize(Mycalc['SpecOutLow'][1]))
        if line3:
            line3.set_xdata(Mycalc['SpecOutHigh'][0])
            line3.set_ydata(Normalize(Mycalc['SpecOutHigh'][1]))

        if do_sticks and not do_Polplot:
            x_new = Mycalc['reduced_sticks'][0]
            y_new = Mycalc['reduced_sticks'][1]

            vstem.set_paths([np.array([[xx, 0],
                                       [xx, yy]]) for (xx, yy) in zip(x_new, y_new)])

            vmarks.set_xdata(x_new)
            vmarks.set_ydata(y_new)

        elif do_sticks and do_Polplot:
            if do_other:
                for (vlstem, vlmarks, props, key, colorfmt, markers) in [
                    (vstempi, vmarkpi, propspi, "pi_sticks", "C0", 'o'),
                    (vstemsp, vmarksp, propssp, "sigmaplus_sticks", "C3", '+'),
                    (vstemsm, vmarksm, propssm, "sigmaminus_sticks", "C2", '_'),
                    (vstemo, vmarko, propso, "other_sticks", "C4", 'D'),
                ]:
                    
                    x_new = Mycalc[key][0]
                    y_new = Mycalc[key][1]
    
                    vlstem.set_paths([np.array([[xx, 0],
                                                [xx, yy]]) for (xx, yy) in zip(x_new, y_new)])
    
                    vlstem.set_color(colorfmt)
                    vlmarks.set_xdata(x_new)
                    vlmarks.set_ydata(y_new)
                    vlmarks.set_marker(markers)
                    vlmarks.set_color(colorfmt)
                    # vlmarks.update(props)
            else:
                for (vlstem, vlmarks, props, key, colorfmt, markers) in [
                    (vstempi, vmarkpi, propspi, "pi_sticks", "C0", 'o'),
                    (vstemsp, vmarksp, propssp, "sigmaplus_sticks", "C3", '+'),
                    (vstemsm, vmarksm, propssm, "sigmaminus_sticks", "C2", '_'),
                  #  (vstemo, vmarko, propso, "other_sticks", "C4", 'D'),
                ]:
                    
                    x_new = Mycalc[key][0]
                    y_new = Mycalc[key][1]
    
                    vlstem.set_paths([np.array([[xx, 0],
                                                [xx, yy]]) for (xx, yy) in zip(x_new, y_new)])
    
                    vlstem.set_color(colorfmt)
                    vlmarks.set_xdata(x_new)
                    vlmarks.set_ydata(y_new)
                    vlmarks.set_marker(markers)
                    vlmarks.set_color(colorfmt)

        figs.canvas.blit(axs.bbox)
    axs.set_ylim(-.1, 1.1)

    # slider connections
    B_Slider.on_changed(update)
    Te_Slider.on_changed(update)
    Res_Slider.on_changed(update)
    if banglevary:
        Bang_Slider.on_changed(update)
    if polvary:
        Polang_Slider.on_changed(update)

    #No clue why the reset button doesn't work anymore. =\
    # def reset(event): 
    #     B_Slider.reset()
    #     Te_Slider.reset()
    #     Res_Slider.reset()
        # if polvary:
        #     Bang_Slider.reset()
        #     Polang_Slider.reset()
    # button.on_clicked(reset)
    enable_manual_input(figs, B_Slider)
    enable_manual_input(figs, Te_Slider)
    enable_manual_input(figs, Res_Slider)
    plt.suptitle(f'{Inputdict["plottitle"]}')

    axs.set_ylabel("Normalized Intensity (Arb units)", weight='semibold')
    axs.set_xlabel('Wavelength (nm)', weight='semibold')

    plt.show()



# %%Lineshapes
# =============================================================================
# The lineshapes here are all built the same way, so a custom lineshape can be defined.
# The code uses these functions to take a "stick" with a given height and perform the convolution.
# The resulting function is scaled such that the peak of the function matches the stick height.
# All of the inputs here expect wavelengths in meters and temperature in Kelvin
# Atomic mass is in mass number. ie Rb = 85, Ar = 40, etc
# The inputs should be func(wave_range, peak, Temp, atomic mass, pixel_size = 0.005, ion_velocity=0, bidir=False)
# This would allow the function to be placed in the Convol_Sticks function, even if things like pixel_size, mass, ion_velocity, or bidirectional doppler shift are not used.
# Ensure that wave_range is a np.array of possible wavelengths for the convolution to be evaluated at.
# =============================================================================

def Gaussian(x, x0, Temp_in, amu, v_ion=0, bidir=False):
    # Dickhauer, but added the ion velocity shift to the mean calculation
    # Gradic2022 including velocity addition+-
    tdop = x0*np.sqrt(kboltz*Temp_in/(amu*m_prot*c_light*c_light))
    if bidir == True:  # First try having v_ion be an array, so that the counter-streaming velocities can be different.
        try:
            mean1 = x0*(1 + v_ion[0]/c_light)
            mean2 = x0*(1 - v_ion[1]/c_light)

            xx1 = (x-mean1)/tdop
            xx2 = (x-mean2)/tdop
            G1 = (1/(np.sqrt(2*np.pi)*tdop))*np.exp(-0.5*xx1**2)
            G2 = (1/(np.sqrt(2*np.pi)*tdop))*np.exp(-0.5*xx2**2)
            G = (G1+G2)/2

        except:  # If not an array, just do plus and minus the same velocity
            mean1 = x0*(1 + v_ion/c_light)
            mean2 = x0*(1 - v_ion/c_light)

            xx1 = (x-mean1)/tdop
            xx2 = (x-mean2)/tdop
            G1 = (1/(np.sqrt(2*np.pi)*tdop))*np.exp(-0.5*xx1**2)
            G2 = (1/(np.sqrt(2*np.pi)*tdop))*np.exp(-0.5*xx2**2)
            G = (G1+G2)/2
    else:  # If bidir is False, do normal doppler shifted.
        mean = x0*(1 + v_ion/c_light)
        xx = (x-mean)/tdop
        G = (1/(np.sqrt(2*np.pi)*tdop))*np.exp(-0.5*xx**2)
    return G


def GaussianInstrum(x, x0, Temp_in, amu, pixel_size=0, v_ion=0, bidir=False):
    # Hey1996 for including instrum as W, but added the ion velocity shift to the mean calculation
    # Gradic2022 including velocity addition+-

    tdop = (x0*np.sqrt(kboltz*Temp_in/(amu*m_prot*c_light*c_light)))**2
    W = np.sqrt((pixel_size)**2 + 8*np.log(2)*tdop)
    if bidir == True:
        try:  # First try having v_ion be an array, so that the counter-streaming velocities can be different.
            mean1 = x0*(1 + v_ion[0]/c_light)
            mean2 = x0*(1 - v_ion[1]/c_light)
            xx1 = (x-mean1)/W
            xx2 = (x-mean2)/W

            G1 = (2*np.sqrt(np.log(2))/(np.sqrt(np.pi)*W))*np.exp(-4*np.log(2)*xx1**2)
            G2 = (2*np.sqrt(np.log(2))/(np.sqrt(np.pi)*W))*np.exp(-4*np.log(2)*xx2**2)
            G = (G1 + G2)/2
        except:  # If not an array, just do plus and minus the same velocity
            mean1 = x0*(1 + v_ion/c_light)
            mean2 = x0*(1 - v_ion/c_light)
            xx1 = (x-mean1)/W
            xx2 = (x-mean2)/W

            G1 = (2*np.sqrt(np.log(2))/(np.sqrt(np.pi)*W))*np.exp(-4*np.log(2)*xx1**2)
            G2 = (2*np.sqrt(np.log(2))/(np.sqrt(np.pi)*W))*np.exp(-4*np.log(2)*xx2**2)
            G = (G1 + G2)/2
    else:  # If bidir is False, do normal doppler shifted.
        mean = x0*(1 + v_ion/c_light)
        xx = (x-mean)/W
        G = (2*np.sqrt(np.log(2))/(np.sqrt(np.pi)*W))*np.exp(-4*np.log(2)*xx**2)
    return G


def Voigt(x, x0, Temp_in, amu, pixel_size=1, v_ion=0, bidir=False):
    tdop = (x0*np.sqrt(kboltz*Temp_in/(amu*m_prot*c_light*c_light)))

    # Calculate a Voigt profile using the dopler broadening as the sigma and pixel size as the gamma for the lorrentzian
    if bidir == True:
        try:  # First try having v_ion be an array, so that the counter-streaming velocities can be different.
            mean1 = x0*(1 + v_ion[0]/c_light)
            mean2 = x0*(1 - v_ion[1]/c_light)
        except:  # If not an array, just do plus and minus the same velocity
            mean1 = x0*(1 + v_ion/c_light)
            mean2 = x0*(1 - v_ion/c_light)
        V1 = scp.special.voigt_profile(x-mean1, tdop, pixel_size)
        V2 = scp.special.voigt_profile(x-mean2, tdop, pixel_size)
        V = (V1 + V2)/2
    else:  # If bidir is False, do normal doppler shifted.
        mean = x0*(1 + v_ion/c_light)
        V = scp.special.voigt_profile(x-mean, tdop, pixel_size)
    return Normalize(V)


def Lorrentzian(xrange, peak, pixel_size=0.002222, Skewness=[1, 1], v_ion=0, bidir=False):
    '''
    Calculates a Lorrentzian function (Can be skewed/stepped with uneven left and right sides)

    Parameters
    ----------
    xrange : Array
        Range over which to calculate the  lorrenztian
    peak : TYPE
        DESCRIPTION.
    pixel_size : TYPE, optional
        DESCRIPTION. The default is 0.002222.
    Skewness : TYPE, List [Left,Right]
        Controls how skewed the left and right tapers are, respectively. Based on the drop off ratio and pixel_size. The default is [1,1].

    Returns
    -------
    sl_func : TYPE
        DESCRIPTION.

    '''
    # The Skewness is the left and right constants [left,right], and change the skewness of the lorrentzian. If both are 1, it is a normal lorrentzian. These values can be modified to fit the instrument profile.

    l_w = Skewness[0]*pixel_size
    r_w = Skewness[1]*pixel_size

    if bidir == True:
        try:  # First try having v_ion be an array, so that the counter-streaming velocities can be different.
            mean1 = peak*(1 + v_ion[0]/c_light)
            mean2 = peak*(1 - v_ion[1]/c_light)
            SL1 = (((l_w)**2 / (((xrange-mean1)**2 + (l_w)**2))) * np.heaviside(mean1-xrange, 1) + ((r_w)**2 / (((xrange-mean1)**2 + (r_w)**2))) * np.heaviside(xrange-mean1, 1))
            SL2 = (((l_w)**2 / (((xrange-mean2)**2 + (l_w)**2))) * np.heaviside(mean2-xrange, 1) + ((r_w)**2 / (((xrange-mean2)**2 + (r_w)**2))) * np.heaviside(xrange-mean2, 1))
            SL = (SL1+SL2)/2
        except:  # If not an array, just do plus and minus the same velocity
            mean1 = peak*(1 + v_ion/c_light)
            mean2 = peak*(1 - v_ion/c_light)
            SL1 = (((l_w)**2 / (((xrange-mean1)**2 + (l_w)**2))) * np.heaviside(mean1-xrange, 1) + ((r_w)**2 / (((xrange-mean1)**2 + (r_w)**2))) * np.heaviside(xrange-mean1, 1))
            SL2 = (((l_w)**2 / (((xrange-mean2)**2 + (l_w)**2))) * np.heaviside(mean2-xrange, 1) + ((r_w)**2 / (((xrange-mean2)**2 + (r_w)**2))) * np.heaviside(xrange-mean2, 1))
            SL = (SL1+SL2)/2
    else:  # If bidir is False, do normal doppler shifted.
        mean = peak*(1 + v_ion/c_light)  # Doppler shift to peak
        SL = (((l_w)**2 / (((xrange-mean)**2 + (l_w)**2))) * np.heaviside(mean-xrange, 1) + ((r_w)**2 / (((xrange-mean)**2 + (r_w)**2))) * np.heaviside(xrange-mean, 1))
    return SL



def Normalize(list_in, scaling=1.0):
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
    return np.array([scaling*x/np.max(list_in) for x in list_in])

# %%Utils


def make_Jmj(J_am):
    # This function takes an input general angular momentum J and creates a list of [J,mJ] for all corresponding mJ's.
    # J_am can be an array itself and the output will be a flat list of all possible J,mJ values, starting with the lowest [J,mj] and
    # cycling through the mJ's for each J.
    # =============================================================================
    #
    #     'Starting with the lowest J is likely to cause some issues and there should, instead, be a generalized way for the user to input the J-values ordering'
    #     'The issue would be because the eigval function sorts in ascending order of energies, so I should do likewise, as that may correct for potential issues'
    #     'with energy ordering and states as tracking past diagonalization is not possible.'
    # =============================================================================
    ###
    if type(J_am) == int or type(J_am) == float:
        Jtmp = [[float(J_am), float(x)] for x in np.linspace(-J_am, J_am, int(2*J_am+1))]
    else:
        Jtmp = [[float(jvs), float(x)] for jvs in J_am for x in np.linspace(-jvs, jvs, int(2*jvs+1))]
    return Jtmp


def polarization(mupper, mlower, Btheta, gamma=0, Pol_Filter=False, T_Alpha=1, T_Beta=0):
    '''
    This polarization calculation is generally based on Isler1997 http://dx.doi.org/10.1063/1.872095
    A correction is made with respect to a negative sign in the definition of epsilon plus/minus to be more in line with typical notation.
    This correction removes the issue of an extra negative sign for epsilon_+1-1 when calculated out and the author having the magnitude only.
    There is also a correction in the calculation to properly account for the relative pi/sigma intensity. See documentation for full details.
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

    if Pol_Filter == False:
        # A toggle is needed here because an aligned polarized filter is not the same as not having a filter.
        pol1 = 1
        pol2 = 1
        ealph = 1
        ebeta = 1
    else:
        pol1 = T_Alpha
        pol2 = T_Beta
        gammaangle = np.deg2rad(gamma)
        ealph = pol1*np.cos(gammaangle) + pol2*np.sin(gammaangle)
        ebeta = pol1*np.sin(gammaangle) + pol2*np.cos(gammaangle)
    theta_rad = np.deg2rad(Btheta)  # convert deg to radians

    epsz = ealph*np.sin(theta_rad)
    e0 = epsz

    # eps0 is for a change in am of 0 (pi)
    eps0 = np.sqrt(e0*e0)

    # eps1 is the for a change in am of 1 (sigma)
    eps1 = 0.5*np.sqrt((np.cos(theta_rad)**2) * (ealph**2) + (ebeta**2))  # While this looks like Isler appendix, we haven't squared it for the intensity calc yet.
    # eps1 = np.sqrt(0.5*(np.cos(theta_rad)**2) * (ealph**2) + (ebeta**2 ))

    # Calculate and check the change in magnetic angular momentum. circularly (sigma) light has a change of +/-1, linear (pi) is a change of 0.
    delta_am = np.array(mupper) - np.array(mlower)
    rad_tmp = [0, 3]
    # Using match case instead of a bunch of if statements.
    match delta_am:
        case 1:
            rad_tmp = [eps1, 1]
        case -1:
            rad_tmp = [eps1, -1]
        case 0:
            rad_tmp = [eps0, 0]
        case _:
            rad_tmp = [0, 3]
    return rad_tmp


def CGCoeff(j1, m1, j2, m2, j3, m3):
    # Calculates Clebsch-Gordan Coefficient and outputs a float instead of sympy result
    thisCG = float(CG(j1, m1, j2, m2, j3, m3).doit())
    return thisCG


def Zeemanlow(gf, mf, B):
    # Lowfield Zeeman approximation. Valid for J or F total AM.
    return (muBcm*gf*mf*B)


def ZeemanHigh(ms, ml, B):
    # Highfield Zeeman approximation, excluding Hyperfine splitting. Energy is still diagonal.
    return muBcm*B*(2*ms + ml)


def ZeemanhighHFS(state_in, Ahfs, g_I, Bmag):
    # This is the high field approximation for zeeman splitting (with HFS). E_high = A*m_i*m_j + uB(gjmj + gimi)B
    # For high field, can't use J, and must use s and L. Keep limited to diagonal though.
    EPB_temp = Ahfs*state_in[3]*state_in[5] + mu_B*(state_in[5] + 2*state_in[7] + g_I*state_in[5])*Bmag
    return EPB_temp


def CGw3(j1, m1, j2, m2, j3, m3):
    # Clebsh Gordon Coeff using wigner3j
    return ((-1)**(-j1 + j2 - m3))*np.sqrt(2*j3 + 1)*float(w3j(float(j1), float(j2), float(j3), float(m1), float(m2), -float(m3)))


# %%Zeeman
def dipolestr(LG, LE):
    '''
    Wigner3j = (J'    J    1  )
               (-m_j' m_j  q  )   
    Wigner6j = {S    L'  J'  }
               {1    J    L  }                   

    Parameters
    ----------
    LG, LE are arrays of quantum numbers. Formatted as LG=[J,mJ,L,mL,S,mS]

    Returns
    -------
    tdp3 : float
        Returns the relative dipole strength factor as a multiple of the reduced dipole matrix. Metcalf 4.32

    '''
    thisq = LE[1] - LG[1]
    # if thisq == -1 or thisq== 0 or thisq== 1:

    tw3 = float(w3j(LG[0], 1, LE[0], LG[1], thisq, -LE[1]))
    tdp3 = ((-1)**(LE[2] + LG[4] - LE[1])) * np.sqrt((2*LE[0]+1)*(2*LG[0]+1))*tw3*float(w6j(LE[2], LE[0], LG[4], LG[0], LG[2], 1))

    # else:
    #     tdp3 = 0
    return tdp3


def stateprocess(L_state, S_state, E_level, thisB, sortE=False):
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

    # Create an array that contains [L,ml,S,ms] for the given L,S values. Starts with highest ml,mS values and itereates through ms first.
    LS0 = [[float(L_state), float(x), float(S_state), float(y)] for x in np.linspace(L_state, -L_state, int(2*L_state+1), dtype=float) for y in np.linspace(S_state, -S_state, int(2*S_state+1), dtype=float)]
    # LS0 = [[L_state, x, S_state, y] for x in np.linspace(L_state, -L_state, int(2*L_state+1),dtype=float) for y in np.linspace(S_state, -S_state, int(2*S_state+1),dtype=float) ]
    # print(LS0)
    # Now create the possible J values
    Jvals = np.arange(np.abs(L_state-S_state), (L_state+S_state+.1))  # ,dtype=float)

    # Use those J values to make all possible mJ values and starts with Lowest J, and lowest mJ
    Jmjs0 = make_Jmj(Jvals)

    # Combine both lists together to be a single large [J,mJ,L,mL,S,mS] list. The dimensions should match for LS0 and Jmjs0.
    bigstatelist = [[x[0], x[1], y[0], y[1], y[2], y[3]] for x, y in zip(Jmjs0, LS0)]

    # Expand the input energies which were given for each J to correspond to the number of mJs
    E_cm0 = [[E_level[x]]*int(2*Jvals[x] + 1) for x in range(len(E_level))]
    # Flatten the results out
    E_cm0_flat = [item for row in E_cm0 for item in row]
    E_cm0_Mat = np.diagflat(E_cm0_flat)
    # This step uses a low field approximation in order to sort the state-list based on the given magnetic field such that the relative state function order is known.
    # This is used for signal strength calculation as otherwise we may not know what order the functions are after they get diagonalized. This step might not be needed.
    LowE = []

    if sortE == True:
        biglist = []
        for i, x in enumerate(E_cm0_flat):
            thislandegJ = lande_gj(bigstatelist[i][0], bigstatelist[i][2], bigstatelist[i][4])  # caluclate lande gj factor for a given state
            if thisB == 0:
                thisZeemanEnergy = x  # This is to prevent issues with attempting to sort for B=0
            else:
                thisZeemanEnergy = Zeemanlow(thislandegJ, bigstatelist[i][1], thisB) + x  # This is the low zeeman approximation total energy
            biglist.append([thisZeemanEnergy, x, Jmjs0[i], LS0[i], bigstatelist[i]])

        sortedlist = sorted(biglist)  # Soted in ascending energy
        LowE = [x[0] for x in sortedlist]
        sortedstates = [x[-1] for x in sortedlist]
        E_cm_flat = [x[1] for x in sortedlist]  # This is the sorted energies, and is in cm^-1
        # Jmjs = [x[2] for x in sortedlist]
        # LS = [x[3] for x in sortedlist]
        E_SO_cm_Mat = np.diagflat([item for item in E_cm_flat])  # Despite the name, these energies are actually in cm-1

    # We now use wigner 3j to get an equivalent decoupling as per how CG coefficients work.
    # The result is a 2d where the columns are the |J,mJ> states and the rows determine what the corresponding |L,mL,S,mS> state was.
    # Eventually would like a way to convert this data to a more human readable form
        dim = len(Jmjs0)
        w3mat = np.zeros((dim, dim))
        # [J,mJ,L,mL,S,mS]
        for i, states1 in enumerate(sortedstates):

            for j, states2 in enumerate(sortedstates):
                # This conditional is saying that J != J', as otherwise the different J levels will cause duplicate CG's.
                if states1[2] != states2[2]:
                    CG1 = 0
                else:
                    # The CG coefficient couples J,mJ to L,mL, S mS
                    CG1 = CGw3(states2[2], states2[3], states2[4], states2[5], states1[0], states1[1])
                w3mat[i, j] = CG1
        return [w3mat, E_SO_cm_Mat, LowE, E_cm_flat, sortedstates]

    else:
        dim = len(Jmjs0)
        w3mat = np.zeros((dim, dim))
        # [J,mJ,L,mL,S,mS]
        for i, states1 in enumerate(bigstatelist):

            for j, states2 in enumerate(bigstatelist):
                # This conditional is saying that J != J', as otherwise the different J levels will cause duplicate CG's.
                if states1[2] != states2[2]:
                    CG1 = 0
                else:
                    # The CG coefficient couples J,mJ to L,mL, S mS
                    CG1 = CGw3(states2[2], states2[3], states2[4], states2[5], states1[0], states1[1])
                w3mat[i, j] = CG1
        return [w3mat, E_cm0_Mat, LowE, E_cm0_flat, bigstatelist]

def Zeeman_signal(LevelG, LZG, LevelE, LZE, Bangle=90, gamma=0, Filter=False):

    # Need to make a list of all possible combinations such that each upper state is repeated by the dim of the lower level, and a second list where the opposite is true.
    # ie, each lower level is repeated by the count of upper level states. This makes two lower x upper length arrays that will be useful for wigner calcs.
    bjlist = [[x, y] for y in LevelE[-1] for x in LevelG[-1]]  # This list will be [[bigstatesG], [bigstatesE]]

    # need to define the polarization effect from the filter.
    rad = []
    for i in range(len(bjlist)):


        rad.append(polarization(bjlist[i][1][1], bjlist[i][0][1], Bangle, gamma, Pol_Filter=Filter))

    # Calculate the signal strength coefficients as the product of radiation polarization term and dipole strength.
    # This is essentially a row vector.
    Signalcoeff = []
    for i in range(len(bjlist)):
        if rad[i][1] == 3:
            Signalcoeff.append(0)
        else:
            Signalcoeff.append(rad[i][0]*dipolestr(bjlist[i][1], bjlist[i][0]))
    Signalcoeffa = np.array(Signalcoeff)

    afcf = np.kron(np.array(LZE[1]), np.array(LZG[1]))  # This uses the tensor product after the magnetic field is taken into account.

    signal_out = []
    # for i in range(len(bjlist)): #This removes any components who would have delta mj greater than 1, though it makes the spectra match worse....
    #     if rad[i][1] ==3:
    #         signal_out.append(0)
    #     else:
    #         signal_out.append((np.sum(afcf[:,i]*Signalcoeffa)**2))

    for i in range(len(bjlist)):
        signal_out.append((np.sum(afcf[:, i]*Signalcoeffa)**2))
    tsigout = np.array(signal_out)
    sigsum = np.sum(tsigout)

    sigout_raw = tsigout
    # sigout=sigout_raw
    sigout = np.array([x/sigsum for x in tsigout])  # Make it so that sum of intensities is 1.
    # Conver the output to wavelength in nm then account for shift from vacuum to air.
    wave_vac = np.array([1e7/(y-x) for y in LZE[0] for x in LZG[0]])
    wave_air = Vac_to_air(wave_vac)

    return [sigout, wave_air, rad, bjlist, afcf, Signalcoeffa, wave_vac]

def Zeeman_sig_norad(LevelG, LZG, LevelE, LZE):

    # Need to make a list of all possible combinations such that each upper state is repeated by the dim of the lower level, and a second list where the opposite is true.
    # ie, each lower level is repeated by the count of upper level states. This makes two lower x upper length arrays that will be useful for wigner calcs.
    bjlist = [[x, y] for y in LevelE[-1] for x in LevelG[-1]]  # This list will be [[bigstatesG], [bigstatesE]]

    # need to define the polarization effect from the filter.

    # Calculate the signal strength coefficients as the product of radiation polarization term and dipole strength.
    # This is essentially a row vector.
    Signalcoeff = []
    for i in range(len(bjlist)):
        Signalcoeff.append(dipolestr(bjlist[i][0], bjlist[i][1]))
    Signalcoeffa = np.array(Signalcoeff)

    afcf = np.kron(np.array(LZE[1]), np.array(LZG[1]))  # This uses the tensor product after the magnetic field is taken into account.

    # Calculate the signal which is the sum of the square of the product of afcf, which si reduced dipole signal (including field effects) and the signal coeff. Isler1997
    signal_out = []
    for i in range(len(bjlist)):
        signal_out.append((np.sum(afcf[:, i]*Signalcoeffa)**2))
    tsigout = np.array(signal_out)
    sigsum = np.sum(tsigout)
    sigout = np.array([x/sigsum for x in tsigout])  # Make it so taht sum of intensities is 1.

    return sigout


def Zeeman_func(Level_in, Binput, Lowfield=False, Highfield=False, Eterm=None):
    '''
    Calculates the Zeeman shift hamiltonian and combines it with the Eo energy before diagonalizing to get the eig vec/vals.

    Parameters
    ----------
    Level_in : TYPE
        DESCRIPTION.
    Binput : TYPE
        DESCRIPTION.

    Returns
    -------
    list
        DESCRIPTION.

    '''
    # Level=[CGmat, Eo_mat,E0_flat,bigstatelist]
    Ldim = len(Level_in[-1])  # Dimension of the atomic level: Length of the # of possible states.
    # Spi is wavefunction basis matrix. Comprised of CG coefficients.
    Psi = Level_in[0]

    PsiL = np.zeros_like(Psi)  # Matrix for storing mL contribution of H_Zeeman
    PsiS = np.zeros_like(Psi)  # Matrix for storing mS contribution of H_Zeeman

    # bigstatelist [J,mJ,L,mL,S,mS]
    # Create diagonal elements for the m values for each qm. Doing this is forcing that <n|m> = delta_nm, so 1 if n=m, but 0 else.
    tmLMat = np.diagflat([x[3] for x in Level_in[-1]])
    tmSMat = np.diagflat([x[5] for x in Level_in[-1]])

    # This is equivalent to the <n|Jz|m> calculation to make the pertubation matrix of zeeman energies
    # The += doesn't seem to change anything, but is kept as my original code used it
    for i in range(Ldim):
        tPsiL = Psi[i]
        for j in range(Ldim):
            PsiL[i, j] = np.linalg.multi_dot([tPsiL, tmLMat, Psi[j]])
            PsiS[i, j] = np.linalg.multi_dot([tPsiL, tmSMat, Psi[j]])

    # Combine the energy matricies for mL and mS

    # Diagonalize and calculate the eigenvalues and eigenvectors
    # Check for High and Low field approximations
    if Lowfield == False and Highfield == False:  # If both flags False - Do intermediate regime
        H_Zeeman = g_S*PsiS + g_L*PsiL

        scaledMat = np.asmatrix(Level_in[1] + Binput*H_Zeeman*(muBcm))

    elif Highfield == True and not Lowfield:  # Check if only High field
        H_ZeemanH = np.array([ZeemanHigh(x[5], x[3], Binput) for x in Level_in[-1]])
        EtermMat = Eterm*np.eye(Ldim)
        if Eterm:
            scaledMat = np.asmatrix(EtermMat + H_ZeemanH)
        else:
            scaledMat = np.asmatrix(Level_in[1] + H_ZeemanH)

    elif Lowfield == True and not Highfield:  # Check if only Low field
        gjmat = np.array([lande_gj(x[0], x[2], x[4]) for x in Level_in[-1]])
        tmJMat = np.diagflat([x[1] for x in Level_in[-1]])

        H_ZeemanL = gjmat*tmJMat

        scaledMat = np.asmatrix(Level_in[1] + Binput*H_ZeemanL*(muBcm))
    else:
        print('Something went wrong in Zeeman_func flag checking')
        # scaledMat = np.asmatrix(Level_in[1] )

    eigZ, eigvecZ = np.linalg.eigh(scaledMat)  # , driver = "evx") #Using driver="evx" seems to change some of the signs of the resulting matricies. In the 0-field case, the eigenvectors should all be +1 I believe.

    return [eigZ, eigvecZ]


def HZeeman(Level_in):
    '''
    Calculates the Zeeman Energy from the Hamiltonian using perturbation theory of <psi' |H'|psi> for H' = H_Zeeman for each m (ml, ms)
    This is a scaled value, equivalent to B=1T, and so scalar multpilication of our desired magnetic field will yeild the correct Energy.



    Parameters
    ----------
    Level_in : TYPE
        DESCRIPTION.
    Binput : TYPE
        DESCRIPTION.

    Returns
    -------
    list
        DESCRIPTION.

    '''
    # Level=[CGmat, Eo_mat,E0_flat,bigstatelist]
    Ldim = len(Level_in[-1])  # Dimension of the atomic level: Length of the # of possible states.
    # Spi is wavefunction basis matrix. COmprised of CG coefficients.
    Psi = Level_in[0]

    PsiR = Level_in[0].T

    PsiL = np.zeros_like(Psi)
    PsiS = np.zeros_like(Psi)

    # bigstatelist [J,mJ,L,mL,S,mS]
    # Create diagonal elements for the m values for each qm. Doing this is forcing that <n|m> = delta_nm, so 1 if n=m, but 0 else.
    tmLMat = np.diagflat([x[3] for x in Level_in[-1]])
    tmSMat = np.diagflat([x[5] for x in Level_in[-1]])

    # This is equivalent to the <n|Jz|m> calculation to make the pertubation matrix of zeeman energies
    # The += doesn't seem to change anything, but is kept as my original code used it
    # for i in range(Ldim):
    #     for j  in range(Ldim):
    #         PsiL[i,j] = Psi[i]@tmLMat@PsiR[i]
    #         PsiS[i,j] = Psi[i]@tmSMat@PsiR[i]
    # PsiL = np.zeros_like(Psi)
    # PsiS = np.zeros_like(Psi)

    # #bigstatelist [J,mJ,L,mL,S,mS]
    # #Create diagonal elements for the m values for each qm. Doing this is forcing that <n|m> = delta_nm, so 1 if n=m, but 0 else.
    # #The resulting PsiL,PsiS are NOT diagonal due to the J-> LS coupling done. However, the result is accurate for applying
    # # the Lz or Sz operator onto the state.
    # tmLMat = np.diagflat([x[3] for x in Level_in[-1]])
    # tmSMat = np.diagflat([x[5] for x in Level_in[-1]])

    # This is equivalent to the <n|Jz|m> calculation to make the pertubation matrix of zeeman energies
    # The += doesn't seem to change anything, but is kept as my original code used it
    for i in range(Ldim):
        tPsiL = Psi[i]
        for j in range(Ldim):
            PsiL[i, j] = np.linalg.multi_dot([tPsiL, tmLMat, Psi[j]])
            PsiS[i, j] = np.linalg.multi_dot([tPsiL, tmSMat, Psi[j]])
    # #Combine the energy matricies for mL and mS
    H_Zeeman = g_S*PsiS + g_L*PsiL  # Combine all components

    E0 = Level_in[1]

    return [E0, H_Zeeman]


# %%ZeemanHFS

def stateprocess_HFS(L_state, S_state, I_state, HFS_Const, g_I, E_level, thisB, sortE=False):

    # Create an array that contains [L,ml,S,ms] for the given L,S values. Starts with lowest ml,mS values and itereates through ms first.
    # LS0 = [[L_state, x, S_state, y] for x in np.linspace(-L_state, L_state, int(2*L_state+1)) for y in np.linspace(-S_state, S_state, int(2*S_state+1)) ]
    # Create an array that contains [L,ml,S,ms,I,mI] for the given L,S,I values. Starts with Highest mL,mS values and itereates through mI first, then mS, then mL
    # Extra type-casting is required due to new numpy.float64 objects interaction with sympy wigner functions
    LSI0 = [[float(L_state), float(x), float(S_state), float(y), float(I_state), float(z)] for x in np.linspace(-L_state, L_state, int(2*L_state+1)) for y in np.linspace(-S_state, S_state, int(2*S_state+1))for z in np.linspace(I_state, -I_state, int(2*I_state+1))]

    # Now create the V values, starting with the lowest J state first.

    # Make all Js corresponding to L,S relation
    Jvals = np.arange(np.abs(L_state-S_state), (L_state+S_state+.1))

    # Make all mJs for the Js
    Jmjs0 = make_Jmj(Jvals)

    # Make a list of J,mJ,I,mI states
    JI0 = [[float(x[0]), float(x[1]), float(I_state), float(y)] for x in Jmjs0 for y in np.linspace(-I_state, I_state, int(2*I_state+1))]

    # Combine J and Is to make possible Fs
    Fvals = [np.arange(np.abs(-x + I_state), (x + I_state+.1)) for x in Jvals]
    # Flatten the list to be one big list
    Fvals = [float(item) for row in Fvals for item in row]
    # Make the possible mFs for all Fs
    Fmfs0 = make_Jmj(Fvals)

    E_cm0 = [[E_level[x]]*int(2*Jvals[x] + 1)*int(2*I_state + 1) for x in range(len(E_level))]
    E0_cm_flat = [item for row in E_cm0 for item in row]

    HFS_A_In = HFS_Const['HFSA']
    try:
        HFSConst = [[HFS_A_In[x]]*int(2*Jvals[x] + 1)*int(2*I_state + 1) for x in range(len(E_level))]
        HFS_Flat = [item for row in HFSConst for item in row]
    except:
        # This handles the case of single length HFS inputs
        HFSConst = [[HFS_A_In]*int(2*Jvals[x] + 1)*int(2*I_state + 1) for x in range(len(E_level))]
        HFS_Flat = [item for row in HFSConst for item in row]

        # HFS_Flat = HFSConst
    # print(HFSConst)

    if 'HFSB' in HFS_Const:
        try:
            HFSConstB = [list(HFS_Const['HFSB'][x])*int(2*Jvals[x] + 1)*int(2*I_state + 1) for x in range(len(E_level))]
            HFS_FlatB = [item for row in HFSConstB for item in row]
        except:
            HFSConstB = [list(HFS_Const['HFSB'])*int(2*Jvals[x] + 1)*int(2*I_state + 1) for x in range(len(E_level))]
            HFS_FlatB = [item for row in HFSConstB for item in row]
    else:
        HFS_FlatB = np.zeros_like(E0_cm_flat)
    # Duplicate energies based on degeneracy and flatten the results

    # Duplicate HFS constants to match F levels based on degeneracy and flatten the results
    # HFSConst = [[HFS_Const[0][x]]*int(2*Jvals[x] +1)*int(2*I_state + 1) for x in range(len(E_level))]
    # Flatten the shape of the Energy array

    # Make a large list of [F,mF,J,mJ,I,mI,L,mL,S,mS]
    bigstatelist = [[x[0], x[1], y[0], y[1], y[2], y[3], z[0], z[1], z[2], z[3]] for x, y, z, in zip(Fmfs0, JI0, LSI0)]

    # This step uses a low field approximation in order to sort the state-list based on the given magnetic field such that the relative state function order is known.
    # This is used for signal strength calculation as otherwise we may not know what order the functions are after they get diagonalized. This step might not be needed.

    if sortE:
        biglist = []
        for i, x in enumerate(E0_cm_flat):
            # thislandegF = lande_gf(bigstatelist[i][0],I_state, bigstatelist[i][2], L_state,S_state,g_I) #Calculate lande gf factor for the given transition
            thislandegF = lande_gf(bigstatelist[i][0], bigstatelist[i][4], bigstatelist[i][2], bigstatelist[i][6], bigstatelist[i][8], g_I)  # Calculate lande gf factor for the given transition

            # Lowfield zeeman approximation with hfs, which should be mf and gf.
            thisZeemanEnergy = Zeemanlow(thislandegF, bigstatelist[i][1], thisB) + x + Delta_E_HFS(I_state, bigstatelist[i][2], bigstatelist[i][0], HFS_Flat[i], HFS_FlatB[i])

            # thisZeemanEnergy = Zeemanlow(thislandegJ,bigstatelist[i][3],thisB) + muBcm*g_I*thisB*bigstatelist[i][5] + x + Delta_E_HFS(I_state,bigstatelist[i][2],bigstatelist[i][0],HFS_Flat[i],HFS_FlatB[i])
            biglist.append([thisZeemanEnergy, x, bigstatelist[i], HFS_Flat[i], HFS_FlatB[i]])
        # Sort the results and pull sorted values
        sortedlist = sorted(biglist)
        sortedstates = [x[2] for x in sortedlist]
        sortedHFC = [x[3] for x in sortedlist]  # Sorted hyperfine constants
        sortedHFCB = [x[4] for x in sortedlist]  # Sorted hyperfine constants

        E_cm = [x[1] for x in sortedlist]
        E_low = [x[0] for x in sortedlist]

        # Calculate Hyperfine energies with sorted states and values.
        E_HFSFlat = [Delta_E_HFS(states1[4], states1[2], states1[0], sortedHFC[i], sortedHFCB[i]) for i, states1 in enumerate(sortedstates)]  # Hyperfine + E0 (spin orbital, Nist)

        # E_HFSFlat = [E_cm[i] + Delta_E_HFS(states1[4],states1[2],states1[0], sortedHFC[i], sortedHFCB[i])    for i,states1 in enumerate(sortedstates)] #Hyperfine + E0 (spin orbital, Nist)
        # E_HFSA = [Delta_E_HFS(states1[4],states1[2],states1[0], sortedHFC[i],0)    for i,states1 in enumerate(sortedstates)] #Just Hyperfine energy
        # E_HFSB = [Delta_E_HFS(states1[4],states1[2],states1[0], 0, B=sortedHFCB[i])    for i,states1 in enumerate(sortedstates)] #Just Hyperfine energy

        # E0_cm_Mat = np.diagflat([item for item in E_HFSFlat]) #Make the diagonal energy matrix
        E0_cm_Mat = np.diagflat([item for item in E_cm])  # Make the diagonal energy matrix

        # We now use wigner 3j to get an equivalent decoupling as per how CG coefficients work.
        # There are two coefficients as we go from |F,mF> -> |J,mJ>|I,mI> then another set for |J,mJ>|I,mI> -> |L,mL>|S,mS>|I,mI>. Keeping track of constants is slightly tricky.
        # There can be issues with overla
        # Eventually would like a way to convert this data to a more human readable form
        dimhfs = len(Fmfs0)
        w3cg1 = np.zeros((dimhfs, dimhfs))
        w3cg2 = np.zeros((dimhfs, dimhfs))
        # Done in a single set of iterations. The check for states1[2] != states2[2] is saying that J != J' to calculate, otherwise it's zero.
        for i, states1 in enumerate(sortedstates):

            for j, states2 in enumerate(sortedstates):
                # This conditional is saying that J != J', as otherwise the different J levels will cause duplicate CG's.
                if states1[2] != states2[2]:
                    CG1 = 0
                else:
                    # Make a large list of [F,mF,J,mJ,I,mI,L,mL,S,mS]

                    # The CG coefficient couples F,mF to J,mJ;I,mI for CG1
                    CG1 = CGw3(states2[2], states2[3], states2[4], states2[5], states1[0], states1[1])
                # This conditional ensures we don't allow nuclear spin to flip or change between states. states[5] is mI.
                if states1[5] != states2[5]:
                    CG2 = 0
                    # CG2 = CGw3(states2[6],states2[7],states2[8],states2[9],states1[2],states1[3])
                    # w3cg2[i,j] = CG2

                else:
                    CG2 = CGw3(states2[6], states2[7], states2[8], states2[9], states1[2], states1[3])

                    w3cg2[i, j] = CG2
                w3cg1[i, j] = CG1

        w3cgprod = np.matmul(w3cg1, w3cg2)

        return [w3cgprod, E0_cm_Mat, E_HFSFlat, E_low, sortedstates]
    else:

        # E_cm = [x[1] for x in biglist]

        # Calculate Hyperfine energies with sorted states and values.
        E_HFSFlat = [Delta_E_HFS(states1[4], states1[2], states1[0], HFS_Flat[i], HFS_FlatB[i]) for i, states1 in enumerate(bigstatelist)]  # Hyperfine + E0 (spin orbital, Nist)

        # E0_HFSFlat = [E_cm[i] + Delta_E_HFS(states1[4],states1[2],states1[0], HFS_Flat[i], HFS_FlatB[i])    for i,states1 in enumerate(bigstatelist)] #Hyperfine + E0 (spin orbital, Nist)
        # E_HFSA = [Delta_E_HFS(states1[4],states1[2],states1[0], sortedHFC[i],0)    for i,states1 in enumerate(sortedstates)] #Just Hyperfine energy
        # E_HFSB = [Delta_E_HFS(states1[4],states1[2],states1[0], 0, B=sortedHFCB[i])    for i,states1 in enumerate(sortedstates)] #Just Hyperfine energy

        E0_cm_Mat = np.diagflat([item for item in E0_cm_flat])  # Make the diagonal energy matrix

        # We now use wigner 3j to get an equivalent decoupling as per how CG coefficients work.
        # There are two coefficients as we go from |F,mF> -> |J,mJ>|I,mI> then another set for |J,mJ>|I,mI> -> |L,mL>|S,mS>|I,mI>. Keeping track of constants is slightly tricky.
        # There can be issues with overla
        # Eventually would like a way to convert this data to a more human readable form
        dimhfs = len(Fmfs0)
        w3cg1 = np.zeros((dimhfs, dimhfs))
        w3cg2 = np.zeros((dimhfs, dimhfs))
        # Done in a single set of iterations. The check for states1[2] != states2[2] is saying that J != J' to calculate, otherwise it's zero.
        for i, states1 in enumerate(bigstatelist):

            for j, states2 in enumerate(bigstatelist):
                # This conditional is saying that J != J', as otherwise the different J levels will cause duplicate CG's.
                if states1[2] != states2[2]:
                    CG1 = 0

                else:
                    # Make a large list of [F,mF,J,mJ,I,mI,L,mL,S,mS]

                    # The CG coefficient couples F,mF to J,mJ;I,mI for CG1
                    CG1 = CGw3(states2[2], states2[3], states2[4], states2[5], states1[0], states1[1])
                # This conditional ensures we don't allow nuclear spin to flip or change between states. states[5] is mI.
                if states1[5] != states2[5]:
                    CG2 = 0
                    # CG2 = CGw3(states2[6],states2[7],states2[8],states2[9],states1[2],states1[3])
                    # w3cg2[i,j] = CG2

                else:
                    CG2 = CGw3(states2[6], states2[7], states2[8], states2[9], states1[2], states1[3])

                    w3cg2[i, j] = CG2
                w3cg1[i, j] = CG1

        w3cgprod = np.matmul(w3cg1, w3cg2)
        return [w3cgprod, E0_cm_Mat, E_HFSFlat, HFS_Flat, bigstatelist]


def dipolestr_HFS(LG, LE):
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
    # States are from bigsortedlist = [F,mF,J,mJ,I,mI,L,mL,S,mS]
    # [0,1 ,2,3 ,4,5 ,6,7 ,8,9 ]

    thisq = LE[1] - LG[1]

    # This is Metcalf 4.33
    # Calculate the wigner 3j and two wigner 6j.
    tw3 = float(w3j(LG[0], 1, LE[0], LG[1], thisq, -LE[1]))
    tw61 = float(w6j(LE[6], LE[2], LG[8], LG[2], LG[6], 1))
    tw62 = float(w6j(LE[2], LE[0], LG[4], LG[0], LG[2], 1))

    tphase = (-1)**(1 + LE[6] + LG[8] + LG[2] + LE[2] + LG[4] - LE[1])
    troot = np.sqrt((2*LG[0] + 1)*(2*LE[0]+1) * (2*LG[2] + 1)*(2*LE[2] + 1))

    tdp3 = tphase * troot*tw3*tw61*tw62

    return tdp3


def Zeeman_signal_HFS(LevelG, LZG, LevelE, LZE, Bangle=90, gamma=0, Filter=False):

    # Need to make a list of all possible combinations such that each upper state is repeated by the dim of the lower level, and a second list where the opposite is true.
    # ie, each lower level is repeated by the count of upper level states. This makes two lower x upper length arrays that will be useful for wigner calcs.
    bflist = [[x, y] for y in LevelE[-1] for x in LevelG[-1]]  # This list will be [[bigstatesG], [bigstatesE]]

    # need to define the polarization effect from the filter.
    rad = []
    for i in range(len(bflist)):
        rad.append(polarization(bflist[i][1][1], bflist[i][0][1], Bangle, gamma, Pol_Filter=Filter))
        # Rad output is [polfactor,'polarization'] with polarization being 'sigma-','pi','sigma+' accordingly
    # Calculate signal strength coefficient as product of dipole strength and polarization component
    Signalcoeff = []
    for i in range(len(bflist)):
        Signalcoeff.append(rad[i][0]*dipolestr_HFS(bflist[i][0], bflist[i][1]))
    Signalcoeffa = np.array(Signalcoeff)

    afcf = np.kron(np.array(LZE[1]), np.array(LZG[1]))  # This uses the tensor product after the magnetic field is taken into account.

    # Calculate the signal. This is like the dot product of a column of the afcf reduced dipole state matrix with the row vector comprised of dipole strengths(with polarization)
    signal_out = []
    for i in range(len(bflist)):
        signal_out.append((np.sum(afcf[:, i]*Signalcoeffa)**2))
    tsigout = np.array(signal_out)
    sigsum = np.sum(tsigout)
    sigout = np.array([x/sigsum for x in tsigout])  # Make it so taht sum of intensities is 1.

    # Convert the result from cm^-1 to nm and include shift to air from vacuum
    wave_vac = np.array([1e7/(y-x) for y in LZE[0] for x in LZG[0]])
    wave_air = Vac_to_air(wave_vac)

    return [sigout, wave_air, rad, bflist, afcf, Signalcoeffa, wave_vac]


def Zeeman_func_HFS(Level_in, Binput, gI, PsiHFS=True):
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
    list containing the eigenvalues and eigenvectores after matrix diagonalization
        [eigZ,eigvecZ]


    '''
    Ldim = len(Level_in[-1])  # Dimension of the atomic level: Length of the # of possible states.
    Psi = Level_in[0]
    PsiR = Level_in[0].T
    PsiL = np.zeros_like(Psi)
    PsiS = np.zeros_like(Psi)
    PsiI = np.zeros_like(Psi)
    PsiH = np.zeros_like(Psi)

    # Create diagonal elements for the m values for each qm.
    tmLMat = np.diagflat([x[7] for x in Level_in[-1]])
    tmSMat = np.diagflat([x[9] for x in Level_in[-1]])
    tmIMat = np.diagflat([x[5] for x in Level_in[-1]])

    # H_HFS = AIJcos(IJ) = AIJ*(F(F+1) - I(I+1) - J(J+1))/(IJ)

    # tH_HFS = np.diagflat([x for x in Level_in[3]])

    tH_HFS = np.diagflat([x for x in Level_in[2]])
    # tH_HFS = np.diagflat([x for x in Level_in[3]])

    # This is equivalent to the <n|Jz|m> calculation to make the pertubation matrix of zeeman energies
    # The += doesn't seem to change anything, but is kept as my original code used it
    # HFS should be treated in the same way as Zeeman, with the matrix multiplication of the wave fxns.

    for i in range(Ldim):
        tPsiL = Psi[i]
        for j in range(Ldim):

            PsiL[i, j] = np.linalg.multi_dot([tPsiL, tmLMat, Psi[j]])
            PsiS[i, j] = np.linalg.multi_dot([tPsiL, tmSMat, Psi[j]])
            PsiI[i, j] = np.linalg.multi_dot([tPsiL, tmIMat, Psi[j]])
            # PsiH[i,j] = np.linalg.multi_dot([tPsiL, tH_HFS, Psi[j]])

            # PsiH[i,j] = PsiR[i]@(tH_HFS@PsiR[j])

            PsiH[i, j] = np.linalg.multi_dot([PsiR[i], tH_HFS, PsiR[j]])

    H_Zeeman = g_S*PsiS + g_L*PsiL - gI*PsiI  # Combine all components

    # Combine the energy matrices
    if PsiHFS == 'True':
        # print('PsiH')
        scaledMat = np.asmatrix(Level_in[1] + PsiH + Binput*H_Zeeman*(muBcm))
    else:
        scaledMat = np.asmatrix(Level_in[1] + Level_in[2] + Binput*H_Zeeman*(muBcm))

    # scaledMat = np.asmatrix(Level_in[1] +Level_in[2]+Binput*H_Zeeman*(muBcm) )
    # Diagonalize and calculate eigenvalues and eigenvectors.
    eigZ, eigvecZ = np.linalg.eigh(scaledMat)  # , driver = "evx")
    # return [eigZ,eigvecZ,[tmLMat,tmSMat,tmIMat,g_S*PsiS, g_L*PsiL,  gI*PsiI]]

    return [eigZ, eigvecZ]


def HZeeman_HFS(Level_in, gI):
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
    list containing the eigenvalues and eigenvectores after matrix diagonalization
        [eigZ,eigvecZ]


    '''
    # stateprocess outputs [w3mat, E_0, E_HFS,HHFS,statelist]
    Ldim = len(Level_in[-1])  # Dimension of the atomic level: Length of the # of possible states.
    Psi = Level_in[0]

    PsiR = Level_in[0].T
    PsiL = np.zeros_like(Psi)
    PsiS = np.zeros_like(Psi)
    PsiI = np.zeros_like(Psi)
    PsiH = np.zeros_like(Psi)

    # Create diagonal elements for the m values for each qm.
    tmLMat = np.diagflat([x[7] for x in Level_in[-1]])
    tmSMat = np.diagflat([x[9] for x in Level_in[-1]])
    tmIMat = np.diagflat([x[5] for x in Level_in[-1]])

    # H_HFS = AIJcos(IJ) = AIJ*(F(F+1) - I(I+1) - J(J+1))/(IJ)

    # tH_HFS = np.diagflat([x for x in Level_in[3]])

    # tH_HFS = np.diagflat([2*x for x in Level_in[2]])
    tH_HFS = np.diagflat([x for x in Level_in[2]])

    # This is equivalent to the <n|Jz|m> calculation to make the pertubation matrix of zeeman energies
    # The += doesn't seem to change anything, but is kept as my original code used it
    # HFS should be treated in the same way as Zeeman, with the matrix multiplication of the wave fxns.
    for i in range(Ldim):
        tPsiL = Psi[i]

        for j in range(Ldim):
            PsiL[i, j] = np.linalg.multi_dot([tPsiL, tmLMat, Psi[j]])
            PsiS[i, j] = np.linalg.multi_dot([tPsiL, tmSMat, Psi[j]])
            PsiI[i, j] = np.linalg.multi_dot([tPsiL, tmIMat, Psi[j]])
            # PsiH[i,j] = np.linalg.multi_dot([tPsiL, tH_HFS, Psi[j]])

            # PsiH[i,j] = np.linalg.multi_dot([PsiR[i], tH_HFS, PsiR[j]])

    H_Zeeman = g_S*PsiS + g_L*PsiL - gI*PsiI  # Combine all components
    PsiH = tH_HFS

    E0 = Level_in[1]
    AHFS0 = np.diagflat([x for x in Level_in[3]])

    EHF0 = Level_in[2]

    out = {}
    out['E0'] = E0
    out['H_Z'] = H_Zeeman
    out['PsiHFS'] = PsiH
    out['AHFS0'] = AHFS0
    out['EHF0'] = EHF0
    return out


# %%ZeemanFan

#        stateprocces HFS return [w3cgprod,E0_cm_Mat,E_HFSFlat,E_low,bigstatelist]


def FastZfan(States_in, Bfield, Eterm_=0, g_I=0, A_HFS=0):
    #

    # Make a large list of [F,mF,J,mJ,I,mI,L,mL,S,mS] HFS or [J,mJ,L,mL,S,mS] no HFS
    ZfansLow = []
    ZfansHigh = []
    for i, states in enumerate(States_in[-1]):
        if g_I == 0:
            thislandegF = lande_gj(states[0], states[2], states[4])
            thisZeemanEnergyL = Zeemanlow(thislandegF, states[1], Bfield) + States_in[1][i][i]
            thisZeemanEnergyH = ZeemanHigh(states[5], states[3], Bfield) + Eterm_

        else:

            thislandegF = lande_gf(states[0], states[4], states[2], states[6], states[8], g_I)  # Calculate lande gf factor for the given transition
            thisZeemanEnergyL = Zeemanlow(thislandegF, states[1], Bfield) + States_in[2][i] + States_in[1][i][i]
            thisZeemanEnergyH = ZeemanhighHFS(states, A_HFS, g_I, Bfield) + States_in[2][i] + States_in[1][i][i]

        ZfansLow.append(thisZeemanEnergyL)
        ZfansHigh.append(thisZeemanEnergyH)

    # print(Zfans[-1])
    return [ZfansLow, ZfansHigh]


def ZeemanFan(Inputdeck, flags=''):
    # Take the Inputdeck, generate states, then calculate energy for each value. Transpose and put it out.
    # Take the required values from the input deck.
    '''
    flags : String
        Can be any combination of G,E,L,H. The default is ''.
        G: Fan plot Ground
        E: Fan Plot of Excited
        L: Include Lowfield Approx
        H: Include Highfield Approx (Only with HFS)
    '''

    S_ground = Inputdeck['s_ground']
    S_excited = Inputdeck['s_excited']
    L_ground = Inputdeck['l_ground']
    L_excited = Inputdeck['l_excited']
    E_ground = Inputdeck['E_ground']
    E_excited = Inputdeck['E_excited']
    Bfield_ext = Inputdeck['Bmag']
    # Bext_angle = Inputdeck['b_angle']
    try:
        EtermE = Inputdeck['EtermE']
        EtermG = Inputdeck['EtermG']
    except:
        EtermE = None
        EtermG = None
    # Generate empty lists for storing values in later. Ground, excited, low and exact. Look to add high approximation too.
    GZ_fan = []
    EZ_fan = []

    # Check if HFS is being included.
    if 'I_spin' in Inputdeck:
        HFSConstsG = {}
        HFSConstsG['HFSA'] = Inputdeck['HFS_G']
        HFSConstsE = {}
        HFSConstsE['HFSA'] = Inputdeck['HFS_E']

        if 'HFSB_G' in Inputdeck:  # Check if either of the B HFS constants are in the input deck, otherwise make a copy of the A coefficients as zeros.
            HFSConstsG['HFSB'] = Inputdeck['HFSB_G']
        else:
            HFSConstsG['HFSB'] = np.zeros_like(E_ground)

        if 'HFSB_E' in Inputdeck:  # Check if either of the B HFS constants are in the input deck, otherwise make a copy of the A coefficients as zeros.
            HFSConstsE['HFSB'] = Inputdeck['HFSB_E']
        else:
            HFSConstsE['HFSB'] = np.zeros_like(E_excited)

        I_spin = Inputdeck['I_spin']
        mu_I = Inputdeck['mu_I']
        thisg_I = lande_gi(I_spin, mu_I)

        output = {}
        if 'G' in flags:

            Ground_Level_HFS = stateprocess_HFS(L_ground, S_ground, I_spin, HFSConstsG, thisg_I, E_ground, 0)
            Ground_Zeeman_HFS = HZeeman_HFS(Ground_Level_HFS, thisg_I)

            EigG = []
            for Bvals in Bfield_ext:

                # tmat = np.asmatrix(Ground_Zeeman_HFS[0] +  Ground_Zeeman_HFS[2] +  Bvals*Ground_Zeeman_HFS[1]*muBcm)
                if 'PsiHFS' in Inputdeck:
                    tmat = np.asmatrix(Ground_Zeeman_HFS['E0'] + Ground_Zeeman_HFS['PsiHFS'] + Bvals*Ground_Zeeman_HFS['H_Z']*muBcm)
                else:
                    G_HFS0 = np.diagflat([x for x in Ground_Zeeman_HFS['EHF0']])
                    tmat = np.asmatrix(Ground_Zeeman_HFS['E0'] + G_HFS0 + Bvals*Ground_Zeeman_HFS['H_Z']*muBcm)

                EigG.append(np.linalg.eigh(tmat)[0])
            output['G'] = np.array(EigG).T

            if 'L' in flags:
                LowG = []
                for Bvals in Bfield_ext:
                    LowG.append(FastZfan(Ground_Level_HFS, Bvals, g_I=thisg_I, )[0])

                output['LowG'] = np.array(LowG).T
            if 'H' in flags:
                HighG = []
                for Bvals in Bfield_ext:
                    HighG.append(FastZfan(Ground_Level_HFS, Bvals)[1])

                output['HighG'] = np.array(HighG).T
        if 'E' in flags:
            Excited_Level_HFS = stateprocess_HFS(L_excited, S_excited,  I_spin, HFSConstsE, thisg_I, E_excited, 0)
            Excited_Zeeman_HFS = HZeeman_HFS(Excited_Level_HFS, thisg_I)

            EigE = []
            for Bvals in Bfield_ext:

                if 'PsiHFS' in Inputdeck:
                    tmat = np.asmatrix(Excited_Zeeman_HFS['E0'] + Excited_Zeeman_HFS['PsiHFS'] + Bvals*Excited_Zeeman_HFS['H_Z']*muBcm)
                else:
                    E_HFS0 = np.diagflat([x for x in Excited_Zeeman_HFS['EHF0']])
                    tmat = np.asmatrix(Excited_Zeeman_HFS['E0'] + E_HFS0 + Bvals*Excited_Zeeman_HFS['H_Z']*muBcm)

                EigE.append(np.linalg.eigh(tmat)[0])
            output['E'] = np.array(EigE).T
            if 'L' in flags:
                LowE = []
                for Bvals in Bfield_ext:
                    LowE.append(FastZfan(Excited_Level_HFS, Bvals, g_I=thisg_I, )[0])

                output['LowE'] = np.array(LowE).T
            if 'H' in flags:
                HighE = []
                for Bvals in Bfield_ext:
                    HighE.append(FastZfan(Excited_Level_HFS, Bvals)[1])

                output['HighE'] = np.array(HighE).T

    else:
        output = {}
        if 'G' in flags:
            Ground_Level = stateprocess(L_ground, S_ground, E_ground, 0)
            Ground_Zeeman = HZeeman(Ground_Level)

            EigG = []
            for Bvals in Bfield_ext:

                tmat = np.asmatrix(Ground_Zeeman[0] + Bvals*Ground_Zeeman[1]*muBcm)

                EigG.append(np.linalg.eigh(tmat)[0])
            output['G'] = np.array(EigG).T

            if 'L' in flags:
                LowG = []
                for Bvals in Bfield_ext:
                    LowG.append(FastZfan(Ground_Level, Bvals)[0])

                output['LowG'] = np.array(LowG).T
            if 'H' in flags:
                HighG = []
                for Bvals in Bfield_ext:
                    HighG.append(FastZfan(Ground_Level, Bvals,Eterm_=EtermG)[1])

                output['HighG'] = np.array(HighG).T
        if 'E' in flags:
            Excited_Level = stateprocess(L_excited, S_excited, E_excited, 0)
            Excited_Zeeman = HZeeman(Excited_Level)

            EigE = []
            for Bvals in Bfield_ext:

                tmat = np.asmatrix(Excited_Zeeman[0] + Bvals*Excited_Zeeman[1]*muBcm)

                EigE.append(np.linalg.eigh(tmat)[0])
            output['E'] = np.array(EigE).T

            if 'L' in flags:
                LowE = []
                for Bvals in Bfield_ext:
                    LowE.append(FastZfan(Excited_Level, Bvals)[0])

                output['LowE'] = np.array(LowE).T
            if 'H' in flags:
                HighE = []
                for Bvals in Bfield_ext:
                    HighE.append(FastZfan(Excited_Level, Bvals,Eterm_=EtermE)[1])

                output['HighE'] = np.array(HighE).T
    return output


def plotZfan(Inputdeck, flags='', savefig=False, figname='placeholder',markercount = 15,
             markersize_ = 6,markersize_H = 6,markersize_L = 6):
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
    flags : String
        Can be any combination of G,E,L,H. The default is ''.
        G: Fan plot Ground
        E: Fan Plot of Excited
        L: Include Lowfield Approx
        H: Include Highfield Approx (Only with HFS)

    figname : String, optional
        DESCRIPTION. The default is 'placeholder'.

    Returns
    -------
    Returns exact and lowfield zeeman fan incase different plotting is required.
    [Gfan,Efan,Gfan_low,Efan_low]

    '''
    # Takes the inputdeck and a range of magnetic fields to calculate a typical zeeman fan.
    # Exact and lowfield are plotted on the same figure if that option is requested.

    Zfantest = ZeemanFan(Inputdeck, flags=flags)  # Calculate
    # Pull energy arrays from Zfancalculation

    Brange = Inputdeck['Bmag']
    if markercount == 0:
        markersize_ = 0
        markersize_H = 0
        markersize_L = 0
    if 'G' in flags:
        # figfg, axsfg = plt.subplots(figsize = (16,8))
        figs, axs = plt.subplots(figsize=(12, 10))
        plt.rcParams.update({'font.size': 18})

        axs.plot(Brange[0], Zfantest['G'][0][0], color='blue',  linewidth=.75, label="Exact", 
                 marker='D', markevery=markercount,markersize=markersize_, alpha=0.8)  # This makes a single point for the purpose of the legend.
        for things in (Zfantest['G']):
            axs.plot(Brange, things, color='blue',  linewidth=.75,
                     marker='D', markevery=markercount,markersize=markersize_, alpha=0.8)  # , label = "Ground")

        if 'L' in flags:
            for thingsL in (Zfantest['LowG']):
                plt.plot(Brange, thingsL, color='maroon', linestyle=':', linewidth=.75,
                         marker='x', markevery=markercount-2,markersize=markersize_L, alpha=0.9)    # , label = "Ground_Low")
            axs.plot(Brange[0], Zfantest['LowG'][0][0], color='maroon', linestyle=':', label="Ground Low",
                     marker='x', markevery=markercount-2,markersize=markersize_L, alpha=0.8)  # This makes a single point for the purpose of the legend.
        if 'H' in flags:
            for thingsH in (Zfantest['HighG']):
                plt.plot(Brange, thingsH, color='green', linestyle='-.', linewidth=.75,
                         marker='o', markevery=markercount+2,markersize=markersize_H, alpha=0.8)   # , label = "Ground_Low")
            axs.plot(Brange[0], Zfantest['HighG'][0][0], color='green', linestyle='-.', label="Ground High",linewidth=.75,
                     marker='o', markevery=markercount+2,markersize=markersize_H, alpha=0.8)  # This makes a single point for the purpose of the legend.
        # plt.ylim(np.min(things),np.max(things))
        axs.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
        plt.tight_layout()
    if 'E' in flags:
        # figfe, axsfe = plt.subplots(figsize = (16,8))

        figs, axs = plt.subplots(figsize=(12, 10))
        plt.rcParams.update({'font.size': 18})

        axs.plot(Brange[0], Zfantest['E'][0][0], color='blue',  linewidth=.75, label="Exact", 
                 marker='D', markevery=markercount,markersize=markersize_, alpha=0.8)  # This makes a single point for the purpose of the legend.
        for things2 in Zfantest['E']:
            axs.plot(Brange, things2, color='blue', linewidth=.75,
                     marker='D', markevery=markercount,markersize=markersize_, alpha=0.8)  # , label = "Ground")
        if 'L' in flags:
            for thingsL in (Zfantest['LowE']):
                axs.plot(Brange, thingsL, color='maroon', linestyle=':', linewidth=.75,
                         marker='x', markevery=markercount-2,markersize=markersize_L, alpha=0.8) # , label = "Ground_Low")
            axs.plot(Brange[0], Zfantest['LowE'][0][0], color='maroon', linestyle=':', label="Excited Low",linewidth=.75,
                     marker='x', markevery=markercount-2,markersize=markersize_L, alpha=0.8) # This makes a single point for the purpose of the legend.
        if 'H' in flags:
            for thingsH in (Zfantest['HighE']):
                axs.plot(Brange, thingsH, color='green', linestyle='-.',linewidth=.75,
                         marker='o', markevery=markercount+2,markersize=markersize_H, alpha=0.8)   # , label = "Ground_Low")
            axs.plot(Brange[0], Zfantest['HighE'][0][0], color='green', linestyle='-.', label="Excited High",linewidth=.75,
                     marker='o', markevery=markercount+2,markersize=markersize_H, alpha=0.8)   # This makes a single point for the purpose of the legend.
        # plt.ylim(np.min(things2),np.max(things2))
        axs.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
        plt.tight_layout()
    # Setting plot parameters
    plt.xlabel('Magnetic Field (T)', weight='semibold')
    plt.ylabel('Energy (cm^-1)', weight='semibold')
    plt.legend(loc='best')
    # plt.title(Inputdeck['plottitle'])
    plt.tight_layout()
    if savefig==True:
        plt.savefig(f'{figname}_Excited.png',dpi=2000) #Save the figure if requested
    return Zfantest


# %%Main


def Zeeman_Main(Inputdeck, **kwargs):
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
                     'Convfxn': 'Gaussian' , #Optional: 'Gaussian', 'GaussianInstrum'', 'Skewed' , 'Lorrentzian'. Future work will allow custom functions or modifying the skewed lorrentzian
                     'Skewness' : [0.62,0.99] #Optional, [Left,Right] constant for the skewness of the lorrentzian function. Between 0 and 1. If both are 1, it is a normal lorrentzian.
                     'Temp' : 300, #Temperature in K. Used for Gaussian convolution
                     'specstep' : 0.002222 , #stepsize of linefunction [nm]
                     'specres' : 10 , #How many points are evaluated per spectrometer pixel window. Higher means smoother curve.
                     'fxnwindow': 0.5 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                     'plottitle': 'Excample Title' , #Title for plotting (optional).                      
                     'spec_window' : [425,435], #Min and max for convolution window (nm)
                     'plot_window' : [425,435], #Min and max for plot window range (nm)
                     'HFS_G' : [505.5e6] , #OPtional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                     'HFS_E': [0, 496.2e6, 440.5e6] , #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                     'I_spin' : 0.5, #Optional: nuclear spin, I. Tabulated by Stone et al.
                     'mu_I' :  .11778476 , #Optional: Nuclear dipole moment as tabulated by Stone.              
                     'HFSB_G' : [505.5e6] , #OPtional: Hyperfine B constant for the lower level, need one entry per J level, lowest J first [Hz]
                     'HFSB_E': [0, 496.2e6, 440.5e6] , #Optional: Hyperfine B constants for upper level, need one entery per J level, lowest J first [Hz]
                     'Calc_HFS_Vals': [False,False], #Optional. Passed regarding wether HFS constants A and B should be approximated for values of 0 above. NOT YET IMPLEMENTED
                     'ion_vel'
                     'PsiHFS'
                     'DoLowSig' : True , #Bool for whether Lowfield approximation and signal are included and plotted
                     'DoHighSig' : True , #Bool for whether Highfield approximation and signal are included and plotted

                     }

    Returns
    -------
    data : Dictionary
        DESCRIPTION.

    '''

    # Define local variables from the inputdeck
    S_ground = Inputdeck['s_ground']
    S_excited = Inputdeck['s_excited']
    L_ground = Inputdeck['l_ground']
    L_excited = Inputdeck['l_excited']
    E_ground = Inputdeck['E_ground']
    E_excited = Inputdeck['E_excited']
    Bfield_ext = Inputdeck['Bmag']
    Bext_angle = Inputdeck['b_angle']

    for k, val in kwargs.items():
        print('**kwargs not yet implemented')

    # Set polarization filter angle to 0 and define no polarization filter.
    Polangle = 0
    Polfilter = False
    if 'Pol_angle' in Inputdeck:  # Check if the inputdeck has a pol_angle and adjust parameters accordingly.
        Polangle = Inputdeck['Pol_angle']
        Polfilter = True

    if 'bidir' in Inputdeck:
        bidirect = Inputdeck['bidir']
    else:
        bidirect = False

    if 'sortE' in Inputdeck:
        sortE = Inputdeck['sortE']
    else:
        sortE = False

    if 'amu' in Inputdeck:
        pass
    else:
        Inputdeck['amu'] = 1  # Setting a default value so it's not required to be input.

    # 1. Establish the "Spectrometer" fallback range
    spec_data = Inputdeck.get('SpectrumData')
    # If data exists, use its bounds; otherwise, use the global default
    data_range = [np.min(spec_data[0]), np.max(spec_data[0])] if spec_data is not None else [1, 1000]
    
    # 2. Resolve windows: Dictionary Value > Spectrometer range > default
    specwindow = Inputdeck.get('spec_window', data_range)
    plotwindow = Inputdeck.get('plot_window', data_range)
    
    
    if 'EtermG' in Inputdeck:
        EtermG = Inputdeck['EtermG']
    else:
        EtermG = None
        
    if 'EtermE' in Inputdeck:
        EtermE = Inputdeck['EtermE']
    else:
        EtermE = None    
        
        
    if 'I_spin' in Inputdeck:  # Check if the input deck has nuclear spin as to calculate HFS splitting and zeeman effect or not.
        HFS_ConstsG = {}
        HFS_ConstsG['HFSA'] = Inputdeck['HFS_G']
        HFS_ConstsE = {}
        HFS_ConstsE['HFSA'] = Inputdeck['HFS_E']
        if 'HFSB_G' in Inputdeck:  # Check if either of the B HFS constants are in the input deck, otherwise make a copy of the A coefficients as zeros.
            HFS_ConstsG['HFSB'] = Inputdeck['HFSB_G']
        else:
            HFS_ConstsG['HFSB'] = np.zeros_like(E_ground)

        if 'HFSB_E' in Inputdeck:  # Check if either of the B HFS constants are in the input deck, otherwise make a copy of the A coefficients as zeros.
            HFS_ConstsE['HFSB'] = Inputdeck['HFSB_E']
        else:
            HFS_ConstsE['HFSB'] = np.zeros_like(E_excited)

        if 'PsiHFS' in Inputdeck:
            PsiHFS = Inputdeck['PsiHFS']
        else:
            PsiHFS = False
            

        # Use spin, spin parity, and nuclear dipole moment to calculate lande g_i nuclear factor. Might be off by a factor of 1836?
        I_spin = Inputdeck['I_spin']
        mu_I = Inputdeck['mu_I']
        thisg_I = lande_gi(I_spin, mu_I)

        # Process the wavefunction matrix and energies. Preliminary sorting is done so that the final state order should be roughly as expected.
        Ground_Level_HFS = stateprocess_HFS(L_ground, S_ground, I_spin, HFS_ConstsG, thisg_I, E_ground, Bfield_ext)
        Excited_Level_HFS = stateprocess_HFS(L_excited, S_excited, I_spin, HFS_ConstsE, thisg_I, E_excited, Bfield_ext)

        # Take the sorted states, the state matrix, and energy matrix and calculate the zeeman splitting (with HFS)
        Ground_Zeeman_HFS = Zeeman_func_HFS(Ground_Level_HFS, Bfield_ext, thisg_I, PsiHFS)
        Excited_Zeeman_HFS = Zeeman_func_HFS(Excited_Level_HFS, Bfield_ext, thisg_I, PsiHFS)
        # Use the eigenvectors and values from the diagonalization of the hamiltonian to calculate the relative signal strength. This function uses a combination of Isler1997 and Metcalf 4.32/4.33
        Z_sig_HFS = Zeeman_signal_HFS(Ground_Level_HFS, Ground_Zeeman_HFS, Excited_Level_HFS, Excited_Zeeman_HFS, Bangle=Bext_angle, gamma=Polangle, Filter=Polfilter)

        # Output a bunch of data
        data = {}
        data['signal'] = Z_sig_HFS[0]
        data['wave_vac'] = Z_sig_HFS[6]
        data['wave_air'] = Z_sig_HFS[1]
        data['HFSG_Energy'] = Ground_Level_HFS[2]
        data['HFSE_Energy'] = Excited_Level_HFS[2]
        data['statelist_Ground'] = Ground_Level_HFS[-1]
        data['statelist_Excited'] = Excited_Level_HFS[-1]
        data['RedDipoleOp'] = Z_sig_HFS[-3]  # This is the tensor/kronecker product of the upper and lower state matrices after the field has been included. afcf in calculations. Might not be what is expected.
        data['w3coeffs_ground'] = Ground_Level_HFS[0]
        data['w3coeffs_excited'] = Excited_Level_HFS[0]
        data['ZeemanSplits_G'] = Ground_Zeeman_HFS[-1]
        data['rad'] = Z_sig_HFS[2]  # Given out so that it can be used to plot specific polarizations.

    # If there is no nuclear spin in the input deck, do a similar calculation but no hpyerfine splitting ins included nor nuclear zeeman splitting.
    else:
        # Process the wavefunction matrix and energies. Preliminary sorting is done so that the final state order should be roughly as expected.
        Ground_Level = stateprocess(L_ground, S_ground, E_ground, Bfield_ext)
        Excited_Level = stateprocess(L_excited, S_excited, E_excited, Bfield_ext)
        data = {}
        # Preliminary workflow for attempting stark shift calculation. Presently doesn't provide reasonable results.

        # Use ground and excited states to calculate zeeman splitting and diangonalize the energy hamiltonian to obtain eigenvalues and vectors.
        Ground_Zeeman = Zeeman_func(Ground_Level, Bfield_ext)
        Excited_Zeeman = Zeeman_func(Excited_Level, Bfield_ext)
        # Use the new eigenvectors as the basis to calculate the relative signal strengths for the transitions
        Z_sig = Zeeman_signal(Ground_Level, Ground_Zeeman, Excited_Level, Excited_Zeeman, Bangle=Bext_angle, gamma=Polangle, Filter=Polfilter)
        Z_sig_norad = Zeeman_sig_norad(Ground_Level, Ground_Zeeman, Excited_Level, Excited_Zeeman)

        if 'DoLowSig' in Inputdeck:
            if Inputdeck['DoLowSig'] == 'N' or Inputdeck['DoLowSig'] == 'No':
                pass
            else:
                Ground_Level = stateprocess(L_ground, S_ground, E_ground, Bfield_ext, sortE=sortE)
                Excited_Level = stateprocess(L_excited, S_excited, E_excited, Bfield_ext, sortE=sortE)
                Ground_Zeeman_Low = Zeeman_func(Ground_Level, Bfield_ext, Lowfield=True)
                Excited_Zeeman_Low = Zeeman_func(Excited_Level, Bfield_ext, Lowfield=True)
                Z_sig_Low = Zeeman_signal(Ground_Level, Ground_Zeeman_Low, Excited_Level, Excited_Zeeman_Low, Bangle=Bext_angle, gamma=Polangle, Filter=Polfilter)
                # Z_sig_Low = Zeemansig_Low(Ground_Level,Excited_Level,Bangle=Bext_angle,gamma=Polangle, Filter=Polfilter)
                data['signal_low'] = Z_sig_Low[0]
                data['wave_vac_low'] = Z_sig_Low[6]
                data['wave_air_low'] = Z_sig_Low[1]
                data['EigGround_low'] = Ground_Zeeman_Low
                data['EigExcited_low'] = Excited_Zeeman_Low
        if 'DoHighSig' in Inputdeck:
            if Inputdeck['DoHighSig'] == 'N' or Inputdeck['DoHighSig'] == 'No':
                pass
            else:
                Ground_Level = stateprocess(L_ground, S_ground, E_ground, Bfield_ext, sortE=sortE)
                Excited_Level = stateprocess(L_excited, S_excited, E_excited, Bfield_ext, sortE=sortE)
                Ground_Zeeman_High = Zeeman_func(Ground_Level, Bfield_ext, Highfield=True, Eterm=EtermG)
                Excited_Zeeman_High = Zeeman_func(Excited_Level, Bfield_ext, Highfield=True,Eterm=EtermE)
                Z_sig_High = Zeeman_signal(Ground_Level, Ground_Zeeman_High, Excited_Level, Excited_Zeeman_High, Bangle=Bext_angle, gamma=Polangle, Filter=Polfilter)
                # Z_sig_Low = Zeemansig_Low(Ground_Level,Excited_Level,Bangle=Bext_angle,gamma=Polangle, Filter=Polfilter)
                data['signal_high'] = Z_sig_High[0]
                data['wave_vac_high'] = Z_sig_High[6]
                data['wave_air_high'] = Z_sig_High[1]
                data['EigGround_high'] = Ground_Zeeman_High
                data['EigExcited_high'] = Excited_Zeeman_High
        data['EigGround'] = Ground_Zeeman
        data['EigExcited'] = Excited_Zeeman
        data['ZeemanSplits_G'] = Ground_Zeeman[-1]

        data['bigstatelist'] = Z_sig[3]
        data['signal'] = Z_sig[0]
        data['wave_vac'] = Z_sig[6]
        data['wave_air'] = Z_sig[1]
        data['statelist_Ground'] = Ground_Level[-1]
        data['statelist_Excited'] = Excited_Level[-1]

        # Energy level output with lowfield approximation of zeeman splitting output.
        data['E0_Ground'] = Ground_Level[3]
        data['E0_Excited'] = Excited_Level[3]

        data['RedDipoleOp'] = Z_sig[-3]
        data['rad'] = Z_sig[2]
        data['Zsig_norad'] = Z_sig_norad
        # A_Ein isn't polarization normalized, A_Ein2 is the typical polarization normalized expression.
        data['A_Ein'] = [A_Ein(data['wave_air'][i], dips)/(2*data['bigstatelist'][i][1][0] + 1) for i, dips in enumerate(data['signal'])]
        data['A_Ein_polavg'] = [2*A_Ein(data['wave_air'][i], dips)/(3*(2*data['bigstatelist'][i][1][0] + 1)) for i, dips in enumerate(data['Zsig_norad'])]

        data['afcf'] = Z_sig[4]

    # There are several instrument functions available for convolution. Lorrentzian, Skewed Lorrentzian, Doppler, and DopplerInstrum, which is a combination of a doppler and instrument (not a Pseudo Voigt though.)

    # These conditionals can very much be cleaned up and simplified in the future.

    if 'ion_vel' in Inputdeck:
        ion_vel = Inputdeck['ion_vel']

    else:
        ion_vel = 0
    if 'Convfxn' in Inputdeck:
        Convtype = Inputdeck['Convfxn']
        convtemp = 300
        SkewedConstants = [1, 1]
        if 'specstep' in Inputdeck:
            stepsize = Inputdeck['specstep']
        else:
            stepsize = 0.002
        if Inputdeck['Convfxn'] == 'Gaussian' or Inputdeck['Convfxn'] == 'GaussianInstrum':
            convtemp = Inputdeck['Temp']
            SkewedConstants = [1, 1]

        if Inputdeck['Convfxn'] == 'Skewed' or Inputdeck['Convfxn'] == 'Lorrentzian':

            if 'Skewness' in Inputdeck:
                SkewedConstants = Inputdeck['Skewness']
            else:
                SkewedConstants = [1, 1]
        if Inputdeck['Convfxn'] == 'Voigt':
            convtemp = Inputdeck['Temp']

        if 'fxnwindow' in Inputdeck:
            fxnwindow = Inputdeck['fxnwindow']
        else:
            fxnwindow = 0.1  # This is how far the bin window will extend beyond the edge of the window. May need to be extended if large broadening.
    else:
        # Setting the convolution window for the function

        # This is the default configuration using a normal Lorrentzian lineshape.
        Convtype = 'Lorrentzian'
        SkewedConstants = [1, 1]
        convtemp = 300
        ion_vel = 0
        if 'fxnwindow' in Inputdeck:
            fxnwindow = Inputdeck['fxnwindow']
        else:
            fxnwindow = 0.1  # This is how far the bin window will extend beyond the edge of the window. May need to be extended if large broadening.
        if 'specstep' in Inputdeck:

            stepsize = Inputdeck['specstep']
        else:
            stepsize = 0.002222  # This is the resolution of the spectrometer pixels in nm.

    if 'specres' in Inputdeck:
        specres = Inputdeck['specres']
    else:
        specres = 10
    ConvolvedSpect, redsticks, binwinds = Convol_Spect(data['wave_air'], data['signal'],
                                                       specwindow, plotwindow, Inputdeck['amu'], stepsize,
                                                       Temperature_in=convtemp, functiontype=Convtype,
                                                       wind_size=fxnwindow, Skewness_consts=SkewedConstants, ionvel=ion_vel, specres=specres, bidir=bidirect)
    
    data['SpecOut'] = ConvolvedSpect
    data['reduced_sticks'] = redsticks
    data['bin_windows'] = binwinds

    pilist = np.array([x[1] for x in data['rad']]) == 0
    sigmaplus = np.array([x[1] for x in data['rad']]) == -1  # The sign is flipped.
    sigmaminus = np.array([x[1] for x in data['rad']]) == 1  # If you start from the bottom, you have the change the sign here
    otherlists = np.array([x[1] for x in data['rad']]) == 3

    # data['radlists'] = np.array([pilist,sigmaplus,sigmaminus])

    # data['radlists'] = np.array([pilist,sigmaplus,sigmaminus])
    data['pi_sticks'] = [data['wave_air'][pilist], data['signal'][pilist]]
    data['sigmaplus_sticks'] = [data['wave_air'][sigmaplus], data['signal'][sigmaplus]]
    data['sigmaminus_sticks'] = [data['wave_air'][sigmaminus], data['signal'][sigmaminus]]
    data['other_sticks'] = [data['wave_air'][otherlists], data['signal'][otherlists]]

    if 'DoLowSig' in Inputdeck:
        ConvolvedSpectLow, redsticksLow, binwindsLow = Convol_Spect(data['wave_air_low'], data['signal_low'],
                                                                    specwindow, plotwindow,  Inputdeck['amu'], stepsize,
                                                                    Temperature_in=convtemp, functiontype=Convtype,
                                                                    wind_size=fxnwindow, Skewness_consts=SkewedConstants, ionvel=ion_vel, specres=specres, bidir=bidirect)

        data['SpecOutLow'] = ConvolvedSpectLow
        data['reduced_sticksLow'] = redsticksLow
        # data['bin_windows'] = binwindsLow

    if 'DoHighSig' in Inputdeck:
        ConvolvedSpectHigh, redsticksHigh, binwindsHigh = Convol_Spect(data['wave_air_high'], data['signal_high'],
                                                                    specwindow, plotwindow,  Inputdeck['amu'], stepsize,
                                                                    Temperature_in=convtemp, functiontype=Convtype,
                                                                    wind_size=fxnwindow, Skewness_consts=SkewedConstants, ionvel=ion_vel, specres=specres, bidir=bidirect)

        data['SpecOutHigh'] = ConvolvedSpectHigh
        data['reduced_sticksHigh'] = redsticksHigh
        # data['bin_windows'] = binwindsLow
    data['Input'] = Inputdeck
    return data

# %%Convolution
def MultiSpec(appendlist, scalelist):
    """
    Combines several Zeeman outputs into single arrays for wavelength, signal,
    temperature, and atomic mass.
    """
    waves = np.array([])
    signals = np.array([])
    temps = np.array([])
    amus = np.array([])
    
    for scaler, obj in zip(scalelist, appendlist):
        waves = np.append(waves, obj['wave_air'])
        signals = np.append(signals, scaler * np.array(obj['signal']))
        
        # Capture metadata for each component to allow different T/amu during convolution
        # Defaults to 300K and 1 amu if not found 
        t_val = obj['Input'].get('Temp', 300)
        a_val = obj['Input'].get('amu', 1)
        
        temps = np.append(temps, np.full_like(obj['wave_air'], t_val))
        amus = np.append(amus, np.full_like(obj['wave_air'], a_val))
        
    return [waves, signals, temps, amus]

def Convol_Sticks(x0, xwind, Temp_in, amu_in, stepsize_nm, function="", Skew_consts=[1, 1], vel_ion=0, specres=10, bidirect=False):
    '''
    Assume that the inputs are all in units of nm.
    The difference in this function is that I'm assuming that this subroutine will be called after the sticks have already been limited to a larger window.
    This means that the full lineshape for a given convolution may not always be plotted, but then ensures that any given larger convolution will always end up the same length of array.
    '''
    # hstep = stepsize_nm # step size in nm!! (Want 10x the steps per pixe for better shape?)
    hstep = (1/specres)*stepsize_nm  # step size in nm!! (Want 1/specres times the steps per pixe for better shape)

    xmin = np.min(xwind)
    xmax = np.max(xwind)

    numwavesteps = int(round((xmax-xmin)/hstep))
    waveout = np.linspace(xmin, xmax, numwavesteps)

    match function:
        case 'Skewed':
            # This one just does convolution over the entire range, should be fine because the range is fed by the binned stick windows.
            tlorry = Lorrentzian(1e-9*waveout, 1e-9*x0, pixel_size=1e-9*stepsize_nm, Skewness=Skew_consts, v_ion=vel_ion, bidir=bidirect)
            tlorry[tlorry > 1] = 1  # Sometimes an odd range for the lorrentzian causes the central peak to be twice as high for a single point, this should correct for that.
            specout = tlorry
        case 'Lorrentzian':
            # This uses the skewed lorrentzian but sets to skewness such that the result is a normal lorrentzian function.
            tlorry = Lorrentzian(1e-9*waveout, 1e-9*x0, pixel_size=1e-9*stepsize_nm, Skewness=[1, 1], v_ion=vel_ion, bidir=bidirect)
            tlorry[tlorry > 1] = 1  # Sometimes an odd range for the lorrentzian causes the central peak to be twice as high for a single point, this should correct for that.
            specout = tlorry
        case 'Gaussian':
            tlorry = Gaussian(1e-9*waveout, 1e-9*x0, Temp_in, amu_in, v_ion=vel_ion, bidir=bidirect)
            # tlorry[tlorry>1] =1 #Sometimes an odd range for the lorrentzian causes the central peak to be twice as high for a single point, this should correct for that.
            specout = tlorry
        case 'GaussianInstrum':
            tlorry = GaussianInstrum(1e-9*waveout, 1e-9*x0, Temp_in, amu_in, pixel_size=1e-9*stepsize_nm, v_ion=vel_ion, bidir=bidirect)
            # tlorry[tlorry>1] =1 #Sometimes an odd range for the lorrentzian causes the central peak to be twice as high for a single point, this should correct for that.
            specout = tlorry
        case 'Voigt':

            tlorry = Voigt(1e-9*waveout, 1e-9*x0, Temp_in, amu_in, pixel_size=1e-9*stepsize_nm, v_ion=vel_ion, bidir=bidirect)
            # tlorry[tlorry>1] =1 #Sometimes an odd range for the lorrentzian causes the central peak to be twice as high for a single point, this should correct for that.
            specout = tlorry

        case _:
            # If no function is otherwise given, default to Lorrentzian
            # This uses the skewed lorrentzian but sets to skewness such that the result is a normal lorrentzian function.
            tlorry = Lorrentzian(1e-9*waveout, 1e-9*x0, pixel_size=stepsize_nm, Skewness=[1, 1], v_ion=vel_ion, bidir=bidirect)
            tlorry[tlorry > 1] = 1  # Sometimes an odd range for the lorrentzian causes the central peak to be twice as high for a single point, this should correct for that.
            specout = tlorry

    return [waveout, np.array(specout)]

def make_stickbins(sticks_in, *data_arrays, windsize=.1):
    wavered = sticks_in
    min0, max0 = np.min(wavered), np.max(wavered)
    tlen = len(wavered)
    tmax = tmin = min0
    wavewind = []

    while tmax != max0:
        if tlen == 1:
            wavewind.append([tmin - windsize, tmax + windsize])
            tred1 = wavered > tmax
            tmin = np.min(wavered[tred1])
            tred2 = wavered < tmin + windsize
            ttop = wavered[tred1 & tred2]
            tmax = np.max(ttop)
            tlen = len(ttop)
        else:
            tred1 = wavered >= tmax
            tred2 = wavered <= (tmax + windsize)
            ttop = wavered[tred1 & tred2]
            tmax = np.max(ttop)
            tlen = 1 if tmax == tmin else len(ttop)
    
    wavewind.append([tmin - windsize, tmax + windsize])
    
    wavebins = []
    binned_data = [[] for _ in data_arrays]
    for windows in wavewind:
        mask = (wavered > windows[0]) & (wavered < windows[1])
        wavebins.append(wavered[mask])
        for i, arr in enumerate(data_arrays):
            binned_data[i].append(np.array(arr)[mask])
            
    return {'binned_waves': wavebins, 'binned_data': binned_data, 'bin_windows': wavewind}

def Convol_Spect(waves_in, signal_in, spec_wind, plot_wind, atomic_weight_amu, resolution_nm, 
                 Temperature_in=300, functiontype="", wind_size=.1, 
                 Skewness_consts=[1, 1], ionvel=0, specres=10, bidir=False):

    # 1. Resolve inputs into arrays
    amus_arr = np.full_like(waves_in, atomic_weight_amu) if not hasattr(atomic_weight_amu, "__len__") else atomic_weight_amu
    temps_arr = np.full_like(waves_in, Temperature_in) if not hasattr(Temperature_in, "__len__") else Temperature_in

    # 2. Physics Reduction (spec_window) - used for the actual math
    phys_reducer = (np.array(signal_in) > (0.05 * np.max(signal_in)/len(signal_in))) & \
                   (waves_in > np.min(spec_wind)) & (waves_in < np.max(spec_wind))
    
    t_wav_phys, t_sig_phys = waves_in[phys_reducer], signal_in[phys_reducer]
    t_temp_phys, t_amu_phys = temps_arr[phys_reducer], amus_arr[phys_reducer]

    # 3. Stick Reduction (plot_window) - THIS IS THE FIX
    # This ensures redsticks contains ONLY what is visible on the plot
    plot_mask = (t_wav_phys >= np.min(plot_wind)) & (t_wav_phys <= np.max(plot_wind))
    redsticks = [t_wav_phys[plot_mask], t_sig_phys[plot_mask]]

    # 4. Binning and Convolution
    tbins = make_stickbins(t_wav_phys, t_sig_phys, t_temp_phys, t_amu_phys, windsize=wind_size)

    tempwavs, tempsigs = [], []
    for sticks, data_list, winds in zip(tbins['binned_waves'], zip(*tbins['binned_data']), tbins['bin_windows']):
        # Optimization: Skip bins entirely outside the plot window
        if winds[1] < np.min(plot_wind) or winds[0] > np.max(plot_wind):
            continue
            
        sigs, temps, amus = data_list
        twav = np.linspace(winds[0], winds[1], int(round(specres*(winds[1] - winds[0])/resolution_nm)))
        tsig = np.zeros_like(twav)
        
        for w, s, t, a in zip(sticks, sigs, temps, amus):
            tspec = Convol_Sticks(w, [winds[0], winds[1]], t, a, resolution_nm, 
                                 function=functiontype, Skew_consts=Skewness_consts, 
                                 vel_ion=ionvel, specres=specres, bidirect=bidir)
            tsig += s * tspec[1]

        tempwavs.extend(twav)
        tempsigs.extend(tsig)

    # 5. Final Spectrum Clipping
    out_w, out_s = np.array(tempwavs), np.array(tempsigs)
    mask = (out_w >= np.min(plot_wind)) & (out_w <= np.max(plot_wind))
    
    final_spectrum = [np.concatenate([[np.min(plot_wind)], out_w[mask], [np.max(plot_wind)]]),
                      np.concatenate([[0], out_s[mask], [0]])]

    return final_spectrum, redsticks, tbins['bin_windows']


# %%Read Inputs





def read_Spectra(filename, scalewave=1, headercount=1, keepheaders=False):
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

        # Read through the rows and scale the wavelength by the scalefactor.
        for rows in reader:
            specdat.append([float(x) for x in rows])

        specdat = np.array(specdat).T  # Transpose so that the ith entry is all wavelengths or intensities.
        specwav = specdat[0]

        Spectraout = []
        [Spectraout.append(x) for i, x in enumerate(specdat[1:])]

        # This small bit gives magnetic field values based on the headers info
        # Bvals = [float(x) for x in headers[0][1:]]

        specdata_out = {}
        specdata_out['wavelength'] = [scalewave*x for x in specwav]  # output wavelength, assumed in nm.
        # specdata_out['Bvals'] = Bvals
        if len(Spectraout) == 1:
            specdata_out['signals'] = Spectraout[0]
        else:
            specdata_out['signals'] = Spectraout

    if keepheaders == True:
        return [specdata_out, headers]
    else:
        return specdata_out


def read_dict(filename):

    with open(f'{filename}', mode='r', newline='') as readme:
        reader = csv.DictReader(readme, delimiter=',', quoting=csv.QUOTE_NONE)
        # headers = []
        # for i in range(headercount):
        #     headers.append(next(reader))
        datraw = []

        # Read through the rows and scale the wavelength by the scalefactor.
        for rows in reader:
            trow = {k: v for k, v in rows.items() if v}

            for keys in trow.keys():
                try:
                    trow[keys] = float(trow[keys])
                except:
                    pass

            datraw.append(trow)

    return datraw

def read_HFSInput(filename):

    with open(f'{filename}', mode='r', newline='') as readme:
        reader = csv.DictReader(readme, delimiter=',', quoting=csv.QUOTE_NONE)
        # headers = []
        # for i in range(headercount):
        #     headers.append(next(reader))
        datraw = []

        # Read through the rows and scale the wavelength by the scalefactor.
        for rows in reader:
            trow = {k: v for k, v in rows.items() if v}

            for keys in trow.keys():
                try:
                    trow[keys] = float(trow[keys])
                except:
                    pass

            datraw.append(trow)

    return datraw
#%% Save and Load
def Savenpy_Dict(filename, Dict_in,spectra = None, Normalize_sig = False):
    savedat = {}

    if spectra:
        specdat = {}
        try: #Try assuming spectra is a dictionary with keys['wavelength' and 'signals'
           specdat['wavelength'] = np.array(spectra['wavelength'] )
        except: #Else, spectra should be [wavelength,signals] , with each as arrays. 
            specdat['wavelength'] = np.array(spectra[0] )
        #Allow for signals to be multiple lists
        try :
            len(np.shape(spectra['signals']))
            if len(np.shape(spectra['signals']))==2:
                if Normalize_sig:
                    specdat['signals'] = [Normalize(x) for x in spectra['signals']] 
                else:
                    specdat['signals'] = [x for x in spectra['signals']] 

            else:
                if Normalize_sig:

                    specdat['signals'] = Normalize(spectra['signals'])
                else:
                    specdat['signals'] = spectra['signals']

        except:
            if len(np.shape(spectra[1]))==2:
                if Normalize_sig:
                    specdat['signals'] = [Normalize(x) for x in spectra[1]] 
                else:
                    specdat['signals'] = [x for x in spectra[1]] 

            else:
                if Normalize_sig:
                    specdat['signals'] = Normalize(spectra[1])
                else:
                    specdat['signals'] = spectra[1]

  
    
        savedat['Specdata'] = specdat
    savedat['data'] = Dict_in
    np.save(f'{filename}',savedat,allow_pickle=True)
    return

def Loadnpy_Dict(filepath):
    #This function takes the numpy.save dictionary and loads it, including the spectrometer data, for future use.
    #Need to test if figures can be saved too...
    thisload = np.load(f'{filepath}.npy', allow_pickle=True)
    loadeddict = thisload.item().get('data')
    try:
        loadedspec = thisload.item().get('Specdata')
        return loadeddict, loadedspec
    except:
        pass
        # print('Failed to load Spectrometer data - May not exist')
    return loadeddict

