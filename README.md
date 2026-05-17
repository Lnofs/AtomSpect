# What is AtomSpect?
AtomSpect is a newly developed Python package designed for simplify the spectroscopic analysis for diagnostic purposes. It has been developed to be lightweight and easy to implement into any diagnostic pipepline - aided by being fully open source with GPL-3.0 licensing.

AtomSpect is NOT intended to be a replacement for the robust ab-initio atomic structure calculations for energy levels. 

# Installation
To use: Unzip the repository and open one of the example scripts in the examples folder. Run the example script, not AtomSpect.py.


# Using AtomSpect
AtomSpect uses an Input Dictionary with several keyword arguments for the calculation. This includes a vast amount of optional keywords that can be passed. A full list and description of the keywords is provided below and in the function description for AtomSpect_Main.

The InputDict can be passed to either: MakeSlider, plotZfan, or AtomSpect_Main. MakeSlider and plotZfan will directly generate plots from the input file.

<img width="631" height="345" alt="image" src="https://github.com/user-attachments/assets/bdba6c53-4c9c-4437-a047-57d94de96b71" />

AtomSpect_Main outputs a dictionary with many results and intermediary steps. This output is passed into the PlotFunction to generate the plot. There are additional optional plot flags as well.
PlotFunction requires a PlotObject and plot_vars. The PlotObject can be either the AtomSpect_Main output or a list of [wavelength,intensity].
The PlotFunction is robust with a lot of optional arguments. It is an array of plot axes which can be addressed via index. This gives the flexability of a generalized 2D plot array. You can feed this axis, which is the output of PlotFunction into PlotFunction.

In addition to the PlotObject, several plotting parameters are set with the plot_vars = [color, linewidth,linestyle, legendlabel, markerstyle(optional)]

Image showing a very simple example of operation. Additional spectra can be added but ensure to set makefig=False and axisin=thisaxs (per the Argon example shown here).
<img width="910" height="398" alt="image" src="https://github.com/user-attachments/assets/afda7519-e135-4097-90c4-ac3717eff447" />

The versatility of the plotting function is showcased across several of the included examples.

## Slider GUI!
An input dictionary can be fed into MakeSlider. Boolean flags exist for:
do_sticks - Plotting component "stick" lines. Default True
do_Polplot - Labels the polarization of each component stick. Default True
do_other - whether components with delta M > 1 are plotted (purple star sticks). Default False
banglevary - Whether angle of magnetic field is allowed to vary
polvary - Whether a polarization filter is along the line of sight and if the filter angle can vary.
Tmax - Maximum temperature [eV]. Default 30
Bmax - Maxium field strength [T]. Default 5

The Slider always lets: Bfield strength (in T), Temperature (eV), and Wi (~Spectrometer pixel resolution in nm) vary.
If you double click on the numerical value under the slider, you can enter a value directly to be calculated. The value entered here CAN be outside the normal slider ranges!

Each change in the slider currently (April 2026) causes a full recalculation, diagonalization, convolution, and re-plotting of the spectra.
Fortunately, this is optimized enough to be real-time viable.

Default ranges:
B: [0,5] T, dB: 0.01
T: [0.01,30] , dT: 0.01
Wi: [0.002,0.1] , dWi: 0.001
Bangle: , dBangle: 0.5
PolAngle: , dPolAngle: 0.5

Max Value 
<img width="1549" height="948" alt="image" src="https://github.com/user-attachments/assets/566f9ece-7355-41d1-9492-0816b87c07d2" />

## Input Dictionary Keys

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
                     'Convfxn': 'Gaussian' , #Optional: 'Gaussian', 'GaussianInstrum'', 'Skewed' , 'Lorrentzian', 'Voigt'.
                     'Skewness' : [0.62,0.99] #Optional, [Left,Right] constant for the skewness of the lorrentzian function. Between 0 and 1. If both are 1, it is a normal lorrentzian.
                     'Temp' : 300, #Temperature in K. Used for Gaussian, Gaussian Instrum, Voigt convolutions
                     'TempeV' : 20 , #Temperature in eV. Can be passed instead of Temp in K. Converted as needed internally.
                     'specstep' : 0.002222 , #stepsize of linefunction [nm]
                     'specres' : 10 , #How many points are evaluated per spectrometer pixel window. Higher means smoother curve.
                     'fxnwindow': 0.5 , #How far from the central peak the convolution will be calculated. Also related to how stick binning works. 
                     'plottitle': 'Excample Title' , #Title for plotting (optional).                      
                     'SpectrumData' : [[Wavelength],[Intensity]] , #List of arrays for wavelength (nm) and intensity (arb)
                     'spec_window' : [425,435], #Min and max for convolution window (nm)
                     'plot_window' : [425,435], #Min and max for plot window range (nm)
                     'HFS_G' : [505.5e6] , #Optional: Hyperfine A constant for the lower level, need one entry per J level, lowest J first [Hz]
                     'HFS_E': [0, 496.2e6, 440.5e6] , #Optional: Hyperfine A constants for upper level, need one entery per J level, lowest J first [Hz]
                     'I_spin' : 0.5, #Optional: nuclear spin, I. Tabulated by Stone et al.
                     'mu_I' :  .11778476 , #Optional: Nuclear dipole moment as tabulated by Stone.              
                     'HFSB_G' : [505.5e6] , #OPtional: Hyperfine B constant for the lower level, need one entry per J level, lowest J first [Hz]
                     'HFSB_E': [0, 496.2e6, 440.5e6] , #Optional: Hyperfine B constants for upper level, need one entery per J level, lowest J first [Hz]
                     'Calc_HFS_Vals': [False,False], #Optional. Passed regarding wether HFS constants A and B should be approximated for values of 0 above. NOT YET IMPLEMENTED
                     'ion_vel' : Float or [Float1,Float2] , #Velocity of particles, mono and bidirectional allowed. A single value can be passed and used for bidirectional splitting if bidir flag is set.
                     'bidir' : True , #Bool as to whether ion motion generates bidirection doppler shift
                     'PsiHFS' : True , #Bool, testing flag left in for completeness
                     'DoLowSig' : True , #Bool for whether Lowfield approximation and signal are included and plotted
                     'DoHighSig' : True , #REQUIRES TERM ENERGY!Bool for whether Highfield approximation and signal are included and plotted
                     'EtermG': 169086.9076, #Term energy of lower state in cm^-1 - Required for DoHighSig
                     'EtermE' : 191444.47832914, #Term energy of upper state in cm^-1 - Required for DoHighSig
                     'sortE' : True , #Bool as to whether the wavefunctions will attempt to be sorted based on expected resulting energy before diagonalization. This significantly slows down operations. 
                      }
                      


# Plotting Zfans
AtomSpect can do typical Zeeman Fan plots with a dedicated set of functions which has been optimized for significant speed.
Compared to a normal InputDictionary, 'Bmag' needs to be an array of magnetic field values. The flags refer to:
"G" - Ground
"E" - Excited
"H" - Add High-field approximation. This requires Term Energy!
"L" - Add Low-field approximation (no diagonalization)

At least "G" or "E" are required, otherwise, any combination of the flags, in any order, may be passed into plotZfan.

<img width="888" height="465" alt="image" src="https://github.com/user-attachments/assets/26f8bbbd-3fc0-449f-9250-d2a3c377518c" />

<img width="443" height="361" alt="image" src="https://github.com/user-attachments/assets/cc4a162e-5524-4922-82e0-025ab29ecf32" />

# Future Developement
Stark Shift. Work is already well underway for implementing the DC stark shift for electric fields into the package.
Hyperfine Constants. Preliminary work is done and results are sound, just need to simplify and streamline the calculation into AtomSpect_Main
Importing/saving data. We plan to allow for importing of atomic energy level data such as from then Open ADAS work. This also includes building up a database of atoms and energies for quick use such that the user doesn't need to type in the details for every line they are interested. This goes great with the next feature...
Fancy GUI! The Slider GUI is okay, but is unable to handle multiple species/lines. A robust GUI is in the initial stages of design which would transform everything about using AtomSpect.

# Data Use
All data provided here is open source and free to use. Additional spectrometer data is always welcome.


I kindly request that any works that utilize this code reference the AtomSpect publication:(doi goes here)

Contact me at: lmn0022@auburn.edu

