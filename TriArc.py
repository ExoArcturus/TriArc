#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 11:39:42 2022

@author: arcturus
"""

# import the usual...
import numpy as np
import pylab as plt


# import constants from pRT
from petitRADTRANS import nat_cst as nc
# import main package from pRT
from petitRADTRANS import Radtrans

# import molecular weights of species
import molecular_weights as mw


# ----- CLASSES -----

# define exoPlanet class, with all inputs for exoplanetary system

class exoPlanet:
    def __init__(self, R_pl, R_st, gravity, P0, temp, 
                 composition, lin_spec, ray_spec, 
                 cont_spec, MMW, mixing_ratios, species_l):
        
        self.R_pl = R_pl                            # radius of planet (surface defined by pressure, P0, for gaseous planets)
        self.R_st = R_st                            # radius of parent star
        self.gravity = gravity                      # surface gravity of planet
        self.P0 = P0                                # surface pressure (pressure at R_pl)
        self.temp = temp                            # temperature (currently isothermal temperature)
        self.composition = composition              # atmospheric composition dictionary (as mass fractions, currently with uniform mixing)
        self.lin_spec = lin_spec                    # list of species contributing to line opacity to consider
        self.ray_spec = ray_spec                    # list of species contributing to Rayleigh opacity to consider
        self.cont_spec = cont_spec                  # list of collisionally-induced-absorption opacities to consider
        self.MMW = MMW                              # mean molecular weight of atmosphere (vertically constant)
        self.mixing_ratios = mixing_ratios          # mixing ratios dictionary of background (not currently utilised)
        self.species_l = species_l                  # complete list of species in atmosphere to consider
        

# define Atmosphere_Model class to keep calculated atmospheric properties packaged together (composition - mass fractions, planet - exoPlanet object, temps - temperature profile, model - opacities)

class Atmosphere_Model:
    def __init__(self, comp, planet, temps, model):
        
        self.comp = comp                            # vertical composition dictionary profile of mass fractions
        self.planet = planet                        # exoPlanet object to keep simple inputs stored
        self.temps = temps                          # vertical temperature profile with pressure
        self.model = model                          # Radtrans object, containing all opacity data

# define Spectral_band class, with all necessary parameters to describe the spectral bands used to target the retrieval

class Spectral_Band:
    def __init__(self, band_label, min_wlen, max_wlen, precision, 
                 species_1, species_2=None,species_3=None,species_4=None,species_5=None,species_6=None):
        
        self.label = band_label                     # used to label the band, SHOULD BE SAME AS KEY OF DICTIONARY containing all spectral bands
        self.min_wlen = min_wlen                    # upper wavelength limit of band in microns
        self.max_wlen = max_wlen                    # lower wavelength limit of band in microns
        self.precision = precision                  # precision of instrument (noise to be added) at the band
        species = []                                
        species.append(species_1)                   # add all inputted species (up to 5) into a list
        if species_2 != None:                       
            species.append(species_2)
        if species_3 != None:
            species.append(species_3)
        if species_4 != None:
            species.append(species_4)
        if species_5 != None:
            species.append(species_5)
        if species_6 != None:
            species.append(species_6)
        self.species = species                      # initialise list of species
        
   
    

    
# Plot parameters for reasonably nice figures
plt.rcParams['font.family'] = "serif"
plt.rcParams['figure.figsize'] = (10,6)
plt.rcParams['font.size'] = 18
    


def logsumexp(arr):
    max_val = np.max(arr)
    arr_shifted = arr - max_val
    arr_exp = np.exp(arr_shifted)
    sum_exp = np.sum(arr_exp)
    result = np.log(sum_exp) + max_val
    return result

# ----- INITIALISATION FUNCTIONS FOR IMPORTING INPUT FILES -----

""" planet input files must include the following:
    
    planet_radius - loaded in as R_pl (in cgs units)
    star_radius - loaded in as R_st (in cgs units)
    surface_gravity - loaded in as gravity (in cgs units)
    surface_pressure - loaded in as P0 (in bars)
    iso_temperature - loaded in as temp (in kelvin), currently isothermal only
    
    line_species - list of all species contributing to line opacity to consider
    rayleigh_species - list of all species contributing to rayleigh opacity to consider
    cont_species - list of all species contributing to continuum CIA opacity
    species_list - list of all species to consider, including not contributing to opacity (affects MMW)
    
    abundances - dictionary of mixing ratios, using molecule name as a key, e.g. abundances['CH4'] = 0.01 (1% methane)
"""

# load in planet input file (e.g. Hycean.py) which should be imported in script utilising TriArc, and creates exoPlanet object.
# also performs conversion to mass fraction and calculates mean molecular weight.


def load_planet(atm):                       # argument is an imported file (i.e. something.py after import something)
    print('Loading in planet input file: ' +str(atm)+'.py')

    # import atmospheric composition (in mixing ratios)
    n = atm.abundances
    
    # calculate mean molecular weight (MMW)
    MMW = 0
    species_list=[]
    species_list = atm.species_list
    
    #remove duplicates from species list using dictionary
    species_list = list(dict.fromkeys(species_list))
    for species in species_list:
        MMW = MMW + mw.mass[species] * n[species]
    print('Atmosphere mean molecular weight: '+str(MMW))
    
    # convert from abundances to mass fractions
    mass_fraction={}
    for species in species_list:
        mf =  n[species] * ( mw.mass[species] / MMW)
        mass_fraction[species]=mf
    print('Atmosphere mass fractions:')
    print(mass_fraction)
    
    # creates exoPlanet object based on inputs loaded in from input file
    generic_planet = exoPlanet(atm.planet_radius,
                               atm.star_radius,
                               atm.surface_gravity,
                               atm.surface_pressure,
                               atm.iso_temperature,
                               mass_fraction,
                               atm.line_species,
                               atm.rayleigh_species,
                               atm.cont_species,
                               MMW,
                               n,
                               species_list)
    print(str(atm)+'.py successfully loaded.')
    return generic_planet

""" band data input files must import the Spectral_Band class from TriArc:
    from TriArc import Spectral_Band
    
    the each band should be specified with:
        band_list={}
        band_list['label'] = Spectral_Band('label', min bound on wavelength, max bound on wavelength, precision, 'species 1', 'species 2' (etc.))
    
    note additional species (beyond one) are optional. the precision should be a default or placeholder value, if not fetching JWST noise from text files.
"""

# load in band list input file (e.g. prebiosignatures.py) which should be imported in script utilising TriArc,
# and creates dictionary of Spectral_Band objects.
# will also fetch JWST noise if JWST_nosie = True, otherwise will use default values which need to be specified in input file

def compile_band_list(band_input_data,
                      JWST_noise=True,
                      spectral_resolution=1000, 
                      use_prism=True, 
                      auto_assign=True):
    
    # create dictionary of Spectral Band objects from band_input_data input file (e.g. prebiosignature.py)
    print('Loading in band data input file: '+str(band_input_data)+'.py')
    all_bands = {}
    all_bands = band_input_data.band_list       

    # overrides default noise values with JWST noises created using PandArc or otherwise if JWST_Noise == True.    
    if JWST_noise == True:
        noise = implement_JWST_noises(all_bands,auto_assign=auto_assign,use_prism=use_prism,spectral_resolution=spectral_resolution)
        for band in all_bands.keys():
            all_bands[band].precision = noise[band]
            
    # otherwise uses default values which must be specified in band input data
    else:
        print('Using default noise specified in band_input_data file.')
    return all_bands


# loads in JWST noise text file labelled with instrument_R_spectral_resolution e.g. NIRSpec G395M_R_1000
# calculates and returns the average noise within the wavelength range specified
# see PandArc.py for how to create JWST noise files (which all automatically be created with the appropriate name)

def fetch_JWST_noise(wlen_min,                          # minimum wavelength 
                     wlen_max,
                     instrument,
                     spectral_resolution=1000):
       
    # load in text file, first column should be wavelengths (in microns), second column should be precision in ppm
    instrument_data = np.loadtxt(instrument+'_R_'+str(spectral_resolution)+'.txt')
    wavelengths = instrument_data[:,0]
    precision = instrument_data[:,1]
    
    # discard values out of range of specified wavelengths
    precision_in_range = precision[np.where(wavelengths < wlen_max)]
    new_wavelengths = wavelengths[np.where(wavelengths < wlen_max)]
    precision_in_range_in_range = precision_in_range[np.where(new_wavelengths > wlen_min)]
    
    # calculate mean precision
    mean_precision = np.mean(precision_in_range_in_range)
    return mean_precision


# automatically assigns instruments based on approach described in Claringbold et al. 2022
# for 'R<=100, use_prism=True' uses NIRSpec Prism, NIRSpec G395M and MIRI LRS.
# for 'R<=100, use_prism=False' uses NIRSpec G140M, NIRSpec G235M, NIRSpec G395M and MIRI LRS.
# for 'R>100' uses NIRSpec 140H, NIRSpec G235H and NIRSpec  
# (note no option for MIRI, so do not include bands > 5 microns at high resolution)!!!

# if not using auto assign, must specify instruments[band] in script

# be careful that each band is entirely within the wavelength range of a single instrument

def auto_assign_instruments(wlen_max,spectral_resolution,use_prism=True):
        
    # high resolution (<5 microns only), R > 100
    if spectral_resolution > 100:
        if wlen_max > 5:
            print('MIRI LRS (>5 micron band) does not operate at this spectral resolution. Use spectral resolution of R=100 or less, or remove the band.')
            instrument = None
        elif 5 >= wlen_max > 2.9:
            instrument = 'NIRSpec G395H'
        elif 2.9 >= wlen_max > 1.8:
            instrument = 'NIRSpec G235H'
        elif wlen_max <= 1.8:
            instrument = 'NIRSpec G140H'
        else:
            print('Invalid wavelength range on band.')
            instrument = None
            
    # low resolution, R<= 100  
    else:
        if wlen_max > 5:
            instrument = 'MIRI LRS'
        elif 5 >= wlen_max > 2.9:
            instrument = 'NIRSpec G395M'
        elif 2.9 >= wlen_max:
            
            # gives option to exclude NIRSpec Prism and instead use G140M and G235M filters
            # this may be required if the target is too bright for the prism (saturation issues)
            if use_prism == True:
                instrument = 'NIRSpec Prism'
            else:
                instrument = 'NIRISS SOSS'
    return instrument
            

# Function to assign JWST data to the bands specified in the band_data file, overriding default values in file.

def implement_JWST_noises(band_data,auto_assign=True,use_prism=True,spectral_resolution=1000,instruments=None):
    precision = {}
    
    # auto-assigns instruments, if manually specified must specify all bands with instruments[band] in script.
    if auto_assign == True:
        instruments = {}
        print('Automatically assigning JWST instruments to each band based on wavelength...')
        for band in band_data.keys():
            instruments[band] = auto_assign_instruments(band_data[band].max_wlen,spectral_resolution)
    
    # repeat: if not using auto assign, must specify instruments[band] in script!!!
         
    print('Loading in JWST noise data...')
    
    # using fetch_JWST_noise for each band to enter in noise data
    for band in band_data.keys():
        precision[band] = fetch_JWST_noise(band_data[band].min_wlen,band_data[band].max_wlen,instruments[band],spectral_resolution=spectral_resolution)
    return precision

# ----- Intermediate functions, used by the high-level functions (e.g. perform retrieval, access detectability) -----


#Function to calculate the transmission spectrum for a given composition (comp), exoPlanet object (planet), temperature profile (temps), Radtrans object (model). 

def transmission_spectrum(comp, planet, temps, model, **kwargs):
    #input uniform mixing ratio in column
    mean_mol_weight = planet.MMW*np.ones_like(temps)
    
    #use pRT function to calculate transmission spectrum
    model.calc_transm(temps, comp, planet.gravity, mean_mol_weight, R_pl=planet.R_pl, P0_bar=planet.P0, **kwargs)
    
    #return transmission spectrum in ( transmission radius / star radius ) squared, in ppm
    return ((model.transm_rad/planet.R_st)**2)*10**6, model.freq

#Function to add noise to a transmission spectrum.

def add_noise(model_spec, noise_sd):
    # add simple gaussian noise to all data points in spectrum
    noisy_signal = np.random.normal(model_spec,noise_sd)
    return noisy_signal

#Function to estimate the goodness of fit between a signal (noisy_signal) and a test signal (signal) using a radial basis likelihood function.

def goodness_of_fit(noisy_signal, signal, sigma):
    # calculate value of normalised gaussian least squares at each wavelength
    gauss_least_squares = np.log(((2 * np.pi * sigma**2)**(-1/2)) * np.exp(-(1/(2*(sigma**2))*(noisy_signal-signal)**2)))
    # find product of normalised gaussian least squares to calculate likelihood function
    log_likelihood = np.sum(gauss_least_squares)
    return log_likelihood

#Function to recalculate the MMW with one species removed from the background atmosphere (not-functional ATM).

def recalc_MMW(background_atmosphere,omitted_species):
    MMW = 0
    species_list = background_atmosphere.species_l
    if omitted_species in species_list:
        species_list.remove(omitted_species)
    for species in species_list:
        MMW = MMW + mw.mass[species] * background_atmosphere.n[species]
    return MMW

#Function to correctly suffix line species such that rebinned opacities are used instead when specifying a spectral resolution other than default (1000).

def select_spectral_resolution(test_species,prebio_species,background_atmosphere,spectral_resolution=1000):
    #use atm for brevity
    atm = background_atmosphere

    # add '_R_100' (example) suffix to line species names, so pRT calls on the correct opacities
    # (requires .h5 files at that resolution to have been made and named approiately)
    new_line_spec = [s + '_R_' + str(spectral_resolution) for s in atm.lin_spec]
    
    # add '_R_100' (example) suffix to keys of composition dictionary
    new_res_composition = {k+'_R_'+str(spectral_resolution): v for k, v in atm.composition.items() if k in atm.lin_spec}
    for m in atm.lin_spec:
        if m in atm.composition.keys():
            atm.composition.pop(m)
 
    atm.composition.update(new_res_composition)
    
    # correct the species list with the appropriately suffixed species
    for count, species in enumerate(atm.species_l):
        if species in atm.lin_spec:
            atm.species_l[count] = species+'_R_'+str(spectral_resolution)
        else:
            atm.species_l[count] = species
    atm.lin_spec = new_line_spec
    
    # suffix the test and prebiosignature species as well
    prebio_species = prebio_species + '_R_'+str(spectral_resolution)
    test_species = test_species + '_R_'+str(spectral_resolution)
    
    # return all the species list (note continuum/rayleigh unchanged, as desired) with now correct suffixes
    return test_species, prebio_species, atm.species_l, atm.lin_spec, atm.composition


#Function to setup the pressure-temperature profile and opacity data to create Atmosphere_Model object.

def setup_atmosphere(background_atmosphere,wlen_min,wlen_max,spectral_resolution=1000,test_species='CO',prebio_species='CO', **spectrum_kwargs):
    #Import background atmosphere object, with arguments R_pl (radius of planet), R_st (radius of star), P0 (surface pressure), gravity, temp, composition (mass fractions), lin_spec, ray_spec, cont_spec, MMW (mean molecular weight), n (vol fractions)
    atm = background_atmosphere
    
    print('Setting up atmosphere...')
    
    #Instructs the use of opacities if rebinned to lower spectral resolution (must be prepared beforehand)
    if spectral_resolution != 1000:
        test_species, prebio_species, atm.species_l, atm.lin_spec, atm.composition = select_spectral_resolution(test_species,prebio_species,atm,spectral_resolution)
    
        
    #Reads in opacity sources in Radtrans object
    model_atmosphere = Radtrans(line_species = atm.lin_spec, rayleigh_species = atm.ray_spec, continuum_opacities = atm.cont_spec, wlen_bords_micron = [wlen_min,wlen_max])

    #Set up pressure-temperature structure of atmosphere (assumes isothermal)
    pressures = np.logspace(-6,np.log10(atm.P0),100)
    model_atmosphere.setup_opa_structure(pressures)
    temperature = atm.temp * np.ones_like(pressures)
    
    #Set up vertical mass_fractions in atmosphere (assumes uniform mixing)
    mass_fractions={}
    for species in atm.species_l:
        mass_fractions[species] = atm.composition[species] * np.ones_like(temperature)
    
    #Stores all atmospheric properties in object
    back = Atmosphere_Model(mass_fractions,atm,temperature,model_atmosphere)
    
    #Calculates a transmission spectrum for background atmosphere
    background_spectrum, frequencies = transmission_spectrum(mass_fractions, atm, temperature, model_atmosphere,**spectrum_kwargs)
    
    print('    Atmosphere setup successfully.')
    
    return back, background_spectrum, test_species, prebio_species
    

#Function to generate test transmission spectra (the 'hypotheses' in Bayes' theorem) for different abundances to test against.

def generate_test_spectra(back_atm,background_spec,test_species,signal_samples=111,**spectrum_kwargs):
    
    print('Creating hypothesis spectra...')
    
    #Choose uniform prior in log space to apply sampling over
    test_abundances = np.logspace(-11,0,num=signal_samples)
    
    #Calculate transmission spectrum for each test abundance sampled, remove background, and add to list.
    signal_list=[]
    for abundance_species in test_abundances:
        back_atm.comp[test_species] = abundance_species * np.ones_like(back_atm.temps)
        signal_bg, frequencies = transmission_spectrum(back_atm.comp, back_atm.planet, back_atm.temps, back_atm.model, **spectrum_kwargs)
        signal = signal_bg - background_spec
        signal_list.append(signal)
    
    print('    Hypothesis spectra successfully created.')
    
    return signal_list, test_abundances

#Primary function, injects prebiosignature into background, performs Bayesian analysis against set of test spectra, to output posterior PDF.
#Uses add_noise and goodness_of_fit functions

def bayesian_analysis(prebio_signal,back_atm,background_spec,prebio_species,noise,signal_list,signal_samples=111,**spectrum_kwargs):
    
    #Adds specified amount of prebiosignature molecule
    back_atm.comp[prebio_species] = prebio_signal * np.ones_like(back_atm.temps)
    
    #Calculate a model transmission spectrum
    model_bg, frequencies = transmission_spectrum(back_atm.comp, back_atm.planet, back_atm.temps, back_atm.model, **spectrum_kwargs)
    
    #Subtract the background
    model =  model_bg - background_spec
    
    #Adds noise to the transmission spectrum, and calculates goodness of fit against the test spectra.
    
    repeats = 900
    fits_list=np.empty([repeats,signal_samples])
    sigma=noise
    
    for i in range(0,repeats):
        all_fits=[]
        noisy_signal = add_noise(model, noise)
        for test_signal in signal_list:
            log_fit = goodness_of_fit(noisy_signal, test_signal, sigma)
            #fit = 10**log_fit
            #print(log_fit)
            all_fits.append(log_fit)
            
        normalisation=logsumexp(all_fits)
        all_fits = all_fits - normalisation
        
        
        fits_list[i]=10**all_fits
    #print(fits_list)
        
    #Calculates the mean fit from all the fits for each noise profile, and plots it.
    
    sum_fits = np.sum(fits_list,axis=0)
    
    #This is the posterior PDF of the bayesian analysis
    fits_mean = sum_fits/repeats
    #print(fits_mean)
    return fits_mean

# ----- HIGH-LEVEL FUNCTIONS FOR RETRIEVALS -----

# ----- Function to plot the transmission spectrum/spectra -----
# Plot up to three spectra on the same plot: background, with a prebiosignature added, with a prebiosignature and noise added
# By default only the prebiosignature added is plotted (set prebio_signal = 0 for this to be the background)
# Use plot_background = True to plot the background
# Use plot_noisy = True to plot the noisy signal in addition to the noiseless spectrum
# background_atmosphere - exoPlanet object input data (after using .load_planet(input file))
# wlen_min - minimum bound on wavelength to plot
# wlen_max - maximum bound on wavelength to plot
# prebio_spec - species to add to background
# plot_label, plot_title - to adorn the plot
# prebiosignal - mass fraction of prebio_spec to add to the background
# spectral resolution - spectral_resolution to plot at
# noise - noise to add to noisy spectrum to be plotted

def plot_spectrum(background_atmosphere,
                  filename='Figure.pdf',
                  wlen_min=1,
                  wlen_max=10,
                  prebio_spec='no_species',
                  plot_label="", plot_title="",
                  prebio_signal=0,
                  spectral_resolution=1000,
                  plot_background=False,
                  noise=0,
                  plot_noisy=False,
                  clear_plot=True,
                  custom_xticks=False,
                  **spectrum_kwargs):
    
    # Set up atmosphere
    atmosphere, background_spectrum, test_species, prebio_species = setup_atmosphere(background_atmosphere,
                                                                                     wlen_min, wlen_max,
                                                                                     spectral_resolution=spectral_resolution,
                                                                                     prebio_species=prebio_spec,
                                                                                     **spectrum_kwargs)
    
    # Inject prebio species into background atmosphere
    atmosphere.comp[prebio_species] = prebio_signal * np.ones_like(atmosphere.temps)
    
    # Compute transmission spectrum
    print('Computing transmission spectrum')
    spectrum, frequencies = transmission_spectrum(atmosphere.comp, atmosphere.planet, atmosphere.temps, atmosphere.model, **spectrum_kwargs)
    
    # Create x ticks for the plot
    x_ticks = np.linspace(wlen_min,wlen_max,10)
    x_ticks_int=[]
    for i in range(len(x_ticks)):
        x_ticks_int.append(int(x_ticks[i]))

    # prepare the plot
    plt.title(plot_title)
    plt.xscale('log')
    plt.xlabel('Wavelength ($\mu$m)')
    plt.ylabel('Transit depth $(R_t/R_*)^2$ (ppm)')
    
    if custom_xticks == False:
        plt.xticks(x_ticks,x_ticks_int)
    
    # adds noise to new spectrum, noisy_spectrum
    if noise != 0:
        noisy_spectrum = add_noise(spectrum,noise)
    
    # plots background if True
    if plot_background == True:
        plt.plot(nc.c/frequencies/1e-4, background_spectrum, label = 'Background', linewidth=1)
        
    # plots noisy signal if True
    if plot_noisy == True:
        plt.errorbar(nc.c/frequencies/1e-4, noisy_spectrum, yerr=40, label = 'With noisy signal', linewidth=1)
        
    # plots signal
    plt.plot(nc.c/frequencies/1e-4, spectrum, label = plot_label, linewidth=1)
    
    plt.legend()
    # saves plot
    plt.savefig(filename)
    if clear_plot == True:
        plt.clf()

# ----- Function to retrieve specific quantities of prebiosignature molecule -----
# band name - band name for identification see Claringbold et al. 2022
# background_atmosphere - background atmospheric and planetary properties in form of exoPlanet object (see top)
# prebio_species - species added to the background
# test_species - species being retrieved for
# wlen_mix, wlen_max - wavelength range of band (in microns)
# noise - noise of instrument in wavelength range (in ppm)
# prebio_abundance - mass fraction of prebio_species added to the background
# spectral_resolution - spectral resolution to evaluate at (currently choice of 100 to 1000)
# plot - determine whether to save a plot of the posterior PDF

def perform_retrieval(band_name,
                      background_atmosphere,
                      prebio_spec, test_spec,
                      wlen_min, wlen_max, 
                      noise, 
                      prebio_abundance, 
                      spectral_resolution=1000, 
                      plot=True):
    
    back, background_spectrum, test_species, prebio_species = setup_atmosphere(background_atmosphere,wlen_min,wlen_max,spectral_resolution,test_spec,prebio_spec)
    
    #generate test spectra
    signal_samples=111
    signal_list, test_abundances = generate_test_spectra(back,background_spectrum,test_species,signal_samples)
    
    
    #Set retrieved species to 0, and recalculate MMW.
    back.comp[test_species]= 0 * np.ones_like(back.temps)
    
    print('Retrieving '+test_species+'...')   
    
    #perform Bayesian analysis between model and test spectra
    unnorm_fits_mean = bayesian_analysis(prebio_abundance,back,background_spectrum,prebio_species,noise,signal_list,signal_samples)
    fits_mean = unnorm_fits_mean/sum(unnorm_fits_mean)
    
    #plot the posterior pdf
    if plot==True:
        plt.plot(test_abundances, fits_mean)
        plt.xscale('log')
        plt.xlabel('Mass Fraction of '+prebio_species)
        plt.ylabel('Bayesian likelihood')
        plt.savefig('Retrieval.pdf')
    
    mean_abundance=0
    std_abundance=0
    var_abundance=0
    log_test_abundances=np.log10(test_abundances)
    
    
    #calculate the mean retrieved abundance
    for i in range(0,signal_samples-1):
        mean_abundance = mean_abundance + (fits_mean[i] * log_test_abundances[i])
    mean_abundance_lin = 10**mean_abundance
    print("    Retrieval complete.")
    
    print("mean retrieved abundance = " + str(mean_abundance_lin))
    print("mean retrieved abundance (log10) = " + str(mean_abundance))
    

    #std of retrieved abundance
    for i in range(0,signal_samples-1):
        var_abundance = var_abundance + (fits_mean[i] * (log_test_abundances[i] - mean_abundance)**2)
    std_abundance = np.sqrt(var_abundance)
    print("std of retrieved abundance (log10) = " + str(std_abundance))

    
    return test_abundances, fits_mean

# ----- Function to calculate the detectability threshold -----
# band_name - band name for identification, see Claringbold et al. 2022
# background_atmosphere - background atmospheric and planetary properties in form of exoPlanet object (see top)
# prebio_species - species added to the background in increasing amounts until detected
# test_species - species being retrieved for
# wlen_mix, wlen_max - wavelength range of band (in microns)
# noise - noise of instrument in wavelength range (in ppm)
# min_abundance - minimum log10(mass_fraction) to start testing at (must be below detection threshold). Will run quicker if this value is higher..
# spectral_resolution - spectral resolution to evaluate at (currently choice of 100 and 1000)
# output_file - .txt file to open and write detection threshold in form: band_name, test_species, prebio_species, detection threshold mass fraction, mean retrieved mass fraction of prebio_species, MMW


def assess_detectability(band_name,
                         background_atmosphere,
                         prebio_spec, test_spec,
                         wlen_min, wlen_max,
                         noise,
                         min_abundance=-6,
                         spectral_resolution=1000,
                         output_file='Results.txt',
                         **spectrum_kwargs):
    
    back_atm, background_spectrum, test_species, prebio_species = setup_atmosphere(background_atmosphere,wlen_min,wlen_max,spectral_resolution,test_spec,prebio_spec,**spectrum_kwargs)
    back = back_atm
    #Generate the hypothesis spectra
    signal_samples=111
    signal_list, test_abundances = generate_test_spectra(back,background_spectrum,test_species,signal_samples,**spectrum_kwargs) 
    
    #Set retrieved species to 0
    back.comp[test_species] = 0 * np.ones_like(back.temps)
    
    #recalculate mean molecular weight (MMW)
    #back.planet.MMW = recalc_MMW(background_atmosphere,test_species)
      
    print('Finding detection threshold of '+prebio_spec)
    
    detection_significance_sum=0
    
    #iterate the Bayesian analysis, performing goodness of fit between the test spectra and a model spectrum, calculating the significance of the detection, and increasing the amount of prebiosignature until a significant detection
    j=0
    no_detection=False
    #iterate until the detection significance is three sigma
    while detection_significance_sum < 0.997:
        
        #define end condition if never a significant detection
        max_j=(-min_abundance)*10 + 1
        
        #choose list in log space to test abundances of the prebiosignature in the model spectrum
        exploration_list = np.logspace(min_abundance,0,(-min_abundance)*10 + 1)
        prebio_signal = exploration_list[j]
        j = j + 1
        
        #perform bayesian analysis on the model spectrum and return posterior PDF
        fits_mean = bayesian_analysis(prebio_signal, back, background_spectrum, prebio_species, noise, signal_list, signal_samples, **spectrum_kwargs)
        
        #create empty variables
        new_fits_mean = []
        new_test_abundances = []
        mean_abundance=0
        
        
        log_test_abundances=np.log10(test_abundances)
        plt.plot(test_abundances,fits_mean/sum(fits_mean))
        plt.xscale('log')
        
        #calculate the mean retrieved abundance
        for i in range(0,signal_samples-1):
            mean_abundance = mean_abundance + (fits_mean[i] * log_test_abundances[i])
        mean_abundance_lin = 10**mean_abundance
        
        fits_mean_size = fits_mean.size
        
        #sum the posterior probability for every abundance which retrieves 1% of the input or more to obtain  "detection significance"
        
        #select only the signals which retrieve 1% of the input
        for i in range(0,signal_samples-1):
            if test_abundances[i] > prebio_signal/100:
                new_fits_mean.append(fits_mean[i])
                new_test_abundances.append(test_abundances[i])
        
        #sum all of the successful retrievals
        detection_significance_sum = sum(new_fits_mean)/sum(fits_mean)
        #print(detection_significance_sum)
        
        sum_elements=0
        
        #ends the loop if the detection significance is greater than 99.7% (3 sigma)
        for i in range(fits_mean_size):
            sum_elements  = sum_elements + fits_mean[i]
            if sum_elements > 0.997:
                break
        
        print(prebio_signal)
        
        
        #ends the loop is a pure atmosphere is tested and failed to retrieve
        if j == max_j:
            print("No detection.")
            no_detection=True
            break
    
    #prints the results of the detectability assessment.
    
    if no_detection == True:
        print(band_name + "\t" + test_species + "\t failed to detect")
        
    
    else:
        print(band_name + "\t" + test_species + "\t complete")
        output = [band_name,prebio_species,test_species,prebio_signal,mean_abundance_lin,back.planet.MMW]
        with open(output_file, 'a') as f:
            f.write('\n')
            for item in output:
                f.write('%s\t' % item)
    
    return prebio_signal


# ----- Function to find detection threshold for a spectral band -----
# band_name - label of band to run (should be key of band_list dictionary in band_input_data)
# atmosphere_input_data - background atmosphere input data file (e.g. Hycean.py) should be imported beforehand
# band_input_data - band input data file (e.g. prebiosignatures.py) should be imported beforehand
# spectral_resolution - specify spectral_resolution (opacity data must be prepared beforehand)
# retrieved_species - specify which species to find detection tresholds for, or 'All' to retrieve all species with features in the band
# JWST_noise - True to use JWST noise from .txt files, False to use defaults specified in band_input_data
# output_file - .txt file to save results to
# auto_assign - whether to automatically assign each band an instrument based on wavelength or spectral resolution
# (if False, must specify instrument[band] beforehand)
# min_abundance - minimum mass fraction to test for detection threshold (see assess_detectability)

def run_band(band_name,
             atmosphere_input_data,
             band_input_data,
             spectral_resolution=1000,
             retrieved_species='All',
             JWST_noise=True,
             output_file='Results.txt',
             auto_assign=True,
             min_abundance=-6,
             use_prism=True,
             **spectrum_kwargs):
    
    # reassign to shorter variable name
    atm = atmosphere_input_data
    
    # load in band_list_data
    band_list = compile_band_list(band_input_data,
                                  JWST_noise=JWST_noise,
                                  auto_assign=auto_assign,
                                  spectral_resolution=spectral_resolution,
                                  use_prism=use_prism)
    
    # pick out specified band
    print('Finding detection thresholds in band '+band_name+'...')
    band = band_list[band_name]
    
    # retireves all species (default)
    if retrieved_species == 'All':
        print('    Retrieving all species with features in band')
        
        # loop over all species in band
        for spec in band.species:   
            
            # load in atmosphere input data
            atmos = load_planet(atm)
            
            # assess detectability using data from band
            assess_detectability(band.label,
                                 atmos,
                                 spec, spec,
                                 band.min_wlen, band.max_wlen,
                                 band.precision,
                                 spectral_resolution=spectral_resolution,
                                 output_file=output_file,
                                 min_abundance=min_abundance,
                                 **spectrum_kwargs)
    
    # retrieve single specified species
    else:
        print('    Retrieving only ' + retrieved_species)
        atmos = load_planet(atm)
        assess_detectability(band.label,
                             atmos,
                             retrieved_species, retrieved_species,
                             band.min_wlen, band.max_wlen,
                             band.precision,
                             spectral_resolution=spectral_resolution,
                             output_file=output_file,
                             min_abundance=min_abundance,
                             **spectrum_kwargs)

# ----- Function to find detection threshold for all spectral band -----
# atmosphere_input_data - background atmosphere input data file (e.g. Hycean.py) should be imported beforehand
# band_input_data - band input data file (e.g. prebiosignatures.py) should be imported beforehand
# spectral_resolution - specify spectral_resolution (opacity data must be prepared beforehand)
# JWST_noise - True to use JWST noise from .txt files, False to use defaults specified in band_input_data
# output_file - .txt file to save results to
# auto_assign - whether to automatically assign each band an instrument based on wavelength or spectral resolution
# (if False, must specify instrument[band] beforehand)
# min_abundance - minimum mass fraction to test for detection threshold (see assess_detectability)

def run_atmosphere(atmosphere_input_data,
                   band_input_data,
                   spectral_resolution=1000,
                   JWST_noise=True,
                   output_file='Results.txt',
                   auto_assign=True,
                   min_abundance=-6,
                   use_prism=True,
                   **spectrum_kwargs):
    
    print('Finding detection thresholds in: '+str(atmosphere_input_data)+'.py')
    print('using band list data from: '+str(band_input_data)+'.py')
    
    # reassign to shorter variable name    
    atm = atmosphere_input_data
    
    # load in band_list input data
    band_list = compile_band_list(band_input_data,
                                  JWST_noise=JWST_noise,
                                  auto_assign=auto_assign,
                                  spectral_resolution=spectral_resolution,
                                  use_prism=use_prism)
    print('Bands to run:')
    print(band_list.keys())
    
    # iterate over bands
    for band in band_list.values():
        print('Finding detection thresholds in band '+band.label+'...')
        
        # iterate over species in bands
        for spec in band.species:
            print('    Retrieving ' + spec)
            
            # load in atmosphere
            atmos = load_planet(atm)
            
            # assess detectability using data from band
            assess_detectability(band.label,
                                 atmos,
                                 spec, spec,
                                 band.min_wlen, band.max_wlen,
                                 band.precision,
                                 spectral_resolution=spectral_resolution,
                                 output_file=output_file,
                                 min_abundance=min_abundance,
                                 **spectrum_kwargs)
    
    print(str(atmosphere_input_data)+' detection thresholds computed successfully!')
    print('Find results in '+output_file)
    

    
    
    
    
