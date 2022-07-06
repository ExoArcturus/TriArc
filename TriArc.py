#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 11:39:42 2022

@author: arcturus
"""


import numpy as np
import pylab as plt

#import constants from pRT (not used currently)
from petitRADTRANS import nat_cst as nc
#import main package from pRT
from petitRADTRANS import Radtrans

#import molecular weights of species
import molecular_weights as mw

# ----- CLASSES -----

# define exoPlanet class, with all inputs for exoplanetary system

class exoPlanet:
    def __init__(self, R_pl, R_st, gravity, P0, temp, composition, lin_spec, ray_spec, cont_spec, MMW, mixing_ratios,species_l):
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
  
    
# ----- INITIALISATION FUNCTIONS FOR IMPORTING INPUT FILES -----
# load in input file (e.g. Hycean.py) which should be imported in script utilising TriArc, and creates exoPlanet object.
# also performs conversion to mass fraction and calculates mean molecular weight.

""" input files must include the following:
    
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

def load_planet(atm):
    
    #import atmospheric composition (in mixing ratios)
    n = atm.abundances
    
    #calculate mean molecular weight (MMW)
    MMW = 0
    species_list=[]
    species_list = atm.species_list
    #remove duplicates from species list using dictionary
    species_list = list(dict.fromkeys(species_list))
    for species in species_list:
        MMW = MMW + mw.mass[species] * n[species]
    print(MMW)
    
    #convert from abundances to mass fractions
    mass_fraction={}
    for species in species_list:
        mf =  n[species] * ( mw.mass[species] / MMW)
        mass_fraction[species]=mf
    #print(mass_fraction)
    generic_planet = exoPlanet(atm.planet_radius,atm.star_radius,atm.surface_gravity,atm.surface_pressure,atm.iso_temperature,mass_fraction,atm.line_species,atm.rayleigh_species,atm.cont_species,MMW,n,species_list)
    return generic_planet

# ----- Intermediate functions, used by the high-level functions (e.g. perform retrieval, access detectability) -----


#Function to calculate the transmission spectrum for a given composition (comp), exoPlanet object (planet), temperature profile (temps), Radtrans object (model). 

def transmission_spectrum(comp, planet, temps, model):
    #input uniform mixing ratio in column
    mean_mol_weight = planet.MMW*np.ones_like(temps)
    
    #use pRT function to calculate transmission spectrum
    model.calc_transm(temps, comp, planet.gravity, mean_mol_weight, R_pl=planet.R_pl, P0_bar=planet.P0)
    
    #return transmission spectrum in ( transmission radius / star radius ) squared, in ppm
    return ((model.transm_rad/planet.R_st)**2)*10**6

#Function to add noise to a transmission spectrum.

def add_noise(model_spec, noise_sd):
    # add simple gaussian noise to all data points in spectrum
    noisy_signal = np.random.normal(model_spec,noise_sd)
    return noisy_signal

#Function to estimate the goodness of fit between a signal (noisy_signal) and a test signal (signal) using a radial basis likelihood function.

def goodness_of_fit(noisy_signal, signal, sigma):
    # calculate value of normalised gaussian least squares at each wavelength
    gauss_least_squares = ((2 * np.pi * sigma**2)**(-1/2)) * np.exp(-(1/(2*(sigma**2))*(noisy_signal-signal)**2))
    # find product of normalised gaussian least squares to calculate likelihood function
    likelihood = np.prod(gauss_least_squares)
    return likelihood

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

def setup_atmosphere(background_atmosphere,wlen_min,wlen_max,spectral_resolution=1000,test_species='CO',prebio_species='CO'):
    #Import background atmosphere object, with arguments R_pl (radius of planet), R_st (radius of star), P0 (surface pressure), gravity, temp, composition (mass fractions), lin_spec, ray_spec, cont_spec, MMW (mean molecular weight), n (vol fractions)
    atm = background_atmosphere
    
    
    #Instructs the use of opacities if rebinned to lower spectral resolution (must be prepared beforehand)
    if spectral_resolution != 1000:
        test_species, prebio_species, atm.species_l, atm.lin_spec, atm.composition = select_spectral_resolution(test_species,prebio_species,background_atmosphere,spectral_resolution)
    
        
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
    background_spectrum = transmission_spectrum(mass_fractions, atm, temperature, model_atmosphere)
    
    return back, background_spectrum, test_species, prebio_species
    

#Function to generate test transmission spectra (the 'hypotheses' in Bayes' theorem) for different abundances to test against.

def generate_test_spectra(back_atm,background_spec,test_species,signal_samples=111):
    
    #Choose uniform prior in log space to apply sampling over
    test_abundances = np.logspace(-11,0,num=signal_samples)
    
    #Calculate transmission spectrum for each test abundance sampled, remove background, and add to list.
    signal_list=[]
    for abundance_species in test_abundances:
        back_atm.comp[test_species] = abundance_species * np.ones_like(back_atm.temps)
        signal_bg = transmission_spectrum(back_atm.comp, back_atm.planet, back_atm.temps, back_atm.model)
        signal = signal_bg - background_spec
        signal_list.append(signal)
    return signal_list, test_abundances

#Primary function, injects prebiosignature into background, performs Bayesian analysis against set of test spectra, to output posterior PDF.
#Uses add_noise and goodness_of_fit functions

def bayesian_analysis(prebio_signal,back_atm,background_spec,prebio_species,noise,signal_list,signal_samples=111):
    
    #Adds specified amount of prebiosignature molecule
    back_atm.comp[prebio_species] = prebio_signal * np.ones_like(back_atm.temps)
    
    #Calculate a model transmission spectrum
    model_bg = transmission_spectrum(back_atm.comp, back_atm.planet, back_atm.temps, back_atm.model)
    
    #Subtract the background
    model =  model_bg - background_spec
    
    #Adds noise to the transmission spectrum, and calculates goodness of fit against the test spectra.
    
    repeats = 9000
    fits_list=np.empty([repeats,signal_samples])
    sigma=noise
    
    for i in range(0,repeats):
        all_fits=[]
        noisy_signal = add_noise(model, noise)
        for test_signal in signal_list:
            fit = goodness_of_fit(noisy_signal, test_signal, sigma)
            all_fits.append(fit)
        normalisation=sum(all_fits)
        all_fits = all_fits / normalisation
        fits_list[i]=all_fits
        
    #Calculates the mean fit from all the fits for each noise profile, and plots it.
    
    sum_fits = np.sum(fits_list,axis=0)
    
    #This is the posterior PDF of the bayesian analysis
    fits_mean = sum_fits/repeats
    return fits_mean

# ----- HIGH-LEVEL FUNCTIONS FOR RETRIEVALS -----

# ----- Function to retrieve specific quantities of prebiosignature molecule -----
# band name - band name for identification see Claringbold et al. 2022
# background_atmosphere - background atmospheric and planetary properties in form of exoPlanet object (see top)
# prebio_species - species added to the background
# test_species - species being retrieved for
# wlen_mix, wlen_max - wavelength range of band (in microns)
# noise - noise of instrument in wavelength range (in ppm)
# prebio_abundance - mass fraction of prebio_species added to the background
# spectral_resolution - spectral resolution to evaluate at (currently choice of 100 to 1000)


def perform_retrieval(band_name,background_atmosphere,prebio_spec,test_spec,wlen_min,wlen_max,noise,prebio_abundance,spectral_resolution=1000):
    
    back, background_spectrum, test_species, prebio_species = setup_atmosphere(background_atmosphere,wlen_min,wlen_max,spectral_resolution,test_spec,prebio_spec)
    #generate test spectra
    signal_samples=111
    signal_list, test_abundances = generate_test_spectra(back,background_spectrum,test_species,signal_samples)
    
    
    #Set retrieved species to 0, and recalculate MMW.
    back.comp[test_species]= 0 * np.ones_like(back.temps)
    
    #recalculate mean molecular weight (MMW)
    #back.planet.MMW = recalc_MMW(background_atmosphere,test_species)
    
    #perform Bayesian analysis between model and test spectra
    fits_mean = bayesian_analysis(prebio_abundance,back,background_spectrum,prebio_species,noise,signal_list,signal_samples)
    
    #plot the posterior pdf
    plt.plot(test_abundances, fits_mean)
    plt.xscale('log')
    plt.xlabel('Abundance of HCN')
    plt.ylabel('Bayesian likelihood')
    
    mean_abundance=0
    std_abundance=0
    var_abundance=0
    log_test_abundances=np.log10(test_abundances)
    
    
    #calculate the mean retrieved abundance
    for i in range(0,signal_samples-1):
        mean_abundance = mean_abundance + (fits_mean[i] * log_test_abundances[i])
    mean_abundance_lin = 10**mean_abundance
    print("mean retrieved abundance = " + str(mean_abundance_lin))
    print("mean retrieved abundance (log10) = " + str(mean_abundance))
    

    #std of retrieved abundance
    for i in range(0,signal_samples-1):
        var_abundance = var_abundance + (fits_mean[i] * (log_test_abundances[i] - mean_abundance)**2)
    std_abundance = np.sqrt(var_abundance)
    print("std of retrieved abundance (log10) = " + str(std_abundance))

    
    return fits_mean





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


def assess_detectability(band_name,background_atmosphere,prebio_spec,test_spec,wlen_min,wlen_max,noise,min_abundance=-6,spectral_resolution=1000,output_file='Results.txt'):
    
    back, background_spectrum, test_species, prebio_species = setup_atmosphere(background_atmosphere,wlen_min,wlen_max,spectral_resolution,test_spec,prebio_spec)
    
    #Generate the hypothesis spectra
    signal_samples=111
    signal_list, test_abundances = generate_test_spectra(back,background_spectrum,test_species,signal_samples) 
    
    #Set retrieved species to 0
    back.comp[test_species] = 0 * np.ones_like(back.temps)
    
    #recalculate mean molecular weight (MMW)
    #back.planet.MMW = recalc_MMW(background_atmosphere,test_species)
      
    
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
        fits_mean = bayesian_analysis(prebio_signal, back, background_spectrum, prebio_species, noise, signal_list, signal_samples)
        
        #create empty variables
        new_fits_mean = []
        new_test_abundances = []
        mean_abundance=0
        
        
        log_test_abundances=np.log10(test_abundances)
        
        
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
        detection_significance_sum = sum(new_fits_mean)
        
        
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
        with open('Results.txt', 'a') as f:
            f.write('\n')
            for item in output:
                f.write('%s\t' % item)

    

    
    
    
    
