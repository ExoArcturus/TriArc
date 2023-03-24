#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 22:18:34 2022

@author: alastair
"""

# MaTriArc (Modelling and reTRIeval: ARCturus) is the primary engine behind TriArc
# Defines the primary object, ExoPlanetAtmosphere, and the setup, forward modelling and retrieval methods
# Powered by petitRADTRANS for forward modelling (including parametric p-T profiles and equilibrium chemistry)
# Uses custom retrieval method for single-parameter retrievals


# import the usual
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#  ---- Using peititRADTRANS (pRT) for forward modelling ----

# import physics module from petitRADTRANS for analytic p-T profiles
import petitRADTRANS.physics as physics

# import chemistry module from pRT for equilibrium chemsitry
import petitRADTRANS.poor_mans_nonequ_chem as chemistry


# import main object from pRT for radiative transfer calculations
from petitRADTRANS import Radtrans

# import constants from pRT to convert from astronomical to cgs units
from petitRADTRANS import nat_cst as nc

# import custom-defined molecular weights for calculating MMW and mass fractions
import molecular_weights as mw


# ---- GENERIC CALCULATION FUNCTIONS ----

# Function to simply calculate the equilibrium temperature of a planet given:
# orbital distance in au
# stellar radius in solar radii
# stellar mass in solar masses
# stellar temperature in kelvin
# bond albedo
# assuming perfect day-night energy redistribution and blackbody radiation


def compute_equilibrium_temperature(orbital_distance, stellar_radius, stellar_mass, stellar_temperature, albedo):
    
    # use Kepler's law to calculate orbital distance from period - DEPRECATED
    
    # orbital_distance = (((period/365)**2) / stellar_mass)**(1/3)
    
    # calculate equilibrium temperature assuming blackbody radiation
    
    equilibrium_temperature = stellar_temperature * (( (1-albedo) * (stellar_radius / (2 * orbital_distance * 215))**2 ))**(1/4)  
    
    return equilibrium_temperature

# Function to add uncorrelated Gaussian noise to a numpy array.
# Conventional use adds noise (standard deviation) to a transmission spectrum (model_spectrum) in ppm

def add_noise(model_spectrum, noise_sd):    
    noisy_model = np.random.normal(model_spectrum, noise_sd)
    return noisy_model
        
        
# Function to estimate the goodness of fit between a signal (noisy_signal) and a test signal (signal) using a radial basis likelihood function.
# sigma should be set to the noise.

def goodness_of_fit(noisy_signal, signal, sigma):
    
    # calculate value of normalised gaussian least squares at each wavelength
    gauss_least_squares = ((2 * np.pi * sigma**2)**(-1/2)) * np.exp(-(1/(2)*((noisy_signal-signal)/sigma)**2))
    # find product of normalised gaussian least squares to calculate likelihood function
    likelihood = np.prod(gauss_least_squares)
    
    return likelihood

# Function to relablel a species with the appropriate suffix for a pRT species at rebinned sepctral resolution 

def label_species_with_resolution(species,spectral_resolution):
    species = species + '_R_' + str(spectral_resolution)
    return species    
    
# Function to obtain Bayesian likelihood for each of a set of hypothesis spectra (test_spectra) given a noisy model spectrum as evidence, and noise in ppm.

def retrieval_fitting(test_spectra, model_spectrum, noise):
    
    # initalise empty list for mean likelihood per spectrum
    
    fit_mean_list=[]
    
    # iterate over each hypothesis spectrum
    
    for test in test_spectra:
        
        # initialise empty list for the likelihood at each wavelength
        
        fit_list=[]
        
        # iterate over each wavelength point in spectrum
        
        for i, wlen in enumerate(test):
            #fit = np.exp( -( (test/noise[i]) - (model_spectrum[i]/noise[i]) )**2) # deprecated goodness of fit calculation
            
            # calculate likelihood for each wavelength point
            fit = goodness_of_fit(model_spectrum[i], test, noise[i])
            # add likelihood to list to compute mean
            fit_list.append(fit)
            
        # compute mean likelihood per spectrum
        fit_mean = np.mean(fit_list)
        
        # add mean likelihood to list
        fit_mean_list.append(fit_mean)
    
    # normalise list of mean likelihood
    normalisation = np.sum(fit_mean_list)
    posterior_pdf = fit_mean_list / normalisation
    
    # return posterior pdfs as mean likelihood as a function of parameter used to generate test spectra
    return posterior_pdf

# ---- CLASSES ----

# ExoPlanetAtmosphere class is for the main object, storing all parameters relevant to the planet atmosphere, including pRT Radtrans object
# Parameters are initialised sequentially, with dedicated methods to define and safely redfine parameters
# Key methods act on this object, including calculating transmission spectra and performing retrievals

class ExoPlanetAtmosphere:
    
    # ExoPlanetAtmosphere object is initialised with key planetary parameters, in cgs units
    
    def __init__(self, R_planet, R_star, gravity, T_eq):
        self.R_planet = R_planet            # radius of planet surface (defined as radius of pressure layer equal to P0, by default 10 bar)
        self.R_star = R_star                # radius of star
        self.gravity = gravity              # gravity at planet surface
        self.T_eq = T_eq                    # equilibrium temperature of planet
        self.abundances = {}                # initialise empty dictionary for abundances
        self.mass_fractions = {}            # initialise empty dictionary for mass fractions
    
    # ----- MACRO SETUP FUNCTION -----
    # macro function to perform all necessary setup steps at once
    # steps are performed in proper order to simplify creation of ExoPlanetAtmosphere object
    # should be performed immediately after initialising ExoPlanetAtmosphere object
    # all necessary inputs are collected

    # inputs:
        # opacity sources (line species, rayleigh species, continuum species)
        # wavelength bounds (wlen_min and wlen_max) in microns
        # spectral resolution, default R=1000
        
        # P0, surface presure in bar (default 10 bar)
        # pressure profile in bar:
            # 'Default' - 100 layers equally log-spaced from the surface (P0) to 1 microbar, in bars
            # user-defined log-spaced array
        
        # temperature profile (one of the following):
            # 'Isothermal' - default, no further input
            # 'Guillot' - with user-defined or default T_int, gamma, and kappa_IR
            # 'Functional'- also requires T(P) function defined in temperature_function, scaled with equilibrium temperature
            # User-defined array, should be same dimensions as pressure profile
            
        # composition (one of the following):
            # chemical_equilibrium=True, with user-defined or default COratio and metallicity, and optionally quench pressure in bar
            
            # abundances is a user-defined dictionary with opacity species as keys (all opacity species need a defined abundance)
            # the vertical profile also needs to be specified:
                # 'Uniform' - default, no further input
                # User-defined dictionary of arrays, with same dimensions as pressure profile, with species as keys

        # **kwargs - used to add cloud or haze particles to transmission spectrum (see Raddtrans.calc_transm kwargs)
        
        # plot - if True will plot transmission spectrum

    def setup_atmosphere(self,
                         line_species,
                         ray_species,
                         cont_species,
                         wlen_min,
                         wlen_max,
                         P0=10,
                         spectral_resolution=1000,
                         abundances=None,
                         kappa_IR=0.01,
                         gamma=0.5,
                         T_int=200,
                         temperature_function=None,
                         other_species=None,
                         COratio=0.55,
                         metallicity=0,
                         Pquench=None,
                         pressure_profile='Automatic', 
                         temperature_profile='Isothermal', 
                         chemical_equilibrium=True, 
                         vertical_mixing_ratios='Uniform',
                         plot=True,
                         **kwargs):
        
        # first determine all opacity sources
        # line_species (corr-k line opacities), ray_species (rayleigh continuum opacities), cont_species (CIA continuum opacities)
        
        self.set_opacity_sources(line_species,ray_species,cont_species, other_species=None)
        
        # obtain necessary inputs for Guillot or Functional temperature profiles
        
        if temperature_profile=='Guillot':
            self.set_guillot_temp_profile(kappa_IR,gamma,T_int)
        elif temperature_profile=='Functional':
            self.set_functional_temp_profile(temperature_function)
            
        # create atmospheric structure profile (with defined surface pressure, pressure profile, temperature profile)
        
        self.generate_atmospheric_structure(P0=P0,pressure_profile=pressure_profile,temperature_profile=temperature_profile)
        
        # obtain necessary inputs for composition or chemical equilibrium
        
        if chemical_equilibrium==False:
            self.set_abundances(abundances)
        else:
            self.set_equilibrium_chemistry(COratio,metallicity,Pquench=Pquench)
        
        # create atmospheric composition
        
        self.generate_composition(equilibrium_chemistry=chemical_equilibrium,vertical_mixing_ratios=vertical_mixing_ratios)

        # create opacities and radiative transfer object

        self.generate_opacities(wlen_min,wlen_max,spectral_resolution=spectral_resolution)

        # creates a transmission spectrum to use as a background, and plots it if plot=True

        self.background_spectrum, self.wavelengths = self.calculate_transmission_spectrum(plot=plot,**kwargs)
    
    
    # ---- INDIVIDUAL SETUP FUNCTIONS ----
    
    # these functions are used to perform initial setup and provides inputs to ExoPlanetAtmosphere to allow modelling and retrieval applications
    # setup functions are order-sensitive, if in doubt, use the macro function setup_atmosphere() above
    # always use these functions to perform initial inputs rather than manually define variables
    # use edit object functions (below) instead of using these to edit inputs, or reinitialise the object
    
    # method to input opacity sources to ExoPlanetAtmosphere object as lists of species strings
    
    def set_opacity_sources(self, line_species, ray_species, cont_species, other_species=None):
        self.line_species =  line_species           # list of species contributing to (corr-k) line opacity
        self.ray_species = ray_species              # list of species contributing to rayleigh opacity
        self.cont_species = cont_species            # list of species pairs contributing to CIA (collisionally-induced-absorption) opacities
        
        # create list of opacity species that require abundances
        species_list = line_species + ray_species   # combine line species and rayleigh species
        if other_species is not None: 
            species_list += other_species           # opportunity to add additional species not contriubting to opacity
        self.species_list = [*set(species_list)]    # remove duplicates from species list
        #print(self.species_list)            
        
    # method to input vertically-averaged mixing ratios (abundances) if not using chemical equilibrium to calculate composition
    # abundances should be a dictionary with relevant opacity species as keys
    # this method should be called upon to input all major atmospheric components, as this is when mean molecular weight is calculated
    # abundances are converted into mass fractions
    
    def set_abundances(self, abundances):
        # iterate over species list (defined in set_opacity_sources)
        for species in self.species_list:
            self.abundances[species] = abundances[species]  # input relevant abundance
            
        # calculate mean molecular weight (MMW)
        self.calculate_MMW()
        
        # convert abundances to mass fractions
        self.compute_mass_fractions()
        
    # method to recalculate mean molecular weight
    def calculate_MMW(self):
        
        # resets MMW to 0
        self.MMW = 0
        
        # iterate over species list
        for species in self.species_list:
            
            # avoids recounting rebinned opacities
            if '_R_' not in species:
                # adds species contribution to MMW
                self.MMW += mw.mass[species] * self.abundances[species]
            
    # method to convert abundances into mass fractions
    # by default will apply to all species, 
    # unless single_species argument is specified, in which case it only apply to argument of single_species
    
    def compute_mass_fractions(self, single_species=False):
        
        # default setting, convet all abundances into mass fractions
        if single_species == False:
            # iterate over species in species list
            for species in self.species_list:
                # compute mass fraction for species as save to dictionary
                mf = self.abundances[species] * (mw.mass[species] / self.MMW)
                self.mass_fractions[species] = mf
        
        # if single_species specified, calculates mass fraction of the species taking the single_species argument 
        else:
            # computes mass fraction as saves to dictionary
            mf = self.abundances[single_species] * (mw.mass[single_species] / self.MMW)
            self.mass_fractions[single_species] = mf
    
    # method to take inputs for calculation of analytic p-T following the two-grey-stream approach of Guillot (2010)
    
    def set_guillot_temp_profile(self, kappa_IR, gamma, T_int):
        self.kappa_IR = kappa_IR    # opacity in the infrared stream
        self.gamma = gamma          # ratio of optical opacity to infrared opacity
        self.T_int = T_int          # planet internal temperature (from internal flux boundary condition)
    
    # method to take input of a function to calculate temperature profile from pressure profile, 
    # i.e. T(P) function scaled with equilibrium temperature 
    
    def set_functional_temp_profile(self, temperature_function):
        self.temperature_function = temperature_function
    
    # method to compute p-T profile, given a surface pressure
    # pressure profile can be automatically calculated ('Automatic') or user-defined (a log-spaced array)
    # temperature profile can be 'Isothermal', Analytic ('Guillot'), 
    # Functional (a user-defined function T(P), or user-defined (array with dimension equal to pressure profile)
    
    def generate_atmospheric_structure(self, P0=10, pressure_profile='Automatic', temperature_profile='Isothermal'):
        
        # define pressure in bar at radius of planet surface
        self.P0 = P0
        
        # uses default pressure_profile of one hundred layers equally spaced in log space between 1 microbar and P0
        if pressure_profile == 'Automatic':
            self.pressures = np.logspace(-6,P0,100)
        # otherwise uses user defined pressure profile, using argument of pressure_profile
        else:
            self.pressures = pressure_profile
        
        # if temperature profile set to 'Isothermal' (default), creates uniform array of equilibrium temperature with same dimensions as pressure profile
        if temperature_profile == 'Isothermal':
            self.temperatures = self.T_eq * np.ones_like(self.pressures)
            
        # if temperature profile set to 'Guillot', analytically computes temperature profile given inputs set in set_guillot_temp_profile
        elif temperature_profile == 'Guillot':
            self.temperatures = physics.guillot_global(self.pressures, self.kappa_IR, self.gamma, self.gravity, self.T_int, self.T_eq)
        
        # if temperature profile set to 'Functional', calculates temperature at each pressure layer with user-defined P_T set in set_func_temp_profile
        elif temperature_profile == 'Functional':
            
            # iterate over each pressure layer
            for i, pressure_layer in enumerate(self.pressures):
                # calculate temperature from T(P) profile at each layer
                self.temperatures[i] = self.T_eq * self.temperature_function(pressure_layer)
        
        # otherwise argument of temperature profile is used (must be same dimension as pressure profile)
        else:
            self.temperatures = temperature_profile
    
    # method to take inputs for equilibrium chemistry calculation
    
    def set_equilibrium_chemistry(self, COratio, metallicity, Pquench=None):
        self.COratio = COratio * np.ones_like(self.pressures)           # C:O ratio of atmosphere (assumed vertically constant)
        self.metallicity = metallicity * np.ones_like(self.pressures)   # [Fe/H] relative to solar of atmosphere (assumed vertically constant)
        self.Pquench = Pquench                                          # pressure at which we assume quenching to occur (vetical mixing timescale equals chemical timescale)
        # quench pressure can be ignored (default)
        
    # method to generate complete atmospheric composition (vertical mixing profiles)
    # if equilibrium_chemistry is set to True, will calculate the vertical mixing profiles will equilibrium chemistry
    # (in this case ignore vertical_mixing_ratios argument)
    
    # otherwise need to use set_abundances, and by default this will establish uniform vertical mixing ratios
    # vertical mixing profiles can be defined (in terms of the pressure_profile) manually as the input to the vertical_mixing_ratios argument
    
    def generate_composition(self, equilibrium_chemistry=False,vertical_mixing_ratios='Uniform'):
        
        # initialise empty dictionary for mass fraction vertical profile
        self.mf_profile={}
        
        # if equilibrium chemistry is true, uses pRT chemistry package to interpolate equilibrium abundances at each pressure layer
        if equilibrium_chemistry == True:
            if self.Pquench is None:
                self.mf_profile = chemistry.interpol_abundances(self.COratio, self.metallicity, self.temperatures, self.pressures)
            else:
                self.mf_profile = chemistry.interpol_abundances(self.COratio, self.metallicity, self.temperatures, self.pressures, Pquench_carbon=self.Pquench)
            
            # calculates average MMW from vertical MMW profile
            self.MMW = np.mean(self.mf_profile['MMW'])
        
        # otherwise uses inputs to calculate abundances
        else:
            # by default the vertical mixing profile is assumed to be uniform
            if vertical_mixing_ratios=='Uniform':
                
                # iterate over all species
                for species in self.species_list:
                    
                    # mass fraction profile is uniform with dimensions of pressure_profile
                    self.mf_profile[species] = self.mass_fractions[species] * np.ones_like(self.pressures)
                    # MMW profile is therefore also uniform
                    self.mf_profile['MMW'] = self.MMW * np.ones_like(self.pressures)
            
            # if vertical mixing ratios is given an argument, it will be used to calculate vertical mass fraction profile 
            else:
                
                # iterate over all species
                for species in self.species_list:
                    # create MMW vertical mixing profile by calculating at each layer
                    MMW_profile = np.ones_like(self.pressures)
                    for i, layer in enumerate(self.pressures):
                        self.calculate_MMW()
                        MMW_profile[i] = self.MMW
                        self.mf_profile['MMW'] = MMW_profile
                    
                    # convert mixing ratios to mass fraction profiles
                    mf = vertical_mixing_ratios[species] * (mw.mass[species] / self.MMW)
                    self.mf_profile[species] = mf
                    
                 
    # method to relabel all references to line species appropriately for rebinned spectral resolutions i.e. adding _R_x suffix
    # only relevant if using spectral resolutions other than R=1000
    # needs to be done after opacity species and mass fraction profiles are defined, but before Radtrans object is created
    # abundances modified with modify_abundance should be automatically adapted, but otherwise they will need to properly relabelled
    
    def set_spectral_resolution(self, spectral_resolution):
        
        # initialise empty list to store unrelabelled species
        self.original_line_species = []
        
        # create new line species list with proper suffixes
        new_line_species = [s + '_R_' + str(spectral_resolution) for s in self.line_species]
        
        # relabel species within keys of mass fraction profiles dictionary
        new_res_abundances = {k+'_R_'+str(spectral_resolution): v for k, v in self.mf_profile.items() if k in self.line_species}
        
        """ # loop to eliminate unused dictionary elements, not recommended
        for m in self.line_species:
            if m in self.mf_profile.keys():
                if m not in self.ray_species:
                    self.mf_profile.pop(m)
        """
                    
        # saves mass fraction profiles with relabelled keys
        self.mf_profile.update(new_res_abundances)
        
        # initialise empty list to record new species
        new_species=[]
        
        # iterate over species in species list
        for count, species in enumerate(self.species_list):
            
            # if the species is a line_species it is relabelled replaces the original
            # if it also a rayleigh species it is both relabelled and added alongside the original
            # if is just a rayleigh species it is left unchanged
            # this is to fit with pRT, where lin opacities are rebinned but rayleigh opacities are not
            if species in self.line_species:
                if species in self.ray_species:
                    new_species.append(species+'_R_'+str(spectral_resolution))
                else:
                    self.species_list[count] = species+'_R_'+str(spectral_resolution)
            else:
                self.species_list[count] = species
                
        # new species (i.e. ones contributing to both line and rayleigh opacity) are added to the species list
        self.species_list.append(new_species)
        
        # original line opacity species are saved
        self.original_line_species = self.line_species
        
        # line species opacity list is updated
        self.line_species = new_line_species
    
    # method to setup opacities structure and radiative transfer object (pRT Radtrans object)
    
    def generate_opacities(self, wlen_min, wlen_max, spectral_resolution=1000):
        self.spectral_res = spectral_resolution                     # R, spectral resolution (delta lambda/lambda)
        if spectral_resolution != 1000:                             # if not default (R=1000), need to relabel appropriately
            self.set_spectral_resolution(spectral_resolution)
        self.wlen_min = wlen_min                                    # minimum wavelengths to compute opacities for (in microns)
        self.wlen_max = wlen_max                                    # maximum wavelengths to compute opacities for (in microns)
        
        # initialise radiative transfer object, reads opacities from input_data, somewhat time-intensive (depending on wavelength range and spectral resolution)
        self.atmosphere = Radtrans(line_species=self.line_species, 
                                   rayleigh_species = self.ray_species, 
                                   continuum_opacities=self.cont_species, 
                                   wlen_bords_micron=[wlen_min,wlen_max])
        
        # setup opacity structure of atmosphere based on pressure profile
        self.atmosphere.setup_opa_structure(self.pressures)

    # ---- EDIT OBJECT FUNCTIONS ----
    
    # use the following functions to safely modify the ExoPlanetAtmosphere object, including abundances and temperature profile

    # method to modify the abundance of a particular species, like adding a trace species
    # for species to contribute to opacity it must have been listed when setting up ExoPlanetAtmosphere object
    # this can be done at any point after initial setup, and automatically adjusts for spectral resolution
    # mean molecular weights and global mass fractions are by default not recalculated
    # currently advised to reinitialise object to recalculate MMW and mass fractions correctly in non-simple cases
    # can input as abundance or mass fraction depending on argument of input_type
    # vertical profile can be 'Uniform' or user-defined array with same dimensions as pressure profile

    def modify_abundance(self, species, abundance, input_type='Abundance',vertical_mixing_profile='Uniform'):
        
        # updates abundance of species and converts to mass fraction
        if input_type == 'Abundance':
            self.abundances[species] = abundance
            self.compute_mass_fractions(single_species=species)
        
        # updates mass fraction of species
        elif input_type == 'Mass Fraction':
            self.mass_fractions[species] = abundance
        
        # adds uniform mixing profile of species
        if vertical_mixing_profile == 'Uniform':
            self.mf_profile[species] = self.mass_fractions[species] * np.ones_like(self.pressures)
            
        # adds user-defined mixing profile of species
        else:
            self.mf_profile[species] = self.vertical_mixing_profile
        
        # renames relevant keys in mass fraction profile dictionary to correct spectral resolution
        if self.spectral_res != 1000:
            profile={}
            if species in self.original_line_species:
                f = self.mf_profile[species]
                species = label_species_with_resolution(species,self.spectral_res)
                profile[species] = f
            
            self.mf_profile.update(profile)
        
    # method to recalculate the equilibrium using a new C:O ratio, metallicity, or quench pressure
    
    def recalculate_equilibrium(self, COratio=None, metallicity=None, Pquench='Unchanged'):
        
        # new C:O ratio (vertically uniform)
        if COratio is not None:
            self.COratio = COratio * np.ones_like(self.pressures)
            
        # new metallicity (vertically uniform)
        if metallicity is not None:
            self.metallicity = metallicity * np.ones_like(self.pressures)
        
        # new quench pressure
        if Pquench != 'Unchanged':
            self.Pquench = Pquench
        
        # recalculates abundances using equilibrium chemistry
        self.generate_composition(equilibrium_chemistry=True)
        
        # relables keys of dictionary with appropriate suffixes if rebinned to different spectral resolution
        if self.spectral_res != 1000:
            profile={}
            for species in self.mf_profile.keys():
                if species in self.original_line_species:
                    f = self.mf_profile[species]
                    species = label_species_with_resolution(species,self.spectral_res)
                    profile[species] = f
            
            self.mf_profile.update(profile)
    
    # recalculates the temperature profile with a new equilibrium temperature, internal temrpeature, gamma, kappa_IR,
    # or can use different temperature profile method (in combination with set_func_temp_profile() if necessary)
    # note this will change the equilibrium chemistry, so you may want to recalculate_equilibrium after changing temperature profile 
    
    def recalculate_temperatures(self, T_eq=None, temperature_profile='Isothermal',T_int=None, gamma=None, kappa_IR=None):
        
        # update equilibrium temperature
        if T_eq is not None:
            self.T_eq = T_eq
        
        # update internal temperature
        if T_int is not None:
            self.T_int = T_int
        
        # update gamma (optical/IR opacities)
        if gamma is not None:
            self.gamma = gamma
        
        # update IR opacity
        if kappa_IR is not None:
            self.kappa_IR = kappa_IR
        
        # redetermine temperature profile
        if temperature_profile == 'Isothermal':
            self.temperatures = self.T_eq * np.ones_like(self.pressures)
        elif temperature_profile == 'Guillot':
            self.temperatures = physics.guillot_global(self.pressures, self.kappa_IR, self.gamma, self.gravity, self.T_int, self.T_eq)
        elif temperature_profile == 'Functional':
            for i, pressure_layer in enumerate(self.pressures):
                self.temperatures[i] = self.T_eq * self.temperature_function(pressure_layer)
        else:
            self.temperatures = temperature_profile

    # ---- SPECTRAL FUNCTIONS ----
    
    # TRANSMISSION SPECTRUM
    # method to calculate a transmission spectrum, 
    # returning the (transit radius / stellar radius)**2 in ppm as first argument
    # and the corresponding wavelengths in microns as the second argument
    # requires temperature, mass fraction and MMW profiles to be defined, as well as planet gravity, surface radius and pressure, and stellar radius
    # spectral kwargs can be passed on to Radtrans.calc_transm, to input cloud and haze parameters
    # option to plot spectrum using plot=True

    def calculate_transmission_spectrum(self, plot=False, **kwargs):
        
        # pRT calculation of transmission spectrum
        
        self.atmosphere.calc_transm(self.temperatures,self.mf_profile,self.gravity,self.mf_profile['MMW'],R_pl=self.R_planet,P0_bar=self.P0, **kwargs)
        
        # convert to appropraite units (ppm and microns)
        
        spectrum = ((self.atmosphere.transm_rad/self.R_star)**2)*10**6
        wavelengths = nc.c/self.atmosphere.freq/1e-4
        
        # plot transmission spectrum if plot=True
        
        if plot==True:
            plt.plot(wavelengths, spectrum)
        
        # returns (transit radius/stellar radius)**2 in ppm and wavelength in microns
        return spectrum, wavelengths
    
    def calculate_emission_spectrum(self, plot=False, **kwargs):
        
        # pRT calculation of emission spectrum 
        
        self.atmosphere.calc_flux(self.temperatures, self.mf_profile, self.gravity, self.mf_profile['MMW'], R_pl=self.R_planet, **kwargs)
        
        spectrum = (self.atmosphere.flux/1e-6)
        wavelengths = nc.c/self.atmosphere.freq/1e-4
        
        # plot transmission spectrum if plot=True
        
        if plot==True:
            plt.plot(wavelengths, spectrum)
        
        # returns (transit radius/stellar radius)**2 in ppm and wavelength in microns
        return spectrum, wavelengths
        
    
    # ---- UTILITY FUNCTIONS -----
    
    # pair of methods to save and reload mass fraction profiles

    def save_atmosphere_abundances(self):
        self.saved_mf_profile = {}
        for species in self.mf_profile.keys():
            self.saved_mf_profile[species] = self.mf_profile[species]
    
    def reload_atmosphere_abundances(self):
        self.mf_profile={}
        for species in self.saved_mf_profile.keys():
            self.mf_profile[species] = self.saved_mf_profile[species]
    
    
    # method to calculate a set of test spectra, with background removed, with varying amounts of a certain species, with uniform mixing profile
    # the abundance of test_species is varied over a range in log space
    # min_abundance is minimum amount of species tested (log10)
    # max_abundance is maximum amount of species tested (log10)
    # signal_samples is number of spectra generated between min and max abundance
    
    def generate_test_spectra(self, test_species, signal_samples=111, min_abundance=-11, max_abundance=0):
        print('Creating hypothesis spectra...')
        
        self.save_atmosphere_abundances()
        
        #Choose uniform prior in log space to apply sampling over
        test_abundances = np.logspace(min_abundance,max_abundance,num=signal_samples)
        """
        if self.spectral_res != 1000:
            test_species = label_species_with_resolution(test_species, self.spectral_res)
        """
        #Calculate transmission spectrum for each test abundance sampled, remove background, and add to list.
        model_list=[]
        for abundance in test_abundances:
            
            self.modify_abundance(test_species, abundance)
            model_with_background = self.calculate_transmission_spectrum()
            
            model = model_with_background[0] - self.background_spectrum
            
            model_list.append(model)
            
        
        print('    Hypothesis spectra successfully created.')
        
        self.reload_atmosphere_abundances()
        
        return model_list, test_abundances
    
    # method to generate a model spectrum, with background removed
    
    def generate_model_spectrum(self):
        
        model_with_background = self.calculate_transmission_spectrum()
        model = model_with_background[0] - self.background_spectrum
     
        return model
    
    # method to interpolate noise, if a function of wavelength, to same sampling as wavelengths of background
    # noise is noise at wavelength given by value at corresponding position in noise_wavelengths (arrays should be same size)
    # if noise is instead specified as a float or int, it will create an appropriately sized constant array
    
    def noise_interpolation(self, noise, noise_wavelengths):
        if type(noise) == int:
            noise_profile = noise * np.ones_like(self.wavelengths)
        elif type(noise) == float:
            noise_profile = noise * np.ones_like(self.wavelengths)
        else:
            noise_profile = np.interp(self.wavelengths, noise_wavelengths, noise)
            
        #print(noise_profile)
        return noise_profile
    
    
    # ---- RETRIEVAL FUNCTIONS ----
    
    # method to perform a Bayesian retrieval, using 'model_spectrum' as evidence, and test_spectra as hypothesis spectra
    # test spectra should be a list of spectrum arrays (each with same dimension as model spectrum)
    # test_parameters is set of parameters to describe test_spectra (e.g. temperatures, abundances), and should have same dimension as the list size of test_spectra
    # noise is a constant, or can be an array as a function of noise_wavelengths
    # repeats is the number of different noise profiled generated which are averaged together
    
   
    def perform_defined_retrieval(self, model_spectrum, test_spectra, test_parameters, noise, plot=True, repeats = 90, noise_wavelengths=None):
        
        noise = self.noise_interpolation(noise, noise_wavelengths)
        
        
        self.background_spectrum, self.wavelengths = self.calculate_transmission_spectrum()
        
        pdf = np.zeros_like(len(test_parameters))
        for i in range(0,repeats-1):
            noisy_model = add_noise(model_spectrum, noise)
            
            pdf = pdf + retrieval_fitting(test_spectra, noisy_model, noise)
        
        posterior_pdf = pdf/repeats
        
        
        print('    Retrieval complete.')
        
        if plot==True:
            plt.plot(test_parameters, posterior_pdf)
            plt.xscale('log')
            plt.ylabel('Bayesian likelihood')
            plt.savefig('Retrieval.pdf')
        
        return posterior_pdf
    
    # method to perform Bayesian retrieval of the abundance of a single species ('test_species'), with a given noise profile
    # noise can be constant (float or int), or a function of wavelength (array, with x set by noise_wavelengths)
    # repeats is the number of different noise profiled generated which are averaged together
    # signal_samples is the number of test_spectra with abundances given by equidistance in log-space between min_abundances and a complete atmosphere
    
    def perform_abundance_retrieval(self, test_species, noise, plot=True, signal_samples=111, repeats = 90, noise_wavelengths=None, min_abundance=-11):
            
        noise = self.noise_interpolation(noise, noise_wavelengths)
        
        self.background_spectrum, self.wavelengths = self.calculate_transmission_spectrum()
        
        model_spectrum = self.generate_model_spectrum()
        
        test_spectra, test_abundances = self.generate_test_spectra(test_species, signal_samples=signal_samples, min_abundance=min_abundance)
        
        pdf = np.zeros_like(signal_samples)
        for i in range(0,repeats-1):
            noisy_model = add_noise(model_spectrum, noise)
            
            pdf = pdf + retrieval_fitting(test_spectra, noisy_model, noise)
        
        posterior_pdf = pdf/repeats
        
        
        print('    Retrieval complete.')
        
        if plot==True:
            plt.plot(test_abundances, posterior_pdf)
            plt.xscale('log')
            plt.xlabel('Abundance of '+test_species)
            plt.ylabel('Bayesian likelihood')
            plt.savefig('Retrieval.pdf')
            
        mean_abundance=0
        std_abundance=0
        var_abundance=0
        log_test_abundances=np.log10(test_abundances)
        
        
        #calculate the mean retrieved abundance
        for i in range(0,signal_samples-1):
            mean_abundance = mean_abundance + (posterior_pdf[i] * log_test_abundances[i])
        mean_abundance_lin = 10**mean_abundance
        
        print("mean retrieved abundance = " + str(mean_abundance_lin))
        print("mean retrieved abundance (log10) = " + str(mean_abundance))
        

        #std of retrieved abundance
        for i in range(0,signal_samples-1):
            var_abundance = var_abundance + (posterior_pdf[i] * (log_test_abundances[i] - mean_abundance)**2)
        std_abundance = np.sqrt(var_abundance)
        print("std of retrieved abundance (log10) = " + str(std_abundance))
        
        return test_abundances, posterior_pdf
    
     
    # ---- PLOTTING FUNCTIONS ----
    
        
    def plot_mixing_ratios(self, plot_species='Opacity species',label=None):
        
        if plot_species == 'Opacity species':
            plot_species = self.species_list
            
        print('Plotting mixing ratios of ' + str(plot_species))
        for species in  plot_species:
            if label is None:
                
                plt.plot(self.mf_profile[species],self.pressures,label=species)
                
            else:
                plt.plot(self.mf_profile[species],self.pressures, label=label)
                
            plt.xlim(1e-12,1)
            plt.xscale('log')
            plt.yscale('log')
            plt.ylim(1e2,1e-6)
            plt.xlabel('Abundance')
            plt.ylabel('Pressure (bar)')
            plt.legend()
            
    def plot_PT_profile(self, **kwargs):
        plt.plot(self.temperatures, self.pressures)
        plt.yscale('log')
        plt.ylim(1e2,1e-6)
        plt.xlim(200,1000)
        plt.xlabel('Temperature (K)')
        plt.ylabel('Pressure (bar)')

# FUNCTION TO LOAD IN A TEXT FILE SPECIFYING NOISE (2nd column) AS A FUNCTION OF WAVELENGTH (1st column)
# Designed to easily incorporate the text files generated by PandArc

def load_noise_profile(instrument_file):
    noise_profile = np.loadtxt(instrument_file)

    noise_profile_wlens = noise_profile[:,0]
    noise_profile_noises = noise_profile[:,1]
    
    return noise_profile_wlens, noise_profile_noises

# FUNCTION TO LOAD IN INPUTS FROM A CHEMISTRY CODE, such as ARGO or VULCAN
# reads in pressure, temperature, and mixing ratios, and creates a corresponding ExoPlanetAtmosphere object

def load_chemistry_input(file_path, species_list, eq_temp, planet_radius, star_radius, gravity,
                         input_type = 'ARGO',
                         wlen_min = 0.5, wlen_max =5,
                         rayleigh_species=['H2','He'],
                         continuum_species=['H2-H2','H2-He'],
                         spectral_resolution=100,
                         plot=True,
                         P0=1,
                         **kwargs):
    
    # read ARGO .dat file (remove intro columns first)
    data = pd.read_csv(file_path, delim_whitespace=True,header=0)
    data = data.to_dict(orient='list')
    
    # read pressure-temperature profile
    
    # ARGO uses p(bar) and T(K) as headers, and records pressure in bar
    if input_type == 'ARGO':
        pressure_profile = np.array(data['p(bar)'])
        temperature_profile = np.array(data['T(K)'])
        
    # VULCAN uses Pressure and Temp as headers, and records pressure in cgs units
    elif input_type == 'VULCAN':
        pressure_profile = np.array(data['Pressure'])*1e-6
        temperature_profile = np.array(data['Temp'])
    
    # read vertical mixing profiles and calculate mean column abundances
    mixing_ratios={}
    abundances={}
    for species in species_list:
        mixing_ratios[species] = np.flip(np.array(data[species]))
        abundances[species] = np.mean(mixing_ratios[species])
        
    species_l=species_list.copy()
        
    if input_type == 'ARGO':
        # deal with the fact ARGO uses H3N instead of NH3
        species_l=species_list.copy()
        if 'H3N' in species_list:
            species_l.remove('H3N')
            species_l.append('NH3')
            mixing_ratios['NH3'] = mixing_ratios['H3N']
            abundances['NH3'] = abundances['H3N']
    
    # remove non-absorbing species from list of species with line opacities
    line_species=species_l
    line_species.remove('H2')
    line_species.remove('He')
    
    if 'N2' in line_species:
        line_species.remove('N2')

    # print vertically averaged abundances
    print(abundances)

    # flip pressure and temperature profiles for each in petitRADTRANS
    pressure_profile = np.flip(pressure_profile)
    temperature_profile = np.flip(temperature_profile)
    
    # create ExoPlanetAtmosphere object with planetary parameters
    Planet = ExoPlanetAtmosphere(planet_radius,star_radius,gravity,eq_temp)
    
    # setups object with all inputs from
    Planet.setup_atmosphere(line_species,
                            rayleigh_species,
                            continuum_species,
                            wlen_min,wlen_max,
                            abundances=abundances,
                            temperature_profile=temperature_profile,
                            pressure_profile=pressure_profile,
                            P0=P0,
                            spectral_resolution=spectral_resolution,
                            vertical_mixing_ratios=mixing_ratios,
                            plot=plot,
                            **kwargs)
    
    return Planet



"""
# --------------------------------------
# EXAMPLE MODELLING AND RETRIEVAL SCRIPT 
#   (for purposes of MPAGS submission)
# --------------------------------------

# Calculate the equilibrium temperature for temperate Jupiter, TOI-1899b, assumming a bond albedo of 0.2
Teq = compute_equilibrium_temperature(0.1587, 0.61, 0.63, 3841, 0.2)
print(Teq)

# Initialise ExoPlanetAtmosphere object with planetary parameters for TOI-189bb
Planet = ExoPlanetAtmosphere(1.15*nc.r_jup,0.607*nc.r_sun,890,Teq)

# Assume Jovian chemical abundances in atmosphere 
abundances={}
abundances['NH3'] = 2e-5
abundances['CH4'] = 2e-4
abundances['H2O'] = 0
abundances['CO'] = 0
abundances['CO2'] = 0
abundances['H2'] = 0.9
abundances['He'] = 0.1

# Load in noise profile generated using PandArc package to simulate noise from JWST

# The first profile we load in uses 2 groups per integration (no saturation)
noise_wlens, noises = load_noise_profile('NIRSpec Prism_R_100.txt')

# The second profile we load in uses 5 groups per integration (partial saturation around 1-1.5 microns)
noise_wlens_2, noises_2 = load_noise_profile('NIRSpec Prism5groups.txt')

# Define the lower and upper wavelength bounds for 
wlen_min = 2.0
wlen_max = 4.0

# Perform setup function to load in all variables in ExoPlanetAtmosphere object
# Here we assume an analytical (Guillot) pressure-temperature profile with reasonable values assumed for opacity
# We also assume a quench pressure of 1 bar and an internal temperature of 50K (evolved gas giant)
# We also assume chemical equilibrium (overwriting the previous assumption of Jovian abundances), with Solar metallicity and C/O ratio
# We use a spectral resolution of 100, which is the same as the noise profile, but this isn't essential (as it will interpolate)

Planet.setup_atmosphere(['CH4','NH3','CO','H2O','CO2','H2S'],
                        ['H2','He'],
                        ['H2-H2','H2-He'],
                        wlen_min, wlen_max,
                        abundances=abundances,
                        chemical_equilibrium=True,
                        temperature_profile='Guillot',
                        gamma=0.4,
                        kappa_IR=0.05,
                        T_int=50,
                        COratio=0.55,
                        metallicity=0,
                        plot=False,
                        spectral_resolution=100,
                        Pquench = 1)

# We can plot the vertical mixing ratios of the major chemical species
#Planet.plot_mixing_ratios(['CH4', 'H2O', 'NH3'])

# We can also plot the pressure-temperture profile
#Planet.plot_PT_profile()

# We can recalcualte the chemical equilibrium with an enhanced C/O ratio
#Planet.recalculate_equilibrium(0.8, 0)

# And plot the mixing ratios
#Planet.plot_mixing_ratios(['CH4', 'H2O', 'NH3'])

# We can also calculate and plot the planet's transmissio spectrum
#Planet.calculate_transmission_spectrum(plot=True)

# Finally, we retrieval the abundance of H2O in our forward model for our two different 
Planet.perform_abundance_retrieval('H2O', noises, noise_wavelengths=noise_wlens)

Planet.perform_abundance_retrieval('H2O', noises_2, noise_wavelengths=noise_wlens_2)




wlen_min = 1.0
wlen_max = 2.0

Planet.setup_atmosphere(['CH4','NH3','CO','H2O','CO2','H2S'],
                        ['H2','He'],
                        ['H2-H2','H2-He'],
                        wlen_min, wlen_max,
                        abundances=abundances,
                        chemical_equilibrium=True,
                        temperature_profile='Guillot',
                        gamma=0.4,
                        kappa_IR=0.05,
                        T_int=50,
                        COratio=0.55,
                        metallicity=0,
                        plot=False,
                        spectral_resolution=100,
                        Pquench = 1)


Planet.perform_abundance_retrieval('H2O', noises, noise_wavelengths=noise_wlens)

#Planet.perform_abundance_retrieval('H2O', 150)

"""

"""
Planet.modify_abundance('H2O',1e-2)
Planet.perform_abundance_retrieval('H2O',50, plot=True)

#Planet.recalculate_equilibrium(2, 0)
#Planet.calculate_transmission_spectrum(plot=True)
#Planet.plot_mixing_ratios(['H2O'])
Planet.modify_abundance('H2O',1e-4)
Planet.perform_abundance_retrieval('H2O',50, plot=True)
"""



"""
test_int_temps = np.linspace(50,250,5)
color_list=['C0','C1','C2','C3','C4', 'C5', 'C6']

for i, temp in enumerate(test_int_temps):
    Planet.set_guillot_temp_profile(0.01, 0.4, 200)
    Planet.recalculate_temperatures(271, temperature_profile='Guillot',T_int=temp)
    Planet.recalculate_equilibrium(0.55,metallicity=0,Pquench=1)
    Planet.plot_mixing_ratios(plot_species=['CO'],label='CO ' + str(temp)+' K', color=color_list[i])
    Planet.plot_mixing_ratios(plot_species=['CH4'],label='CH4 ' + str(temp)+' K',linestyle='--', color=color_list[i])
"""

"""
test_COs = np.linspace(0.2,1.4,5)
Planet.recalculate_temperatures(350)
    
for i, COs in enumerate(test_COs):
    Planet.recalculate_equilibrium(COs, 1, Pquench=1)
    #Planet.modify_abundance('NH3',0)
    Planet.plot_mixing_ratios(plot_species=['H2O'],label='H2O ' + str(COs)+'', color=color_list[i])
    Planet.plot_mixing_ratios(plot_species=['CH4'],label='CH4 ' + str(COs)+'',linestyle='--', color=color_list[i])
    #Planet.plot_mixing_ratios(plot_species=['H2'],label='H2 ' + str(COs)+'',linestyle=':', color=color_list[i])
    #Planet.calculate_transmission_spectrum(plot=True,Pcloud=0.01)



"""
"""

Planet.set_opacity_sources(['CH4','NH3','CO','H2O','CO2'], ['H2','He'], ['H2-H2','H2-He'])
Planet.generate_atmospheric_structure()

Planet.set_equilibrium_chemistry(0.5, 0)
Planet.generate_composition(equilibrium_chemistry=True)

Planet.generate_opacities(1, 5)

Planet.calculate_transmission_spectrum(plot=True)

temperatures = Planet.temperatures
print(Planet.temperatures)


Planet.set_equilibrium_chemistry(1,0)
Planet.generate_composition(equilibrium_chemistry=True)
Planet.calculate_transmission_spectrum(plot=True)
"""

"""
print(Planet.temperatures)
Planet.recalculate_temperatures(400)
Planet.recalculate_equilibrium(0.5,metallicity=0)
Planet.calculate_transmission_spectrum(plot=True)

Planet.recalculate_equilibrium(1.0,metallicity=0)
Planet.calculate_transmission_spectrum(plot=True)

Planet.recalculate_equilibrium(1.5,metallicity=0)
Planet.calculate_transmission_spectrum(plot=True)




#TOI.perform_retrieval('CH4', 10, signal_samples=111)
"""



        
    


    

        
            

    
