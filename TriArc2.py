#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 22:18:34 2022

@author: alastair
"""

import numpy as np
import matplotlib.pyplot as plt


import petitRADTRANS as pRT



import petitRADTRANS.physics as physics
import petitRADTRANS.poor_mans_nonequ_chem as chemistry

from petitRADTRANS import Radtrans
from petitRADTRANS import nat_cst as nc

import molecular_weights as mw


def compute_equilibrium_temperature(period, stellar_radius, stellar_mass, stellar_temperature, albedo):
    
    orbital_distance = (((period/365)**2) / stellar_mass)**(1/3)
    
    
    equilibrium_temperature = stellar_temperature * ( (1-albedo) * (stellar_radius / (2 * orbital_distance * 215))**2 )**(1/4)  
    
    return equilibrium_temperature


class ExoPlanet:
    def __init__(self, R_planet, R_star, gravity, T_eq):
        self.R_planet = R_planet
        self.R_star = R_star
        self.gravity = gravity
        self.T_eq = T_eq
        self.abundances = {}
        self.mass_fractions = {}
    
    def set_opacity_sources(self, line_species, ray_species, cont_species, other_species=None):
        self.line_species =  line_species
        self.ray_species = ray_species
        self.cont_species = cont_species
        species_list = line_species + ray_species #+ other_species
        if other_species is not None:
            species_list += other_species
        self.species_list = [*set(species_list)]
        print(self.species_list)
        
    def set_abundances(self, abundances):
        for species in self.species_list:
            self.abundances[species] = abundances[species]
        self.calculate_MMW()
        self.compute_mass_fractions()
        
    def calculate_MMW(self):
        self.MMW = 0
        for species in self.species_list:
            if '_R_' not in species:
                self.MMW += mw.mass[species] * self.abundances[species]
            
    def compute_mass_fractions(self, single_species=False):
        if single_species == False:
            for species in self.species_list:
                mf = self.abundances[species] * (mw.mass[species] / self.MMW)
                self.mass_fractions[species] = mf
        else:
            mf = self.abundances[single_species] * (mw.mass[single_species] / self.MMW)
            self.mass_fractions[single_species] = mf
    
    def set_guillot_temp_profile(self, kappa_IR, gamma, T_int):
        self.kappa_IR = kappa_IR
        self.gamma = gamma
        self.T_int = T_int
    
    def set_functional_temp_profile(self, temperature_function):
        self.temperature_function = temperature_function
    
    def generate_atmospheric_structure(self, P0=10, pressure_profile='Automatic', temperature_profile='Isothermal'):
        self.P0 = P0
        if pressure_profile == 'Automatic':
            self.pressures = np.logspace(-6,P0,100)
        else:
            self.pressures = pressure_profile
        
        if temperature_profile == 'Isothermal':
            self.temperatures = self.T_eq * np.ones_like(self.pressures)
        elif temperature_profile == 'Guillot':
            self.temperatures = physics.guillot_global(self.pressures, self.kappa_IR, self.gamma, self.gravity, self.T_int, self.T_eq)
        elif temperature_profile == 'Functional':
            for i, pressure_layer in enumerate(self.pressures):
                self.temperatures[i] = self.T_eq * self.temperature_function(pressure_layer)
        else:
            self.temperatures = temperature_profile
        
    def set_equilibrium_chemistry(self, COratio, metallicity, Pquench=None):
        self.Pquench = None
        self.COratio = COratio * np.ones_like(self.pressures)
        self.metallicity = metallicity * np.ones_like(self.pressures)
        if Pquench is not None:
            self.Pquench = Pquench
        
        
    def generate_composition(self, equilibrium_chemistry=False,vertical_mixing_ratios='Uniform'):
        self.mf_profile={}
        if equilibrium_chemistry == True:
            if self.Pquench is None:
                self.mf_profile = chemistry.interpol_abundances(self.COratio, self.metallicity, self.temperatures, self.pressures)
            else:
                self.mf_profile = chemistry.interpol_abundances(self.COratio, self.metallicity, self.temperatures, self.pressures, Pquench_carbon=self.Pquench)
            
            self.MMW = np.mean(self.mf_profile['MMW'])
        
        else:
            if vertical_mixing_ratios=='Uniform':
                for species in self.species_list:
                    self.mf_profile[species] = self.mass_fractions[species] * np.ones_like(self.pressures)
                    self.mf_profile['MMW'] = self.MMW * np.ones_like(self.pressures)
            
            else:
                for species in self.species_list:
                    for layer in vertical_mixing_ratios[species]:
                        self.calculate_MMW()
                        vertical_mixing_ratios['MMW'][layer] = self.MMW
                    mf = vertical_mixing_ratios[species] * (mw.mass[species] / self.MMW)
                    self.mf_profile[species] = mf
                    
                    
    def set_spectral_resolution(self, spectral_resolution):
        self.original_line_species = []
        
        new_line_species = [s + '_R_' + str(spectral_resolution) for s in self.line_species]
        new_res_abundances = {k+'_R_'+str(spectral_resolution): v for k, v in self.mf_profile.items() if k in self.line_species}
        """
        for m in self.line_species:
            if m in self.mf_profile.keys():
                if m not in self.ray_species:
                    self.mf_profile.pop(m)
        """
                    
        self.mf_profile.update(new_res_abundances)
        
        
        new_species=[]
        for count, species in enumerate(self.species_list):
            if species in self.line_species:
                if species in self.ray_species:
                    new_species.append(species+'_R_'+str(spectral_resolution))
                else:
                    self.species_list[count] = species+'_R_'+str(spectral_resolution)
            else:
                self.species_list[count] = species
        self.species_list.append(new_species)
        self.original_line_species = self.line_species
        self.line_species = new_line_species
        
    
    def generate_opacities(self, wlen_min, wlen_max, spectral_resolution=1000):
        self.spectral_res = spectral_resolution
        if spectral_resolution != 1000:
            self.set_spectral_resolution(spectral_resolution)
        self.wlen_min = wlen_min
        self.wlen_max = wlen_max
        self.atmosphere = Radtrans(line_species=self.line_species, rayleigh_species = self.ray_species, continuum_opacities=self.cont_species, wlen_bords_micron=[wlen_min,wlen_max])
        self.atmosphere.setup_opa_structure(self.pressures)

    
    def calculate_transmission_spectrum(self, plot=False, **kwargs):
        
        self.atmosphere.calc_transm(self.temperatures,self.mf_profile,self.gravity,self.mf_profile['MMW'],R_pl=self.R_planet,P0_bar=self.P0, **kwargs)
        
        
        spectrum = ((self.atmosphere.transm_rad/self.R_star)**2)*10**6
        
        wavelengths = nc.c/self.atmosphere.freq/1e-4
        
        if plot==True:
            plt.plot(wavelengths, spectrum)
        
        
        return spectrum, wavelengths
        
    
    
    # Macro function to perform all necessary setup steps at once
    
    def setup_atmosphere(self,
                         line_species,
                         ray_species,
                         cont_species,
                         wlen_min,
                         wlen_max,
                         P0=10,
                         spectral_resolution=1000,
                         abundances=None,
                         kappa_IR=None,
                         gamma=None,
                         T_int=None,
                         temperature_function=None,
                         other_species=None,
                         COratio=0.55,
                         metallicity=0,
                         Pquench=None,
                         pressure_profile='Automatic', 
                         temperature_profile='Isothermal', 
                         chemical_equilibrium=False, 
                         vertical_mixing_ratios='Uniform',
                         plot=True,
                         **kwargs):
    
        self.set_opacity_sources(line_species,ray_species,cont_species, other_species=None)
        
        if temperature_profile=='Guillot':
            self.set_guillot_temp_profile(kappa_IR,gamma,T_int)
        elif temperature_profile=='Functional':
            self.set_functional_temp_profile(temperature_function)
            
        self.generate_atmospheric_structure(P0=P0,pressure_profile=pressure_profile,temperature_profile=temperature_profile)
        
        if chemical_equilibrium==False:
            self.set_abundances(abundances)
        else:
            self.set_equilibrium_chemistry(COratio,metallicity,Pquench=Pquench)
        
        self.generate_composition(equilibrium_chemistry=chemical_equilibrium,vertical_mixing_ratios=vertical_mixing_ratios)

        self.generate_opacities(wlen_min,wlen_max,spectral_resolution=spectral_resolution)

        self.background_spectrum, self.wavelengths = self.calculate_transmission_spectrum(plot=plot,**kwargs)
        
    def generate_test_spectra(self, test_species, signal_samples=111):
        print('Creating hypothesis spectra...')
        
        #Choose uniform prior in log space to apply sampling over
        test_abundances = np.logspace(-11,0,num=signal_samples)
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
        
        return model_list, test_abundances
    
    
    
    def generate_model_spectrum(self):
        
        model_with_background = self.calculate_transmission_spectrum()
        model = model_with_background[0] - self.background_spectrum
     
        return model

    def perform_abundance_retrieval(self, test_species, noise, plot=True, signal_samples=111, repeats = 90):
            
        model_spectrum = self.generate_model_spectrum()
        
        test_spectra, test_abundances = self.generate_test_spectra(test_species, signal_samples=signal_samples)
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
        
    
    def modify_abundance(self, species, abundance, input_type='Abundance',vertical_mixing_profile='Uniform'):
        
        
        if input_type == 'Abundance':
            self.abundances[species] = abundance
            self.compute_mass_fractions(single_species=species)
            
        elif input_type == 'Mass Fraction':
            self.mass_fractions[species] = abundance
        
        if vertical_mixing_profile == 'Uniform':
            self.mf_profile[species] = self.mass_fractions[species] * np.ones_like(self.pressures)
        else:
            self.mf_profile[species] = self.vertical_mixing_profile
            
        if self.spectral_res != 1000:
            profile={}
            if species in self.original_line_species:
                f = self.mf_profile[species]
                species = label_species_with_resolution(species,self.spectral_res)
                profile[species] = f
            
            self.mf_profile.update(profile)
    
    def recalculate_equilibrium(self, COratio, metallicity, Pquench=None):
        self.COratio = COratio * np.ones_like(self.pressures)
        self.metallicity = metallicity * np.ones_like(self.pressures)
        self.Pquench = Pquench
        self.generate_composition(equilibrium_chemistry=True)
        
        if self.spectral_res != 1000:
            profile={}
            for species in self.mf_profile.keys():
                if species in self.original_line_species:
                    f = self.mf_profile[species]
                    species = label_species_with_resolution(species,self.spectral_res)
                    profile[species] = f
            
            self.mf_profile.update(profile)
        
    def recalculate_temperatures(self, temp, temperature_profile='Isothermal',T_int=None):
        if T_int is not None:
            self.T_int = T_int
        
        if temperature_profile == 'Isothermal':
            self.temperatures = temp * np.ones_like(self.pressures)
        elif temperature_profile == 'Guillot':
            self.temperatures = physics.guillot_global(self.pressures, self.kappa_IR, self.gamma, self.gravity, self.T_int, temp)
        elif temperature_profile == 'Functional':
            for i, pressure_layer in enumerate(self.pressures):
                self.temperatures[i] = self.T_eq * self.temperature_function(pressure_layer)
        else:
            self.temperatures = temperature_profile
        
    def plot_mixing_ratios(self, plot_species='Opacity species',**kwargs):
        
        if plot_species == 'Opacity species':
            plot_species = self.species_list
            
        print('Plotting mixing ratios of ' + str(plot_species))
        for species in  plot_species:
            plt.plot(self.mf_profile[species],self.pressures,**kwargs)
            plt.xlim(1e-12,1)
            plt.xscale('log')
            plt.yscale('log')
            plt.ylim(1e2,1e-6)
            plt.legend()
            
    def plot_PT_profile(self, **kwargs):
        plt.plot(self.temperatures, self.pressures)
        plt.yscale('log')
        plt.ylim(1e2,1e-6)
        plt.xlim(200,1000)
    
    
        
    
    
    
    
def add_noise(model_spectrum, noise_sd):
    noisy_model = np.random.normal(model_spectrum, noise_sd)
    return noisy_model
        
        
#Function to estimate the goodness of fit between a signal (noisy_signal) and a test signal (signal) using a radial basis likelihood function.

def goodness_of_fit(noisy_signal, signal, sigma):
    
    # calculate value of normalised gaussian least squares at each wavelength
    gauss_least_squares = ((2 * np.pi * sigma**2)**(-1/2)) * np.exp(-(1/(2)*((noisy_signal-signal)/sigma)**2))
    # find product of normalised gaussian least squares to calculate likelihood function
    likelihood = np.prod(gauss_least_squares)
    
    return likelihood


def label_species_with_resolution(species,spectral_resolution):
    species = species + '_R_' + str(spectral_resolution)
    return species    
    

def retrieval_fitting(test_spectra, model_spectrum, noise):
    fit_mean_list=[]
    
    for test in test_spectra:
        fit_list=[]
        
        for i, wlen in enumerate(test):
            #fit = np.exp( -( (test/noise) - (model_spectrum[i]/noise) )**2)
            fit = goodness_of_fit(model_spectrum[i], test, noise)
            fit_list.append(fit)
            
        fit_mean = np.mean(fit_list)
        
        fit_mean_list.append(fit_mean)
    
    normalisation = np.sum(fit_mean_list)
    #print(fit_mean_list)
    
    posterior_pdf = fit_mean_list / normalisation
    
    return posterior_pdf
        
Teq = compute_equilibrium_temperature(30, 0.61, 0.63, 3841, 0.6)
print(Teq)
        
Planet = ExoPlanet(1.3*nc.r_jup,0.6*nc.r_sun,1000,Teq)
abundances={}
abundances['NH3'] = 2e-5
abundances['CH4'] = 2e-4
abundances['H2O'] = 0
abundances['CO'] = 0
abundances['CO2'] = 0
abundances['H2'] = 0.9
abundances['He'] = 0.1

Planet.setup_atmosphere(['CH4','NH3','CO','H2O','CO2'],
                        ['H2','He'],
                        ['H2-H2','H2-He'],
                        1, 4,
                        abundances=abundances,
                        chemical_equilibrium=True,
                        temperature_profile='Isothermal',
                        gamma=0.9,
                        kappa_IR=0.01,
                        T_int=100,
                        COratio=0.5,
                        metallicity=0,
                        plot=False,
                        spectral_resolution=100)


#Planet.plot_PT_profile()
Planet.recalculate_equilibrium(0.2, 1)

Planet.perform_abundance_retrieval('H2O',60, plot=True)

Planet.recalculate_equilibrium(2, 1)

Planet.perform_abundance_retrieval('H2O',60, plot=True)




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



        
    


    

        
            

    
