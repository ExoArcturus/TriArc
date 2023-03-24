#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 15:06:28 2022

@author: alastair
"""

# DetArc (DETection thresholds: Arcturus) contains the functions to obtain detection thresholds for a given atmosphere
# Uses the ExoPlanetAtmosphere class from MaTriArc engine, powered by petitRADTRANS

# import the usual
import numpy as np
import matplotlib.pyplot as plt


# import ExoPlanetAtmosphere class, quick contains all setup, modelling and retrieval methods for an ExoPlanet atmosphere
from MaTriArc import ExoPlanetAtmosphere



def assess_detectability(atmosphere, 
                         noise_profile,
                         injected_species,
                         retrieved_species, 
                         spectral_resolution=1000,
                         output_file='Results.txt',
                         min_abundance = -6,
                         noise_wavelengths=None,
                         **spectrum_kwargs):
    
    atmosphere.save_atmosphere_abundances()
    
    test_spectra, test_abundances = atmosphere.generate_test_spectra(retrieved_species)
    
    atmosphere.mf_profile[retrieved_species] = 0 * np.ones_like(atmosphere.pressures)
    
    print('Finding detection threshold of ' + retrieved_species)
    
    detection_significance_sum=0
    
    
    
    
    #iterate the Bayesian analysis, performing goodness of fit between the test spectra and a model spectrum, calculating the significance of the detection, and increasing the amount of prebiosignature until a significant detection
    
    # FIRST PASS (COURSE)
    j=0
    no_detection=False
    #iterate until the detection significance is three sigma
    while detection_significance_sum < 0.997:
        
        #define end condition if never a significant detection
        max_j = (-min_abundance) + 1
        
        #choose list in log space to test abundances of the prebiosignature in the model spectrum
        exploration_list = np.logspace(min_abundance,0,(-min_abundance) + 1)
        
        injected_abundance = exploration_list[j]
        
        j = j + 1
        
        atmosphere.modify_abundance(injected_species, injected_abundance)
        model_spectrum = atmosphere.generate_model_spectrum()
        
        #perform bayesian analysis on the model spectrum and return posterior PDF
        posterior_pdf = atmosphere.perform_defined_retrieval(model_spectrum, test_spectra, test_abundances, noise_profile, noise_wavelengths=noise_wavelengths)

        #create empty variables
        new_fits_mean = []
        new_test_abundances = []
        mean_abundance=0
        
        
        log_test_abundances=np.log10(test_abundances)
        
        
        #calculate the mean retrieved abundance
        for i in range(0,len(test_abundances)):
            mean_abundance = mean_abundance + (posterior_pdf[i] * log_test_abundances[i])
        mean_abundance_lin = 10**mean_abundance
        
        posterior_pdf_size = posterior_pdf.size
        
        #sum the posterior probability for every abundance which retrieves 1% of the input or more to obtain  "detection significance"
        
        #select only the signals which retrieve 1% of the input
        for i in range(0,len(test_abundances)):
            if test_abundances[i] > mean_abundance_lin/100:
                new_fits_mean.append(posterior_pdf[i])
                new_test_abundances.append(test_abundances[i])
        
        #sum all of the successful retrievals
        detection_significance_sum = sum(new_fits_mean)
        
        
        sum_elements=0
        
        #ends the loop if the detection significance is greater than 99.7% (3 sigma)
        for i in range(posterior_pdf_size):
            sum_elements  = sum_elements + posterior_pdf[i]
            if sum_elements > 0.997:
                break
        
        print(injected_abundance)
        
        
        #ends the loop if a pure atmosphere is tested and failed to retrieve
        if j == max_j:
            print("No detection.")
            no_detection=True
            break
    
    log_injected_abundance = np.log10(injected_abundance)
    
    # SECOND PASS (FINE)
    k=1
    no_detection=False
    
    #iterate until the detection significance is three sigma
    while detection_significance_sum < 0.997:
        
        #define end condition if never a significant detection
        max_j = 10
        
        #choose list in log space to test abundances of the prebiosignature in the model spectrum
        exploration_list = np.logspace(log_injected_abundance-1,log_injected_abundance,10)
        
        injected_abundance = exploration_list[k]
        
        k = k + 1
        
        atmosphere.modify_abundance(injected_species, injected_abundance)
        model_spectrum = atmosphere.generate_model_spectrum()
        
        #perform bayesian analysis on the model spectrum and return posterior PDF
        posterior_pdf = atmosphere.perform_defined_retrieval(model_spectrum, test_spectra, test_abundances, noise_profile, noise_wavelengths=noise_wavelengths)

        #create empty variables
        new_fits_mean = []
        new_test_abundances = []
        mean_abundance=0
        
        
        log_test_abundances=np.log10(test_abundances)
        
        
        #calculate the mean retrieved abundance
        for i in range(0,len(test_abundances)):
            mean_abundance = mean_abundance + (posterior_pdf[i] * log_test_abundances[i])
        mean_abundance_lin = 10**mean_abundance
        
        posterior_pdf_size = posterior_pdf.size
        
        #sum the posterior probability for every abundance which retrieves 1% of the input or more to obtain  "detection significance"
        
        #select only the signals which retrieve 1% of the input
        for i in range(0,len(test_abundances)):
            if test_abundances[i] > mean_abundance_lin/100:
                new_fits_mean.append(posterior_pdf[i])
                new_test_abundances.append(test_abundances[i])
        
        #sum all of the successful retrievals
        detection_significance_sum = sum(new_fits_mean)
        
        
        sum_elements=0
        
        #ends the loop if the detection significance is greater than 99.7% (3 sigma)
        for i in range(posterior_pdf_size):
            sum_elements  = sum_elements + posterior_pdf[i]
            if sum_elements > 0.997:
                break
        
        print(injected_abundance)
        
        
        #ends the loop if a pure atmosphere is tested and failed to retrieve
        if j == max_j:
            print("No detection.")
            no_detection=True
            break
        
    #prints the results of the detectability assessment.
    
    if no_detection == True:
        print(retrieved_species + "\t failed to detect")
        
    
    else:
        print(retrieved_species + "\t detection threshold found")
        output = [injected_species,retrieved_species,injected_abundance,mean_abundance_lin,atmosphere.MMW]
        with open(output_file, 'a') as f:
            f.write('\n')
            for item in output:
                f.write('%s\t' % item)
                
    
    
    return injected_abundance
    
    
    
    

