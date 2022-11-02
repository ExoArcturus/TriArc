#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 13:23:53 2022

@author: arcturus
"""

# import the usual
import numpy as np

# supress warnings (random warnings about Vega pop up anyway, just ignore them...)
import warnings
warnings.filterwarnings('ignore')

# import PandExo essential functions
import pandexo.engine.justdoit as jdi
import pandexo.engine.justplotit as jpi

# Intermediate function to unpack data from input file into exo_dict, which stores all inputs of the run
# obs_target should be an input file e.g. GJ1132b.py which should be already imported in the script

""" observation target input file should contain:
    
    number_of_transits - number of transits to observe (integer)
    transit_duration_time - duration of each transit (in seconds)
    out_of_transit_time - duration of observation out of transit i.e. baseline (in seconds)
    noise_floor - in ppm, by default leave as 0
    
    star_magnitude - magnitude of star at a reference wavelength
    ref_wavelength - reference wavelength for magnitude, e.g. J=1.25, H=1.6, K=2.22
    eff_tempeature - effective temperature of star
    metallicity - metallicity of star relative to solar as log Fe/H 
    log_gravity - log of surface gravity in cgs
    
    planet_radius - in Earth radii
    star_radius - in solar radii

"""

def create_exo_dict(obs_target):
    
    # Use PandExo to create exo_dict, and start filling in required values
    exo_dict = jdi.load_exo_dict()
    
    # Assume 80% saturation limit for now.
    exo_dict['observation']['sat_level'] = 80
    exo_dict['observation']['sat_unit'] = '%'
    
    # Use input for number of transits
    exo_dict['observation']['noccultations'] = obs_target.number_of_transits
    # fixed binning
    exo_dict['observation']['R'] = None
    
    # use total observing time for baseline (out of transit)
    exo_dict['observation']['baseline_unit'] = 'total'
    exo_dict['observation']['baseline'] = obs_target.out_of_transit_time
    
    exo_dict['observation']['noise_floor'] = obs_target.noise_floor
    
    # use phoenix to simulate stellar profile (can instead use 'user' for own stellar spectrum, this functionality hasn't been added to PandArc)
    exo_dict['star']['type'] = 'phoenix'
    exo_dict['star']['mag'] = obs_target.star_magnitude
    exo_dict['star']['ref_wave'] = obs_target.ref_wavelength
    exo_dict['star']['temp'] = obs_target.eff_temperature
    exo_dict['star']['metal'] = obs_target.metallicity
    exo_dict['star']['logg'] = obs_target.log_gravity
    
    # assume a constant transit depth so noise can be used for multiple different spectra.
    exo_dict['planet']['type'] = 'constant'                  
    exo_dict['planet']['transit_duration'] = obs_target.transit_duration
    exo_dict['planet']['td_unit'] = 's'
    exo_dict['planet']['radius'] = obs_target.planet_radius
    exo_dict['planet']['r_unit'] = 'R_earth'            # in earth radii
    exo_dict['star']['radius'] = obs_target.star_radius
    exo_dict['star']['r_unit'] = 'R_sun'                # in solar radii
    exo_dict['planet']['f_unit'] = 'rp^2/r*^2'          # standard units for transit depth
    return exo_dict
    

# ----- HIGHER LEVEL FUNCTION -----
# Use this to create .txt files with noise values, to be read by TriArc
# obs_target should point to imported file with necessary inputs (see above function)
# spectral_resolution can be specified, but can introduce errors with the binning, so test spectral resolutions by plotting PandExo output first
# spectral_resolution = 100,250,1000 are all tested and work well
# you can add instruments to compute using instrument_list argument. This isn't recommended, test with PandExo first.
# for use in code, perform before TriArc functions, and run again when examining new targets/spectral resolutions (will overwrite the .txt files)

def simulate_JWST_noise(obs_target,spectral_resolution=1000,instrument_list='Default'):
    
    # Error state (shouldn't occur)
    if obs_target == None:
        print('No observation target specified. Make sure to include exo_dictionary file if using JWST noise.')
        
    # create exo_dict from input file
    exo_dict = create_exo_dict(obs_target)
    
    # determines instrument list to run (in addition to manually added ones) based on spectral resolution
    
    # for high resolution, use high-res NIRSpec instruments
    if spectral_resolution > 100:
        default_instrument_list = ['NIRSpec G140H','NIRSpec G235H','NIRSpec G395H']
        
    # otherwise run medium-spec NIRSpec instruments, and low-res NIRSpec Prism and MIRI LRS.
    else:
        default_instrument_list = ['NIRSpec Prism','NIRISS SOSS','NIRSpec G395M','MIRI LRS']
        
    # if no instruments manually added, just use default list (above)
    if instrument_list == 'Default':
        instrument_list = default_instrument_list
        
    # otherwise use manually specified list, adding in default instruments as well
    else:
        for instrument in default_instrument_list not in instrument_list:
            instrument_list.append(instrument)
            
    # run PandExo to simulate observations
    print('Simulating JWST observations to generate noise...')
    results = jdi.run_pandexo(exo_dict,instrument_list,save_file=False)
    
    # high-res mode, rebins instruments to spectral_resolution, add outputs wavelength (x) and precision at that wavelength.
    if spectral_resolution > 100:
        for i, inst in enumerate(instrument_list):
            
            # extracts values from simulated transmission spectrum
            x,y,e = jpi.jwst_1d_spec(results[i][inst],plot=False,model=False,R=spectral_resolution)
            x = x[0]
            e = e[0]*1e6 #PPM UNITS!
            
            # saves file as instrument_R_spectral_resolution.txt for use by TriArc
            
            # compile appropriate file name
            file_name = inst+'_R_'+str(spectral_resolution)+'.txt'
            
            # output wavelengths and precision as two column-file
            output = np.column_stack((x,e))
            print('Saving '+file_name)
            np.savetxt(file_name,output)
            
    # low-res mode, rebins instrument to specified resolution, except MIRI and Prism which are kept at native
    else:
        for i, inst in enumerate(instrument_list):
            if 'MIRI' not in inst:
                if 'Prism' not in inst:
                    x,y,e = jpi.jwst_1d_spec(results[i][inst],plot=False,model=False,R=spectral_resolution)
                    x = x[0]
                    e = e[0]*1e6 #PPM UNITS!
                else:
                    # extracts values directly from results, so no need for jpi rebinning
                    x = results[i][inst]['FinalSpectrum']['wave']
                    e = results[i][inst]['FinalSpectrum']['error_w_floor']*1e6
            else:
                # extracts values directly from results, so no need for jpi rebinning
                x = results[i][inst]['FinalSpectrum']['wave']
                e = results[i][inst]['FinalSpectrum']['error_w_floor']*1e6
                
            # saves file as instrument_R_spectral_resolution.txt for use by TriArc
            
            # compile appropriate file name
            file_name = inst+'_R_'+str(spectral_resolution)+'.txt'
            
            # output wavelengths and precision as two column-file
            output = np.column_stack((x,e))
            print('Saving '+file_name)
            np.savetxt(file_name,output)
    
    print('    Instrument precision profiles successfully generated.')
    
def run_exo_dict(exo_dict,spectral_resolution=1000,instrument_list='Default'):
    
    # determines instrument list to run (in addition to manually added ones) based on spectral resolution
    
    # for high resolution, use high-res NIRSpec instruments
    if spectral_resolution > 100:
        default_instrument_list = ['NIRSpec G140H','NIRSpec G235H','NIRSpec G395H']
        
    # otherwise run medium-spec NIRSpec instruments, and low-res NIRSpec Prism and MIRI LRS.
    else:
        default_instrument_list = ['NIRSpec Prism','NIRISS SOSS','NIRSpec G395M','MIRI LRS']
        
    # if no instruments manually added, just use default list (above)
    if instrument_list == 'Default':
        instrument_list = default_instrument_list
        
    # otherwise use manually specified list, adding in default instruments as well
    else:
        for instrument in default_instrument_list not in instrument_list:
            instrument_list.append(instrument)
            
    # run PandExo to simulate observations
    print('Simulating JWST observations to generate noise...')
    results = jdi.run_pandexo(exo_dict,instrument_list,save_file=False)
    
    # high-res mode, rebins instruments to spectral_resolution, add outputs wavelength (x) and precision at that wavelength.
    if spectral_resolution > 100:
        for i, inst in enumerate(instrument_list):
            
            # extracts values from simulated transmission spectrum
            x,y,e = jpi.jwst_1d_spec(results[i][inst],plot=False,model=False,R=spectral_resolution)
            x = x[0]
            e = e[0]*1e6 #PPM UNITS!
            
            # saves file as instrument_R_spectral_resolution.txt for use by TriArc
            
            # compile appropriate file name
            file_name = inst+'_R_'+str(spectral_resolution)+'.txt'
            
            # output wavelengths and precision as two column-file
            output = np.column_stack((x,e))
            print('Saving '+file_name)
            np.savetxt(file_name,output)
            
    # low-res mode, rebins instrument to specified resolution, except MIRI and Prism which are kept at native
    else:
        for i, inst in enumerate(instrument_list):
            if 'MIRI' not in inst:
                if 'Prism' not in inst:
                    x,y,e = jpi.jwst_1d_spec(results[i][inst],plot=False,model=False,R=spectral_resolution)
                    x = x[0]
                    e = e[0]*1e6 #PPM UNITS!
                else:
                    # extracts values directly from results, so no need for jpi rebinning
                    x = results[i][inst]['FinalSpectrum']['wave']
                    e = results[i][inst]['FinalSpectrum']['error_w_floor']*1e6
            else:
                # extracts values directly from results, so no need for jpi rebinning
                x = results[i][inst]['FinalSpectrum']['wave']
                e = results[i][inst]['FinalSpectrum']['error_w_floor']*1e6
                
            # saves file as instrument_R_spectral_resolution.txt for use by TriArc
            
            # compile appropriate file name
            file_name = inst+'_R_'+str(spectral_resolution)+'.txt'
            
            # output wavelengths and precision as two column-file
            output = np.column_stack((x,e))
            print('Saving '+file_name)
            np.savetxt(file_name,output)
    
    print('    Instrument precision profiles successfully generated.')



    
    
    