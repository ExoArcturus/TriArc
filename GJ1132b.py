#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 15:01:12 2022

@author: arcturus
"""

number_of_transits = 5                  
transit_duration = 0.4*60*60            # in seconds
out_of_transit_time = 2*60*60           # in seconds
noise_floor = 0                         # in ppm, leave at zero

star_magnitude = 9.425                  # at reference wavelength
ref_wavelength = 1.25                   # reference wavelength in microns, J = 1.25, H = 1.6, K = 2.22
eff_temperature = 3270                  # effective temperature of star
metallicity = -0.12                     # metallicity relative to solar as log Fe/H
log_gravity = 4.881                     # log of surface gravity in cgs

planet_radius = 1.2                     # in Earth radii
star_radius = 0.21                      # in solar radii
