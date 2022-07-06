#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 16:50:44 2022

@author: arcturus
"""

# EXAMPLE INPUT FILE

from petitRADTRANS import nat_cst as nc

# in cgs units
planet_radius = 2.51*nc.r_earth
star_radius = 0.21*nc.r_sun
surface_gravity = 1356

# in bar
surface_pressure = 100

# in Kelvin (isothermal atmosphere)
iso_temperature = 300

# include all species that you may want to retrieve
line_species = ['CH4','CO','CO2','H2O','NH3','HCN','C2H2','H2S','NO','SO2']
rayleigh_species = ['H2','He']
cont_species = ['H2-H2','H2-He']
species_list = ['CH4','CO','CO2','H2O','NH3','H2','He','HCN','C2H2','H2S','NO','SO2','N2']

# abundances in mixing ratios (uniform vertical mixing profiles)
abundances = {}
abundances['H2'] = 0.90
abundances['He'] = 0.0894
abundances['N2'] =  0.0 
abundances['CH4'] = 0.0005
abundances['CO'] = 0.0
abundances['CO2'] = 0.0 
abundances['H2O'] = 0.01
abundances['HCN'] = 0.0 
abundances['NH3'] = 0.0001 

abundances['H2S'] = 0.0 
abundances['O2'] = 0.0 
abundances['O3'] = 0.0 
abundances['PH3'] = 0.0 
abundances['C2H2'] = 0.0 
abundances['SO2'] = 0.0 
abundances['NO'] =  0.0 



