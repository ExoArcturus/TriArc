#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 11:42:38 2022

@author: arcturus
"""

from petitRADTRANS import nat_cst as nc

planet_radius = 1.00*nc.r_earth
star_radius = 0.21*nc.r_sun
surface_gravity = 981
surface_pressure = 45
iso_temperature = 424

line_species = ['CH4','CO','CO2','H2O','NH3','HCN','C2H2','H2S','NO','SO2']
rayleigh_species = ['H2','CO']
cont_species = ['H2-H2']
species_list = ['CH4','CO','CO2','H2O','NH3','H2','He','HCN','C2H2','H2S','NO','SO2','N2']

abundances = {}
abundances['H2'] = 0.989
abundances['He'] = 0.0
abundances['N2'] =  0.0041
abundances['CH4'] = 0.0
abundances['CO'] = 0.00575
abundances['CO2'] = 0.0 
abundances['H2O'] = 0.0000001
abundances['HCN'] = 0.0
abundances['NH3'] = 0.0 

abundances['H2S'] = 0.0 
abundances['O2'] = 0.0 
abundances['O3'] = 0.0 
abundances['PH3'] = 0.0 
abundances['C2H2'] = 0.0
abundances['SO2'] = 0.0 
abundances['NO'] =  0.0 