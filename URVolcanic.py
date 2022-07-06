#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 11:15:00 2022

@author: arcturus
"""

from petitRADTRANS import nat_cst as nc

planet_radius = 1.20*nc.r_earth
star_radius = 0.21*nc.r_sun
surface_gravity = 1089
surface_pressure = 1
iso_temperature = 480

line_species = ['CH4','CO','CO2','H2O','NH3','HCN','C2H2','H2S','NO','SO2']
rayleigh_species = ['H2','N2']
cont_species = ['H2-H2']
species_list = ['CH4','CO','CO2','H2O','NH3','H2','He','HCN','C2H2','H2S','NO','SO2','N2']

abundances = {}
abundances['H2'] = 0.90
abundances['He'] = 0.001
abundances['N2'] =  0.089
abundances['CH4'] = 0.003
abundances['CO'] = 0.003
abundances['CO2'] = 0.0 
abundances['H2O'] = 0.00002
abundances['HCN'] = 0.003
abundances['NH3'] = 0.0 

abundances['H2S'] = 0.0 
abundances['O2'] = 0.0 
abundances['O3'] = 0.0 
abundances['PH3'] = 0.0 
abundances['C2H2'] = 0.0
abundances['SO2'] = 0.0 
abundances['NO'] =  0.0 