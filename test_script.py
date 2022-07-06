#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 13:29:13 2022

@author: arcturus
"""

# import TriArc package for functions load_planet, perform_retrieval and assess_detectability
import TriArc

# import background atmosphere input files (Hycean.py and URVolcanic.py)
import Hycean
import URVolcanic

# load in input files using load_planet, with the imported input file (e.g. Hycean.py) as the argument of the load_planet function
Hycean_atm = TriArc.load_planet(Hycean)
URVolcanic_atm = TriArc.load_planet(URVolcanic)

# retrieve HCN from a Hycean background atmosphere with 0.0001 mass fraction of HCN added.
TriArc.perform_retrieval('2d',Hycean_atm,'HCN','HCN',2.9,3.1,17,0.0001,spectral_resolution=100)

# calculate the detection threshold of acetylene (C2H2) in a URVolcanic background atmosphere with a noise of 17 ppm, between 2.9 and 3.1 microns
TriArc.assess_detectability('2d',URVolcanic_atm,'C2H2','C2H2',2.9,3.1,17,spectral_resolution=1000,output_file='URVolcanic.txt')

