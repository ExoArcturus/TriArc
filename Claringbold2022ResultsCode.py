#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 13:29:13 2022

@author: arcturus
"""

# import TriArc package for functions load_planet, perform_retrieval and assess_detectability
import TriArc

import PandArc

# import background atmosphere input files (Hycean.py and URVolcanic.py)
import Hycean
import URVolcanic
import PostImpact_100kyr
import PostImpact_10Myr

# import band list
import prebiosignatures

#import observing parameters
import GJ1132b_R250

PandArc.simulate_JWST_noise(GJ1132b_R250,spectral_resolution=100)

TriArc.run_atmosphere(Hycean,
                      prebiosignatures,
                      JWST_noise = True,
                      spectral_resolution = 100,
                      output_file = 'Results_Hycean_R250.txt')


TriArc.run_atmosphere(URVolcanic,
                      prebiosignatures,
                      JWST_noise = True,
                      spectral_resolution = 100,
                      output_file = 'Results_URVolcanic.txt')

TriArc.run_atmosphere(PostImpact_100kyr,
                      prebiosignatures,
                      JWST_noise = True,
                      spectral_resolution = 100,
                      output_file = 'Results_PI100kyr.txt',
                      min_abundance = -7)

TriArc.run_atmosphere(PostImpact_10Myr,
                      prebiosignatures,
                      JWST_noise = True,
                      spectral_resolution = 100,
                      output_file = 'Results_PI10Myr.txt',
                      min_abundance = -7)
