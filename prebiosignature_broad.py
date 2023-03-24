#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 11:23:28 2022

@author: arcturus
"""

from TriArc import Spectral_Band

band_list={}


band_list['NIRISS'] = Spectral_Band('NIRISS', 1, 2.3, 30, 'NH3', 'H2S')
band_list['G395M_1'] = Spectral_Band('G395M_1', 2.9, 4.9, 30, 'HCN', 'NH3', 'CO', 'CH4', 'C2H2')
band_list['G395M_2'] = Spectral_Band('G395M_2', 2.9, 4.9, 30, 'H2CO', 'HC3N', 'SO2', 'H2S')
band_list['MIRI'] = Spectral_Band('MIRI', 5.1, 6.9, 30, 'NO', 'NH3', 'HC3N', 'H2CO')
