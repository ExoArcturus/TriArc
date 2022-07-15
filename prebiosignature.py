#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 11:23:28 2022

@author: arcturus
"""

from TriArc import Spectral_Band

band_list={}

band_list['1d'] = Spectral_Band('1d', 1.55, 1.65, 30, 'H2S','NH3')
band_list['1f'] = Spectral_Band('1f', 1.90, 2.00, 30, 'H2S')
band_list['2a'] = Spectral_Band('2a', 2.1 , 2.3,  30, 'NH3','CH4')
band_list['2c'] = Spectral_Band('2c', 2.65, 2.75, 30, 'NO','H2S')
band_list['2d'] = Spectral_Band('2d', 2.9 , 3.1,  30, 'HCN','C2H2','NH3')
band_list['3a'] = Spectral_Band('3a', 3.3 , 3.5,  30, 'CH4')
band_list['3b'] = Spectral_Band('3b', 3.5 , 3.6,  30, 'HCN','CH4')
band_list['3c'] = Spectral_Band('3c', 3.65, 3.8,  30, 'C2H2','CH4')
band_list['3d'] = Spectral_Band('3d', 3.9 , 4.0,  30, 'HCN','SO2')

band_list['4a'] = Spectral_Band('4a', 4.25, 4.45, 30, 'SO2')
band_list['4b'] = Spectral_Band('4b', 4.55, 4.65, 30, 'C2H2')
band_list['4c'] = Spectral_Band('4c', 4.7 , 4.8,  30, 'HCN')
band_list['4d'] = Spectral_Band('4d', 4.8 , 4.9,  30, 'CO')

# bands that require use of MIRI instrument, only work in low-resolution mode (spectral resolutoin <= 100)
band_list['5b'] = Spectral_Band('5b', 5.3 , 5.6,  30, 'NO')
band_list['6a'] = Spectral_Band('6a', 6.2 , 6.9,  30, 'NH3')
band_list['7a'] = Spectral_Band('7a', 7.1 , 7.6,  30,' SO2')
