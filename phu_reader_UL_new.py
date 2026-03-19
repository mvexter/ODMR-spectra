# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 12:49:27 2025
@author: exter

Program read .phu file, like time-binning file produced by TimeHarp 300

Program inspired by: https://github.com/Photon-HDF5/phconvert/blob/master/notebooks/Example%20reading%20PicoQuant%20PHU%20files.ipynb
Program uses phconvert package from https://pypi.org/project/phconvert/
All routines stripped from this package and some are modified !! 
(to correct errors related to dimension of arrays)
"""

# import phconvert as phc
import matplotlib.pyplot as plt
import numpy as np
import pqreader_UL # package contains all necessary subroutines
from scipy.ndimage import gaussian_filter1d

# From pqreader_from_github
#
# phconvert - Reference library to read and save Photon-HDF5 files
#
# Copyright (C) 2014-2015 Antonino Ingargiola <tritemio@gmail.com>
#
"""
This module contains functions to load and decode files from PicoQuant
hardware.

The main functions to decode PicoQuant files (PTU, HT3, PT3, T3R) are respectively:

- :func:`load_ptu`
- :func:`load_ht3`
- :func:`load_pt3`
- :func:`load_t3r`

These functions return the arrays timestamps (also called macro-time or timetag), 
detectors (or channel), nanotimes (also called micro-time or TCSPC time) and an
additional metadata dict.

Other lower level functions are:

- :func:`ptu_reader` to load metadata and raw t3 records from PTU files
- :func:`ht3_reader` to load metadata and raw t3 records from HT3 files
- :func:`pt3_reader` to load metadata and raw t3 records from PT3 files
- :func:`process_t3records` to decode the t3 records and return
  timestamps (after overflow correction), detectors and TCSPC nanotimes.
- :func:`process_t3records_t3rfile` to decode the t3 records for t3r files.
- :func:`process_t2records` to decode the t2 records and return
  timestamps (after overflow correction) and detectors.

The functions performing overflow/rollover correction
can take advantage of numba, if installed, to significanly speed-up
the processing.
"""

# ------------------------------
filename = "Y:/Mark Mathot/Data/20251127/Rabi-20dBm-Bperp2.phu"
filename = "Y:/Mark Mathot/Data/20251126rabi/rabi1.phu"
filename = "Y:/Mark Mathot/Data/20251126rabi/rabi2.phu"
filename = "Y:/Mark Mathot/Data/20251126rabi/rabi3.phu"
filename = "Y:/Mark Mathot/Data/20251126rabi/rabi4.phu"
filename = "Y:/Mark Mathot/Data/20251126rabi/rabi5.phu"
filename = "Y:/Mark Mathot/Data/20251126rabi/rabi6.phu"
filename = "Y:/Mark Mathot/Data/20251127/rabi-longrun.phu"
filename = "Y:/Mark Mathot/Data/20251127/rabi8.phu"
filename = "Y:/Mark Mathot/Data/20251127/Rabi-20dBm-Bperp.phu"
filename = "Y:/Mark Mathot/Data/20251127/Rabi-12dBm-Bperp-shifted.phu"
filename = "Y:/Mark Mathot/Data/20251127/Rabi-16dBm-Bperp-strongshift.phu"
filename = "Y:/Mark Mathot/Data/20251127/rabi-20dbm-b02.phu"
filename = "Y:/Mark Mathot/Data/20251127/Rabi-16dBm-LowPower-LongRun.phu"
filename = "Y:/Mark Mathot/Data/20251201/Rabi-20dBm-Vref165mV-pol-45.phu"
filename2 = "Y:/Mark Mathot/Data/20251201/Rabi-20dBm-Vref165mV-pol-45.phu"
filename = "Y:/Mark Mathot/Data/20251218/Rabi-14dBm-2.772GHz-162mV-shortrun1.phu"
filename = "Y:/Mark Mathot/Data/20260112/Rabi-25dBm-2.773GHz-220mV-12+28mus-40000sec.phu"

# Step 1: readphu file
hist, bin_size, meta = pqreader_UL.load_phu(filename)
time = np.arange(hist.shape[1]) * bin_size[0] 
data_raw = hist[0]
micros = 1e6
# Step 2: Smooth data with Gaussian filter
sigma_smooth = 3 # was 3
data = gaussian_filter1d(data_raw, sigma=sigma_smooth)
## WAS: hist, bin_size, meta = phc.pqreader.load_phu(filename)
# hist.shape  # the array containing all the histograms
# bin_size  # one bin size per histogram
## REMOVED
# assert all(bin_size == bin_size[0])  # all curves must have the same bin size
# meta.keys()
# meta['acquisition_duration']  # in seconds
# n = 600
# ns = 1e9
# micros = 1e6
# time = np.arange(hist.shape[1]) * bin_size[0]  # in s

plt.figure()
plt.plot(time * micros, data)
plt.plot(time * micros, data_raw, 'b-')
# Add second curve
# hist, bin_size, meta = pqreader_UL.load_phu(filename2)
# data_raw2 = hist[0]
# data2 = gaussian_filter1d(data_raw2, sigma=sigma_smooth)
# plt.plot(time * micros, data2,'r')
plt.xlabel('Time ($\mu$s)');
plt.legend()
#plt.xlim(48, 51) 
#plt.ylim(230000,330000)
ymax = 1.005*np.max(hist[0])
ymin = 0.93*ymax
# plt.xlim(0,12)
# plt.ylim(ymin,ymax)
# plt.ylim(1.15e7,1.165e7)
plt.title(filename)
plt.grid();
## _ptu_print_tags(meta['tags'])
pqreader_UL._ptu_print_tags(meta['tags'])
## WAS: phc.pqreader._ptu_print_tags(meta['tags'])
np.equal([1,2,3], [1,2,3])
np.alltrue
plt.show()