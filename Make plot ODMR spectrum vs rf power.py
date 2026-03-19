# -*- coding: utf-8 -*-
"""
Created on Wed Dec 10 11:39:22 2025
@author: exter

Program combines measured data (width & amplitudes of ODMR spectra) into plots
"""

import numpy as np
import matplotlib.pyplot as plt

# Create data array: 
# X = rf power (in dBm)
# A = depth (%)
# W = width (MHz))
# Data at low laser power (Vref = 48 mV), labjournal p. 30 + p. 17 + 24
X1 = (-5, 0, 0, 5, 10, 12, 14, 16, 18, 20, 20, 26)
A1 = (0.51, 1.19, 0.96, 1.12, 1.2, 1.55, 1.53, 1.61, 1.52, 1.29, 1.59, 1.78)
W1 = (0.48, 0.63, 0.60, 0.68, 0.91, 1.03, 1.18, 1.56, 2.24, 3.87, 4.0, 7.1)
# Data at high laser power (Vref = 1.32 V)
# 0 dBm data = 20251121-12-15-06
X2 = (0, 12, 14, 16, 18, 20)
A2 = (0.18, 0.65, 0.64, 0.69, 0.94, 0.97)
W2 = (0.94, 1.01, 1.01, 1.36, 1.53, 1.85)
# NEW data at low laser power, with laser focus at surface (Vref = 71 mV)!!
X3 = (-17, -14, -11, -8, -5, -2, 2, 5, 8, 11, 14, 17, 20, 23, 26)
A3 = (0.37, 0.47, 0.63, 0.85, 1.01, 1.16, 1.56, 1.56, 1.75, 1.85, 1.70, 1.59, 1.64, 1.69, 1.78)
W3 = (0.38, 0.30, 0.36, 0.34, 0.41, 0.43, 0.58, 0.73, 0.89, 1.21, 1.86, 3.26, 3.98, 5.53, 6.25)
# NEW data at high laser power (vref = 1.28 V), with laser focus at surface !!
X4 = (-6, -3, 0, 2, 2, 5, 8, 10, 12, 14, 16, 16, 18, 20, 23, 26)
A4 = (0.38, 0.48, 0.67, 0.81, 0.76, 0.96, 1.17, 1.39, 1.28, 1.46, 1.57, 1.41, 1.53, 1.56, 1.68, 1.60)
W4 = (0.54, 0.52, 0.64, 0.74, 0.75, 0.95, 1.66, 1.12, 2.44, 1.98, 2.48, 2.99, 3.80, 4.55, 5.06, 6.87)
# EXTRA data at low laser power, with laser focus at surface !!
X5 = (0, 1, 3, 6, 16, 19)
A5 = (1.31, 1.33, 1.44, 1.71, 1.56, 1.66)
W5 = (0.53, 0.55, 0.64, 0.86, 2.91, 3.45)
# EXTRA data at high laser power, with laser focus at surface !!
X6 = (8, 10, 12, 12)
A6 = (1.28, 1.24, 1.42, 1.28)
W6 = (1.18, 1.55, 1.74, 2.13)
# Rabi frequency deduced from time domain measurements (page 62, spot with 1.95 MHz Rabi at 20 dBm)
x = np.arange(-30,27,0.5)
# List of reasonable parameters (added by hand) to describe the observations
fref = 1.6  # Rabi frequency in MHz at 20 dBm pump power
# fref2 = 1.65 # Same quantity measured on different day (at different x,y?)
fhom = 0.05 # Homogeneous linewidth in MHz
finh = 0.15 # Inhomogeneous linewidth in MHz
Gammap = 0.071/(2*np.pi) # Spin pump rate (about 1 MHz/Vref[in V])
Gammap2 = 1.28/(2*np.pi) # Spin pump rate (about 1 MHz/Vref[in V])
c = 4 # c=1/X with X = Gammap/Gammac (see ODMR notes)
f2e = fhom + c*Gammap # Effective damping rate in MHz
f2e2 = fhom + c*Gammap2 # Effective damping rate in MHz
# Values for Rabi frequency at different rf power
y = fref * np.exp(np.log(10)*(x-20)/20)
# y2 = fref2 * np.exp(np.log(10)*(x-20)/20)
# Values for ODMR width at two laser power, or Gammap
width = finh + np.sqrt(f2e*f2e + f2e*y*y/Gammap) 
width2 = finh + np.sqrt(f2e2*f2e2 + f2e2*y*y/Gammap) 
# Values for ODMR amplitude
Amax = 1.75 # Maximum amplitude in percent
GammaR = y*y/Gammap
A = Amax * GammaR/(GammaR+Gammap)
GammaR2 = y*y/Gammap2
A2 = Amax * GammaR2/(GammaR2+Gammap2)
 
# Create plot
fig, ax = plt.subplots()
ax.set_title('ODMR linewidth (HWHM single N line) verus rf power')
ax.semilogy(X3,W3, '*', label='low laser power at interface')
ax.semilogy(X4,W4, 'o', label='high laser power at interface')
ax.semilogy(X5,W5, 'c*', label='low laser power at interface')
ax.semilogy(X6,W6, 'ro', label='high laser power at interface')
# ax.semilogy(X1,W1, 'x', label='low laser power (Vref = 0.048 V)')
# ax.semilogy(X2,W2, '.', label='high laser power (Vref = 1.32 V)')
ax.plot(x,y,'g-', label='Rabi frequency')
# ax.plot(x,y2,'g-', label='Rabi frequency (diff. calibration)')
ax.plot(x,width, 'b-', linestyle = 'dashed') # Theoretical prediction
ax.plot(x,width2, 'r-', linestyle = 'dashed')
ax.set_xlabel("rf power in dBm")
ax.set_xlim(-30,27)
ax.set_ylim(0.1,10)
ax.set_ylabel("Linewidth $\gamma/2\pi$ [MHz]")
ax.grid(True)
ax.legend(fontsize='small')
# ax.tight_layout()

# Create plot
fig, ax = plt.subplots()
ax.set_title('ODMR linewidth (HWHM single N line) verus rf power')
ax.plot(X3,W3, '*', label='low laser power at interface')
ax.plot(X4,W4, 'o', label='high laser power at interface')
ax.plot(X5,W5, 'c*', label='low laser power at interface')
ax.plot(X6,W6, 'ro', label='high laser power at interface')
# ax.semilogy(X1,W1, 'x', label='low laser power (Vref = 0.048 V)')
# ax.semilogy(X2,W2, '.', label='high laser power (Vref = 1.32 V)')
ax.plot(x,width, 'b-', linestyle = 'dashed')
# ax.plot(x,width2, 'r-', linestyle = 'dashed')
ax.plot(x,y,'-', label='Rabi frequency')
ax.set_xlabel("rf power in dBm")
ax.set_xlim(-30,27)
ax.set_ylim(0,7.5)
ax.set_ylabel("Linewidth $\gamma/2\pi$ [MHz]")
ax.grid(True)
ax.legend(fontsize='small')
# ax.tight_layout()

# Create plot
fig, ax = plt.subplots()
ax.set_title('ODMR amplitude (in %) verus rf power')
ax.plot(X3,A3, '*', label='low laser power at interface')
ax.plot(X4,A4, 'o', label='high laser power at interface')
ax.plot(X5,A5, 'c*', label='low laser power at interface')
ax.plot(X6,A6, 'ro', label='high laser power at interface')
ax.plot(x,A, 'b-', linestyle = 'dashed')
ax.plot(x,A2, 'r-', linestyle = 'dashed')
# ax.plot(x,width2, 'r-', linestyle = 'dashed')
# ax.plot(X1,A1, 'x', label='low laser power (Vref = 0.048 V)')
# ax.plot(X2,A2, '.', label='high laser power (Vref = 1.32 V)')
ax.set_xlabel("rf power in dBm")
ax.set_xlim(-30,27)
ax.set_ylim(0,2)
ax.set_ylabel("ODMR amplitude [in %]")
ax.grid(True)
ax.legend(fontsize='small')
