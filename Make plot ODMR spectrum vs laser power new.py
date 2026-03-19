# -*- coding: utf-8 -*-
"""
Created on Wed Dec 10 11:39:22 2025
@author: exter

Program combines measured data (width & amplitudes of ODMR spectra) into plots
"""

import numpy as np
import matplotlib.pyplot as plt

# Create data array: 
# X = optical power(V))
# A = depth (%)
# W = width (MHz))
# Data from page 26 of labjournal Martin 
# Data at 0 dBm rf power
# Also data from 21 November 2025
X1 = (0.048,0.198,0.65,1.32)
A1 = (0.96,0.60,0.38,0.18)
W1 = (0.60, 0.62, 0.74, 0.94)
# Data at 10 dBm rf power
X2 = (0.4,1.47,3.0,4.6,8.7)
A2 = (1.48,0.66,0.63,0.36,0.15)
W2 = (0.94,1.35,1.41,1.59,1.66)

# NEW data (28 jan. 2026) taken at 8 dBm (X3) and 10 dBm (X4)
X3 = (0.142, 0.33, 0.61, 0.91, 1.29, 1.76, 2.33, 3.01, 3.06, 4.45)
A3 = (1.95, 1.69, 1.73, 1.41, 1.21, 1.04, 0.94, 0.89, 0.77, 0.63)
W3 = (0.82, 0.73, 0.82, 0.81, 0.84, 0.86, 0.92, 0.89, 0.82, 0.88)
# Data at 10 dBm rf power
A4 = (1.98, 1.90, 1.81, 1.40, 1.36, 1.26, 1.14, 0.99, 0.85, 0.71)
W4 = (0.91, 0.90, 0.91, 0.99, 1.05, 0.90, 1.04, 1.06, 0.93, 0.96)
# Data at 0 dBm rf power (30 jan. 2026)
X5 = (0.148,0.986, 0.025, 0.049, 0.082, 0.26, 0.348, 0.63, 1.85, 3.03, 5.85)
A5 = (1.27, 0.53, 1.34, 1.32, 1.34, 1.19, 1.06, 0.75, 0.37, 0.28, 0.21)
W5 = (0.447, 0.545, 0.453, 0.469, 0.474, 0.483, 0.484, 0.548, 0.552, 0.575, 0.566)

# Calculated curve based on Eq. (9) in notes
Amax = 1.8
OmegaR = 2*np.pi*1.56/3.16 # 3.16 = \sqrt(10) for 10 dB less rf power
OmegaR2 = OmegaR * OmegaR
x = np.arange(0.01,10,0.01)
Gammap = x # Spin pump rate = x (Vref of laser power)
c = 2
y = Amax/(1+c*Gammap*Gammap/OmegaR2)
Amax2 = 1.35
OmegaR3 = OmegaR2/10 # 0 dBm instead of 10 dBm
y2 = Amax2/(1+c*Gammap*Gammap/OmegaR3)

# Create plot
fig, ax = plt.subplots()
ax.set_title('ODMR linewidth (HWHM single N line) verus laser power')
ax.semilogx(X3,W3, 'o', label='rf = 8 dBm')
ax.semilogx(X3,W4, 'x', label='rf = 10 dBm')
ax.semilogx(X5,W5, 'xk', label='rf = 0 dBm')
ax.set_xlabel("Laser power in V (1 V = 0.33 mW)")
ax.set_xlim(0.01,10)
ax.set_ylim(0,1.5)
ax.set_ylabel("Linewidth $\gamma/2\pi$ [MHz]")
ax.grid(True)
ax.legend()
# ax.tight_layout()

# Create plot
fig, ax = plt.subplots()
ax.set_title('ODMR amplitude (in %) verus laser power')
ax.semilogx(X3,A3, 'o', label='rf = 8 dBm')
ax.semilogx(X3,A4, 'x', label='rf = 10 dBm')
ax.semilogx(X5,A5, 'xk', label='rf = 0 dBm')
ax.semilogx(x,y, 'r', label='some calculation', linestyle = 'dashed')
ax.semilogx(x,y2, 'k', label='some calculation', linestyle = 'dashed')
ax.set_xlabel("Laser power in V (1 V = 0.33 mW)")
ax.set_xlim(0.01,10)
ax.set_ylim(0,2.2)
ax.set_ylabel("ODMR amplitude [in %]")
ax.grid(True)
ax.legend()