# -*- coding: utf-8 -*-
"""
Created on Wed Dec 3 2025
@author: exter

Program read .phu file from TimeHarp 300 and performs several fits
Step 1: read .phu file
Step 2: smooth data & normalize
Step 3: fit exponential increase during spin pumping
        NOTE: (xmin, xmax) = fitting range needed as input
              initial guess of fitting parameters needed
Step 4: fit damped harmonic oscillation during ON time of rf source
        NOTE: (xmin, xmax) = fitting range needed as input
              initial guess of fitting parameters needed
Step 5: plot the results & add relevant parameters to figure

Modification 13 January 2026:
    - Removed single Rabi fit without additional exponent
    - Added option for bi-exponential fit of spin pumping
"""

# import phconvert as phc
import matplotlib.pyplot as plt
import numpy as np
import pqreader_UL # package contains all necessary subroutines
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit

# Define exponential curve
def exponent(x, *params):
    # WARNING: xmin is global variable !!
    A, B, tau = params[:3]
    return B - A * np.exp(-(x-xmin)/tau)

def biexponent(x, *params):
    # WARNING: xmin is global variable !!
    A1, B, tau1, A2, tau2 = params[:5]
    return B - A1 * np.exp(-(x-xmin)/tau1) - A2 * np.exp(-(x-xmin)/tau2)

def damped_harmonic(x, * params):
    # WARNING: xleft & xright are global variables !!
    A, B, tdecay, period, phase = params[:5]
    return B + A * np.exp(-x/tdecay) * np.cos(2*np.pi*(x/period)-phase)

def damped_harmonic_plus_decay(x, * params):
    # WARNING: xleft & xright are global variables !!
    A, B, tdecay, period, phase, A1, t1 = params[:7]
    return B + A * np.exp(-x/tdecay) * np.cos(2*np.pi*(x/period)-phase) + A1 * np.exp(-x/t1)

# Step 1: Specify .phu file & read file

filename = "Y:/Mark Mathot/Data/20260113/Rabi-25dBm-2.773GHz-128mV-12+88mus-3600sec.phu"
filename = "Y:/Bas Teisman/Data/20260128 - 20260130/Rabi-2.773GHz-20dBm-148mV-12+38mus-30jan.phu"
filename = "Y:/Bas Teisman/Data/20260223/Rabi-2.6837GHz-20dBm-82mV-20+100mus.phu"

# filename = "Y:/Mark Mathot/Data/20260116/Rabi-20dBm-2.773GHz-1780mV-8+22mus-900sec.phu"

# ON and OFF time of pulses microwave
xon = 20   # Time duration ON period
xoff = 100  # Time duration OFF period
xspin = xoff - 0.4 # Time duration of window used for fitting spin pump 
# xspin = 30 # Manual overrule
xrabi = xon - 0.3  # Time duration of window used for fitting Rabi oscillation
xrabi = 8  # Manaul overrule
xdetail = xon 
bi_exponent = True # False = single exponential fit, True = bi-exponential
plot_fit = True # True = show fits of Rabi oscillations

# Smoothing of data with filterwidth sigma
sigma = 4 # standard deviation of Gauss => 15% contribution at +/- 2

# Read .phy file
hist, bin_size, meta = pqreader_UL.load_phu(filename)
time_raw = np.arange(hist.shape[1]) * bin_size[0] * 1e6 # Time in microseconds
data_raw = hist[0]
time_per_bin = meta['tags']['HistResDscr_MDescResolution']['value'] * 1e6
sync_rate = meta['tags']['HistResDscr_SyncRate']['value'] # Sync rate [s^-1]
xend = 1e6/sync_rate # measurement time per run (in microseconds)
xend = xend + 1 # Add 1 microsecond to run to show the step at the end
# Remove final part of measurement = part beyond repetition time
mask = (time_raw <= xend)
time = time_raw[mask]
data = data_raw[mask]
xstart = 0.12 # Start of rf pulse, deduced from measurements (WHY 0.12 mus?)

# Step 2a: Smooth data with Gaussian filter
data = gaussian_filter1d(data, sigma=sigma)
# Step 2b: Normalize data based on values before switch ON of rf
mask = (time <= xstart)
offset = np.average(data[mask])
data = data/offset # normalization

# Step 3: Exponential decay in final part: preparation and actual fit
## FITTING RANGE & INITIAL GUESSES
xmin = xon + 0.2 # Start of fitting range
xmax = xmin + xspin # End of fitting rage
A_guess = 0.03 # Guess of modulation depth
tau_guess = 1 # Guess of decay time
offset = 1 # Guess of offset (normalized data)
# select data for fitting
mask = (time >= xmin) & (time <= xmax)
x_data = time[mask]
y_data = data[mask]

if bi_exponent == False:  # Perform single exponential fit 
    p0 = [A_guess, offset, tau_guess]
    # Find optimum fit parameters with curve_fit()
    popt, _ = curve_fit(exponent, x_data, y_data, p0=p0, maxfev=80000)
    A, B, tau = popt[:3]
    print('Fit parameters of exponential fit')
    print('Offset B = {:.4f}'.format(B))
    print('Amplitude A = {:.4f}'.format(A))
    print('Contrast A/B = {:.4f}'.format(A/B))
    print('Decay time tau = {:.4f}'.format(tau))
    print(' ')

if bi_exponent == True:  # Perform double exponential fit 
    p0 = [A_guess, offset, tau_guess, A_guess, tau_guess]
    # Find optimum fit parameters with curve_fit()
    popt, _ = curve_fit(biexponent, x_data, y_data, p0=p0, maxfev=80000)
    A1, B, tau1, A2, tau2 = popt[:5]
    print('Fit parameters of exponential fit')
    print('Offset B = {:.4f}'.format(B))
    print('Amplitude A1 = {:.4f}'.format(A1))
    print('Decay time tau1 = {:.4f}'.format(tau1))
    print('Amplitude A2 = {:.4f}'.format(A2))
    print('Decay time tau2 = {:.4f}'.format(tau2))
    print(' ')

# Step 4: Damped oscillation in initial part: preparation and actual fit
## FITTING RANGE & INITIAL GUESSES
xminbegin = 0.15 
xmaxbegin = xminbegin + xrabi
tdecay_guess = 0.5
period_guess = 0.25
A_guess = 0.02
phase_guess = 0
# # select data for fitting
# maskbegin = (time >= xminbegin) & (time <= xmaxbegin)
# x_data = time[maskbegin]
# y_data = data[maskbegin]
# p0begin = [A_guess, 0.99*offset, tdecay_guess, period_guess, phase_guess]
# # Find optimum fit parameters with curve_fit()
# poptbegin, _ = curve_fit(damped_harmonic, x_data, y_data, p0=p0begin, maxfev=80000)
# Abegin, Bbegin, tdecay, period, phase = poptbegin[:5]
# C = 1 - Bbegin
# print('Fit parameters of damped oscillation fit')
# print('Offset B = {:.4f}'.format(Bbegin))
# print('Amplitude A = {:.4f}'.format(Abegin))
# print('Contrast C = 1-B = {:.4f}'.format(C))
# print('Decay time tau = {:.4f}'.format(tdecay))
# print('Period = {:.4f}'.format(period))
# print('')

# Step 4: Damped oscillation in initial part: preparation and actual fit
## EXTRA INITIAL GUESSES
A1_guess = 0.01
t1_guess = 1
# select data for fitting
maskbegin = (time >= xminbegin) & (time <= xmaxbegin)
x_data = time[maskbegin]
y_data = data[maskbegin]
p0begin = [A_guess, 0.99*offset, tdecay_guess, period_guess, phase_guess, A1_guess, t1_guess]
# Find optimum fit parameters with curve_fit()
poptbegin2, _ = curve_fit(damped_harmonic_plus_decay, x_data, y_data, p0=p0begin, maxfev=80000)
Abegin2, Bbegin2, tdecay2, period2, phase2, A1, t1 = poptbegin2[:7]
C = 1 - Bbegin2
print('Fit parameters of damped oscillation fit + exponential decay')
print('Offset B = {:.4f}'.format(Bbegin2))
print('Amplitude A = {:.4f}'.format(Abegin2))
print('Contrast C = 1-B = {:.4f}'.format(C))
# print('Contrast A/B = {:.4f}'.format(Abegin2/Bbegin2))
print('Decay time tau = {:.4f}'.format(tdecay2))
print('Period = {:.4f}'.format(period2))
print('Extra decay, amplitude A1 = {:.4f}'.format(A1))
print('Extra decay, time t1 = {:.4f}'.format(t1))
print('')

# Step 5A: Plot full measurement (time axis in microseconds)
plt.figure(figsize=[9,6])
plt.plot(time, data)
plt.xlim(0,xend) 
plt.ylim(0.96,1.005)

# Plot curve of exponential fit
if bi_exponent == False:
    fit_curve = exponent(time, *popt)
    plt.plot(time, fit_curve, 'r-', lw=1.5, label="Exponential fit")
    plt.text(0.8*xend, 0.98, 'Tau = {:.2f} $\mu$s'.format(tau))
    plt.axvline(xmin, color='gray', linestyle='--', lw=1.2, alpha=0.8)
    plt.axvline(xmax, color='gray', linestyle='--', lw=1.2, alpha=0.8)

if bi_exponent == True:
    fit_curve = biexponent(time, *popt)
    plt.plot(time, fit_curve, 'r-', lw=1.5, label="Exponential fit")
    plt.text(0.8*xend, 0.98, 'Tau1 = {:.2f} $\mu$s'.format(tau1))
    plt.text(0.8*xend, 0.97, 'Tau2 = {:.2f} $\mu$s'.format(tau2))
    plt.axvline(xmin, color='gray', linestyle='--', lw=1.2, alpha=0.8)
    plt.axvline(xmax, color='gray', linestyle='--', lw=1.2, alpha=0.8)

# Plot curve of damped harmonic fit
# fit_curve2 = damped_harmonic(time, *poptbegin)
# plt.plot(time, fit_curve2, 'm-', lw=1.5, label="Damped harmonic fit")
# plt.text(xminbegin, 1.003, 'period = {:.3f} $\mu$s'.format(period), color='magenta')
# plt.text(xminbegin, 1.001, 'tdecay = {:.3f} $\mu$s'.format(tdecay), color='magenta')
# plt.axvline(xminbegin, color='grey', linestyle='--', lw=1.2, alpha=0.8)
# plt.axvline(xmaxbegin, color='grey', linestyle='--', lw=1.2, alpha=0.8)

# Plot curve of damped harmonic fit + extra decay
fit_curve2 = damped_harmonic_plus_decay(time, *poptbegin2)
plt.plot(time, fit_curve2, 'k-', lw=1.5, label="Damped harmonic fit")
plt.text(xmaxbegin, 1.003, 'period = {:.3f} $\mu$s'.format(period2), color='black')
plt.text(xmaxbegin, 1.001, 'tdecay = {:.3f} $\mu$s'.format(tdecay2), color='black')
plt.text(xmaxbegin, 0.999, 'A = {:.3f} '.format(Abegin2), color='black')
C = 1 -Bbegin2
plt.text(xmaxbegin, 0.997, 'C = 1-B = {:.3f} '.format(C), color='black')
plt.text(xmaxbegin, 0.995, 'extra decay = {:.3f} $\mu$s'.format(t1), color='black')
plt.axvline(xminbegin, color='gray', linestyle='--', lw=1.2, alpha=0.8)
plt.axvline(xmaxbegin, color='gray', linestyle='--', lw=1.2, alpha=0.8)

# Finish layout of plot 
plt.xlabel('Time [$\mu$s]') 
plt.ylabel('Counts per time bin [-]')
# Add vertical lines to indicate fit regio
plt.title(filename)
plt.grid();
plt.show()

# Step 5B: Plot detail of Rabi oscillations
plt.figure(figsize=[9,6])
plt.plot(time, data)
plt.xlim(0,xdetail) 
plt.ylim(0.96,1.005)

# Plot curve of damped harmonic fit
# if plot_fit:
#     fit_curve2 = damped_harmonic(time, *poptbegin)
#     plt.plot(time, fit_curve2, 'm', lw=1.5, label="Damped harmonic fit", linestyle = 'dashed')
#     plt.text(xminbegin, 1.003, 'period = {:.3f} $\mu$s'.format(period), color='magenta')
#     plt.text(xminbegin, 1.000, 'tdecay = {:.3f} $\mu$s'.format(tdecay), color='magenta')
#     plt.axvline(xminbegin, color='grey', linestyle='--', lw=1.2, alpha=0.8)
#     plt.axvline(xmaxbegin, color='grey', linestyle='--', lw=1.2, alpha=0.8)

# Plot curve of damped harmonic fit + extra decay
if plot_fit:
    fit_curve2 = damped_harmonic_plus_decay(time, *poptbegin2)
    plt.plot(time, fit_curve2, 'k', lw=1.5, label="Damped harmonic fit", linestyle = 'dashed')
    plt.text(xmaxbegin, 1.003, 'period = {:.3f} $\mu$s'.format(period2), color='black')
    plt.text(xmaxbegin, 1.000, 'tdecay = {:.3f} $\mu$s'.format(tdecay2), color='black')
    plt.text(xmaxbegin, 0.997, 'A = {:.3f} '.format(Abegin2), color='black')
    C = 1 -Bbegin2
    plt.text(xmaxbegin, 0.994, 'C = 1-B = {:.3f} '.format(C), color='black')
    plt.text(xmaxbegin, 0.991, 'extra decay = {:.3f} $\mu$s'.format(t1), color='black')
    plt.axvline(xminbegin, color='gray', linestyle='--', lw=1.2, alpha=0.8)
    plt.axvline(xmaxbegin, color='gray', linestyle='--', lw=1.2, alpha=0.8)

# Finish layout of plot 
plt.xlabel('Time [$\mu$s]') 
plt.ylabel('Counts per time bin [-]')
# Add vertical lines to indicate fit regio
plt.title(filename)
plt.grid();
plt.show()