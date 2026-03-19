# Date: 19 November 2025
# creator: Martin van Exter (modified from program Mark Mathot)
# Variation on program tripleLorentzfit: 
# Select frequency range for fitting and offset
#
# Program perform triple-Lorentz on (average of) ODMR spectrum/spectra
# Step 1: Read data files, with the option to remove bad runs
# Step 2: Smooth data with Gaussian filter
# Step 3: Select fitting range, with the option to exclude some freq. ranges
# Step 4: Perform fit, with option to limit freedom of fit paramaters
#         NOTE: Even option to work with 4 or 6 fit parameters (with L/R shifts)
# Steo 5: Plot results
# NOTE: program uses some global variables, like spacing_MHz

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
from scipy.signal import peak_widths

# Subroutines
def lorentzian(x, A, x0, gamma):
    return -A * (gamma**2) / ((x - x0) ** 2 + gamma**2)

def triple_model(x, *params):
    A, x0, gamma, offset = params[:4]
    if allow_small_shifts:
        dL, dR = params[4], params[5]
    else:
        dL = dR = 0.0
#    spacing_GHz = spacing_MHz / 1000.0
    p = (lorentzian(x, A, x0 - spacing_GHz + dL, gamma) +
         lorentzian(x, A, x0, gamma) +
         lorentzian(x, A, x0 + spacing_GHz + dR, gamma))
    return p + offset

def triple_model_fixed_offset(x, *params):
    A, x0, gamma = params[:3]
    p = (lorentzian(x, A, x0 - spacing_GHz, gamma) +
         lorentzian(x, A, x0, gamma) +
         lorentzian(x, A, x0 + spacing_GHz, gamma))
    return p + offset

# Step 1: Read data files, with the option to remoe bad runs
# filename = "Y:/Mark Mathot/Data/20251031/20251031-13-16-45freqscan.npy"
# filename = "Y:/Mark Mathot/Data/20251118/20251118-15-16-43freqscan.npy"
# filename = "Y:/Mark Mathot/Data/20251118/20251118-13-09-24freqscan.npy"
# filename = "Y:/Mark Mathot/Data/20251118/20251118-14-11-04freqscan.npy"
# filename = "Y:/Mark Mathot/Data/20251118/20251118-14-16-26freqscan.npy"
# filename = "Y:/Mark Mathot/Data/20251029/20251029-13-28-28freqscan.npy"
# filename = "Y:/Mark Mathot/Data/20251029/20251029-14-51-51freqscan.npy"
filename = "Y:/Mark Mathot/Data/20251118/20251118-14-02-58freqscan.npy"
filename = "Y:/Mark Mathot/Data/20251121/20251121-16-02-54freqscan.npy"
# filename = "Y:/Mark Mathot/Data/20251118/20251118-16-54-28freqscan.npy"
# filename = "Y:/Mark Mathot/Data/20251120/20251120-13-59-49freqscan.npy"
# filename = "Y:/Mark Mathot/Data/20251121/20251121-15-32-55freqscan.npy"
filename = "Y:/Mark Mathot/Data/20251118/20251118-13-43-12freqscan.npy"
filename = "Y:/Mark Mathot/Data/20251125/20251125-08-43-24freqscan.npy"
filename = "Y:/Mark Mathot/Data/20251125/20251125-13-21-50freqscan.npy"
filename = "Y:/Mark Mathot/Data/20251202/20251202-11-24-51freqscan.npy"
filename = "Y:/Mark Mathot/Data/20260105/20260106-08-16-52freqscan.npy"
filename = "Y:/Bas Teisman/Data/20260302/20260302-10-22-39freqscan.npy"

data = np.load(filename)
freqs = data[:, 0]              #GHz
runs = data[:, 1:]
n_runs = runs.shape[1]

# Runs niet meenemen waarvan waardes raar zijn
skip_runs = [0] # was [40] for file 20251031-13-16-45
skip_zero_based = sorted([r - 1 for r in skip_runs if 1 <= r <= n_runs])
use_cols = [i for i in range(n_runs) if i not in skip_zero_based]
selected_runs = runs[:, use_cols]
avg_counts = np.mean(selected_runs, axis=1)
dynamic_range = np.max(avg_counts) - np.min(avg_counts)

# Step 2: Smooth data with Gaussian filter
sigma_smooth = 1 # was 3
smoothed = gaussian_filter1d(avg_counts, sigma=sigma_smooth)
inverted = -smoothed
# Two lines removed, because not used
# prominence_factor = 4
# prom = np.std(inverted) / prominence_factor # Standard deviation of smoothed data
# Find frequency corresponding to the minimum of avg_counts
center = 2.6837    # Center frequency used for fitting
window_MHz = 10     # Window used for fitting  
freqlow_offset  = 2.660
freqhigh_offset = 2.670       # Freq. window used to find offset 
mask = (freqs >= freqlow_offset) & (freqs <= freqhigh_offset)
offset = np.mean(smoothed, where = mask)    
# print('Shape of broad_centers = ',broad_centers.shape)

# Step 3: Select fitting range, with the option to exclude some freq. ranges
# FOR NOW REMOVE OPTION TO EXCLUDE REGION FROM FITTING
#exclude_regions = [
#    (2.9325,2.9353),   #linker kant van midden
#    (2.9409,2.945),  #rechter kant van midden
#]

# PREPARATION FOR FITTING
spacing_MHz = 2.16 #2             # ruimte tussen lorentz
num_fine = 2                  # aantal lorentz
allow_small_shifts = False     # kleine shifts side peaks
max_delta_MHz = 0.5           # max allowed shifts of L/R peaks

distance_pts = 100
expected_num_dips = 1
max_amp_factor = 2.0

spacing_GHz = spacing_MHz / 1000.0
window_GHz = window_MHz / 1000.0
max_delta_GHz = max_delta_MHz / 1000.0

results = [] # Create array for storage of results
mask = (freqs >= center - window_GHz/2) & (freqs <= center + window_GHz/2)
x_data = freqs[mask]
y_data = avg_counts[mask]

# FOR NOW REMOVE OPTION TO EXCLUDE PART OF SPECTRUM FROM FITTING
    # de geknipte stukken eruit halen
#    exclude_mask = np.ones_like(x_data, dtype=bool)
#    for lo, hi in exclude_regions:
#        exclude_mask &= ~((x_data >= lo) & (x_data <= hi))
#    x_fit = x_data[exclude_mask]
#    y_fit = y_data[exclude_mask]
x_fit = x_data
y_fit = y_data

    # guesses
depth_guess = offset - np.min(y_fit)
A_guess = max(depth_guess, 0.01 * dynamic_range)
gamma_guess = 0.001
x0_guess = center

p0 = [A_guess, x0_guess, gamma_guess]
# Question: What is the function of the block below?
popt, _ = curve_fit(triple_model_fixed_offset, x_fit, y_fit, p0=p0, maxfev=80000)
results.append((center, popt))
main_offset = offset

# Steo 5: Plot results
plt.figure(figsize=(11,6))
# Plot averaged data
plt.plot(freqs, avg_counts / main_offset, label="data (normalized)")
# Plot fit curves
for center, popt in results:
    fit_curve = triple_model_fixed_offset(freqs, *popt)
    plt.plot(freqs, fit_curve / main_offset, 'r-', lw=1.5, label="Lorentz fit")
    
    # Code to find FWHM, based on maxium and fit_curve: too complicated? 
    # Invert data for easy search of peak with scipy program peak_widths
    inverted_fit = -fit_curve
    peak_index = np.argmax(inverted_fit)
    # Scipy program peak_width returns: width, max, left_value, right_value
    results_half = peak_widths(inverted_fit, peaks=[peak_index], rel_height=0.5)
    # Question: Why three FWHM?
    fwhm = results_half[0][0] * (freqs[1] - freqs[0])
    fwhm_left_freq = freqs[int(results_half[2][0])]
    fwhm_right_freq = freqs[int(results_half[3][0])]
    min_val = fit_curve[peak_index]
    flat_avg = np.max(fit_curve)
    fwhm_height = min_val + (flat_avg - min_val) / 2
    
    contrast = 1 - min_val/flat_avg
    
    # Add horizontal lines to plot
    plt.hlines(y=fwhm_height / main_offset,
               xmin=fwhm_left_freq, xmax=fwhm_right_freq,
               colors='purple', linestyles='--',
               label=f"FWHM ≈ {fwhm*1e3:.3f} MHz")
    
    # Add fit parameters to plot
    A, x0, gamma = popt[:3]
    centers_fit = np.array([x0 - spacing_GHz, x0, x0 + spacing_GHz])
    print(f"Dip op {center:.6f} GHz")
    print(f"  Fitted centers: {[round(x, 5) for x in centers_fit]}")
    print(f"  Amplitude (normalized): {A/main_offset:.5f}")
    print(f"  Contrast (normalized): {contrast:.5f}")
    print(f"  Gamma (HWHM param): {gamma*1e3:.3f} MHz")
    print(f"  FWHM (measured): {fwhm*1e3:.3f} MHz")
    print(f"  Offset: {offset:.0f}")
# Add vertical lines to indicated excluded regions to plot
# for lo, hi in exclude_regions:
#    plt.axvline(lo, color='gray', linestyle='--', lw=1.2, alpha=0.8)
#    plt.axvline(hi, color='gray', linestyle='--', lw=1.2, alpha=0.8)
# Finish layout of plot
# Show frequency window used for fitting
plt.axvline(center - window_GHz/2, color='gray', linestyle='--', lw=1.2, alpha=0.8)
plt.axvline(center + window_GHz/2, color='gray', linestyle='--', lw=1.2, alpha=0.8)

plt.xlabel("Frequentie (GHz)")
plt.ylabel("Normalized counts")
plt.title(f"Triple Lorentz Fit {filename}")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
