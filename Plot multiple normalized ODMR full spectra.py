# Date: 16 january 2026
# creator: Martin van Exter 
# Program plots multiple normalized ODMR spectra, taken at different rf power

import numpy as np
import matplotlib.pyplot as plt


# Import data from files and average
N = 3 # number of data files
power0 = 'pol. = -45 degrees'
filename0 = "Y:/Mark Mathot/Data/20260119/20260119-09-31-28freqscan.npy" 
power1 = 'pol. = 0 degrees'
filename1 = "Y:/Mark Mathot/Data/20260119/20260119-11-26-06freqscan.npy" 
power2 = 'pol. = +45 degrees'
filename2 = "Y:/Mark Mathot/Data/20260119/20260119-13-03-28freqscan.npy" 
power = np.array([power0, power1, power2])
filename = np.array([filename0, filename1, filename2])

# Create plot (initialization)
fig, ax = plt.subplots()

# Read ODMR spectra from all files in a loop
for i in range(0,N,1):
    filenamex = filename[i]
    data = np.load(filenamex)
    freqs = data[:, 0]
    runs = data[:, 1:]
    X = runs.shape[1] # Data file contains information from all separate runs!
    avg_counts = np.mean(runs[:, :X], axis=1)
    freqmin = np.min(freqs)
    freqmax = np.max(freqs)
    freqmin_range = freqmin + 0.1 * (freqmax-freqmin) # Point at 10% of range
    mask = (freqs < freqmin_range)
    front_part = avg_counts[mask]
    offset = np.average(front_part)
    signal = avg_counts/offset
    ax.plot(freqs, signal, '-', label = power[i])

# Finish plot
ax.set_title('Multiple ODMR spectra shot at different polarizations')
ax.set_xlabel("Frequency [GHz]")
ax.set_xlim(2.70,3.01)
ax.set_ylim(0.93,1.01)
ax.set_ylabel("Normalized count rate")
# ax.set_title(f"{filename}")
ax.grid(True)
ax.legend(fontsize='small')
# ax.tight_layout()

