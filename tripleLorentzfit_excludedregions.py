import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
from scipy.signal import peak_widths


# filename = "Y:/Mark Mathot/Data/20251031/20251031-13-16-45freqscan.npy"
# filename = "Y:/Mark Mathot/Data/20251103/20251103-15-46-44freqscan.npy"
# filename = "Y:/Mark Mathot/Data/20251103/20251103-13-35-16freqscan.npy"
# filename = "Y:/Mark Mathot/Data/20251118/20251118-14-02-58freqscan.npy"
filename = "Y:/Mark Mathot/Data/20251212/20251212-11-27-16freqscan.npy"

spacing_MHz = 2.2             #ruimte tussen lorentz
num_fine = 2                  # aantal lorentz
allow_small_shifts = False     #kleine shifts side peaks
max_delta_MHz = 0.5           #max 

#runs niet meenemen waarvan waardes raar zijn
skip_runs = [40]

#wegknippen specifieke stukken
# exclude_regions = [
#    (2.9325,2.9353),   #linker kant van midden
#    (2.9409,2.945),  #rechter kant van midden
#]
exclude_regions = [
    (2.925,2.925),   #linker kant van midden
    (2.925,2.925),  #rechter kant van midden
]

sigma_smooth = 3
prominence_factor = 4
distance_pts = 100
expected_num_dips = 1
window_MHz = 60.0             
max_amp_factor = 2.0


data = np.load(filename)
freqs = data[:, 0]              #GHz
runs = data[:, 1:]
n_runs = runs.shape[1]

#de runs skippen
skip_zero_based = sorted([r - 1 for r in skip_runs if 1 <= r <= n_runs])
use_cols = [i for i in range(n_runs) if i not in skip_zero_based]
selected_runs = runs[:, use_cols]
avg_counts = np.mean(selected_runs, axis=1)
dynamic_range = np.max(avg_counts) - np.min(avg_counts)

#smoothing en bepalen dip
smoothed = gaussian_filter1d(avg_counts, sigma=sigma_smooth)
inverted = -smoothed
prom = np.std(inverted) / prominence_factor
#force single main dip: take minimum of the averaged curve
broad_centers = np.array([freqs[np.argmin(avg_counts)]])



def lorentzian(x, A, x0, gamma):
    return -A * (gamma**2) / ((x - x0) ** 2 + gamma**2)

def triple_model(x, *params):
    A, x0, gamma, offset = params[:4]
    if allow_small_shifts:
        dL, dR = params[4], params[5]
    else:
        dL = dR = 0.0
    spacing_GHz = spacing_MHz / 1000.0
    p = (lorentzian(x, A, x0 - spacing_GHz + dL, gamma) +
         lorentzian(x, A, x0, gamma) +
         lorentzian(x, A, x0 + spacing_GHz + dR, gamma))
    return p + offset


results = []
spacing_GHz = spacing_MHz / 1000.0
window_GHz = window_MHz / 1000.0
max_delta_GHz = max_delta_MHz / 1000.0

for center in broad_centers:
    mask = (freqs >= center - window_GHz/2) & (freqs <= center + window_GHz/2)
    x_data = freqs[mask]
    y_data = avg_counts[mask]

    #de geknipte stukken eruit halen
    exclude_mask = np.ones_like(x_data, dtype=bool)
    for lo, hi in exclude_regions:
        exclude_mask &= ~((x_data >= lo) & (x_data <= hi))
    x_fit = x_data[exclude_mask]
    y_fit = y_data[exclude_mask]

    #guesses
    offset_guess = np.median(y_fit)
    depth_guess = offset_guess - np.min(y_fit)
    A_guess = max(depth_guess, 0.01 * dynamic_range)
    gamma_guess = 0.0005
    x0_guess = center

    if allow_small_shifts:
        p0 = [A_guess, x0_guess, gamma_guess, offset_guess, 0.0, 0.0]
    else:
        p0 = [A_guess, x0_guess, gamma_guess, offset_guess]

    lower = [0.0, x_fit[0], 1e-8, np.min(y_fit)]
    upper = [max_amp_factor * dynamic_range, x_fit[-1], 0.01, np.max(y_fit)]
    if allow_small_shifts:
        lower += [-max_delta_GHz, -max_delta_GHz]
        upper += [ max_delta_GHz,  max_delta_GHz]

    try:
        popt, _ = curve_fit(triple_model, x_fit, y_fit, p0=p0,
                            bounds=(lower, upper), maxfev=80000)
        results.append((center, popt))
    except Exception as e:
        print(f"lukt niet op {center:.6f} GHz: {e}")

#normalize
if results:
    main_offset = results[0][1][3]
else:
    main_offset = 1.0


plt.figure(figsize=(11,6))


plt.plot(freqs, avg_counts / main_offset, label="data (normalized)")

for center, popt in results:

    fit_curve = triple_model(freqs, *popt)
    plt.plot(freqs, fit_curve / main_offset, 'r-', lw=1.5, label="Lorentz fit")


    inverted_fit = -fit_curve
    peak_index = np.argmax(inverted_fit)


    results_half = peak_widths(inverted_fit, peaks=[peak_index], rel_height=0.5)


    fwhm = results_half[0][0] * (freqs[1] - freqs[0])
    fwhm_left_freq = freqs[int(results_half[2][0])]
    fwhm_right_freq = freqs[int(results_half[3][0])]


    min_val = fit_curve[peak_index]
    flat_avg = np.max(fit_curve)
    fwhm_height = min_val + (flat_avg - min_val) / 2


    plt.hlines(y=fwhm_height / main_offset,
               xmin=fwhm_left_freq, xmax=fwhm_right_freq,
               colors='purple', linestyles='--',
               label=f"FWHM ≈ {fwhm*1e3:.3f} MHz")


    A, x0, gamma, offset = popt[:4]
    dL, dR = (popt[4:6] if allow_small_shifts else (0.0, 0.0))
    centers_fit = np.array([x0 - spacing_GHz + dL, x0, x0 + spacing_GHz + dR])

    print(f"Dip op {center:.6f} GHz")
    print(f"  Fitted centers: {[round(x, 6) for x in centers_fit]}")
    print(f"  Amplitude (normalized): {A/main_offset:.5f}")
    print(f"  Gamma (HWHM param): {gamma*1e3:.3f} MHz")
    print(f"  FWHM (measured): {fwhm*1e3:.3f} MHz")
    print(f"  Offset: {offset:.4f}")


for lo, hi in exclude_regions:
    plt.axvline(lo, color='gray', linestyle='--', lw=1.2, alpha=0.8)
    plt.axvline(hi, color='gray', linestyle='--', lw=1.2, alpha=0.8)

plt.xlabel("Frequentie (GHz)")
plt.ylabel("Normalized counts")
plt.title(f"Triple Lorentz Fit {filename}")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
