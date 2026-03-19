# Date: 19 November 2025
# creator: Martin van Exter (modified from program Mark Mathot)
# Program for quick inspection of (avearge of) ODMR spectrum/spectra
# Program adds cursor, or XY-hair, to graph for inspection (class Cursor)
#              cursor only availble in %matplotlib qt

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.backend_bases import MouseEvent


class Cursor:
    """
    A cross hair cursor.
    """
    def __init__(self, ax):
        self.ax = ax
        self.horizontal_line = ax.axhline(color='k', lw=0.8, ls='--')
        self.vertical_line = ax.axvline(color='k', lw=0.8, ls='--')
        # text location in axes coordinates
        self.text = ax.text(0.72, 0.9, '', transform=ax.transAxes)

    def set_cross_hair_visible(self, visible):
        need_redraw = self.horizontal_line.get_visible() != visible
        self.horizontal_line.set_visible(visible)
        self.vertical_line.set_visible(visible)
        self.text.set_visible(visible)
        return need_redraw

    def on_mouse_move(self, event):
        if not event.inaxes:
            need_redraw = self.set_cross_hair_visible(False)
            if need_redraw:
                self.ax.figure.canvas.draw()
        else:
            self.set_cross_hair_visible(True)
            x, y = event.xdata, event.ydata
            # update the line positions
            self.horizontal_line.set_ydata([y])
            self.vertical_line.set_xdata([x])
            self.text.set_text(f'x={x:1.2f}, y={y:1.2f}')
            self.ax.figure.canvas.draw()

# Import data from files and average
# filename = "Y:/Mark Mathot/Data/20251008/20251008-15-16-46freqscan.npy"
# filename = "Y:/Mark Mathot/Data/20251008/20251008-13-38-36freqscan.npy"
# filename = "Y:/Mark Mathot/Data/20251112/20251112-12-31-50freqscan.npy"
# filename = "Y:/Mark Mathot/Data/20251031/20251031-13-16-45freqscan.npy"
filename = "Y:/Mark Mathot/Data/20251120/20251120-16-54-05freqscan.npy"
filename = "Y:/Mark Mathot/Data/20251120/20251120-18-07-57freqscan.npy"
filename = "Y:/Mark Mathot/Data/20251121/20251121-15-32-55freqscan.npy"
filename = "Y:/Mark Mathot/Data/20251121/20251121-16-02-54freqscan.npy"
filename = "Y:/Mark Mathot/Data/20251125/20251126-08-41-05freqscan.npy"
filename = "Y:/Mark Mathot/Data/20251202/20251202-11-20-40freqscan.npy"
# filename = "Y:/Mark Mathot/Data/20251120/20251120-13-59-49freqscan.npy"
filename = "Y:/Mark Mathot/Data/20251217/20251217-11-23-01freqscan.npy"

data = np.load(filename)
print('Size of datafile = ', np.shape(data))

freqs = data[:, 0]
runs = data[:, 1:]
X = runs.shape[1] # Number of runs to average over
# print(X)
# X = 30  
avg_counts = np.mean(runs[:, :X], axis=1)
freqmin = np.min(freqs)
freqmax = np.max(freqs)
countmin = np.min(avg_counts)
countmax = np.max(avg_counts)

# Create plot
fig, ax = plt.subplots()
ax.set_title('Simple cursor')
ax.plot(freqs, avg_counts, '-')
ax.set_xlabel("Frequency [GHz]")
ax.set_xlim(freqmin,freqmax)
ax.set_ylim(countmin,countmax)
ax.set_ylabel("Counts [C/s]")
ax.set_title(f"{filename}")
# ax.grid(True)
# ax.legend()
# ax.tight_layout()

# Add cursor to plot
cursor = Cursor(ax)
fig.canvas.mpl_connect('motion_notify_event', cursor.on_mouse_move)
# Simulate a mouse move to (0.5, 0.5), needed for online docs
t = ax.transData
MouseEvent(
    "motion_notify_event", ax.figure.canvas, *t.transform((0.5, 0.5))
)._process()


