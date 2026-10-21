import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import colour
import math
from cherenkov_functions import *

#import CIE wavelength data (csv file)
cie_data = pd.read_csv('CIE_values.csv')

#retrieve values from csv file
wavelengths = cie_data['wavelength'].values
x_bar = cie_data['x_bar'].values
y_bar = cie_data['y_bar'].values
z_bar = cie_data['z_bar'].values

#interpolate data to create continuous line (instead of a ton of points)
common_wavelengths = np.linspace(380, 780, 1000)  # High-resolution grid
x_interp = interp1d(wavelengths, x_bar, kind='cubic', fill_value="extrapolate")
y_interp = interp1d(wavelengths, y_bar, kind='cubic', fill_value="extrapolate")
z_interp = interp1d(wavelengths, z_bar, kind='cubic', fill_value="extrapolate")

x_values = x_interp(common_wavelengths)
y_values = y_interp(common_wavelengths)
z_values = z_interp(common_wavelengths)

"""Now generate values."""

#constants
C = 299792458
Q = 1.602177E-19
MAGNETIC_PERMEABILITY = 0.00000125663706127

#parameters
wavelength_init = 200
n = 1
beta_init = 0.9

print(coordinate_calc(0.747,200))

x,y = coordinate_calc(0.9, wavelength_init)

fig, (ax, ax_zoom) = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'wspace': 0.3})
fig.subplots_adjust(bottom=0.15, wspace=0.3)

colour.plotting.plot_chromaticity_diagram_CIE1931(show=False, axes=ax)
ax.set_title("CIE 1931 Chromaticity Diagram")
colour.plotting.plot_chromaticity_diagram_CIE1931(show=False, axes=ax_zoom)
ax_zoom.set_title("Zoomed-In View")
ax_zoom.set_xlim(x - 0.01, x + 0.01)  
ax_zoom.set_ylim(y - 0.01, y + 0.01)
ax_zoom.grid(False)
ax_zoom.set_aspect('equal', adjustable='box')
point_zoom, = ax_zoom.plot([x], [y], 'ro', markersize=10)


# Plot the initial point on the chromaticity diagram

point, = ax.plot(x, y, 'ro', markersize=10, label=f"Point ({x:.2f}, {y:.2f})")
legend = ax.legend(loc='upper right')

axbeta = plt.axes([0.10,0.0,0.8,0.05])
beta_slider = Slider(axbeta, 'Beta', 0.6, 1.0, valinit = 0.9, valstep=0.01)

def update(val):
    beta = beta_slider.val
    x,y = coordinate_calc(beta, wavelength_init)
    print(x,y)

    point.set_data([x],[y])
    legend.texts[0].set_text(f"Point ({x:.4f}, {y:.4f})")

    point_zoom.set_data([x], [y])
    zoom = 0.01
    ax_zoom.set_xlim(x - zoom, x + zoom)
    ax_zoom.set_ylim(y - zoom, y + zoom)

    fig.canvas.draw_idle()  
    

beta_slider.on_changed(update)
plt.show()


