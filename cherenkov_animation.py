import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
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

x,y = coordinate_calc(0.9, wavelength_init)


fig, (ax, ax_zoom) = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'wspace':0.3})
"""'width_ratios': [3, 1]}"""

colour.plotting.plot_chromaticity_diagram_CIE1931(show=False, axes=ax)
ax.set_title("CIE 1931 Chromaticity Diagram")
ax.set_xlim(0.00, 0.3)
ax.set_ylim(0.0, 0.3)
point_main, = ax.plot([x], [y], 'ro', markersize=10, label=f"Point ({x:.2f}, {y:.2f})")
legend_main = ax.legend(loc='upper right')
colour.plotting.plot_chromaticity_diagram_CIE1931(show=False, axes=ax_zoom)
ax_zoom.set_title("Zoomed-In View")
ax_zoom.grid(False)
ax_zoom.set_aspect('equal', adjustable='box')
point_zoom, = ax_zoom.plot([x], [y], 'ro', markersize=10)
ax_zoom.set_xlim(x - 0.001, x + 0.001)
ax_zoom.set_ylim(y - 0.001, y + 0.001)


beta_values = np.concatenate([
    np.linspace(0.75, 0.99, 75),  
    np.linspace(0.99, 0.75, 75)   
])

# Function to update the animation for each frame
def update(frame):
    beta = beta_values[frame]
    
    x,y = coordinate_calc(beta, 200)

    point_main.set_data([x], [y])
    legend_main.texts[0].set_text(f"Point ({x:.4f}, {y:.4f}) \n Beta: {beta:.4f}c")

    point_zoom.set_data([x], [y])
    ax_zoom.set_xlim(x - 0.01, x + 0.01)
    ax_zoom.set_ylim(y - 0.01, y + 0.01)
    return point_main, point_zoom, ax_zoom, legend_main

# Create the animation
ani = FuncAnimation(
    fig, update, frames=len(beta_values), interval=50, blit=False
)


print("Saving animation...")
ani.save("cherenkov_radiation_water_1.gif", writer="imagemagick", fps=20, dpi=200)
print("Animation saved.")

plt.show()