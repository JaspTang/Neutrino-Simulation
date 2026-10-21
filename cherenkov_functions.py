import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import colour
import math

#constants
C = 299792458
Q = 1.602177E-19
MAGNETIC_PERMEABILITY = 0.00000125663706127

#parameters
wavelength_init = 200
n = 1
beta_init = 0.9

#Sellmeier equation
def n(w):
    s_wavelength = w * 0.001
    term_1 = ((5.672526 * 0.1) * s_wavelength * s_wavelength) / (s_wavelength * s_wavelength - 5.085550 * 0.001)
    term_1 = round(term_1, 7)
    term_2 = ((1.736581 * 0.1) * s_wavelength * s_wavelength) / (s_wavelength * s_wavelength - 1.814939 * 0.01)
    term_2 = round(term_2, 7)
    term_3 = ((2.121532 * 0.01) * s_wavelength * s_wavelength) / (s_wavelength * s_wavelength - 2.617261 * 0.01)
    term_3 = round(term_3, 7)
    term_4 = ((1.138493 * 0.1) * s_wavelength * s_wavelength) / (s_wavelength * s_wavelength - 1.073889 * 10)
    term_4 = round(term_4, 7)
    n = math.sqrt(1 + term_1 + term_2 + term_3 + term_4)
    n = round(n, 6)
    return n

#Frank-Tamm Formula for energy emitted
def calc_energy(w, beta):
    w_meter = w * math.pow(10, -9)
    d2e = (math.pi * C * C * Q * Q * MAGNETIC_PERMEABILITY * pow(w_meter, -3)) * (1 - math.pow(beta * beta * n(w) * n(w),-1))
    d2e = d2e * math.pow(10, -9)
    return d2e

#XYZ to xy conversion
def xyz_to_xy(X, Y, Z):
    total = X + Y + Z
    if total == 0:
        return 0, 0
    x = X / total
    y = Y / total
    return x, y

#Do calculations
def coordinate_calc(beta, wavelength):

    energy_values = []
    wavelength_values = []

    for x in range(99):
        energy_values.append(calc_energy(wavelength, beta))
        wavelength_values.append(wavelength)
        wavelength = wavelength + 7

    #CIE data
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

    energy_interp = interp1d(wavelength_values, energy_values, kind='cubic', fill_value="extrapolate")
    energy_values = energy_interp(common_wavelengths)

    k = 100 / np.trapezoid(energy_values * y_values, common_wavelengths)

    X = k * np.trapezoid(energy_values * x_values, common_wavelengths)
    Y = k * np.trapezoid(energy_values * y_values, common_wavelengths)
    Z = k * np.trapezoid(energy_values * z_values, common_wavelengths)

    x, y = xyz_to_xy(X, Y, Z)
    return x,y