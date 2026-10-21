# Cherenkov Radiation in Water: Spectrum and Visible Color

An exploration of how the speed of a charged particle changes the spectrum and visible chromaticity of Cherenkov radiation in water. The project combines a wavelength-dependent refractive-index model, the Frank–Tamm radiation formula, and CIE 1931 color matching to place the emitted light on a chromaticity diagram.

## What is being simulated?

The central question is: **what visible color does the Cherenkov emission spectrum predict, and how does that color change with particle speed?**

A charged particle moving through water can travel faster than light's **phase velocity in the water**, while still moving slower than the vacuum speed of light. Above this threshold it emits Cherenkov radiation. The calculation here predicts the spectrum at emission and converts its visible portion into a standardized description of human color perception.

This connects to neutrino detection because neutrinos can interact in water and produce charged particles, such as electrons or muons. Those charged particles emit the light observed by water Cherenkov detectors. See the [Super-Kamiokande detector explanation](https://www-sk.icrr.u-tokyo.ac.jp/en/sk/about/detector/).

The code models the **optical emission and color calculation** in that chain. It does not generate neutrino interactions, simulate neutrino oscillations, propagate particle tracks, or reconstruct detector events. Particle speed is supplied directly through the dimensionless parameter $\beta=v/c$.

### Outputs

- An interactive CIE 1931 diagram with a speed slider and a second, zoomed view of the predicted chromaticity.
- An animated sweep from $\beta=0.75$ to $0.99$ and back.
- A notebook exploring color-matching functions, spreadsheet-based emission spectra, and XYZ/xy conversion.
- Saved diagrams and a GIF from the original exploration.

![Saved comparison of particle speeds on the CIE 1931 diagram](CIE%201931%20Plotted%20Speeds%2C%20Water.png)

The animation is available as [cherenkov_radiation_water_1.gif](cherenkov_radiation_water_1.gif), approximately 22.6 MB. It shows changing chromaticity coordinates, rather than a spatial animation of a particle or its radiation cone.

## Physics and equations

### 1. Particle speed and the Cherenkov threshold

The speed input is

$$
\beta=\frac{v}{c}, \qquad c=299\,792\,458\ \mathrm{m\,s^{-1}}.
$$

At vacuum wavelength $\lambda$, light's phase velocity in the medium is $c/n(\lambda)$. Cherenkov emission requires

$$
\beta n(\lambda)>1, \qquad \beta_{\mathrm{threshold}}(\lambda)=\frac{1}{n(\lambda)}.
$$

For an illustrative visible refractive index of about 1.33, the threshold is roughly $0.75c$. The precise threshold varies with wavelength because water is dispersive.

The associated emission angle satisfies

$$
\cos\theta_C(\lambda)=\frac{1}{\beta n(\lambda)}.
$$

This angle explains the radiation cone and detector rings, but the current code does not compute angles or render ring geometry. See the [Particle Data Group review of particle passage through matter, section on Cherenkov radiation](https://pdg.lbl.gov/2022/reviews/rpp2022-rev-passage-particles-matter.pdf).

### 2. Water dispersion: the Sellmeier equation

The function `n(w)` in `cherenkov_functions.py` computes the refractive index from a four-term Sellmeier model:

$$
n^2(\lambda)=1+\sum_{i=1}^{4}\frac{B_i\lambda_{\mu m}^{2}}{\lambda_{\mu m}^{2}-C_i},
\qquad \lambda_{\mu m}=10^{-3}\lambda_{nm}.
$$

The coefficients below are transcribed from the implementation. The $B_i$ are dimensionless and the $C_i$ have units of $\mu\mathrm{m}^2$.

| Term | $B_i$ | $C_i$ |
| --- | --- | --- |
| 1 | 0.5672526 | 0.005085550 |
| 2 | 0.1736581 | 0.01814939 |
| 3 | 0.02121532 | 0.02617261 |
| 4 | 0.1138493 | 10.73889 |

Each rational term is rounded to seven decimal places before summation; the final square root is rounded to six decimal places. Wavelength is passed to the function in nanometers.

These coefficients are treated as the project's existing water model. Their original source, temperature, and fitted validity interval are not recorded in the source files. Published water dispersion measurements, such as [Daimon and Masumura (2007)](https://doi.org/10.1364/AO.46.003811), provide relevant background, but this README does not assert that the coefficients above match a particular published fit.

### 3. Emission spectrum: the Frank–Tamm formula

For a singly charged particle, `calc_energy(w, beta)` evaluates the wavelength form of the Frank–Tamm energy spectrum:

$$
\frac{d^2E}{d\ell\,d\lambda_m}
=\frac{\pi\mu_0 q^2c^2}{\lambda_m^3}
\left(1-\frac{1}{\beta^2n^2(\lambda)}\right),
\qquad \beta n(\lambda)>1.
$$

Here $\ell$ is particle path length and $\lambda_m$ is the vacuum wavelength in meters. The constants in the implementation are

| Symbol | Code name | Value |
| --- | --- | --- |
| $c$ | `C` | $299792458\ \mathrm{m\,s^{-1}}$ |
| $q$ | `Q` | $1.602177\times10^{-19}\ \mathrm{C}$ |
| $\mu_0$ | `MAGNETIC_PERMEABILITY` | $1.25663706127\times10^{-6}\ \mathrm{H\,m^{-1}}$ |

The code first converts the input wavelength using $\lambda_m=10^{-9}\lambda_{nm}$. It then multiplies the result by $10^{-9}$ to express the spectral density **per nanometer** rather than per meter of wavelength. Its returned quantity therefore has units of $\mathrm{J\,m^{-1}\,nm^{-1}}$: energy per unit particle path length per unit wavelength, rather than total energy or power.

The $\lambda^{-3}$ factor favors short wavelengths in the energy spectrum. The commonly quoted photon-number spectrum has a $\lambda^{-2}$ dependence; multiplying each photon by its energy $hc/\lambda$ gives the extra inverse-wavelength factor. Dispersion and the threshold term modify this basic spectral shape.

**Current threshold limitation:** `calc_energy` evaluates the expression at every speed without setting it to zero below threshold. It can therefore return negative values when $\beta n(\lambda)<1$. Physical emission should instead be zero there:

$$
S_{\mathrm{physical}}(\lambda,\beta)=
\begin{cases}
\dfrac{\pi\mu_0q^2c^2}{\lambda_m^3}\left(1-\dfrac{1}{\beta^2n^2(\lambda)}\right)10^{-9}, & \beta n(\lambda)>1,\\
0, & \beta n(\lambda)\leq1.
\end{cases}
$$

That piecewise expression describes the physical correction; the original source has been preserved. Slider positions below threshold, and spectra that straddle the threshold, should not be interpreted as physical color predictions. Normalizing a negative spectrum can still produce plausible-looking coordinates.

### 4. Sampling and interpolation

`coordinate_calc(beta, wavelength)` builds a spectrum using 99 samples spaced 7 nm apart. With the default starting wavelength of 200 nm, these samples span **200–886 nm**:

$$
\lambda_j=200+7j\ \mathrm{nm},\qquad j=0,\ldots,98.
$$

The function then:

1. Loads the color-matching table from `CIE_values.csv`.
2. Builds cubic interpolators for the spectrum and each matching function using SciPy's `interp1d`.
3. Evaluates them on a common grid of **1,000 points from 380 to 780 nm**.
4. Uses NumPy's trapezoidal integration to compute XYZ values.

Only that visible interval contributes to the color calculation. The CSV itself spans 360–830 nm at 1 nm spacing. All interpolators enable extrapolation, although the default visible grid lies inside both input ranges. Changing the starting wavelength can cause spectral extrapolation. Cubic interpolation can also overshoot near a sharp threshold; a physically corrected calculation would need to account for this.

### 5. From spectrum to CIE XYZ

Let $S(\lambda)$ be the interpolated spectral energy density and $\bar{x}(\lambda)$, $\bar{y}(\lambda)$, and $\bar{z}(\lambda)$ be the CIE 1931 color-matching functions. The code computes

$$
X=k\int_{380}^{780}S(\lambda)\bar{x}(\lambda)\,d\lambda,
$$

$$
Y=k\int_{380}^{780}S(\lambda)\bar{y}(\lambda)\,d\lambda,
\qquad
Z=k\int_{380}^{780}S(\lambda)\bar{z}(\lambda)\,d\lambda,
$$

with normalization

$$
k=\frac{100}{\int_{380}^{780}S(\lambda)\bar{y}(\lambda)\,d\lambda}.
$$

This fixes $Y=100$ whenever the denominator is nonzero. It allows comparison of spectral color but removes absolute intensity: two points on the diagram do not tell you which particle produces more light. A zero spectrum has no defined chromaticity and requires special handling; the current code does not guard the normalization denominator.

The XYZ matching functions describe a standard observer's color matches. They are not direct measurements of the sensitivities of the three cone-cell types. The [official CIE 1931 2° observer dataset](https://cie.co.at/datatable/cie-1931-colour-matching-functions-2-degree-observer) is a reference for these functions; the local CSV's original download provenance is not recorded.

### 6. From XYZ to chromaticity coordinates

`xyz_to_xy(X, Y, Z)` removes the overall scale:

$$
x=\frac{X}{X+Y+Z},\qquad y=\frac{Y}{X+Y+Z}.
$$

These $(x,y)$ values locate the spectrum on the CIE 1931 chromaticity diagram. The `colour-science` library draws the diagram, and Matplotlib overlays the calculated point. The red marker identifies a coordinate; its red fill is not a rendering of the emitted light's predicted color. The function returns `(0, 0)` when the XYZ sum is zero, although physically the chromaticity of zero light is undefined.

## Project files

| File | Purpose |
| --- | --- |
| `cherenkov_functions.py` | Dispersion, spectral emission, XYZ integration, and xy conversion. |
| `cherenkov_simulation.py` | Interactive diagram; slider spans $\beta=0.60$–$1.00$ in steps of 0.01, initially 0.90. |
| `cherenkov_animation.py` | Generates 150 frames sweeping $\beta=0.75$–$0.99$ and back; exports a GIF at 20 fps and 200 dpi. |
| `cherenkov.ipynb` | Original notebook using locally stored CIE data and externally hosted spreadsheet spectra. |
| `CIE_values.csv` | Wavelength and `x_bar`, `y_bar`, `z_bar` columns used by the color calculation. |
| `CIE 1931*.png` | Saved plots from the exploration. |
| `cherenkov_radiation_water_1.gif` | Saved animated chromaticity sweep. |
| `requirements.txt` | Python dependencies for the scripts and notebook. |

## Running the project

Use Python 3.10 or newer with a fresh virtual environment. NumPy 2.0 or newer is required because the code calls [`numpy.trapezoid`](https://numpy.org/doc/stable/reference/generated/numpy.trapezoid.html). Dependencies are specified as minimum versions, rather than a fully locked environment.

From PowerShell in this project directory:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe cherenkov_simulation.py
```

Run commands from the repository root: the scripts read `CIE_values.csv` using a relative path. The interactive window requires a graphical Matplotlib backend. The package providing `import colour` is named **`colour-science`**.

To calculate a coordinate without opening the interactive application:

```powershell
.\.venv\Scripts\python.exe -c "from cherenkov_functions import coordinate_calc; print(coordinate_calc(0.9, 200))"
```

### Export the animation

The existing exporter requests the external **ImageMagick** writer. Install ImageMagick and ensure Matplotlib can locate it before running:

```powershell
.\.venv\Scripts\python.exe cherenkov_animation.py
```

This writes `cherenkov_radiation_water_1.gif`, replacing an existing file with that name. To use a Python-only GIF writer instead, change `writer="imagemagick"` to `writer="pillow"` in the `ani.save(...)` call; Pillow is included in the dependency list.

### Open the notebook

```powershell
.\.venv\Scripts\python.exe -m jupyterlab cherenkov.ipynb
```

Run cells in order. The notebook's second cell downloads an Excel export from the Google Sheet identified in its source and requires internet access, access to that spreadsheet, and `openpyxl`. The standalone scripts calculate their own spectra and do not require that spreadsheet.

The notebook reads a column labeled `Simulation Values (0.7c)`, while a later plot title labels that series as 0.75c. Treat those legacy labels cautiously until the spreadsheet values are checked. The saved figures and GIF are historical outputs, rather than a guarantee of results regenerated with current dependencies.

## Assumptions and interpretation

- A singly charged particle moves at a prescribed, constant speed through a homogeneous medium represented by the existing water dispersion model.
- The calculation describes emitted spectral energy per path length. It does not integrate a finite track or account for slowing down and secondary-particle production.
- Absorption, scattering, propagation distance, optical boundaries, photosensor response, and detector geometry are omitted. Consequently, the plotted chromaticity describes the modeled source spectrum rather than the light recorded by a particular detector.
- The model predicts a standardized chromaticity, not absolute brightness, dark-adapted vision, or a calibrated screen RGB color.
- The slider includes $\beta=1$ as a limiting value; a massive charged particle cannot reach the vacuum speed of light.
- The missing threshold cutoff and normalization guards are the main obstacles to interpreting the entire slider range physically. The default $\beta=0.9$ case is above the visible-water threshold.

## References

- [Super-Kamiokande: detector and neutrino-detection context](https://www-sk.icrr.u-tokyo.ac.jp/en/sk/about/detector/).
- [Particle Data Group: Passage of Particles Through Matter](https://pdg.lbl.gov/2022/reviews/rpp2022-rev-passage-particles-matter.pdf), Cherenkov-radiation section.
- [Daimon and Masumura, Applied Optics 46, 3811–3820 (2007)](https://doi.org/10.1364/AO.46.003811), background on measured water dispersion; exact provenance of this project's coefficients remains unverified.
- [CIE: color-matching functions of the CIE 1931 standard colorimetric observer](https://cie.co.at/datatable/cie-1931-colour-matching-functions-2-degree-observer), DOI 10.25039/CIE.DS.xvudnb9b.
- [NumPy: trapezoidal numerical integration](https://numpy.org/doc/stable/reference/generated/numpy.trapezoid.html).
