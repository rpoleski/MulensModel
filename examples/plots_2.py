"""
Create Figure 2. of Poleski and Yee (2019)
"Modeling microlensing events with MulensModel"
Astronomy and Computing 26, 35
https://ui.adsabs.harvard.edu/abs/2019A&C....26...35P/abstract
https://arxiv.org/abs/1803.01003

Example magnification curves.

"""
from matplotlib import pyplot
import os

from MulensModel import Model, SatelliteSkyCoord, MODULE_PATH


# Define model parameters.
params = {'t_0': 2456900, 'u_0': 0.2, 't_E': 50.}
params_pi_E = {'pi_E_N': 0.35, 'pi_E_E': 0.5}
params_planet = {'rho': 0.002, 's': 1.5, 'q': 0.001, 'alpha': 168.1}
ra_dec = '18:00:00.00 -28:30:00.0'

# Set models and satellite settings.
model_pspl = Model(params)
model_planet = Model({**params, **params_planet})

# Calculate finite source magnification using VBBL method for this
# range of dates:
model_planet.set_magnification_methods([2456937, 'VBBL', 2456945])

# Parallax settings:
model_parallax = Model({**params, **params_pi_E}, coords=ra_dec)
model_parallax.parallax(earth_orbital=True, satellite=True)
satellite = SatelliteSkyCoord(
    os.path.join(
        MODULE_PATH, 'data/ephemeris_files', 'Spitzer_ephemeris_01.dat'))
# This file gives the Spitzer ephemeris and is part of MulensModel package.

# Plot the magnification curves.
plot_kwargs = {'subtract_2450000': True, 'lw': 2.}
pyplot.figure(figsize=(8, 8))
pyplot.axes([0.1, 0.43, 0.85, 0.55])
model_planet.plot_magnification(label='planetary', **plot_kwargs)
model_parallax.plot_magnification(
    label='annual parallax', linestyle='-.', **plot_kwargs)
model_pspl.plot_magnification(label='PSPL', linestyle='--', **plot_kwargs)
model_parallax.plot_magnification(
    label='satellite parallax', satellite_skycoord=satellite, **plot_kwargs)
pyplot.legend(loc='best')

pyplot.axes([0.1, 0.07, 0.85, 0.25])  # Lower panel starts here.
model_planet.plot_trajectory(caustics=True)
pyplot.xlim(-1.52, 1.61)
pyplot.xlabel(r"$\theta_x$")
pyplot.ylabel(r"$\theta_y$")
pyplot.savefig('figure_2.png')
