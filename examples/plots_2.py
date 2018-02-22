from matplotlib import pyplot
import os

from MulensModel import Model, SatelliteSkyCoord, MODULE_PATH


# Define model parameters.
params = {'t_0': 2456900, 'u_0': 0.2, 't_E': 50.}
params_pi_E = {'pi_E_N': 0.35, 'pi_E_E': 0.5}
params_planet = {'rho': 0.002, 's': 1.5, 'q': 0.001, 'alpha': 348.1}
ra_dec = '18:00:00.00 -28:30:00.0'

# Set models and satellite settings.
model_pspl = Model(params)
model_planet = Model({**params, **params_planet})
model_planet.set_magnification_methods([2456937, 'VBBL', 2456945]) # Calculate
# finite source magnification using VBBL method for this range of dates.
model_parallax = Model({**params, **params_pi_E}, coords=ra_dec)
satellite = SatelliteSkyCoord(os.path.join(MODULE_PATH, 'data',
    'Spitzer_ephemeris_01.dat')) # This file gives the Spitzer ephemeris 
    # and is part of MulensModel package.

# Plot the magnification curves.
plot_kwargs = {'subtract_2450000': True, 'lw': 2.}
model_planet.plot_magnification(label='planetary', **plot_kwargs)
model_pspl.plot_magnification(label='PSPL', linestyle='--', **plot_kwargs)
model_parallax.parallax(earth_orbital=True)
model_parallax.plot_magnification(label='annual parallax', linestyle='-.', 
    **plot_kwargs)
model_parallax.parallax(satellite=True)
model_parallax.plot_magnification(label='satellite parallax', 
    satellite_skycoord=satellite, **plot_kwargs)

pyplot.legend(loc='best')
pyplot.savefig('figure_2.png')
