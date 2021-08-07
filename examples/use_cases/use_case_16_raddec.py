"""
Use cases for passing RA, DEC to MulensData, Model, and Event.
Also plots ground-based and satellite-based data, models, and trajectories.

Based on OGLE-2014-BLG-0939 from Yee et al. 2015 ApJ 802, 76
(satellite parallax measurement)
"""
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt
import os

import MulensModel as mm


data_dir = os.path.join(mm.DATA_PATH, 'photometry_files', 'OB140939')
ephemeris_dir = os.path.join(mm.DATA_PATH, 'ephemeris_files')

ra = '17:47:12.25'
dec = '-21:22:58.2'
ra_dec = ra + " " + dec

# Specifying coordinates to calculate HJD from JD
data_1 = mm.MulensData(
    file_name=os.path.join(data_dir, 'ob140939_OGLE.dat'))

data_2 = mm.MulensData(
    file_name=os.path.join(data_dir, 'ob140939_OGLE.dat'))

coords = SkyCoord(ra_dec, unit=(u.hourangle, u.deg))
data_3 = mm.MulensData(
    file_name=os.path.join(data_dir, 'ob140939_OGLE.dat'))

# Specifying coordinates to calculate a model with parallax
t_0 = 2456836.22
u_0 = 0.922
t_E = 22.87
pi_E_N = -0.248
pi_E_E = 0.234

ground_model = mm.Model(
    {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'pi_E': [pi_E_N, pi_E_E]})
ground_model.coords = '17:47:12.25 -21:22:58.2'
space_model = mm.Model(
    {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'pi_E': [pi_E_N, pi_E_E]},
    ra=ra, dec=dec,
    ephemerides_file=os.path.join(
        ephemeris_dir, 'Spitzer_ephemeris_01.dat'))

# Access Galactic and ecliptic coordinates:
print('l {0}'.format(ground_model.coords.galactic_l))
print('b {0}'.format(ground_model.coords.galactic_b))
print('ecliptic lat {0}'.format(ground_model.coords.ecliptic_lat))
print('ecliptic lon {0}'.format(ground_model.coords.ecliptic_lon))

plt.figure()
ground_model.plot_magnification(label='ground', color='black')
space_model.plot_magnification(label='space', color='red')
plt.title('OB140939 Models with Parallax')
plt.legend()

# Specifying coordinates for an event
ground_data = mm.MulensData(
    file_name=os.path.join(data_dir, 'ob140939_OGLE.dat'),
    plot_properties={'label': 'OGLE', 'color': 'black'})
space_data = mm.MulensData(
    file_name=os.path.join(data_dir, 'ob140939_Spitzer.dat'),
    plot_properties={'label': 'Spitzer', 'color': 'red'},
    ephemerides_file=os.path.join(
        ephemeris_dir, 'Spitzer_ephemeris_01.dat'))

model_params = mm.ModelParameters(
        {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'pi_E_N': pi_E_N,
         'pi_E_E': pi_E_E})
event = mm.Event(datasets=[ground_data, space_data],
                 model=mm.Model(parameters=model_params),
                 coords=ra_dec)

plt.figure()
event.plot_model(color=ground_data.plot_properties['color'])
event.plot_data()

(fs_ogle, fb_ogle) = event.get_ref_fluxes()
# JCY: IMO, data_ref should be a fixed property of event. Otherwise, the user
# could forget which dataset was set as the reference and request values for
# the wrong dataset. Now that there are ways to access fluxes for
# an arbitrary dataset, ref_fluxes should *only* refer to
# the actual reference dataset.
space_model.plot_lc(
    source_flux=fs_ogle, blend_flux=fb_ogle,
    color=space_data.plot_properties['color'])

plt.title('OB140939 Models with Data')
plt.legend(loc='best')
plt.xlim(2456780., 2456880.)
plt.ylim(15.4, 14.6)

plt.figure()
ground_model.parameters.parameters['rho'] = 0.02
ground_model.plot_trajectory(color=ground_data.plot_properties['color'])
ground_model.plot_source(
    ground_data.time, color=ground_data.plot_properties['color'])

space_model.parameters.parameters['rho'] = 0.015
space_model.plot_trajectory(color=space_data.plot_properties['color'])
space_model.plot_source(
    space_data.time, color=space_data.plot_properties['color'])

plt.title('Trajectory as Seen from Ground and Space')
plt.axis('equal')
plt.xlim(-1., 1.)
plt.show()
