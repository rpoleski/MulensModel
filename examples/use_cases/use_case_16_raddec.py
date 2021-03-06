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
    file_name=os.path.join(data_dir, 'ob140939_OGLE.dat'),
    coords=ra_dec)

data_2 = mm.MulensData(
    file_name=os.path.join(data_dir, 'ob140939_OGLE.dat'),
    ra=ra, dec=dec)

coords = SkyCoord(ra_dec, unit=(u.hourangle, u.deg))
data_3 = mm.MulensData(
    file_name=os.path.join(data_dir, 'ob140939_OGLE.dat'), coords=coords)

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
ground_model.plot_magnification(label='ground')
space_model.plot_magnification(label='space')
plt.title('OB140939 Models with Parallax')
plt.legend()

# Specifying coordinates for an event
ground_data = mm.MulensData(
    file_name=os.path.join(data_dir, 'ob140939_OGLE.dat'))
space_data = mm.MulensData(
    file_name=os.path.join(data_dir, 'ob140939_Spitzer.dat'),
    ephemerides_file=os.path.join(
        ephemeris_dir, 'Spitzer_ephemeris_01.dat'))

model_params = mm.ModelParameters(
        {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'pi_E_N': pi_E_N,
         'pi_E_E': pi_E_E})
event = mm.Event(datasets=[ground_data, space_data],
                 model=mm.Model(parameters=model_params),
                 coords=ra_dec)

plt.figure()
event.plot_model()
event.plot_data(label_list=['OGLE', 'Spitzer'])

(fs_ogle, fb_ogle) = event.get_ref_fluxes(data_ref=event.datasets[0])
space_model.plot_lc(f_source=fs_ogle, f_blend=fb_ogle)

plt.title('OB140939 Models with Data')
plt.legend(loc='best')
plt.xlim(2456780., 2456880.)
plt.ylim(15.4, 14.6)

plt.figure()
ground_model.plot_trajectory()
space_model.plot_trajectory()

ground_model.set_datasets([ground_data])
ground_model.parameters.parameters['rho'] = 0.02
ground_model.plot_source_for_datasets()

space_model.set_datasets([space_data])
space_model.parameters.parameters['rho'] = 0.015
space_model.plot_source_for_datasets()

plt.title('Trajectory as Seen from Ground and Space')
plt.axis('equal')
plt.xlim(-1., 1.)
plt.show()
