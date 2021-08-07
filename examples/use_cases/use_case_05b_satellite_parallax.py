"""
This is a use case for Spitzer + Kepler. The goal is to partially reproduce
Figure 3 from Zhu et al. 2017 ApJL 849 L31. = MOA-2016-BLG-290

"""
import numpy as np
import matplotlib.pyplot as plt

import MulensModel as mm


raise NotImplementedError('satellite keyword for MulensData not supported')
# Needs a new trial event, preferably a point lens from Zhu?.

# Import Data
ogle_data = mm.MulensData(
    file_name='OGLE File Name',
    plot_properties={'label': 'OGLE', 'color': 'black'})
moa_data = mm.MulensData(
    file_name='MOA File Name',
    plot_properties={'label': 'MOA', 'color': 'orange'})
spitzer_data = mm.MulensData(
    file_name='Spitzer File Name',
    satellite='Spitzer',  # this keyword does not work.
    plot_properties={'label': 'Spitzer', 'color': 'red'})
kepler_data = mm.MulensData(
    file_name='Kepler File Name', satellite='Kepler',
    plot_properties={'label': 'Kepler', 'color': 'blue'})

# Ephemrides Files
ephemerides_files = {
    'Spitzer': 'Spitzer Ephemerides File Name',
    'Kepler': 'Kepler Ephemerides File Name'}

# Create Model
model = mm.Model(
    {'t_0': 0., 'u_0': 0.1, 't_E': 1.}, ephemerides_files=ephemerides_files)
print('Default Parallax Settings:')
print(model.parallax())
model.parallax(satellite=True, earth_orbital=False, topocentric=False)
print('User Parallax Settings:')
print(model.parallax())

# Create Event
event = mm.Event(
    datasets=[ogle_data, moa_data, spitzer_data, kepler_data], model=model,
    coords="Coordinates")
print('Event Parallax Settings:')
print(event.model.parallax())

# Make Nice Plots
event.plot_data(subtract_2450000=True)
event.plot_model(color='black', subtract_2450000=True)
event.plot_model(satellite='Spitzer', color='red', subtract_2450000=True)
event.plot_model(satellite='Kepler', color='blue', subtract_2450000=True)
plt.legend()
plt.show()
