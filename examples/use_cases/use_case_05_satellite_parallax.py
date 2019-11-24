"""
This is Spitzer use case.
"""
import numpy as np
import os

import MulensModel as mm


raise NotImplementedError('satellite keyword for MulensData not supported')

# Import Data
file_name = os.path.join(mm.DATA_PATH, "photometry_files", "OB140939",
                         "ob140939_Spitzer.dat")
spitzer_data = mm.MulensData(
    file_name=file_name,
    satellite="Spitzer",  # this keyword does not work.
    ephemerides_file=os.path.join(
        mm.DATA_PATH, "ephemeris_files", "Spitzer_ephemeris_01.dat"))

# Create Model
model = mm.Model({'t_0': 0., 'u_0': 0.1, 't_E': 1.})
print('Default Parallax Settings:')
print(model.parallax())
model.parallax(satellite=True, earth_orbital=False, topocentric=False)
print('User Parallax Settings:')
print(model.parallax())

# Create Event
spitzer_event = mm.Event(
    datasets=spitzer_data, model=model, coords="17:50:00 -29:00:05")
print('Event Parallax Settings:')
print(spitzer_event.model.parallax())
