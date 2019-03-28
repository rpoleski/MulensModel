import numpy as np
import os
#Spitzer Example

import MulensModel


raise NotImplementedError('satellite keyword for MulensData not supported')

#Import Data
file_dir = os.path.join(MulensModel.MODULE_PATH, "data")
file_name = os.path.join(file_dir, "photometry_files", "OB140939",
                         "ob140939_Spitzer.dat")
spitzer_data = MulensModel.MulensData(
    file_name=file_name,
    satellite="Spitzer", #this keyword does not work.
    ephemerides_file=os.path.join(
        file_dir, "ephemeris_files", "Spitzer_ephemeris_01.dat"))

#Create Model
model = MulensModel.Model({'t_0': 0., 'u_0': 0.1, 't_E': 1.})
print('Default Parallax Settings:')
print(model.parallax())
model.parallax(satellite=True, earth_orbital=False, topocentric=False)
print('User Parallax Settings:')
print(model.parallax())

#Create Event
spitzer_event = MulensModel.Event(
    datasets=spitzer_data, model=model, coords="17:50:00 -29:00:05")
print('Event Parallax Settings:')
print(spitzer_event.model.parallax())

