import numpy as np
import os

import MulensModel

#Spitzer Example
raise NotImplementedError('satellite keyword for MulensData not supported')

#Import Data
file_dir = os.path.join(MulensModel.MODULE_PATH, "data")
spitzer_data = MulensModel.MulensData(
    file_name=os.path.join(file_dir, "ob151100_Spitzer_ref_v1.dat"), 
    satellite="Spitzer", #this keyword does not work.
    ephemrides_file=os.path.join(file_dir, "Spitzer_ephemrides_01.dat"))

#Create Model
model = MulensModel.Model(
    {'t_0': 0., 'u_0': 0.1, 't_E': 1.})
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

