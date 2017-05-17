import numpy as np

import MulensModel

#Spitzer Example

#Import Data
spitzer_data = MulensModel.MulensData(
    file_name="../../data/ob151100_Spitzer_ref_v1.dat", satellite="Spitzer", 
    ephemrides_file="../../data/Spitzer_ephemrides_01.dat")

#Create Model
model = MulensModel.Model()
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



