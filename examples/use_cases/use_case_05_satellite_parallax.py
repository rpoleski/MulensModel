import numpy as np

import MulensModel

#Basic Example
example_time = np.arange(7000., 7600., 0.1)
example_magnitude = np.zeros(example_time.size)
example_err = np.zeros(example_time.size.size)

example_dataset = MulensModel.MulensData(
    np.array(example_time, example_magnitude, example_err), satellite='K2')

print(example_dataset.satellite)

example_model = MulensModel.Model()
example_model.parallax(satellite=True, earth_orbital=False, topocentric=False) 
#default is all True

#Spitzer Example
spitzer_data = MulensModel.MulensData(
    file_name="Spitzer_data.dat", date_fmt="HJD", satellite="Spitzer", 
    ephemrides_file="Spitzer_ephemrides.dat")

model = MulensModel.Model()

spitzer_event = MulensModel.Event(
    data=spitzer_data, model=model, coords="17:50:00 -29:00:05")




