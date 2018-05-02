import matplotlib.pyplot as pl
from astropy import units as u 
import os
import numpy as np

import MulensModel

"""
Calculate the chi2 between a model and some data.
"""

# Load the data
data = MulensModel.MulensData(
    file_name=os.path.join(
        MulensModel.MODULE_PATH, 'data/photometry_files', 'phot_ob160023.dat'))

# Define the model
model = MulensModel.Model(
    {'t_0': 2457518.902, 'u_0': 0.590, 't_E': 133.34*u.day})

# Combine the model and the data
event = MulensModel.Event(datasets=data, model=model) 
print(event.get_chi2())

# Get magnifications for selected dates
model_times = np.arange(2457200, 2457800, 100.)
print(model_times)
print(event.model.magnification(model_times)) 

# Get fluxes for all datasets
for (i, dataset) in enumerate(event.datasets):
    (f_source, f_blend) = event.model.get_ref_fluxes(data_ref=dataset)
    print("dataset {:}: F_s = {:.3f} F_b = {:.3f}".format(i, f_source[0], f_blend))

