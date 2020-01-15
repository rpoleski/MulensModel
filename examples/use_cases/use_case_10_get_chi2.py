"""
Calculate the chi2 between a model and some data.
"""
import matplotlib.pyplot as pl
from astropy import units as u
import os
import numpy as np

import MulensModel as mm


# Load the data
data = mm.MulensData(
    file_name=os.path.join(
        mm.DATA_PATH, 'photometry_files', 'OB08092', 'phot_ob08092_O4.dat'))

# Define the model
model = mm.Model(
    {'t_0': 2457518.902, 'u_0': 0.590, 't_E': 133.34*u.day})

# Combine the model and the data
event = mm.Event(datasets=data, model=model)
print(event.get_chi2())

# Get magnifications for selected dates
model_times = np.arange(2457200, 2457800, 100.)
print(model_times)
print(event.model.magnification(model_times))

# Get fluxes for all datasets
fmt = "dataset {:}: F_s = {:.3f} F_b = {:.3f}"
for (i, dataset) in enumerate(event.datasets):
    (f_source, f_blend) = event.fits[dataset].source_flux, event.fits[datasets].blend_flux
    print(fmt.format(i, f_source[0], f_blend))
