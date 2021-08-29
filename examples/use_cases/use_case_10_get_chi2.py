"""
Calculate the chi2 between a model and some data.
"""
import matplotlib.pyplot as plt
from astropy import units as u
import os
import numpy as np

import MulensModel as mm


# Load the data
data = mm.MulensData(
    file_name=os.path.join(
        mm.DATA_PATH, 'photometry_files', 'OB08092', 'phot_ob08092_O4.dat'),
    add_2450000=True)

# Define the model
model = mm.Model(
    {'t_0': 2455379.571, 'u_0': 0.523, 't_E': 17.94})

# Combine the model and the data
event = mm.Event(datasets=data, model=model)
print(event.get_chi2())

# Get magnifications for selected dates
model_times = np.arange(2455300, 2455400, 10.)
print(model_times)
print(event.model.get_magnification(model_times))

event.plot_data()
event.plot_model()

# Get fluxes for all datasets
fmt = "dataset {:}: F_s = {:.3f} F_b = {:.3f}"
for (i, dataset) in enumerate(event.datasets):
    (f_source, f_blend) = event.get_flux_for_dataset(dataset)
    print(fmt.format(i, f_source[0], f_blend))

plt.show()
