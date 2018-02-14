import matplotlib.pyplot as pl
from astropy import units as u 
import os

import MulensModel

"""
Calculate the chi2 between a model and some data.
"""

# Load the data
data = MulensModel.MulensData(
    file_name=os.path.join(
        MulensModel.MODULE_PATH, 'data', 'phot_ob160023.dat'))

# Define the model
model = MulensModel.Model(
    {'t_0': 2457518.902, 'u_0': 0.590, 't_E': 133.34*u.day})

# Combine the model and the data
event = MulensModel.Event(datasets=data, model=model) 
print(event.get_chi2())
