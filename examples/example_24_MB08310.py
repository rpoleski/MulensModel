"""
Output information about datasets, model, and event.
Created from example_05_MB08310.py

From `Janczak et al. 2010, ApJ 711, 731
<https://ui.adsabs.harvard.edu/abs/2010ApJ...711..731J/abstract>`_.

"""
import glob
import os

import MulensModel as mm


# Read in MB08310 data files (see data/MB08310) as MulensData objects.
# Grabbing all data files in the MB08310 folder
files = glob.glob(os.path.join(mm.DATA_PATH, "photometry_files",
                               "MB08310", "*.tbl"))

# Read in the data
datasets_default = []
for file_ in sorted(files):
    data = mm.MulensData(file_name=file_, comments=["\\", "|"])
    datasets_default.append(data)

# Print dataset information
for data in datasets_default:
    print(data)

# Define basic point lens model
t_0 = 2454656.39975
u_0 = 0.00300
t_E = 11.14
t_star = 0.05487
plens_model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 't_star': t_star})
method = 'finite_source_uniform_Gould94'
plens_model.set_magnification_methods([t_0-2.*t_star, method, t_0+2.*t_star])

# Print model information
print('\n')
print(plens_model)
print('\n')

# Combine the data and model into an event
event_default = mm.Event(datasets=datasets_default, model=plens_model)
event_default.data_ref = 6

# Print event information
print(event_default)

