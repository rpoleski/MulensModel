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

coordinates = "17:54:14.53 −34:46:40.99"

filters = {'Auck': 'R', 'Bron': 'R', 'Canopus': 'I', 'CTIO_H': 'H',
           'CTIO_I': 'I', 'Danish': 'I', 'MOA': 'R'}

# Read in the data
datasets_default = []
for file_ in sorted(files):
    file_elements = os.path.basename(file_).split('_')
    name = file_elements[0]
    if name == 'CTIO':
        name += '_'
        name += file_elements[1]

    data = mm.MulensData(
        file_name=file_, comments=["\\", "|"], bandpass=filters[name])
    if name == 'MOA':
        data.scale_errorbars(1.6, 0.001)

    datasets_default.append(data)


print("Printing datasets:")
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
plens_model.set_limb_coeff_u('I', 0.547)
plens_model.set_limb_coeff_u('V', 0.714)
plens_model.set_limb_coeff_u('R', 0.633)
plens_model.set_limb_coeff_u('H', 0.368)

print('\nPrinting model:')
print(plens_model)

# Combine the data and model into an event
event_default = mm.Event(
    datasets=datasets_default, model=plens_model, coords=coordinates)
event_default.data_ref = 6

print("\nPrinting event:")
print(event_default)
