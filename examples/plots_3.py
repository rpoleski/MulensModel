'''
This example shows OGLE-2003-BLG-235/MOA-2003-BLG-53, the first microlensing planet. 
See Bond et al. 2004 (http://adsabs.harvard.edu/abs/2004ApJ...606L.155B). 
The data were downloaded from the NASA Exoplanet Archive:
https://exoplanetarchive.ipac.caltech.edu/cgi-bin/DisplayOverview/
nph-DisplayOverview?objname=OGLE-2003-BLG-235L+b&type=CONFIRMED_PLANET
'''
import matplotlib.pyplot as plt
import os

import MulensModel as MM


# Import data
data_dir = os.path.join(MM.MODULE_PATH, 'data', 'OB03235')
OGLE_data = MM.MulensData(
    file_name=os.path.join(data_dir, 'OB03235_OGLE.tbl.txt'),
    comments=['\\', '|'])
MOA_data = MM.MulensData(
    file_name=os.path.join(data_dir, 'OB03235_MOA.tbl.txt'),
    comments=['\\', '|'], phot_fmt='flux')

# Define a model with 2-body lens:
my_1S2L_model = MM.Model({'t_0': 2452848.06, 'u_0': 0.1317, 't_E': 61.5,
    'rho': 0.00096, 'q': 0.0039, 's': 1.120, 'alpha': 223.72})

# Since rho is set, define a time range and method to apply 
# finite source effects:
my_1S2L_model.set_magnification_methods([2452833., 'VBBL', 2452845.])

# Combine the data and model into an Event:
my_event = MM.Event(datasets=[MOA_data, OGLE_data], model=my_1S2L_model)

# Make the plot:
t_range = [2452810., 2452890.]
my_event.plot_data(subtract_2450000=True, data_ref=1,
    label_list=['MOA', 'OGLE'], s=5)
my_event.plot_model(subtract_2450000=True, data_ref=1, t_range=t_range,
    n_epochs=4000, color='black')

plt.legend(loc='best')
plt.xlim(t_range[0]-2450000., t_range[1]-2450000.)
plt.ylim(19., 16.7)
plt.savefig('figure_3.png')
